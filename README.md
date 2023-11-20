# Melon-supplementary
Instruction on building the NCBI/GTBD database of Melon.
- [Prerequisite](#prerequisite)
   * [Step 1: Install necessary packages](#step-1-install-necessary-packages)
   * [Step 2: Download protein sequences from https://ftp.ncbi.nlm.nih.gov/blast/db/](#step-2-download-protein-sequences-from-httpsftpncbinlmnihgovblastdb)
   * [Step 3: Download taxonomy files from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/](#step-3-download-taxonomy-files-from-httpsftpncbinlmnihgovpubtaxonomy)
   * [Step 4: Download profile HMMs from https://www.genome.jp/ftp/db/kofam/](#step-4-download-profile-hmms-from-httpswwwgenomejpftpdbkofam)
   * [Step 5: Collect NCBI assemblies](#step-5-collect-ncbi-assemblies)
   * [Step 5 (alternative): Collect GTDB assemblies](#step-5-alternative-collect-gtdb-assemblies)
   * [Step 6: Download assemblies and generate an accession2assembly mapping](#step-6-download-assemblies-and-generate-an-accession2assembly-mapping)
- [Construction of the protein database](#construction-of-the-protein-database)
   * [Step 1: Extract protein sequences from BLAST databases](#step-1-extract-protein-sequences-from-blast-databases)
   * [Step 2: Re-annotating protein sequences](#step-2-re-annotating-protein-sequences)
   * [Step 3: Cluster to reduce redundancy](#step-3-cluster-to-reduce-redundancy)
- [Construction of the nucleotide database](#construction-of-the-nucleotide-database)
   * [Step 1: Map assemblies to the protein databases](#step-1-map-assemblies-to-the-protein-databases)
   * [Step 2: Parse output files and extract sequences](#step-2-parse-output-files-and-extract-sequences)
   * [Step 3: Cluster to remove duplicated sequences](#step-3-cluster-to-remove-duplicated-sequences)
- [Compress](#compress)

## Prerequisite
### Step 1: Install necessary packages
> [!NOTE]
> Create a new environment using `mamba create` or `conda create` if necessary.

```bash
mamba install -c bioconda -c conda-forge 'taxonkit>=0.14.3' 'seqkit>=2.5.1' 'hmmer>=3.3.2' 'mmseqs2>=v14.7e284' 'blast>=2.14.0' 'diamond>=2.1.8' 'tqdm>=4.65.0' 'pandas>=2.0.3'
```

### Step 2: Download protein sequences from https://ftp.ncbi.nlm.nih.gov/blast/db/
> [!NOTE]
> Use `ascp` from IBM Aspera or `update_blastdb.pl` from `BLAST+` to speed up if needed.

```bash
mkdir -p proteins
cd proteins

## download all splits in a loop
for folder in env_nr nr # refseq_protein
do
    mkdir -p $folder
    cd $folder
    curl ftp://ftp.ncbi.nlm.nih.gov/blast/db/ \
        | grep '[^ ]*.gz$' -o \
        | grep ^$folder \
        | xargs -P 64 -I {} wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/blast/db/{}

    ## decompress all files
    for file in *.tar.gz; do tar -xvf $file; done
    rm -rf *.tar.gz
    cd ..
done
cd ..
```

### Step 3: Download taxonomy files from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
> [!WARNING]  
> Please backup `$HOME/.taxonkit` before executing the following code block otherwise it will be overwritten.

```bash
mkdir -p taxonomy
cd taxonomy

## download taxonomy dump and protein accessions
wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz

## extract bacteria/archaea taxids using taxonkit
tar -xvf taxdump.tar.gz

mkdir -p $HOME/.taxonkit 
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit

taxonkit list --ids 2 --indent '' > bacteria.id
taxonkit list --ids 2157 --indent '' > archaea.id

cd ..
```

### Step 4: Download profile HMMs from https://www.genome.jp/ftp/db/kofam/
> [!NOTE]
> We validated ver. 2023-04-01 using RefSeq proteins. Using other versions, including the latest version is also OK but some profile HMMs may not function well.

```bash
mkdir -p kegg
cd kegg

wget -q --show-progress https://www.genome.jp/ftp/db/kofam/archives/2023-04-01/ko_list.gz
wget -q --show-progress https://www.genome.jp/ftp/db/kofam/archives/2023-04-01/profiles.tar.gz
gzip -d ko_list.gz
tar -xvf profiles.tar.gz

## collect ribosomal protein ko
python -c "
import requests

## https://www.genome.jp/kegg-bin/get_htext#A1
response = requests.get('https://www.genome.jp/kegg-bin/download_htext?htext=ko03011.keg&format=htext&filedir=')

## https://www.genome.jp/kegg/annotation/br01610.html
with open('prokaryote.subset.id', 'w') as w:
    for line in response.text.rstrip().split('\n'):
        if line:
            ls = line.split()
            if ls[0] == 'B':
                if ls[1] in {'Bacteria', 'Archaea'}:
                    kingdom = ls[1].lower()
                else:
                    kingdom = None

            if ls[0] == 'D' and kingdom is not None:
                if ls[1] not in {'K19032', 'K19033'}: # skip these two not in <Ribosomal protein gene clusters in prokaryotes>
                    w.write(ls[1] + '\t' + kingdom + '\n')
"

## copy profiles to a directory
mkdir -p prokaryote.full prokaryote.subset
cut prokaryote.subset.id -f 1 | uniq | xargs -P 64 -I {} cp profiles/{}.hmm prokaryote.subset/{}.hmm
cut profiles/prokaryote.hal -f 1 | uniq | xargs -P 64 -I {} cp profiles/{} prokaryote.full/{}
cd ..
```

### Step 5: Collect NCBI assemblies
> [!NOTE]
> If you want to build a GTDB database, please jump to [Step 5 (alternative): Collect GTDB assemblies](#step-5-alternative-collect-gtdb-assemblies).

```bash
mkdir -p assemblies/fna
cd assemblies

## download metadata file
wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

python -c "
import pandas as pd
import subprocess


def get_taxonomy(taxid):
    output = subprocess.run([
        'taxonkit', 'reformat',
        '--taxid-field', '1',
        '--show-lineage-taxids',
        '--fill-miss-rank',
        '--miss-taxid-repl', '0',
        '--miss-rank-repl', 'unclassified',
        '--trim',
        '-f', '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}'],
        input='\n'.join(taxid)+'\n', text=True, capture_output=True, check=True)
    print(output.stderr)

    taxonomy = {}
    for line in output.stdout.rstrip().split('\n'):
        ls = line.rstrip().split('\t')
        taxonomy[int(ls[0])] = ';'.join([ls[i+7] + '|' + ls[i] for i in range(1, len(ls)-7)])

    return taxonomy


## read assemblies
assembly = pd.read_table('assembly_summary_refseq.txt', skiprows=1, low_memory=False).rename({'#assembly_accession': 'assembly'}, axis=1)
assembly = assembly[(assembly.group.isin(['archaea', 'bacteria'])) & ((assembly['ftp_path'] != 'na'))]

## drop unknown species
assembly['taxonomy'] = assembly['species_taxid'].map(get_taxonomy(assembly['species_taxid'].astype(str).unique()))
assembly = assembly[assembly['taxonomy'].str.split(';').str.get(-1) != '0|unclassified']

## create an index file for tracing
assembly['index'] = pd.Categorical(assembly['taxonomy'], ordered=False).codes + 1
assembly[['assembly', 'taxonomy', 'index']].to_csv('assembly2species.tsv', sep='\t', index=False, header=None)

## create a list for wget
assembly['fna'] = assembly['ftp_path'] + '/' + assembly['ftp_path'].str.split('/').str.get(-1) + '_genomic.fna.gz'

for kingdom in ['archaea', 'bacteria']:
    fna = assembly[assembly['group'] == kingdom]['fna'].to_list()
    n = 8 if kingdom == 'archaea' else 256 # split into 8 or 256 chunks so that each contains same number of genomes
    chunks = [fna[i::n] for i in range(n)]

    for i, chunk in enumerate(chunks):
        with open('fna/' + kingdom + '.split.' + str(i) + '.id', 'w') as w:
            w.write('\n'.join(chunk) + '\n')
"
```

### Step 5 (alternative): Collect GTDB assemblies
```bash
mkdir -p assemblies/fna
cd assemblies

## download metadata files
wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq_historical.txt
wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank_historical.txt

wget -q --show-progress https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tar.gz
wget -q --show-progress https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz
tar -xvf ar53_metadata.tar.gz
tar -xvf bac120_metadata.tar.gz

python -c "
import pandas as pd

## merge gtdb with ncbi ftp link
ncbi = pd.concat([
    pd.read_table('assembly_summary_refseq.txt', skiprows=1, low_memory=False),
    pd.read_table('assembly_summary_refseq_historical.txt', skiprows=1, low_memory=False),
    pd.read_table('assembly_summary_genbank.txt', skiprows=1, low_memory=False),
    pd.read_table('assembly_summary_genbank_historical.txt', skiprows=1, low_memory=False),
]).rename({'#assembly_accession': 'assembly'}, axis=1)

gtdb = pd.concat([
    pd.read_table('ar53_metadata_r214.tsv'),
    pd.read_table('bac120_metadata_r214.tsv')
])
gtdb['assembly'] = gtdb.accession.str.split('_', n=1).str.get(-1)

## some may no longer be available
assembly = pd.merge(gtdb, ncbi[['assembly', 'ftp_path']], how='left', on='assembly')
assembly = assembly[(assembly.ftp_path != 'na') & (assembly.ftp_path.notnull())]
assembly['taxonomy'] = assembly['gtdb_taxonomy'].str.replace('[a-z]__', '', regex=True)
assembly['group'] = assembly['taxonomy'].str.split(';').str.get(0).str.lower()

## create an index file for tracing
assembly['index'] = pd.Categorical(assembly['taxonomy'], ordered=False).codes + 1
assembly[['assembly', 'taxonomy', 'index']].to_csv('assembly2species.tsv', sep='\t', index=False, header=None)

## create a list for wget
assembly['fna'] = assembly['ftp_path'] + '/' + assembly['ftp_path'].str.split('/').str.get(-1) + '_genomic.fna.gz'

for kingdom in ['archaea', 'bacteria']:
    fna = assembly[assembly['group'] == kingdom]['fna'].to_list()
    n = 8 if kingdom == 'archaea' else 256 # split into 8 or 256 chunks so that each contains same number of genomes
    chunks = [fna[i::n] for i in range(n)]

    for i, chunk in enumerate(chunks):
        with open('fna/' + kingdom + '.split.' + str(i) + '.id', 'w') as w:
            w.write('\n'.join(chunk) + '\n')
"
```

### Step 6: Download assemblies and generate an accession2assembly mapping
```bash
## download assemblies then cat
cd fna
find *.id | xargs -P 64 -I {} bash -c '
    wget -i ${1} -q --show-progress -P ${1%.id}; \
    find ${1%.id} -name "*.fna.gz" | xargs cat > ${1%.id}.fna.gz' - {}
cd ..

## create a dictionary that maps sequence to assembly
python -c "
import glob
import os
import gzip
import pandas as pd
from tqdm.contrib.concurrent import process_map


def parse_file(file):
    lines = []
    with gzip.open(file, 'rt') as f:
        for line in f:
            if line[0] == '>':
                lines.append([line[1:].split()[0], '_'.join(os.path.basename(file).split('_')[:2])])
    return lines


pd.DataFrame([
    x for y in process_map(parse_file, glob.glob('fna/**/*.fna.gz'), max_workers=64, chunksize=1) for x in y
], columns=['accession', 'assembly']).to_csv('accession2assembly.tsv', sep='\t', index=False, header=None)
"
cd ..
```

## Construction of the protein database
### Step 1: Extract protein sequences from BLAST databases
For `env_nr`, the taxonomy information is not available, we therefore create a diamond database using the full nr, and later assign bacteria/archaea labels for `env_nr` using LCA.

> [!NOTE]
> Extracting sequences and building diamond database can take a while. To speed up, consider running each `blastdbcmd` block in parallel.

```bash
## env_nr
blastdbcmd -db proteins/env_nr/env_nr -entry all > proteins/env_nr.full.fa

## make a diamond database using full nr
blastdbcmd -db proteins/nr/nr -entry all > proteins/nr.full.fa
diamond makedb --in proteins/nr.full.fa --db proteins/nr.full --taxonmap taxonomy/prot.accession2taxid.FULL.gz --taxonnodes taxonomy/nodes.dmp --taxonnames taxonomy/names.dmp
rm -rf proteins/nr.full.fa

## extract sequences from blastdb
for kingdom in archaea bacteria
do
    blastdbcmd -db proteins/nr/nr -target_only -taxidlist taxonomy/${kingdom}.id > proteins/nr.${kingdom}.fa
    blastdbcmd -db proteins/nr/nr -target_only -taxidlist taxonomy/${kingdom}.id -outfmt '%a@%T@%o'| tr '@' '\t' > proteins/nr.${kingdom}.id # accession@taxid@oid
done

## find shared sequences by original id (oid)
python -c "
from collections import defaultdict

oid = defaultdict(set)
with open('proteins/nr.archaea.id') as f:
    for line in f:
        ls = line.rstrip().split('\t')
        oid[ls[-1]].add(ls[0])

accession = set()
with open('proteins/nr.bacteria.id') as f:
    for line in f:
        ls = line.rstrip().split('\t')
        if ls[-1] in oid:
            accession.add(ls[0])
            accession.update(oid.get(ls[-1]))

with open('proteins/nr.shared.id', 'w') as w:
    w.write('\n'.join(accession) + '\n')
"
```

### Step 2: Re-annotating protein sequences
Annotating against the full KEGG profile HMMs will take days, so we first use a subset of prokaryotic profile HMMs (91 ribosomal protein profile HMMs) to get some candidate sequences, then use the full profile HMMs to get the complete annotations.

```bash
mkdir -p prot/out/prokaryote.subset

for file in proteins/*.fa
do
    filename=${file%.fa}
    filename=${filename##*/}

    ls kegg/prokaryote.subset \
    | xargs -P 8 -I {} hmmsearch \
        --domtblout prot/out/prokaryote.subset/$filename.{} \
        -E 2147483647 --domE 2147483647 \
        --noali \
        --cpu 8 \
        kegg/prokaryote.subset/{} $file > /dev/null
done
```

Extract candidate sequences from BLAST databases. Discard also cross-kingdom sequences. 
```bash
python -c "
import glob
import os
import subprocess
from tqdm import tqdm
from collections import defaultdict

## kick out cross-kingdom sequences shared by bacteria and archaea
shared_accession = set()
with open('proteins/nr.shared.id') as f:
    for line in f:
        shared_accession.add(line.rstrip())

## read pre-defined threshold score for each ko
ko = defaultdict(dict)
with open('kegg/ko_list') as f:
    next(f)
    for line in f:
        ls = line.rstrip().split('\t')
        ko[ls[0]]['threshold'], ko[ls[0]]['score_type'] = ls[1:3]

## parse domain output files, record accession of sequences
accession = defaultdict(set)
for file in tqdm(glob.glob('prot/out/prokaryote.subset/*.hmm')):
    filename = os.path.basename(file).rsplit('.', 2)[0]
    with open(file) as f:
        for line in f:
            if line[0] != '#':
                ls = line.rstrip().split(maxsplit=22)
                if ls[0] not in shared_accession:
                    ks = ko.get(ls[3])
                    if (
                        ks['score_type'] == 'full' and float(ls[7]) > float(ks['threshold']) * 0.75 or
                        ks['score_type'] == 'domain' and float(ls[13]) > float(ks['threshold']) * 0.75
                    ):
                        accession[filename].add(ls[0])

for key, val in accession.items():
    with open('prot/out/{}.fa'.format(key), 'w') as w:
        subprocess.run([
            'seqkit', 'grep', '-f', '-', 'proteins/{}.fa'.format(key)
        ], check=True, text=True, input='\n'.join(val) + '\n', stdout=w)
"
```

Run diamond LCA then extract sequences with matched kingdom.
```bash
diamond blastp \
    --db proteins/nr.full.dmnd \
    --query prot/out/env_nr.full.fa \
    --outfmt 102 \
    --top 0 \
    --threads 64 \
    | taxonkit reformat \
        --taxid-field 2 \
        --format "{k}" > prot/out/env_nr.full.id

python -c "
import subprocess
from collections import defaultdict

accession = defaultdict(set)
with open('prot/out/env_nr.full.id') as f:
    for line in f:
        ls = line.split()
        if len(ls) == 4 and ls[-1] in {'Bacteria', 'Archaea'}:
            accession['env_nr.' + ls[-1].lower()].add(ls[0])

for key, val in accession.items():
    with open('prot/out/{}.fa'.format(key), 'w') as w:
        subprocess.run([
            'seqkit', 'grep', '-f', '-', 'prot/out/env_nr.full.fa'
        ], check=True, text=True, input='\n'.join(val) + '\n', stdout=w)
"

rm -rf prot/out/env_nr.full.*
```

Rerun `hmmsearch` but this time with the full set of prokaryotic profile HMMs. 
```bash
mkdir -p prot/out/prokaryote.full

for file in prot/out/*.fa
do
    filename=${file%.fa}
    filename=${filename##*/}

    ls kegg/prokaryote.full \
    | xargs -P 64 -I {} hmmsearch \
        --domtblout prot/out/prokaryote.full/$filename.{} \
        -E 2147483647 --domE 2147483647 \
        --noali \
        --cpu 1 \
        kegg/prokaryote.full/{} $file > /dev/null
done
```

Parse domain output files of `hmmsearch`, select only sequences that pass the pre-defined threshold.
```bash
mkdir -p prot/seq

python -c "
import glob
import pandas as pd
from collections import defaultdict
from tqdm.contrib.concurrent import process_map


def parse_file(file):
    lines = []
    with open(file) as f:
        for line in f:
            if line[0] != '#':
                ls = line.rstrip().split(maxsplit=22)
                if (ks := ko.get(ls[3])) is not None:
                    if (
                        ks['score_type'] == 'full' and float(ls[7]) > float(ks['threshold']) * 0.75 or
                        ks['score_type'] == 'domain' and float(ls[13]) > float(ks['threshold']) * 0.75
                    ):
                        lines.append([ls[0], ls[3], ks['definition'], int(ls[17]), int(ls[18]), int(ls[2]), float(ls[21])])
    return lines


## read pre-defined threshold score for each ko
ko_subset = set()
with open('kegg/prokaryote.subset.id') as f:
    for line in f:
        ko_subset.add(line.split()[0])

ko = defaultdict(dict)
with open('kegg/ko_list') as f:
    next(f)
    for line in f:
        ls = line.rstrip().split('\t')
        if ls[2] != '-':
            ko[ls[0]]['threshold'], ko[ls[0]]['score_type'] = ls[1:3]
            ko[ls[0]]['definition'] = ls[-1].split('protein ')[-1].lower() if ls[0] in ko_subset else 'others'

for source in ['env_nr', 'nr']:
    for kingdom in ['bacteria', 'archaea']:
        domain = pd.DataFrame([
            x for y in process_map(parse_file, glob.glob('prot/out/prokaryote.full/{}.{}.*'.format(source, kingdom)), max_workers=64, chunksize=1) for x in y
        ], columns=['accession', 'knum', 'definition', 'tstart', 'tend', 'tlen', 'acc'])
        domain['tcov'] = (domain['tend'] - domain['tstart'] + 1) / domain['tlen']

        ## remove low quality and invalid sequences
        domain = domain[
            (domain.groupby('accession')['knum'].transform('nunique') == 1) &
            (domain['definition'] != 'others') &
            (domain['tcov'] > 0.75)
        ].sort_values(['accession', 'tcov'], ascending=False).groupby('accession', as_index=False).first()

        ## create a accession2description mapping
        domain['description'] = domain['definition'] + '-tcov:' + domain['tcov'].map('{:.2f}'.format) + '-acc:' + domain['acc'].map('{:.2f}'.format) + '-' + source + '-' + kingdom
        accession2description = domain.set_index('accession')['description'].to_dict()

        ## save filtered sequences
        with open('prot/seq/{}.{}.fa'.format(source, kingdom), 'w') as w, open('prot/out/{}.{}.fa'.format(source, kingdom)) as f:
            for line in f:
                if line[0] == '>':
                    accession = line[1:].split()[0]
                    if (description := accession2description.get(accession)) is not None:
                        w.write('>' + description + '-' + accession + '\n')
                        save = True
                        continue
                    else:
                        save = False
                if save:
                    w.write(line)
"
```

### Step 3: Cluster to reduce redundancy
Combine all extracted files, shuffle, then split for each combination of kingdom and ribosomal protein.
```bash
mkdir -p prot/raw
cat prot/seq/*.fa | seqkit sort | seqkit shuffle -s 0 > prot/raw.fa

python -c "
from collections import defaultdict

sequence = defaultdict(list)
with open('prot/raw.fa') as f:
    for line in f:
        if line[0] == '>':
            ls = line[1:].rstrip().split('-')
            subset = ls[4] + '.' + ls[0].replace('/', '_')
        sequence[subset].append(line)

for key, val in sequence.items():
    with open('prot/raw/{}.fa'.format(key), 'w') as w:
        w.write(''.join(val))
"
```

> [!NOTE]
> If your are working on HPC, please consider submitting a SLURM job for each `mmseqs` to make them run in parallel, and increase `--threads` to reduce the computational time!
```bash
mkdir -p prot/clustered

ls prot/raw/*.fa | sort | xargs -P 8 -I {} bash -c '
    filename=${1%.fa*}; \
    filename=${filename##*/}; \
    echo $filename; \
    mmseqs easy-cluster \
        $1 prot/clustered/$filename prot/clustered/$filename \
        -c 0.98 --min-seq-id 0.98 --cov-mode 0 \
        -s 7.5 --cluster-reassign --threads 8 -v 0 > /dev/null' - {}
```

Merge all clustered files to get the final protein database. Split them into bacteria/archaea for later usage.
```bash
cat prot/clustered/*rep_seq.fasta | seqkit sort | seqkit shuffle -s 0 > prot/clustered.fa

python -c "
import subprocess
from collections import defaultdict

archaea = {'l2', 'l11', 'l10e', 'l15e', 'l18e', 's3ae', 's19e', 's28e'}
bacteria = {'l2', 'l11', 'l20', 'l27', 's2', 's7', 's9', 's16'}

accession = defaultdict(set)
with open('prot/prot.fa', 'w') as w, open('prot/clustered.fa') as f:
    for line in f:
        if line[0] == '>':
            if (gene := line[1:].split('-')[0]) in archaea | bacteria:
                save = True
                if (
                    gene in archaea and (kingdom := line.split('-')[-2]) == 'archaea' or
                    gene in bacteria and (kingdom := line.split('-')[-2]) == 'bacteria'
                ):
                    accession[kingdom].add(line[1:].rstrip())
            else:
                save = False
        if save:
            w.write(line)

for key, val in accession.items():
    with open('prot/prot.{}.fa'.format(key), 'w') as w:
        subprocess.run([
            'seqkit', 'grep', '-f', '-', 'prot/prot.fa'
        ], check=True, text=True, input='\n'.join(val) + '\n', stdout=w)
"
```

## Construction of the nucleotide database
### Step 1: Map assemblies to the protein databases
> [!NOTE]
> This step will take a while to finish. If you are working on HPC, please submit a SLURM job for each `*.fna.gz` file and increase `--threads` to reduce computational time. If not, combining all `*.fna.gz` files into a single one then run `diamond` with tuned `-b -c` may help to speed up.

```bash
mkdir -p nucl/out

ls assemblies/fna/*.fna.gz | sort | xargs -P 8 -I {} bash -c '
    filename=${1%.fna*}; \
    filename=${filename##*/}; \
    echo $filename; \
    if [[ ${filename} == "bacteria"* ]]; then \
        db=prot/prot.bacteria.fa;
    elif [[ ${filename} == "archaea"* ]]; then \
        db=prot/prot.archaea.fa;
    fi; \
    diamond blastx \
        --db $db \
        --query $1 \
        --out nucl/out/${filename}.txt \
        --outfmt 6 qseqid sseqid pident length qlen qstart qend slen sstart send evalue bitscore \
        --evalue 1e-15 --subject-cover 75 \
        --range-culling -F 15 --range-cover 25 \
        --max-hsps 0 --max-target-seqs 25 \
        --threads 8 --quiet' - {}
```

### Step 2: Parse output files and extract sequences
Parse output files of `diamond` by ensuring less than 25% query overlap, keep only a single sequence per qseqid-gene combination. Discard also sequences close to the boundaries.

```bash
mkdir -p nucl/seq

python -c "
import os
import glob
import subprocess
import pandas as pd
from math import floor, ceil
from collections import defaultdict
from tqdm.contrib.concurrent import process_map


def sort_coordinate(start, end):
    return (start - 1, end, '+') if start < end else (end - 1, start, '-')


def compute_overlap(coordinates):
    qstart, qend, sstart, send = coordinates
    overlap = min(qend, send) - max(qstart, sstart)
    return max(overlap / (qend - qstart), overlap / (send - sstart))


def extract_sequence(file):
    filename = os.path.basename(file).split('.fna.gz')[0]
    qrange = defaultdict(set)
    lines = []

    with open('nucl/out/' + filename + '.txt') as f:
        for line in f:
            ls = line.split()
            qseqid, sseqid = ls[0], ls[1]
            qstart, qend, strand = sort_coordinate(int(ls[5]), int(ls[6]))
            if (
                qseqid not in qrange or
                all([compute_overlap((qstart, qend, *x)) < 0.25 for x in qrange.get(qseqid)])
            ):
                qrange[qseqid].add((qstart, qend))
                qlen, slen = int(ls[4]), int(ls[7]) * 3
                qcoord = floor((qstart + qend) / 2) if strand == '+' else ceil((qstart + qend) / 2)

                ## make sure not too close to the boundary
                if qend + slen < qlen and qstart - slen > 0:
                    lines.append([qseqid, qcoord - 5000, qcoord + 5000, sseqid, float(ls[2]), strand, sseqid.split('-')[0]])

    ## keep only one sequence per qseqid + gene
    lines = pd.DataFrame(lines, columns=['qseqid', 'qstart', 'qend', 'sseqid', 'pid', 'strand', 'gene'])
    lines = lines.sort_values(['qseqid', 'gene', 'pid', 'strand'], ascending=False).groupby(['qseqid', 'gene'], as_index=False).first()
    lines.drop('gene', axis=1).to_csv('nucl/out/' + filename + '.bed', sep='\t', header=None, index=False)

    ## extract sequneces from original files
    with open('nucl/seq/' + filename + '.fa', 'w') as f:
        subprocess.run([
            'seqkit', 'subseq', file,
            '--bed', 'nucl/out/' + filename + '.bed'
        ], stdout=f, stderr=subprocess.DEVNULL, check=True)


process_map(extract_sequence, glob.glob('assemblies/fna/*.fna.gz'), max_workers=64, chunksize=1)
"
```

### Step 3: Cluster to remove duplicated sequences
Create a file for each species-gene combo.
```bash
mkdir -p nucl/raw
cat nucl/seq/*.fa | seqkit sort | seqkit shuffle -s 0 > nucl/raw.fa

python -c "
import pandas as pd
from collections import defaultdict

assembly2species = {}
with open('assemblies/assembly2species.tsv') as f:
    for line in f:
        ls = line.rstrip().split('\t')
        assembly2species[ls[0]] = ls[1:]

accession2assembly = {}
with open('assemblies/accession2assembly.tsv') as f:
    for line in f:
        ls = line.rstrip().split('\t')
        accession2assembly[ls[0]] = ls[1]

sequence = defaultdict(list)
accession = set()
with open('nucl/raw.fa') as f:
    for line in f:
        if line[0] == '>':
            ls = line.split()
            accession.add(ls[0][1:].rsplit('_', 1)[0])
            subset = ls[1].split('-')[0].replace('/', '_') + '.' + str(assembly2species.get(accession2assembly.get(ls[0][1:].rsplit('_', 1)[0]))[-1])
        sequence[subset].append(line)

## split the sequences into two parts
with open('nucl/nucl_a.fa', 'w') as w:
    for key, val in sequence.items():
        record = ''.join(val)
        if record.count('>') == 1:
            w.write(record)
        else:
            with open('nucl/raw/{}.fa'.format(key), 'w') as ww:
                ww.write(record)

## save the metadata file for later usage
pd.DataFrame([
    [x] + assembly2species.get(accession2assembly.get(x))[0].split(';') for x in accession
], columns=['accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']).sort_values('accession').to_csv('nucl/raw.tsv', index=False, sep='\t')
"
```

> [!NOTE]
> If your are working on HPC, please consider submitting a SLURM job for each `mmseqs` to make them run in parallel, and increase `--threads` to reduce the computational time! Lower `-P` if you encounter memory issues!
```bash
mkdir -p nucl/clustered

ls nucl/raw/*.fa | sort | xargs -P 64 -I {} bash -c '
    filename=${1%.fa*}; \
    filename=${filename##*/}; \
    echo $filename; \
    mmseqs easy-cluster \
        $1 nucl/clustered/$filename nucl/clustered/$filename \
        -c 0.9998 --min-seq-id 0.9998 --cov-mode 1 \
        -s 7.5 --cluster-reassign --threads 1 -v 0 > /dev/null' - {}
```

Combine all clustered files to get the nucleotide database.
```bash
find nucl/clustered/*rep_seq.fasta -exec cat {} \; > nucl/nucl_b.fa
cat nucl/nucl_a.fa nucl/nucl_b.fa | seqkit sort | seqkit shuffle -s 0 > nucl/nucl.fa

python -c "
import pandas as pd
from collections import defaultdict

sequence = defaultdict(list)
accession = set()
with open('nucl/nucl.fa') as f:
    for line in f:
        if line[0] == '>':
            ls = line.split()
            accession.add(ls[0][1:].rsplit('_', 1)[0])
            subset = ls[1].split('-')[-2] + '.' + ls[1].split('-')[0].replace('/', '_')
        sequence[subset].append(line)

for key, val in sequence.items():
    with open('nucl/nucl.{}.fa'.format(key), 'w') as w:
        w.write(''.join(val))

metadata = pd.read_table('nucl/raw.tsv')
metadata[metadata.accession.isin(accession)].to_csv('nucl/metadata.tsv', index=False, sep='\t')
"
```

## Compress
Get all necessary files into the database.
```bash
mkdir -p database
cp prot/prot.fa nucl/metadata.tsv nucl/nucl.bacteria*.fa nucl/nucl.archaea*.fa database
tar -zcvf database.tar.gz database
```
