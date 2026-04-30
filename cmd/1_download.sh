echo "----- 1_download.sh -----"
echo "Collecting assembly metadata ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p assembly/fna
wget -qN https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -P assembly
wget -qN https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq_historical.txt -P assembly
wget -qN https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt -P assembly
wget -qN https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank_historical.txt -P assembly

wget -qN https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv.gz -P assembly
wget -qN https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz -P assembly

gzip --force --keep -d assembly/ar53_metadata.tsv.gz
gzip --force --keep -d assembly/bac120_metadata.tsv.gz

python -c "
import pandas as pd

ncbi = pd.concat([
    pd.read_table('assembly/assembly_summary_refseq.txt', skiprows=1, low_memory=False),
    pd.read_table('assembly/assembly_summary_refseq_historical.txt', skiprows=1, low_memory=False),
    pd.read_table('assembly/assembly_summary_genbank.txt', skiprows=1, low_memory=False),
    pd.read_table('assembly/assembly_summary_genbank_historical.txt', skiprows=1, low_memory=False),
]).rename({'#assembly_accession': 'assembly'}, axis=1)

gtdb = pd.concat([
    pd.read_table('assembly/ar53_metadata.tsv', low_memory=False),
    pd.read_table('assembly/bac120_metadata.tsv', low_memory=False)
])
gtdb['assembly'] = gtdb.accession.str.split('_', n=1).str.get(-1)

assembly = pd.merge(gtdb, ncbi[['assembly', 'ftp_path']], how='left', on='assembly').fillna('na')
assembly[assembly.ftp_path == 'na'].to_csv('assembly/deprecated.tsv', index=False, sep='\t')

assembly = assembly[(assembly.ftp_path != 'na') | (assembly.gtdb_representative == 't')]
assembly['taxonomy'] = assembly['gtdb_taxonomy'].str.replace('[a-z]__', '', regex=True)
assembly['group'] = assembly['taxonomy'].str.split(';').str.get(0).str.lower()

assembly['index'] = pd.Categorical(assembly['taxonomy'], ordered=False).codes + 1
assembly[['assembly', 'taxonomy', 'index']].to_csv('assembly/assembly2species.tsv', sep='\t', index=False, header=None)
assembly['fna'] = assembly['ftp_path'] + assembly['ftp_path'].str.split('/').str.get(-2) + '_genomic.fna.gz'

for kingdom in ['archaea', 'bacteria']:
    fna = assembly[assembly['group'] == kingdom]['fna']
    fna = fna[fna.notnull()].to_list()
    n = 8 if kingdom == 'archaea' else 256 # split into 8 or 256 chunks so that each contains same number of genomes
    chunks = [fna[i::n] for i in range(n)]

    for i, chunk in enumerate(chunks):
        with open('assembly/fna/' + kingdom + '.split.' + str(i) + '.id', 'w') as w:
            w.write('\n'.join(chunk) + '\n')

with open('assembly/info.tsv', 'w') as f:
    f.write(f'species\t{assembly.taxonomy.nunique()}\n')
    f.write(f'assemblies\t{assembly.assembly.nunique()}\n')
    f.write(f'deprecated_reps\t{len(assembly[assembly.ftp_path == 'na'])}\n')

with open('assembly/assertion.tsv', 'w') as f:
    f.write(f'{assembly.assembly.nunique()}\n')
"

echo "Downloading GTDB species representatives ..."
## ----------------------------------------------------------------------------------------------------
wget -qN https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz
tar -xf gtdb_genomes_reps.tar.gz

python -c "
import pandas as pd
import glob
import os
import shutil

dset = set(pd.read_table('assembly/deprecated.tsv').assembly)
files = glob.glob('gtdb_genomes_reps*/**/*.fna.gz', recursive=True)
files = [file for file in files if os.path.basename(file).rsplit('_genomic')[0] in dset]
os.makedirs('assembly/fna/deprecated_reps', exist_ok=True)
for file in files:
    shutil.copy(file, f'assembly/fna/deprecated_reps/{os.path.basename(file)}')
"

cat assembly/fna/deprecated_reps/*.fna.gz > assembly/fna/deprecated_reps.fna.gz
## ----------------------------------------------------------------------------------------------------
