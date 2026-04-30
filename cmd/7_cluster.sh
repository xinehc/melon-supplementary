echo "----- 7_cluster.sh -----"
echo "Splitting marker-containing sequences ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p nucl/raw
find nucl/seq -maxdepth 1 -name '*.fa' | sort | xargs cat > nucl/raw.fa

python -c "
import pandas as pd
import json

from collections import defaultdict

with open('mg.json') as f:
    mg = json.load(f)

assembly2species = {}
with open('assembly/assembly2species.tsv') as f:
    for line in f:
        ls = line.rstrip().split('\t')
        assembly2species[ls[0]] = ls[1:]

accession2assembly = {}
with open('assembly/accession2assembly.tsv') as f:
    for line in f:
        ls = line.rstrip().split('\t')
        accession2assembly[ls[0]] = ls[1]

sequence = defaultdict(list)
accession = set()
with open('nucl/raw.fa') as f:
    for line in f:
        if line[0] == '>':
            ls = line.split()
            gs = ls[1].split('-')
            ac = ls[0][1:].rsplit('_', 1)[0]
            species = assembly2species.get(accession2assembly.get(ac))
            phylum = species[0].split(';')[0].split('|')[-1]
            if (
                phylum == 'Archaea' and gs[4] == 'archaea' and gs[0] in mg['archaea'] or
                phylum == 'Bacteria' and gs[4] == 'bacteria' and gs[0] in mg['bacteria']
            ):
                save = True
                accession.add(ac)
                subset = gs[0].replace('/', '_') + '.' + species[-1]
            else:
                save = False
        if save:
            sequence[subset].append(line)

with open('nucl/nucl_a.fa', 'w') as w:
    for key, val in sequence.items():
        record = ''.join(val)
        if record.count('>') == 1:
            w.write(record)
        else:
            with open('nucl/raw/{}.fa'.format(key), 'w') as ww:
                ww.write(record)

pd.DataFrame([
    [x] + assembly2species.get(accession2assembly.get(x))[0].split(';') for x in accession
], columns=['accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']).sort_values('accession').to_csv('nucl/raw.tsv', index=False, sep='\t')
"
## ----------------------------------------------------------------------------------------------------





echo "Clustering marker-containing sequences ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p nucl/clustered

find nucl/raw -maxdepth 1 -name '*.fa' | xargs -P 8 -I {} bash -c '
    filename=${1%.fa*};
    filename=${filename##*/};
    mkdir -p nucl/clustered/$filename;
    mmseqs easy-cluster \
        $1 nucl/clustered/$filename/$filename nucl/clustered/$filename/$filename \
        -c 0.9995 --min-seq-id 0.9995 --cov-mode 1 \
        -s 7.5 --cluster-reassign --threads 4 -v 0 > /dev/null;
    mv nucl/clustered/$filename/${filename}_rep_seq.fasta nucl/clustered/${filename}_rep_seq.fasta;
    rm -rf nucl/clustered/$filename/;
' - {}

find nucl/clustered -maxdepth 1 -name '*rep_seq.fasta' -exec cat {} \; > nucl/nucl_b.fa
cat nucl/nucl_a.fa nucl/nucl_b.fa | seqkit sort --quiet | seqkit shuffle --quiet -s 0 > nucl/nucl.fa

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
## ----------------------------------------------------------------------------------------------------
