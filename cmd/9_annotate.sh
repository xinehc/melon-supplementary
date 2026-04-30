echo "----- 9_annotate.sh -----"
echo "Annotating PlusPFP sequences ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p plus/out
for file in plus/faa/*.faa
do
    filename=${file%.faa}
    filename=${filename##*/}

    ls profile/prokaryote.subset \
    | xargs -P 8 -I {} hmmsearch \
        --domtblout plus/out/$filename.{} \
        -E 2147483647 --domE 2147483647 \
        --noali \
        --cpu 8 \
        profile/prokaryote.subset/{} $file > /dev/null
done

python -c "
import glob
import pandas as pd
import multiprocessing
multiprocessing.set_start_method('fork', force=True)

from collections import defaultdict
from tqdm.contrib.concurrent import process_map
from tqdm import tqdm

def parse_file(file):
    lines = []
    with open(file) as f:
        for line in f:
            if line[0] != '#':
                ls = line.rstrip().split(maxsplit=22)
                if (ks := ko.get(ls[3])) is not None:
                    if (
                        ks['score_type'] == 'full' and float(ls[7]) > float(ks['threshold']) or
                        ks['score_type'] == 'domain' and float(ls[13]) > float(ks['threshold'])
                    ):
                        lines.append([ls[0], ls[3], ks['definition'], int(ls[17]), int(ls[18]), int(ls[2]), float(ls[21])] + [file.split('/')[-1].split('.')[0]])
    return lines

ko = defaultdict(dict)
with open('profile/ko_list') as f:
    next(f)
    for line in f:
        ls = line.rstrip().split('\t')
        if ls[2] != '-':
            ko[ls[0]]['threshold'], ko[ls[0]]['score_type'] = ls[1:3]
            ko[ls[0]]['definition'] = ls[-1].split('protein ')[-1].lower()

domain = pd.DataFrame([
    x for y in process_map(parse_file, glob.glob(f'plus/out/*.hmm'), max_workers=32, chunksize=1, leave=False) for x in y
], columns=['accession', 'knum', 'definition', 'tstart', 'tend', 'tlen', 'acc', 'taxonomy'])
domain['tcov'] = (domain['tend'] - domain['tstart'] + 1) / domain['tlen']

domain = domain[domain['tcov'] > 0.25].sort_values(['accession', 'tcov'], ascending=False).groupby('accession', as_index=False).first()
domain['description'] = domain['definition'] + '-tcov:' + domain['tcov'].map('{:.2f}'.format) + '-acc:' + domain['acc'].map('{:.2f}'.format) + '-' + 'refseq' + '-' + domain['taxonomy']
accession2description = domain.set_index('accession')['description'].to_dict()

with open('plus/plus.fa', 'w') as w:
    for file in tqdm(sorted(glob.glob('plus/faa/*.faa')), leave=False):
        with open(file) as f:
            for line in f:
                if line[0] == '>':
                    accession = line[1:].split()[0]
                    if (description := accession2description.get(accession)) is not None and ', partial' not in line:
                        w.write('>' + description + '-' + accession + '\n')
                        save = True
                        continue
                    else:
                        save = False
                if save:
                    w.write(line)
"
## ----------------------------------------------------------------------------------------------------






echo "Clustering PlusPFP sequences ..."
## ----------------------------------------------------------------------------------------------------
mmseqs easy-cluster \
    plus/plus.fa plus/plus plus/plus \
    -c 0.95 --min-seq-id 0.95 --cov-mode 0 \
    -s 7.5 --cluster-reassign --threads 32 -v 0 > /dev/null
## ----------------------------------------------------------------------------------------------------





echo "Compressing ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p database
cp nucl/metadata.tsv nucl/nucl.bacteria*.fa nucl/nucl.archaea*.fa database
cat plus/plus_rep_seq.fasta prot/prot.fa | seqkit sort --quiet | seqkit shuffle -s 0 --quiet > database/prot.fa
tar --sort=name -zcf database.tar.gz database
## ----------------------------------------------------------------------------------------------------
