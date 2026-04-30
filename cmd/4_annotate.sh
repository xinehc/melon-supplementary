echo "----- 4_annotate.sh -----"
echo "Annotating protein sequences (subset) ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p prot/out/prokaryote.subset

for file in protein/*.fa
do
    filename=${file%.fa}
    filename=${filename##*/}

    ls profile/prokaryote.subset \
    | xargs -P 8 -I {} hmmsearch \
        --domtblout prot/out/prokaryote.subset/$filename.{} \
        -E 2147483647 --domE 2147483647 \
        --noali \
        --cpu 4 \
        profile/prokaryote.subset/{} $file > /dev/null
done

python -c "
import glob
import os
import subprocess
from tqdm import tqdm
from collections import defaultdict

shared_accession = set()
with open('protein/nr.shared.id') as f:
    for line in f:
        shared_accession.add(line.rstrip())

ko = defaultdict(dict)
with open('profile/ko_list') as f:
    next(f)
    for line in f:
        ls = line.rstrip().split('\t')
        ko[ls[0]]['threshold'], ko[ls[0]]['score_type'] = ls[1:3]

accession = defaultdict(set)
for file in tqdm(glob.glob('prot/out/prokaryote.subset/*.hmm'), leave=False):
    filename = os.path.basename(file).rsplit('.', 2)[0]
    with open(file) as f:
        for line in f:
            if line[0] != '#':
                ls = line.rstrip().split(maxsplit=22)
                if ls[0] not in shared_accession:
                    ks = ko.get(ls[3])
                    if (
                        ks['score_type'] == 'full' and float(ls[7]) > float(ks['threshold']) or
                        ks['score_type'] == 'domain' and float(ls[13]) > float(ks['threshold'])
                    ):
                        accession[filename].add(ls[0])

for key, val in accession.items():
    with open('prot/out/{}.fa'.format(key), 'w') as w:
        subprocess.run([
            'seqkit', 'grep', '--quiet', '-f', '-', 'protein/{}.fa'.format(key)
        ], check=True, text=True, input='\n'.join(val) + '\n', stdout=w)
"
## ----------------------------------------------------------------------------------------------------





echo "Annotating protein sequences (full) ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p prot/out/prokaryote.full prot/seq

for file in prot/out/*.fa
do
    filename=${file%.fa}
    filename=${filename##*/}

    ls profile/prokaryote.full \
    | xargs -P 32 -I {} hmmsearch \
        --domtblout prot/out/prokaryote.full/$filename.{} \
        -E 2147483647 --domE 2147483647 \
        --noali \
        --cpu 1 \
        profile/prokaryote.full/{} $file > /dev/null
done


python -c "
import os
import glob
import pandas as pd
import multiprocessing
multiprocessing.set_start_method('fork', force=True)

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
                        ks['score_type'] == 'full' and float(ls[7]) > float(ks['threshold']) or
                        ks['score_type'] == 'domain' and float(ls[13]) > float(ks['threshold'])
                    ):
                        lines.append([ls[0], ls[3], ks['definition'], int(ls[17]), int(ls[18]), int(ls[2]), float(ls[21])])
    return lines

ko_subset = {os.path.basename(file).split('.hmm')[0] for file in glob.glob('profile/prokaryote.subset/*.hmm')}

ko = defaultdict(dict)
with open('profile/ko_list') as f:
    next(f)
    for line in f:
        ls = line.rstrip().split('\t')
        if ls[2] != '-':
            ko[ls[0]]['threshold'], ko[ls[0]]['score_type'] = ls[1:3]
            ko[ls[0]]['definition'] = ls[-1].split('protein ')[-1].lower() if ls[0] in ko_subset else 'others'

source = 'nr'
for kingdom in ['bacteria', 'archaea']:
    domain = pd.DataFrame([
        x for y in process_map(parse_file, glob.glob('prot/out/prokaryote.full/{}.{}.*'.format(source, kingdom)), max_workers=32, chunksize=1, leave=False) for x in y
    ], columns=['accession', 'knum', 'definition', 'tstart', 'tend', 'tlen', 'acc'])
    domain['tcov'] = (domain['tend'] - domain['tstart'] + 1) / domain['tlen']

    domain = domain[
        (domain.groupby('accession')['knum'].transform('nunique') == 1) &
        (domain['definition'] != 'others') &
        (domain['tcov'] > 0.75)
    ].sort_values(['accession', 'tcov'], ascending=False).groupby('accession', as_index=False).first()

    domain['description'] = domain['definition'] + '-tcov:' + domain['tcov'].map('{:.2f}'.format) + '-acc:' + domain['acc'].map('{:.2f}'.format) + '-' + source + '-' + kingdom
    accession2description = domain.set_index('accession')['description'].to_dict()

    with open('prot/seq/{}.{}.fa'.format(source, kingdom), 'w') as w, open('prot/out/{}.{}.fa'.format(source, kingdom)) as f:
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
