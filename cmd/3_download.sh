echo "----- 3_download.sh -----"
echo "Downloading kofam ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p profile/prokaryote.full profile/prokaryote.subset
wget -qN https://www.genome.jp/ftp/db/kofam/archives/2026-01-01/ko_list.gz -P profile
wget -qN https://www.genome.jp/ftp/db/kofam/archives/2026-01-01/profiles.tar.gz -P profile
gzip --force --keep -d profile/ko_list.gz
tar -xf profile/profiles.tar.gz -C profile

python -c "
import shutil
import requests
import json

with open('mg.json') as f:
    mg = json.load(f)

ko_full = set()
with open('profile/profiles/prokaryote.hal') as f:
    for line in f:
        ko_full.add(line.rstrip().split('.hmm')[0])

ko_subset = set()
with open('profile/ko_list') as f:
    next(f)
    for line in f:
        ls = line.rstrip().split('\t')
        gene = ls[-1].split('ribosomal protein ')[-1].lower()
        if ls[2] != '-' and gene in set(mg['archaea'] + mg['bacteria']):
            ko_subset.add(ls[0])

for ko in ko_subset:
    shutil.copy(f'profile/profiles/{ko}.hmm', f'profile/prokaryote.subset/{ko}.hmm')

for ko in ko_full:
    shutil.copy(f'profile/profiles/{ko}.hmm', f'profile/prokaryote.full/{ko}.hmm')
"
## ----------------------------------------------------------------------------------------------------





echo "Downloading nr ..."
## ----------------------------------------------------------------------------------------------------
source=nr
curl -s ftp://ftp.ncbi.nlm.nih.gov/blast/db/ \
    | grep '[^ ]*.gz$' -o \
    | grep ^$source \
    | xargs -P 8 -I {} wget -qN --retry-on-http-error=503 https://ftp.ncbi.nlm.nih.gov/blast/db/{} -P protein/$source

for file in protein/$source/*.tar.gz; do tar -xf $file -C protein/$source; done
rm -rf protein/$source/*.tar.gz

cd protein/nr
blastdbcmd -db nr -target_only -taxids 2157 > ../nr.archaea.fa
blastdbcmd -db nr -target_only -taxids 2157 -outfmt '%a@%o' > ../nr.archaea.id

blastdbcmd -db nr -target_only -taxids 2 > ../nr.bacteria.fa
blastdbcmd -db nr -target_only -taxids 2 -outfmt '%a@%o' > ../nr.bacteria.id
cd ../..

python -c "
from collections import defaultdict

oid = defaultdict(set)
with open('protein/nr.archaea.id') as f:
    for line in f:
        ls = line.rstrip().split('@')
        oid[ls[-1]].add(ls[0])

accession = set()
with open('protein/nr.bacteria.id') as f:
    for line in f:
        ls = line.rstrip().split('@')
        if ls[-1] in oid:
            accession.add(ls[0])
            accession.update(oid.get(ls[-1]))

with open('protein/nr.shared.id', 'w') as w:
    w.write('\n'.join(accession) + '\n')
"
## ----------------------------------------------------------------------------------------------------
