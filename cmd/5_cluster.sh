echo "----- 5_cluster.sh -----"
echo "Clustering protein sequences ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p prot/raw prot/clustered
find prot/seq -maxdepth 1 -name '*.fa' | sort | xargs cat > prot/raw.fa

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

ls prot/raw/*.fa | sort | xargs -P 4 -I {} bash -c '
    filename=${1%.fa*};
    filename=${filename##*/};
    mmseqs easy-cluster \
        $1 prot/clustered/$filename prot/clustered/$filename \
        -c 0.95 --min-seq-id 0.95 --cov-mode 0 \
        -s 7.5 --cluster-reassign --threads 8 -v 0 > /dev/null' - {}

cat prot/clustered/*rep_seq.fasta | seqkit --quiet sort | seqkit shuffle --quiet -s 0 > prot/prot.fa
## ----------------------------------------------------------------------------------------------------
