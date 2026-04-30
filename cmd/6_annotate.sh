echo "----- 6_annotate.sh -----"
echo "Annotating marker genes ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p nucl/out

find assembly/fna -maxdepth 1 -name '*.fna.gz' | sort | xargs -P 4 -I {} bash -c '
    filename=${1%.fna*};
    filename=${filename##*/};
    diamond blastx \
        --db prot/prot.fa \
        --query $1 \
        --out nucl/out/${filename}.txt \
        --outfmt 6 qseqid sseqid pident length qlen qstart qend slen sstart send evalue bitscore \
        --evalue 1e-25 --subject-cover 75 \
        --range-culling --frameshift 15 --range-cover 25 \
        --max-hsps 0 --max-target-seqs 25 \
        --threads 8 --quiet
' - {}
## ----------------------------------------------------------------------------------------------------





echo "Extracting marker-containing sequences ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p nucl/seq

python -c "
import os
import glob
import subprocess
import pandas as pd
import multiprocessing
multiprocessing.set_start_method('fork', force=True)

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
                pident = float(ls[2])
                qlen, slen = int(ls[4]), int(ls[7]) * 3
                qcoord = floor((qstart + qend) / 2) if strand == '+' else ceil((qstart + qend) / 2)

                ## make sure not too close to the boundary
                if qend + 2500 < qlen and qstart - 2500 > 0:
                    lines.append([qseqid, qcoord - 5000, qcoord + 5000, sseqid, pident, strand, sseqid.split('-')[0]])

    lines = pd.DataFrame(lines, columns=['qseqid', 'qstart', 'qend', 'sseqid', 'pident', 'strand', 'gene'])
    lines = lines.sort_values(['qseqid', 'gene', 'pident', 'strand'], ascending=False).groupby(['qseqid', 'gene'], as_index=False).first()
    lines.drop('gene', axis=1).to_csv('nucl/seq/' + filename + '.bed', sep='\t', header=None, index=False)

    with open('nucl/seq/' + filename + '.fa', 'w') as f:
        subprocess.run([
            'seqkit', 'subseq', file,
            '--quiet', '--bed', 'nucl/seq/' + filename + '.bed'
        ], stdout=f, stderr=subprocess.DEVNULL, check=True)

process_map(extract_sequence, glob.glob('assembly/fna/*.fna.gz'), max_workers=32, chunksize=1, leave=False)
"
## ----------------------------------------------------------------------------------------------------
