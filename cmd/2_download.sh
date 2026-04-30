echo "----- 2_download.sh -----"
echo "Downloading GTDB assemblies ..."
## ----------------------------------------------------------------------------------------------------
find assembly/fna -maxdepth 1 -name '*.id' | sort | xargs -P 8 -I {} bash -c '
    wget -i ${1} -qN --retry-on-http-error=503 -P ${1%.id};
    find ${1%.id} -maxdepth 1 -name "*.fna.gz" | sort | xargs cat > ${1%.id}.fna.gz;
' - {}

expected_count=$(head -n 1 assembly/assertion.tsv)
file_count=$(find assembly/fna -mindepth 2 -name '*.fna.gz' | wc -l | xargs)
if (( file_count != expected_count )); then
    echo "Assertion failed: Found $file_count files, but expected $expected_count. An error may have occurred; please check whether any assemblies are deprecated or rerun this script."
fi
## ----------------------------------------------------------------------------------------------------





echo "Creating accession to assembly mapping ..."
## ----------------------------------------------------------------------------------------------------
python -c "
import glob
import os
import gzip
import pandas as pd
import multiprocessing
multiprocessing.set_start_method('fork', force=True)

from tqdm.contrib.concurrent import process_map

def parse_file(file):
    lines = []
    with gzip.open(file, 'rt') as f:
        for line in f:
            if line[0] == '>':
                lines.append([line[1:].split()[0], '_'.join(os.path.basename(file).split('_')[:2])])
    return lines

pd.DataFrame([
    x for y in process_map(parse_file, glob.glob('assembly/fna/**/*.fna.gz'), max_workers=32, chunksize=1, leave=False) for x in y
], columns=['accession', 'assembly']).to_csv('assembly/accession2assembly.tsv', sep='\t', index=False, header=None)
"
## ----------------------------------------------------------------------------------------------------
