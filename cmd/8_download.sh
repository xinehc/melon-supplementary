echo "----- 8_download.sh -----"
echo "Downloading PlusPFP sequences ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p plus/faa

python -c "
import pandas as pd

assembly = pd.read_table('assembly/assembly_summary_refseq.txt', skiprows=1, low_memory=False).rename({'#assembly_accession': 'assembly'}, axis=1)
assembly = assembly[assembly['ftp_path'] != 'na']

assembly.loc[(assembly.organism_name.str.contains('Homo sapiens')) & (assembly.refseq_category == 'reference genome'), 'group'] = 'human'
assembly = assembly[assembly.group.isin({'fungi', 'protozoa', 'viral', 'plant', 'human'})]
assembly = assembly[assembly.assembly_level.isin({'Chromosome', 'Complete Genome'})]
assembly['faa'] = assembly['ftp_path'] + assembly['ftp_path'].str.split('/').str.get(-2) + '_protein.faa.gz'

for taxonomy in assembly.group.unique():
    faa = assembly[assembly['group'] == taxonomy]['faa'].to_list()
    n = 8
    chunks = [faa[i::n] for i in range(n)]

    for i, chunk in enumerate(chunks):
        if chunk:
            with open('plus/faa/' + taxonomy + '.split.' + str(i) + '.id', 'w') as w:
                w.write('\n'.join(chunk) + '\n')
"

find plus/faa -maxdepth 1 -name '*.id' | sort | xargs -P 8 -I {} bash -c '
    wget -i ${1} -qN --retry-on-http-error=503 -P ${1%.id};
    find ${1%.id} -maxdepth 1 -name "*.faa.gz" | xargs cat | gzip -d > ${1%.id}.faa;
' - {}
## ----------------------------------------------------------------------------------------------------
