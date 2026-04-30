ulimit -n 65536

## download necessary files
bash cmd/1_download.sh
bash cmd/2_download.sh
bash cmd/3_download.sh

## build protein database
bash cmd/4_annotate.sh
bash cmd/5_cluster.sh

## build nucleotide database
bash cmd/6_annotate.sh
bash cmd/7_cluster.sh

## setup PlusPFP
bash cmd/8_download.sh
bash cmd/9_annotate.sh