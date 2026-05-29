# Melon-supplementary
Scripts for building the GTBD database of [Melon](https://github.com/xinehc/melon).

#### Install necessary packages
```bash
conda install -c bioconda -c conda-forge 'taxonkit>=0.20.0' 'seqkit>=2.13.0' 'hmmer>=3.4' 'mmseqs2>=18.8cc5c' 'blast>=2.16.0' 'diamond>=2.1.24' 'tqdm' 'pandas' 'requests' 'tar' 'wget'
```

#### Build the GTDB database
```bash
bash run.sh
```
