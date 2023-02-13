#!/bin/bash
# source ~/miniconda3/etc/profile.d/conda.sh
source ~/anaconda3/etc/profile.d/conda.sh
conda activate charm
mkdir -p ./stat/split

find -L ./Rawdata/ -name "*R1*.gz" | parallel --tag 'zcat {} | wc -l' | sort > ./stat/split/raw.fq.stat
find ./processed -name "*.dna.R1.fq.gz" | parallel --tag 'zcat {} | wc -l' | sort > ./stat/split/dna.fq.stat
find ./processed -name "*.rna.R1.fq.gz" | parallel --tag 'zcat {} | wc -l' | sort > ./stat/split/rna.fq.stat
find ./processed -name "*.ct.R1.fq.gz" | parallel --tag 'zcat {} | wc -l' | sort > ./stat/split/ct.fq.stat
find ./processed -name "*.atac.R1.fq.gz" | parallel --tag 'zcat {} | wc -l' | sort > ./stat/split/atac.fq.stat

