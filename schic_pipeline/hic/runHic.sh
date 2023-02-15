#!/bin/bash
#usage: ./runCHARM.sh

#回到上一级目录
cd ../
mkdir -p slurm_log

snakemake --use-conda --cluster 'sbatch --exclude node03 --qos=medium --output=slurm_log/slurm-%j.out --cpus-per-task={threads} -t 7-00:00 -J 3D' --jobs 100 --resources nodes=100 --rerun-incomplete -s ./hic/hic.smk --keep-going 

mkdir -p ./analysis
cp hic/stat.ipynb ./analysis/stat.ipynb
