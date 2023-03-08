#! /bin/bash
#SBATCH -J python
#SBATCH -o python.out
#SBATCH --partition=compute
#SBATCH --cpus-per-task=40
#SBATCH --qos=medium

python get_pos_data.py