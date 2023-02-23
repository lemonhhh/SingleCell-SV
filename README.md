# SingleCell-SV

## Introduction
Structural variation (SV) is generally defined as a region of DNA approximately 1 kb and larger in size and can include inversions and balanced translocations or genomic imbalances (insertions and deletions).

## Files
- data_process.ipynb

Generate input dataset.

- data_explore.ipynb

Explore data.


## DataSet
The training data we used is from [GSE84920](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84920)

The high-confidence SVs are from *Integrative detection and analysis of structural variation in cancer genomes*

pos: 
- Cancer cell line
neg:
- Cancer cell line (no breakpoint region)
- GM12878
    - loop


## Model

A simple CNN model


## Usage

```shell
conda env create -f sv.yml

conda create -n sv-env python=3.8
conda activate sv-env
pip install -r requirement.txt
```



## Reconstruct 3D model

Algorithm is from hickit.

SNP lists are from *Comprehensive, integrated, and phased whole-genome analysis of the primary ENCODE cell line K562*