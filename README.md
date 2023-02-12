# SingleCell-SV

## Introduction
Structural variation (SV) is generally defined as a region of DNA approximately 1 kb and larger in size and can include inversions and balanced translocations or genomic imbalances (insertions and deletions).

## Files
- rawpair2cool.sh
把从GEO数据库中下载的数据转成cool格式

- data_process.ipynb
生成模型的输入数据,包括pos label的和neg label的。

## DataSet
pos: 高置信SV list中breakpoints附近的小矩阵。
neg:
- 癌症细胞系(K562)
    - 不与breakpoint有overlap的区域
    来自snHiC的数据集
- GM12878
  来自Dip-c的数据集
    - 随机点
    - loop
    - AB区间转换的点

## Model

## Usage

```shell
conda env create -f sv.yml

conda create -n sv-env python=3.8
conda activate sv-env
pip install -r requirement.txt
```


```shell
./rawpair2cool.sh
```

## Reconstruct 3D model

snp文件来自Comprehensive, integrated, and phased whole-genome analysis of the primary ENCODE cell line K562