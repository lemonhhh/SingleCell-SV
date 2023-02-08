# SingleCell-SV

## Introduction
Structural variation (SV) is generally defined as a region of DNA approximately 1 kb and larger in size and can include inversions and balanced translocations or genomic imbalances (insertions and deletions).

## Files
- rawpair2cool.sh
把从GEO数据库中下载的数据转成cool格式

- data_process.ipynb
生成模型的输入数据,包括pos label的和neg label的。

## 数据集
只关注K562细胞系（暂时）

pos: 高置信SV list中breakpoints附近的小矩阵。
neg:
- 癌症细胞系
    - 不与breakpoint有overlap的区域
- GM12878
    - 随机点
    - loop
    - AB区间转换的点


## Usage

'''
conda env create -f sv.yml
'''

./rawpair2cool.sh
