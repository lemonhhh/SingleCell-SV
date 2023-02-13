#!/bin/bash
# source ~/miniconda3/etc/profile.d/conda.sh
source ~/anaconda3/etc/profile.d/conda.sh
conda activate charm

#创建文件夹
# mkdir -f ./stat

#统计原始fastq文件行数
find -L ./Rawdata/ -name "*R1*.gz" | parallel --tag 'zcat {} | grep "^+" | wc -l' | sort > ./stat/raw.fq.stat
#统计DNA行数
find ./processed -name "*.dna.R1.fq.gz" | parallel --tag 'zcat {} | grep "^+" | wc -l' | sort > ./stat/dna.fq.stat
#统计RNA行数
#clean后的read数量
find ./processed -name "*.rna.clean.R1.fq.gz" | parallel --tag 'zcat {} | grep "^+" | wc -l' | sort > ./stat/rna.fq.stat

#统计ct
#read数量
find ./processed -name "*.ct.R1.fq.gz" | parallel --tag 'zcat {} | grep "^+" | wc -l' | sort > ./stat/ct.read.stat
#frag数量
find ./processed/ct_all/*.pairend.sort.bam | parallel --tag 'samtools flagstat {} | grep read1 | sed "s/ + 0 read1//g"'| sort > ./stat/ct.frag.stat

#bulk的，不去malbalc的重
find ./processed -name "AcACEC.ct.R2.rmdup.bam" | parallel --tag 'samtools flagstat {} | grep primary | head -n 1 | sed "s/ + 0 primary//g"' >> ./stat/ct.frag.stat
find ./processed -name "AcCH.ct.R2.rmdup.bam" | parallel --tag 'samtools flagstat {} | grep primary | head -n 1 | sed "s/ + 0 primary//g"' >> ./stat/ct.frag.stat


#统计atac
#read数量
find ./processed -name "*.atac.R1.fq.gz" | parallel --tag 'zcat {} | grep "^+" | wc -l' | sort > ./stat/atac.read.stat
#frag数量
find ./processed/atac_all/*.pairend.sort.bam | parallel --tag 'samtools flagstat {} | grep read1 | sed "s/ + 0 read1//g"'| sort >> ./stat/atac.frag.stat
find ./processed -name "AcACEC.atac.R2.rmdup.bam" | parallel --tag 'samtools flagstat {} | grep primary | head -n 1 | sed "s/ + 0 primary//g"' >> ./stat/atac.frag.stat



#这个文件是2dprocess rules中生成的
find ./ -name "raw.pairs.gz" | parallel --tag 'zcat {} |grep -v "^#"| wc -l' | sort > ./stat/raw.pairs.stat
#统计contact数量
find ./processed -name "contacts.pairs.gz" | parallel --tag 'zcat {} | grep -v "^#" |wc -l' | sort > ./stat/pairs.dedup.stat
#三次clean分别在干什么
find ./result/cleaned_pairs/c1 -name "*.pairs.gz" | parallel --tag 'zcat {} |grep -v "^#" | wc -l' | sort > ./stat/pairs.c1.stat
find ./result/cleaned_pairs/c12 -name "*.pairs.gz" | parallel --tag 'zcat {} | grep -v "^#" |wc -l' | sort > ./stat/pairs.c12.stat
find ./result/cleaned_pairs/c123 -name "*.pairs.gz" | parallel --tag 'zcat {} | grep -v "^#" |wc -l' | sort > ./stat/pairs.c123.stat

#染色体间的
for i in `find ./result/cleaned_pairs/c123 -name "*.pairs.gz" |sort`; do echo -n $i;echo -n -e "\t";zcat $i| grep -v "^#" | awk '$2!=$4 {print $0}' | wc -l; done > ./stat/inter.pairs.c123.stat
for i in `find ./result/cleaned_pairs/c12 -name "*.pairs.gz" |sort`; do echo -n $i;echo -n -e "\t";zcat $i| grep -v "^#" | awk '$2!=$4 {print $0}' | wc -l; done > ./stat/inter.pairs.c12.stat


#结构信息
find processed/ -name "*rms.info" | xargs -I {} grep --with-filename "top3 RMS RMSD" {} | sed -e "s/processed\///g" -e "s/\/3d_info\/50k.align.rms.info:\[M::__main__\] top3 RMS RMSD: /\t/g" | sort > ./stat/rmsd.info

find processed -name "*.yperx.txt" | parallel 'paste <(echo {}) <(cat {})'| sort > ./stat/yperx.stat
