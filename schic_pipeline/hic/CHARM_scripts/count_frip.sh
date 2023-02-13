#!/bin/bash

four_sample=("ACEC001" "ACEC002" "ACEC003" "ACEC004" "ACEC005")
two_sample=("CH001" "CH002" "CH003" "CH004" "CH005")
#four_sample和two_sample相加
ct_sample=("${four_sample[@]}" "${two_sample[@]}")

ct_bulk_sample=("AcACEC" "AcCH")
atac_bulk_sample=("AcACEC")


if [ -f "./stat/ct_frip.txt" ];then
  rm ./stat/ct_frip.txt
  touch ./stat/ct_frip.txt
fi

if [ -f "./stat/atac_frip.txt" ];then
  rm ./stat/atac_frip.txt
  touch ./stat/atac_frip.txt
fi

# ct_peak_file="processed/AcCH/ct/peak/AcCH_peaks.narrowPeak "
ct_peak_file="./stat/encode_h3k27ac.broadPeak"
for ct in ${ct_sample[@]}
do
  #bam行数
  ct_bam_file="./processed/ct_all/${ct}.ct.pairend.sort.bam"
  # ct_bam_file="./processed/ct_all/${ct}.ct.R2.sort.bam"
  ct_bam_line=`samtools view -c ${ct_bam_file}`

  #和peak的交集有多少行
  
  ct_intersect=`bedtools sort -i ${ct_peak_file} | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a ${ct_bam_file} -b stdin -ubam | samtools view -c`
  
  #计算frip
  ct_frip=`echo "scale=4;${ct_intersect}/${ct_bam_line}" | bc`

  #写入./stat/ct_frip.txt文件
  echo -e "${ct}\t${ct_frip}" >> ./stat/ct_frip.txt
done

for ct in ${ct_bulk_sample[@]}
do
  
  #bulk只需要去PCR dup之后的文件
  ct_bam_file="./processed/${ct}/ct/${ct}.ct.R2.rmdup.bam"
  #总共有多少行
  ct_bam_line=`samtools view -c ${ct_bam_file}`

  #ground truth的peak文件
  ct_intersect=`bedtools sort -i ${ct_peak_file} | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a ${ct_bam_file} -b stdin -ubam | samtools view -c`
  
  #计算frip
  ct_frip=`echo "scale=4;${ct_intersect}/${ct_bam_line}" | bc`

  #写入./stat/ct_frip.txt文件
  echo -e "${ct}\t${ct_frip}" >> ./stat/ct_frip.txt
done


# # atac_peak_file="./Ground_truth/double_insertion/atac_output_test/peak/ext250/PD10_XD_omniATAC_ext250_peaks.narrowPeak"
# atac_peak_file="./processed/AcACEC/atac/peak/ext250/AcACEC_ext250_peaks.narrowPeak"
# for atac in ${four_sample[@]}
# do
#   #总共有多少行
#   # atac_bam_file="./processed/atac_all/${atac}.atac.R2.sort.bam"
#   atac_bam_file="./processed/atac_all/${atac}.atac.pairend.sort.bam"
#   atac_bam_line=`samtools view -c ${atac_bam_file}`

#   #和peak的交集有多少行 
#   atac_intersect=`bedtools sort -i ${atac_peak_file} | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a ${atac_bam_file} -b stdin -ubam | samtools view -c`
  
#   #计算frip
#   atac_frip=`echo "scale=4;${atac_intersect}/${atac_bam_line}" | bc`
#   #写入./stat/ct_frip.txt文件
#   echo -e "${atac}\t${atac_frip}" >> ./stat/atac_frip.txt
# done

# for atac in ${atac_bulk_sample[@]}
# do
#   #总共有多少行
#   # atac_bam_file="./processed/atac_all/${atac}.atac.R2.sort.bam"
#   atac_bam_file="./processed/${atac}/atac/${atac}.atac.R2.rmdup.bam"
#   atac_bam_line=`samtools view -c ${atac_bam_file}`
#   #和peak的交集有多少行
#   atac_intersect=`bedtools sort -i ${atac_peak_file} | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a ${atac_bam_file} -b stdin -ubam | samtools view -c`
  
#   #计算frip
#   atac_frip=`echo "scale=4;${atac_intersect}/${atac_bam_line}" | bc`

#   #写入./stat/ct_frip.txt文件
#   echo -e "${atac}\t${atac_frip}" >> ./stat/atac_frip.txt
# done