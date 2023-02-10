####################################
#       CHARM_pipline              #
#@author Z Liu                     #
#@Ver 0.2.0                        #
#@date 2021/8/12                   #
####################################

#############CONFIG#################

import os

#input
# SAMPLES = [i.split(sep='_')[0] for i in os.listdir("./Rawdata")]
#SAMPLES = os.listdir("./Rawdata")

SAMPLES = ['ACEC001','ACEC002','ACEC003','ACEC004','ACEC005','CH001','CH002','CH003','CH004','CH005']
shiftExtendSize = [125, 250]

configfile: "CHARM/config.yaml"

#############RULE_ALL###############

rule all:
    input:
        #RNA
        expand("result/RNA_Res/counts.{type}.{genome}.tsv",type=["gene","exon"],genome=["total","genome1","genome2"] if config["if_RNA_snp_split"] else ["total"]),
        expand("result/RNA_Res/counts.{type}.{genome}.format.tsv",type=["gene","exon"],genome=["total","genome1","genome2"] if config["if_RNA_snp_split"] else ["total"]),
        expand("processed/RNA_all/umibycell.{sample}.rna.R1.fq",sample=SAMPLES),
        #HiC
        expand("result/cleaned_pairs/c12/{sample}.pairs.gz",sample=SAMPLES),
        expand("result/dip_pairs/{sample}.dip.pairs.gz",sample=SAMPLES),
       
        #cuttag 
        expand("processed/{sample}/ct/{sample}.ct.R1.fq.gz", sample=SAMPLES if config["if_cuttag"] else []),
        expand("processed/ct_all/{sample}.ct.pairend.sort.bam", sample=SAMPLES if config["if_cuttag"] else []),
        expand("processed/ct_all/{sample}.ct.R2.sort.bam", sample=SAMPLES if config["if_cuttag"] else []),
  
        
        #atac  part
        #split出来的文件
        # expand("processed/{sample}/atac/{sample}.atac.R1.fq.gz", sample=SAMPLES if config["if_atac"] else []),
        #R1 R2 map的结果
        expand("processed/atac_all/{sample}.atac.pairend.sort.bam", sample=SAMPLES if config["if_atac"] else []),
        #dedup后R2map的结果
        expand("processed/atac_all/{sample}.atac.R2.sort.bam",sample=SAMPLES if config["if_atac"] else []),

        
        
        #3d info
        expand("processed/{sample}/3d_info/{sample}.{res}.align.rms.info",sample=SAMPLES if config["if_structure"] else [],res=["20k","50k","200k","1m"] if config["if_structure"] else []),
        
        expand("processed/{sample}/3d_info/{res}.{rep}.3dg", sample=SAMPLES if config["if_structure"] else [],
            res=["20k","50k","200k","1m"] if config["if_structure"] else [],
            rep=list(range(5)) if config["if_structure"] else []),
        
        expand("result/cif_cpg/{sample}.{res}.{rep}.cpg.cif", sample=SAMPLES if config["if_structure"] else [],
           res=["20k","50k","200k","1m"] if config["if_structure"] else [],
           rep=list(range(5)) if config["if_structure"] else []),
        
        #radial position
        expand("result/radialPos/{res}/{sample}.rp.{res}.{rep}.color", sample=SAMPLES if config["if_structure"] else [],
           res=["20k","50k","200k","1m"] if config["if_structure"] else [],
           rep=list(range(5)) if config["if_structure"] else []),
    
    threads: config["resources"]["generateStat_cpu_threads"] 
    
    #统计结果，行数之类的
    shell:"""
        ./CHARM/CHARM_scripts/generateStat.sh
    """
  
  

include: "rules/CHARM_split.rules"
include: "rules/CHARM_cuttag.rules"
include: "rules/CHARM_atac.rules"
include: "rules/scHiC_2dprocess.rules"
include: "rules/scHiC_3dprocess.rules"
include: "rules/CHARM_RNA.rules"



