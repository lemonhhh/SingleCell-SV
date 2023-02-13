import os


# SAMPLES = ['SRR3344085','SRR3344176']
# SAMPLES = ['SRR3344175']
SAMPLES = ['SRR3344149']

configfile: "hic/config.yaml"

#############RULE_ALL###############

rule all:
    input:
        #HiC
        expand("result/cleaned_pairs/c12/{sample}.pairs.gz",sample=SAMPLES),
        expand("result/dip_pairs/{sample}.dip.pairs.gz",sample=SAMPLES),
        
        #3d info
        expand("processed/{sample}/3d_info/{sample}.{res}.align.rms.info",sample=SAMPLES,res=["20k","50k","200k","1m"]),
        #3dg
        expand("processed/{sample}/3d_info/{res}.{rep}.3dg", sample=SAMPLES,
            res=["20k","50k","200k","1m"],
            rep=list(range(5))),
        #可视化
        expand("result/cif_cpg/{sample}.{res}.{rep}.cpg.cif", sample=SAMPLES,
           res=["20k","50k","200k","1m"],
           rep=list(range(5))),
        
 
    
    threads: config["resources"]["generateStat_cpu_threads"] 
    
    #统计结果，行数之类的
    shell:"""
        echo "Finished"
        # ./hic/CHARM_scripts/generateStat.sh
    """

include: "rules/scHiC_2dprocess.rules"
include: "rules/scHiC_3dprocess.rules"




