# snakemake -s pipeline.smk --cores 5 --use-conda

import os
import argparse
import csv
from pathlib import Path
import pandas as pd

##################################### 定义全局变量
# 样本名 
SAMPLES = [
    "F1_1A", "F1_2A", "F1_3A",
    "F2_1A", "F2_2A", "F2_3A",
    "FG_1A", "FG_2A", "FG_3A",
    "L1_1A", "L1_2A", "L1_3A",
    "L2_1A", "L2_2A", "L2_3A",
    "H_LX", "H_O"
]

# uvigs的父目录，要求文件名的格式为sample_uvigs.fa
virus_dir = "/home/shijiabin/2025_2ME/03_virus_identification/03_pipeline2-2/07_all_sample_results"

# virify 的位置
virify = "/home/shijiabin/opt/emg-viral-pipeline-3.0.2/main.nf"

# 参考数据库的路径
db_dir = "/home/shijiabin/db/virify_db"

# hmm 数据库版本 vpHMM_database_v{1,2,3}.hmm
viphog_version = "v3"

# 数据库分类信息版本 additional_data_vpHMMs_v{1,2,3,4}.tsv
meta_version = "v4"


##################################### 定义流程规则
# 最终要生成的文件
rule all:
    input:
        expand("04_final_taxonomy/{sample}_final_taxonomy.tsv", sample=SAMPLES)

# 聚类
rule cluster:
    input:
        virus = lambda wildcards: os.path.join(virus_dir, f"{wildcards.sample}_uvigs.fa")
    output:
        dir = directory("01_cluster/{sample}"),
        votus = "01_cluster/{sample}/votus.fa"
    log:
        "log/cluster/{sample}.log"
    shell:
        "cluster.sh -i {input.virus} -t {input.virus} -a 95 -o {output.dir} > {log} 2>&1"

# 序列重命名
rule rename:
    input:
        fasta = "01_cluster/{sample}/votus.fa"
    output:
        mapping = "02_rename/{sample}_mapping.tsv",
        original_name = "02_rename/{sample}_original_name.txt",
        renamed_fasta = "02_rename/{sample}_renamed.fa"
    log:
        "log/rename/{sample}.log"
    shell:
        """
        rename.py {input.fasta} -m {output.mapping} -o {output.renamed_fasta} > {log} 2>&1
        
        # 取映射文件第1列
        cut -f 1 {output.mapping} > {output.original_name}
        """

# 执行 virify 预测分类
rule virify:
    input:
        renamed_fasta = "02_rename/{sample}_renamed.fa"
    output:
        dir = directory("03_virify/v3v4_ncbi2022_hex_bex/{sample}"),
        taxonomy = "03_virify/v3v4_ncbi2022_hex_bex/{sample}/original_taxonomy.tsv"
    log:
        "log/virify/{sample}.log"
    shell:
        """
        nextflow run {virify} \
            --fasta  {input.renamed_fasta} \
            --output  {output.dir} \
            --viphog_version {viphog_version} \
            --meta_version {meta_version} \
            --databases {db_dir} \
            --onlyannotate true \
            --publish_all \
            --hmmextend \
            --blastextend \
            --length 0 \
            -profile local,docker > {log} 2>&1

        mv {output.dir}/{wildcards.sample}_renamed/06-taxonomy/*_taxonomy.tsv {output.taxonomy}
        """

# 还原 taxonomy 中的序列名
rule restore:
    input:
        original_name = "02_rename/{sample}_original_name.txt",
        taxonomy = "03_virify/v3v4_ncbi2022_hex_bex/{sample}/original_taxonomy.tsv"
    output:
        final_taxonomy = "04_final_taxonomy/{sample}_final_taxonomy.tsv"
    shell:
        """
        # 排序序列的分类谱系文件
        head -n 1 {input.taxonomy} > {wildcards.sample}_tmp
        tail -n +2 {input.taxonomy} | sort -t$'\t' -k1,1V >> {wildcards.sample}_tmp

        # 合并
        cut -f 2- {wildcards.sample}_tmp | paste {input.original_name} - > {output.final_taxonomy}

        # 删除临时文件
        rm {wildcards.sample}_tmp
        """
        
