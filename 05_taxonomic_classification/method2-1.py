#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# Date: 2026-01-12
#
# Description:
#       Phage taxonomy annotation pipeline:
#       virify (v3 版数据库和 v4 版分类信息2022版)
#       1. 预测ORF (prodigal, -p meta)
#       2. 预测出的蛋白序列输入 vcontact2
#       3. 同一病毒簇中, 若含有参考数据库的序列, 则将该参考序列的分类信息分配给其他成员
#   参考文献: 2022-An integrated detection, annotation and taxonomic classification pipeline using virus-specific protein profile hidden Markov models
# 
# Dependencies:
#       prodigal --version 2.6.3
#       vcontact2 --version 0.11.3
######################################

import os
import subprocess
import pandas as pd
import logging
from pathlib import Path
from Bio import SeqIO

###################################### step1 运行 virify
def run_virify(virify_path, votus_fa, out_dir, viphog_version, meta_version, db_dir, onlyannotate=True, length=0):
    '''
    使用 virify 预测病毒分类注释
    '''
    os.makedirs(out_dir, exist_ok=True)

    cmd = ["nextflow", "run", virify_path, 
           "--fasta", votus_fa, 
           "--output", out_dir, 
           "--viphog_version", str(viphog_version), 
           "--meta_version", str(meta_version), 
           "--databases", str(db_dir), 
           "--onlyannotate", str(onlyannotate), 
           "--length", str(length),
           "--publish_all"]

    subprocess.run(cmd, check=True)





def main():
    """
    完整分析流程
    """
    # votus 路径
    votus = Path("/home/shijiabin/2025_2ME/03_virus_identification/methodA/votus.fna")

    # virify 的位置
    virify_path = Path("/home/shijiabin/opt/emg-viral-pipeline-3.2.0/main.nf")

    # 参考数据库的路径
    db_dir = Path("/home/shijiabin/db/virify_db")

    # hmm 数据库版本 vpHMM_database_v{1,2,3}.hmm
    viphog_version = "v3"

    # 数据库分类信息版本 additional_data_vpHMMs_v{1,2,3,4}.tsv
    meta_version = "v4"

    #### step1 运行 virify
    # 输入
    votus
    # 输出
    virify_dir = Path("01_virify")

    run_virify(virify_path, votus, virify_dir, viphog_version, meta_version, db_dir)

if __name__ == "__main__":
    main()

    
