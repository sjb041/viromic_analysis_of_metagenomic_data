#!/usr/bin/env bash

######################################
# Date: 2026-01-20
# Description:
#     Cluster genomes to Genus-clusters based on AAI
#     https://github.com/snayfach/MGV/tree/master/aai_cluster
# Dependencies:
#     prodigal --version 2.6.3
#     diamond --version 2.0.6
#     python --version 2.7.15 --> numpy biopython
#     mcl --version 14-137
######################################

# 脚本存放目录
scripts="$HOME/2025_2ME/scripts/aai_cluster" 

#### step1 预测 ORF
dir="02_genus_clusters" && mkdir -p $dir && cd $dir

votus="/home/shijiabin/2025_2ME/03_virus_identification/methodA/votus.fna"

prodigal -i $votus -a votu_orfs.faa -p meta

#### step2 DIAMOND 进行成对蛋白质序列比对
mkdir diamond

# 建库
diamond makedb --in votu_orfs.faa --db diamond/db --threads 10

# 执行 all-vs-all 比对（很快）
diamond blastp -q votu_orfs.faa -d diamond/db -e 1e-5 --max-target-seqs 10000  --query-cover 50 --subject-cover 50 --outfmt 6 -o votu_orfs_diamond_results.tsv
rm -r diamond

#### step3 计算每对vOTU之间共享基因的百分比和平均氨基酸相似性（AAI）
micromamba run -n py2.7 python $scripts/amino_acid_identity.py --in_faa votu_orfs.faa --in_blast votu_orfs_diamond_results.tsv --out_tsv votu_aai.tsv
rm votu_orfs_diamond_results.tsv

#### step4 建立votu间的连接网络
micromamba run -n py2.7 python $scripts/filter_aai.py --in_aai votu_aai.tsv --min_percent_shared 20 --min_num_shared 16 --min_aai 40 --out_tsv genus_edges.tsv
rm votu_aai.tsv

#### step5 基于vOTU之间的连接使用MCL进行聚类
mcl genus_edges.tsv -te 8 -I 2.0 --abc -o genus_clusters.txt
rm genus_edges.tsv
