#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# Date: 2026-01-09
#
# Description:
#       Phage taxonomy annotation pipeline:
#       Prodigal + vConTACT2 + geNomad
#       1. 预测ORF (prodigal, -p meta)
#       2. 预测出的蛋白序列输入 vcontact2
#       3. 同一病毒簇中, 若含有参考数据库的序列, 则将该参考序列的分类信息分配给其他成员
#   参考文献: 2025-Metagenomic analysis reveals gut phage diversity across three mammalian models
# 
# Dependencies:
#       geNomad --version 1.11.2
#       prodigal --version 2.6.3
#       vcontact2 --version 0.11.3
######################################

import os
import pandas as pd
import subprocess
from contextlib import redirect_stdout

###################################### step1 genomad 分类注释 
def run_genomad(input_virus, output_dir, db, threads):
    os.makedirs(output_dir, exist_ok=True)
    cmd = f"micromamba run -n genomad genomad annotate {input_virus} {output_dir} {db} --lenient-taxonomy --threads {threads}"
    subprocess.run(cmd, shell=True, check=True)

###################################### step2 预测 ORF
def run_prodigal(input_virus, output_protein):
    os.makedirs(os.path.dirname(output_protein), exist_ok=True)
    cmd = f"prodigal -i {input_virus} -a {output_protein} -p meta"
    subprocess.run(cmd, shell=True, check=True)

###################################### step3 运行 vcontact2, 生成 viral clusters (VCs)
def run_vcontact2(input_protein, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    g2g_file = f"{output_dir}/g2g.csv"
    cmd1 = f"micromamba run -n vc2 vcontact2_gene2genome -s 'Prodigal-FAA' -p {input_protein} -o {g2g_file}"
    cmd2 = f"micromamba run -n vc2 vcontact2 --rel-mode 'Diamond' --db 'None' --pcs-mode MCL --vcs-mode ClusterONE --raw-proteins {input_protein} --proteins-fp {g2g_file} --output-dir {output_dir}"
    subprocess.run(cmd1, shell=True, check=True)
    subprocess.run(cmd2, shell=True, check=True)

###################################### step4 分配分类信息
# 分配分类信息的方法
def assign_classification(genomad_path, vcontact_path, merge_out, final_out):
    """
    输入：
        genomad_path(votus_taxonomy.tsv路径)
        vcontact_path(genome_by_genome_overview.csv路径)
    输出：
        merge_out(合并后文件路径)
        将votus_taxonomy的分类谱系按等级拆分,然后合并到genome_by_genome_overview中
    
        final_out(填充后文件路径)
        首先按VC编号分组,
        对每个VC组,进行如下处理：
        计算组内成员的Genus中非缺失值出现频数, 频数最高的成员,将它的各分类等级列的值赋值给其他Genus为缺失值的成员
        如果成员的Genus中没有非缺失值,则计算组内成员的Family中非缺失值出现频数....
    """
    # 定义分类等级列(从高到低)
    taxonomy_levels = ['Viruses', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

    # 读取输入文件
    df_genomad = pd.read_csv(genomad_path, sep='\t')
    df_vcontact = pd.read_csv(vcontact_path)

    # 将lineage列按分号进行拆分,分配拆分后的分类等级到各列
    split_result = df_genomad['lineage'].str.split(';', expand=True)
    df_genomad[taxonomy_levels] = split_result

    # 提取出 Genome, preVC, VC Status, VC 列
    df_vcontact = df_vcontact[['Genome', 'preVC', 'VC Status', 'VC']]

    # 将df_genomad的各分类等级添加到df_vcontact中
    df_vcontact = pd.merge(
        df_vcontact,
        df_genomad[['seq_name'] + taxonomy_levels],
        left_on='Genome',
        right_on='seq_name',
        how='left'
    )

    # 将空白串替换为缺失值
    df_vcontact = df_vcontact.replace('', pd.NA)

    # 保存合并结果
    df_vcontact.to_csv(merge_out, sep='\t', index=False)

    # 按VC分组
    groups = df_vcontact.groupby('VC')

    # 对每组成员的分类等级进行缺失值填充
    for vc, group in groups:
        #print(f"VC: {vc}")

        # 获取该组在原数据框中的索引
        group_indices = group.index
        
        # 按照分类等级从低到高进行处理(Genus -> Family -> Order ...)
        # 从Genus开始向上循环(不包括Species)
        genus_index = taxonomy_levels.index('Genus')  # 获取Genus在分类等级列表中的索引
        for level_idx in range(genus_index, -1, -1):  # 从Genus开始向上
            current_level = taxonomy_levels[level_idx]
            target_levels = taxonomy_levels[:level_idx+1]  # 包含当前等级及其以上所有等级
            
            # 获取该组中当前分类等级非缺失值
            non_na_values = group[current_level].dropna()

            if len(non_na_values) == 0:
                #print(f"{current_level} 等级没有分类信息,跳过当前等级的填充,向上一等级寻找...")
                continue
            
            if len(non_na_values) > 0:
                # 计算非缺失值出现频数
                value_counts = non_na_values.value_counts()
                
                # 获取频数最高的值
                most_frequent_value = value_counts.index[0]

                # 找到具有最高频数值的成员(取第一个即可)
                reference_member = group[group[current_level] == most_frequent_value].iloc[0]
                
                # 获取需要填充的成员(当前等级为缺失值的成员)
                members_to_fill = group[group[current_level].isna()]
                
                if len(members_to_fill) > 0:
                    # 打印需要填充的成员信息
                    print(f"{vc}")
                    print(f"{current_level} 等级有分类信息 {most_frequent_value}")
                    print(f"需要填充的成员为：")
                    print(members_to_fill[['Genome']].to_string(index=False, header=False))
                    print(f"\n")

                    # 对需要填充的成员进行填充
                    for idx in members_to_fill.index:
                        for level in target_levels:
                            # 只填充缺失值
                            if pd.isna(df_vcontact.loc[idx, level]):
                                df_vcontact.loc[idx, level] = reference_member[level]      
                #else:
                    #print(f"所有成员当前等级都有分类信息,不需要填充,{vc} 处理完毕")

                # 一旦找到可以处理的等级,就跳出循环
                break


    # 保存填充结果
    df_vcontact.to_csv(final_out, sep='\t', index=False)

# 根据 genomad 分类注释 + viral clusters, 分配分类注释
def run_assign_classification(input_genomad_tax, input_overview, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    merge_taxonomy = os.path.join(output_dir, "merge_taxonomy.tsv")
    assign_taxonomy = os.path.join(output_dir, "assign_taxonomy.tsv")
    
    log_file = f"{output_dir}/assign.log"
    
    with open(log_file, "w") as f, redirect_stdout(f):
        assign_classification(input_genomad_tax, input_overview, merge_taxonomy, assign_taxonomy)

# =============================
# Main Workflow
# =============================

def main():
    """
    完整分析流程
    """
    # votus 路径
    votus = Path("/home/shijiabin/2025_2ME/03_virus_identification/methodA/votus.fna")

    # genomad 数据库
    db = Path("/home/shijiabin/db/genomad_db")
    
    # cpus
    threads = 10
    
    #### step1 genomad 分类注释
    # 输入
    votus
    # 输出
    outdir_genomad = Path("01_genomad")
    genomad_tax = output_dir_genomad / "votus_annotate/votus_taxonomy.tsv"

    run_genomad(input_votus, outdir_genomad, db)

    #### step2 预测 ORF
    # 输入
    votus
    # 输出
    protein = Path("02_prodigal/protein.faa")

    run_prodigal(votus, protein)

    #### step3 运行 vcontact2, 生成 viral clusters (VCs)
    # 输入
    protein
    # 输出
    outdir_vcontact2 = Path("03_vcontact2")
    g2g = outdir_vcontact2 / "g2g.csv"
    overview = outdir_vcontact2 / "genome_by_genome_overview.csv"
    
    run_vcontact2(protein, outdir_vcontact2)

    ##### step4 分配分类信息
    # 输入
    genomad_tax
    overview
    # 输出
    outdir_assign = Path("04_assign")
    tax_assign = outdir_assign / "votus_taxonomy.tsv"
    
    run_assign_classification(genomad_tax, overview, outdir_assign)

    ### 提取最终结果
    #final_taxonomy = merge_taxonomy = os.path.join(output_dir_assign, f"{sample}_merge_taxonomy.tsv")
    #final_tax = f"votus_{method}/{sample}_taxonomy.tsv"
    #if os.path.exists(final_taxonomy):
        #subprocess.run(f"cp {final_taxonomy} {final_tax}", shell=True, check=True)
