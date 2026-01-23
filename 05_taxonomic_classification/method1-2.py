#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# Date: 2026-01-10
#
# Description:
#       Phage taxonomy annotation pipeline:
#       Prodigal + vConTACT2 + geNomad
#       1. genomad 分类注释
#       2. 预测ORF (prodigal, -p meta)
#       3. 预测出的蛋白序列输入 vcontact2, 生成 viral clusters
#       4. 根据 genomad 的注释信息和 viral clusters, 分配分类注释
#       5. 优化分类注释表
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
from pathlib import Path
from contextlib import redirect_stdout
from Bio import SeqIO

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
def run_vcontact2(input_protein, output_dir, threads):
    os.makedirs(output_dir, exist_ok=True)
    g2g_file = f"{output_dir}/g2g.csv"
    cmd1 = f"micromamba run -n vc2 vcontact2_gene2genome -s 'Prodigal-FAA' -p {input_protein} -o {g2g_file}"
    cmd2 = f"micromamba run -n vc2 vcontact2 --rel-mode 'Diamond' --db 'None' --pcs-mode MCL --vcs-mode ClusterONE --raw-proteins {input_protein} --proteins-fp {g2g_file} --output-dir {output_dir} -t {threads}"
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

# 按行修改分类注释
def normalize_taxonomy(row, tax_cols):
    """
    逐级规范分类注释：

    - 空值 / NaN / Unassigned → 视为缺失
    - 第一层 Kingdom 缺失 → Unknown
    - 其余层：
        * 若本层缺失：
            - 上一层 Unknown → Unknown
            - 上一层 Unclassified* → 继承
            - 否则 → Unclassified + 上一层
        * 若本层有值 → 保留
    """

    invalid_values = {None, "", "Unknown", "unknown", "Unclassified", "unclassified", "Unassigned", "unassigned"}

    for i, col in enumerate(tax_cols):
        val = row[col]

        # 统一处理缺失
        if pd.isna(val) or val in invalid_values:
            val = None

        if i == 0:
            # 第一层
            row[col] = val if val is not None else "Unknown"
        else:
            prev = row[tax_cols[i - 1]]

            # 保证 prev 是字符串
            if pd.isna(prev) or prev in invalid_values:
                prev = "Unknown"

            if val is not None:
                row[col] = val
            else:
                if prev == "Unknown":
                    row[col] = "Unknown"
                elif isinstance(prev, str) and prev.startswith("Unclassified"):
                    row[col] = prev
                else:
                    row[col] = f"Unclassified{prev}"

    return row

# 优化分类注释表
def optimize_annotation_table(fasta_file, tax_file, output_tax_file):
    """
    优化病毒分类注释表

    以 FASTA 文件中的 Genome ID 为全集，对现有注释表进行补全：
    - FASTA 中存在但注释表中缺失的 Genome 会被补充为新行
    - 新增行的分类信息为空
    - 已有注释信息保持不变

    提取第一列和 Kingdom~Species 之间的列(要找到这两个列名然后提取他们两以及之间的列)
    由高等级开始处理(第二列), 如果第二列为空,那么填充为 Unknown;
    然后处理下一列,如果为空或 Unknown,则判断上一列,若上一列为 Unknown,则填充为 Unknown;
    若上一列不为 Unknown 且不以 Unclassified 开头,则填充为 Unclassified + "上一列的值";
    若上一列以 Unclassified 开头,则直接复制上一列的值;若本列不为空或 Unknown,则保留本列的值. 
    再下一列的处理同第三列
    """

    #### step1 提取原分类注释表的第一列 + Kingdom~Species 列
    # 读取分类注释表
    tax_df = pd.read_csv(tax_file, sep="\t")

    # 找到 Kingdom 和 Species 的列索引
    cols = tax_df.columns.tolist()
    start_col = cols.index("Kingdom")
    # 检查 Species 列是否存在，如果存在则分类等级包含到 Species，否则使用 Genus
    if "Species" in cols:
        end_col = cols.index("Species")
    else:
        end_col = cols.index("Genus")

    # 提取第一列 + Kingdom~Species, 得到数据框用于合并
    keep_cols = [cols[0]] + cols[start_col:end_col + 1]
    tax_df_sub = tax_df[keep_cols].copy()

    # 修改 tax_df_sub 第一列的列名为 vOTU
    tax_df_sub.rename(columns={cols[0]: "vOTU"}, inplace=True)
    
    #### step2 从 FASTA 提取 ID 全集, 用于补充分类注释表的行, 那些补充进去的行注释为空
    # 提取 id 全集, 创建一个数据框用于合并
    fasta_ids = [rec.id for rec in SeqIO.parse(fasta_file, "fasta")]
    fasta_df = pd.DataFrame({"vOTU": fasta_ids})

    # 以 fasta_ids 为全集进行左连接
    merged_df = fasta_df.merge(
        tax_df_sub,
        on="vOTU",
        how="left"
    )

    #### step3 优化那些已有的注释
    # 那些要优化的列
    tax_cols = merged_df.columns[1:].tolist()
    # 如果t ax_cols 最后一个不是 Species 列，则给 merged_df 添加 Species 列，值为空
    if tax_cols[-1] != "Species":
        merged_df["Species"] = ""  # 添加Species列，值为空字符串
        tax_cols = merged_df.columns[1:].tolist()  # 重新获取tax_cols以包含新添加的Species列

    # 行级处理
    merged_df = merged_df.apply(normalize_taxonomy, axis=1, tax_cols=tax_cols)

    # 输出
    merged_df.to_csv(output_tax_file, sep="\t", index=False)

    return merged_df

# =============================
# Main Workflow
# =============================

def main():
    """
    完整分析流程
    """
    # votus 路径
    votus = Path("/home/shijiabin/2025_2ME/03_virus_identification/votus.fna")

    # genomad 数据库
    db = Path("/home/shijiabin/db/genomad_db/genomad_db_v1.9")
    
    # cpus
    threads = 10
    
    #### step1 genomad 分类注释
    # 输入
    votus
    # 输出
    outdir_genomad = Path("01_genomad")
    genomad_tax = outdir_genomad / "votus_annotate/votus_taxonomy.tsv"

    #run_genomad(votus, outdir_genomad, db, threads)

    #### step2 预测 ORF
    # 输入
    votus
    # 输出
    protein = Path("02_prodigal/protein.faa")

    #run_prodigal(votus, protein)

    #### step3 运行 vcontact2, 生成 viral clusters (VCs)
    # 输入
    protein
    # 输出
    outdir_vcontact2 = Path("03_vcontact2")
    g2g = outdir_vcontact2 / "g2g.csv"
    overview = outdir_vcontact2 / "genome_by_genome_overview.csv"
    
    #run_vcontact2(protein, outdir_vcontact2, threads)

    ##### step4 分配分类信息
    # 输入
    genomad_tax
    overview
    # 输出
    outdir_assign = Path("04_assign")
    tax_assign = outdir_assign / "assign_taxonomy.tsv"
    
    run_assign_classification(genomad_tax, overview, outdir_assign)

    ##### step5 优化分类注释表
    # 输入
    votus
    tax_assign
    # 输出
    final_tax = Path("votus_taxonomy.tsv")
    optimize_annotation_table(votus, tax_assign, final_tax)
if __name__ == "__main__":
    main()
