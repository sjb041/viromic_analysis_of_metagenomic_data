#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# Date: 2026-01-09
#
# Description:
#       Phage taxonomy annotation pipeline:
#       Prodigal → vConTACT2 → taxonomy assignmen (ProkaryoticViralRefSeq211-Merged)
#       1. 预测ORF (prodigal, -p meta)
#       2. 预测出的蛋白序列输入 vcontact2
#       3. 同一病毒簇中, 若含有参考数据库的序列, 则将该参考序列的分类信息分配给其他成员
#   参考文献: 2024-The gut ileal mucosal virome is disturbed in patients with Crohn’sdiseaseand exacerbates intestinal inflammation in mice
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

###################################### step1 预测 ORF
def run_prodigal(votus_fa, protein_faa, log):
    '''
    预测 ORF, 使用 meta 模式
    input:
        votus_fa
    output:
        protein_faa, log
    '''
    os.makedirs(os.path.dirname(protein_faa), exist_ok=True)
    os.makedirs(os.path.dirname(log), exist_ok=True)

    cmd = ["prodigal", "-i", votus_fa, "-a", protein_faa, "-p", "meta"]

    with open(log, "w") as lf:
        subprocess.run(cmd, stdout=lf, stderr=lf, check=True)

###################################### step2 运行 vcontact2, 生成 viral clusters (VCs)
def run_vcontact2(protein_faa, db, outdir, g2g, log, threads):
    '''
    运行 vcontact2, 生成 viral clusters (VCs)
    input:
        protein_faa
    output:
        g2g, genome_by_genome_overview, log
    '''
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.dirname(log), exist_ok=True)

    with open(log, "w") as lf:
        lf.write("===== vcontact2_gene2genome =====\n")

        # 生成 g2g 文件
        subprocess.run(
            [
                "micromamba", "run", "-n", "vc2",
                "vcontact2_gene2genome",
                "-s", "Prodigal-FAA",
                "-p", protein_faa,
                "-o", g2g
            ],
            stdout=lf,
            stderr=lf,
            check=True
        )

        lf.write("\n===== vcontact2 clustering =====\n")

        # 生成 viral clusters (VCs)
        subprocess.run(
            [
                "micromamba", "run", "-n", "vc2",
                "vcontact2",
                "--rel-mode", "Diamond",
                "--db", db,
                "--pcs-mode", "MCL",
                "--vcs-mode", "ClusterONE",
                "--raw-proteins", protein_faa,
                "--proteins-fp", g2g,
                "--output-dir", outdir,
                "--threads", str(threads)
            ],
            stdout=lf,
            stderr=lf,
            check=True
        )

###################################### step3 分配分类信息
# 设置日志记录器
def setup_logger(log_file):
    """设置日志记录器"""
    # 清除现有的handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    # 配置新的日志
    logging.basicConfig(
        filename=str(log_file),
        level=logging.INFO,
        format='%(message)s',
        force=True  # 覆盖之前的配置
    )
    return logging.getLogger()

# 根据 vcontact2 的输入 g2g 文件和输出表格, 分配输入序列的分类等级
def assign_taxonomy(g2g_file, overview_file, overview_filtered, final_tax, log_file):
    """
    基于vConTACT2聚类结果为目标序列分配病毒分类信息
    
    该函数通过分析vConTACT2生成的基因组概览文件,为目标序列(来自FASTA文件)分配
    病毒分类信息。主要基于preVC科级别和VC属级别聚类结果,通过统计
    参考序列的分类信息频数来推断目标序列的分类。
    
    参数:
    g2g_file (str): vConTACT2的输入g2g文件,用于获取目标序列的序列名
    overview_file (str): vConTACT2生成的genome_by_genome_overview.csv文件路径
    output_overview_filtered (str): 输出筛选后的完整结果文件路径TSV格式
    output_tax (str): 输出目标序列分类结果文件路径TSV格式,仅前7列
    log_file (str): 日志文件路径,保存警告信息
    
    处理流程:
    1. 从FASTA文件提取目标序列ID列表
    2. 筛选overview文件,只保留同时包含目标序列和参考序列的preVC组
    3. 为每个preVC组的目标序列分配分类信息:
    - Family及以上级别:基于该组参考序列中出现频数最高的分类
    - Genus级别:基于VC小组参考序列中出现频数最高的分类
    """
    os.makedirs(os.path.dirname(overview_filtered), exist_ok=True)
    os.makedirs(os.path.dirname(final_tax), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # 设置日志记录器
    logger = setup_logger(log_file)
    
    # 提取FASTA序列名
    seqid_list = pd.read_csv(g2g_file)['contig_id'].unique().tolist()
    # 读取overview文件
    overview = pd.read_csv(overview_file)

    ###### 筛选overview,只保留那些既有目标序列又有参考序列的preVC组
    # 按照Excel的筛选步骤
    # 1. 筛选那些目标序列preVC列不为空的行
    df1 = overview[(overview['Genome'].isin(seqid_list)) & 
                (overview['preVC'].notna()) & 
                (overview['preVC'] != '')]

    # 2. 使用1的preVC列筛选overview
    preVC_list = list(df1['preVC'].unique())     # preVC列去重,转化为列表
    df2 = overview[overview['preVC'].isin(preVC_list)]

    # 3. 只保留那些既有目标序列又有参考序列的preVC组
    # 找出满足条件的preVC组
    preVC_with_both = []
    for preVC in preVC_list:
        group = df2[df2['preVC'] == preVC]
        has_target = (group['Genome'].isin(seqid_list)).any()  # 组内是否有目标序列
        has_reference = (~group['Genome'].isin(seqid_list)).any()  # 组内是否有参考序列  
        # 只有当组内同时包含目标序列和参考序列时才保留该组
        if has_target and has_reference:
            preVC_with_both.append(preVC)
    # 筛选preVC组
    filtered_df = df2[df2['preVC'].isin(preVC_with_both)]

    ###### 为每个目标序列分配Family和Genus信息
    # Family分配：
    # 对每个preVC组，统计参考序列中各Family的出现频数
    # 排除'Unassigned'值
    # 选择出现频数最高的Family作为该组目标序列的Family
    # 如果有多个Family频数相同则报错
    # Genus分配：
    # 在preVC组内，根据VC进一步分组
    # 排除'Unassigned'值
    # 选择出现频数最高的Genus作为该VC小组目标序列的Genus
    # 如果有多个Genus频数相同则报错
    df = filtered_df.copy() 
    # 按preVC组处理
    for prevc in df['preVC'].unique():
        group_mask = (df['preVC'] == prevc)
        group_data = df[group_mask]
        
        # 获取参考序列用于Family分类
        reference_mask = ~group_data['Genome'].isin(seqid_list)
        reference_family_counts = group_data[reference_mask]['Family'].value_counts()
        
        # 移除'Unassigned'值
        valid_families = reference_family_counts[reference_family_counts.index != 'Unassigned']
        
        if len(valid_families) > 0:
            # 检查是否有多个最高频数
            max_count = valid_families.max()
            top_families = valid_families[valid_families == max_count]
            
            if len(top_families) > 1:
                warning_msg = f"警告: preVC {prevc} 中存在多个相同频数的Family: {list(top_families.index)}, 使用第一个: {top_families.index[0]}"
                print(warning_msg)
                logger.warning(warning_msg)
            
            predicted_family = top_families.index[0]
            
            # 为该组所有目标序列填充Family及更高级分类
            target_mask = group_mask & (df['Genome'].isin(seqid_list))
            
            # 填充Family及更高级分类信息
            reference_taxonomy = group_data[reference_mask].iloc[0]  # 取第一个参考序列的分类信息作为模板
            
            df.loc[target_mask, 'Family'] = predicted_family
            df.loc[target_mask, 'Order'] = reference_taxonomy['Order']
            df.loc[target_mask, 'Class'] = reference_taxonomy['Class']
            df.loc[target_mask, 'Phylum'] = reference_taxonomy['Phylum']
            df.loc[target_mask, 'Kingdom'] = reference_taxonomy['Kingdom']
            
            # 按VC小组处理Genus
            for vc in group_data['VC'].unique():
                if pd.notna(vc):
                    vc_mask = group_mask & (df['VC'] == vc)
                    vc_group_data = df[vc_mask]
                    
                    # 获取参考序列用于Genus分类
                    vc_reference_mask = ~vc_group_data['Genome'].isin(seqid_list)
                    reference_genus_counts = vc_group_data[vc_reference_mask]['Genus'].value_counts()
                    
                    # 移除'Unassigned'值
                    valid_genera = reference_genus_counts[reference_genus_counts.index != 'Unassigned']
                    
                    if len(valid_genera) > 0:
                        # 检查是否有多个最高频数
                        max_count = valid_genera.max()
                        top_genera = valid_genera[valid_genera == max_count]

                        if len(top_genera) > 1:
                            warning_msg = f"警告: VC {vc} 中存在多个相同频数的Genus: {list(top_genera.index)}, 使用第一个: {top_genera.index[0]}"
                            print(warning_msg)
                            logger.warning(warning_msg)
                        
                        predicted_genus = top_genera.index[0]
                        
                        # 为该VC小组的所有目标序列填充Genus
                        vc_target_mask = vc_mask & (df['Genome'].isin(seqid_list))
                        df.loc[vc_target_mask, 'Genus'] = predicted_genus

    # 只保留目标序列的分类结果，并且只保留前7列
    target_classified = df[df['Genome'].isin(seqid_list)].iloc[:, :7]

    # 保存结果
    filtered_df.to_csv(overview_filtered, sep=',', index=False)
    target_classified.to_csv(final_tax, sep=',', index=False)

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

    invalid_values = {None, "", "Unknown", "Unassigned"}

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
    - 新增行的分类信息 Kingdom~Genus 为空
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
    tax_df = pd.read_csv(tax_file, sep=",")

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

    # vcontact2使用的数据库版本
    db = "ProkaryoticViralRefSeq211-Merged"

    # cpus
    threads = 10

    ##### step1 预测 ORF
    # 输入
    votus
    # 输出
    protein = Path("01_prodigal/protein.faa")
    log_prodigal = Path("00_log/prodigal.log")

    #run_prodigal(votus, protein, log_prodigal)

    ##### step2 运行 vcontact2, 生成 viral clusters (VCs)
    # 输入
    protein
    # 输出
    outdir_vcontact2 = Path("02_vcontact2")
    g2g = outdir_vcontact2 / "g2g.csv"
    overview = outdir_vcontact2 / "genome_by_genome_overview.csv"
    log_vcontact2 = Path("00_log/vcontact2.log")

    #run_vcontact2(protein, db, outdir_vcontact2, g2g, log_vcontact2, threads)

    ##### step3 分配分类信息
    # 输入
    g2g
    overview
    # 输出
    overview_filtered = outdir_vcontact2 / "overview_filtered.csv"
    final_tax = Path("03_tax/votus_taxonomy.csv")
    log_assign = Path("00_log/assign.log")

    assign_taxonomy(g2g, overview, overview_filtered, final_tax, log_assign)

    ##### step4 优化分类注释表
    # 输入
    votus
    final_tax
    # 输出
    output_tax_file = Path("votus_taxonomy.tsv")

    optimize_annotation_table(votus, final_tax, output_tax_file)

if __name__ == "__main__":
    main()
