#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# Date: 2026-01-15
#
# Description:
#       Phage taxonomy annotation pipeline:
#       virify (v3 版数据库和 v4 版分类信息2022版)
#       1. votus 输入 virify (--onlyannotate, --length 0, --viphog_version v3, --meta_version v4)
#   参考文献: 2022-An integrated detection, annotation and taxonomic classification pipeline using virus-specific protein profile hidden Markov models
# 
# Dependencies:
#       virify --version 3.2.0
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

###################################### step2 优化分类注释表
# 修改列名
def modify_column_names(tax_file, output_tax_file):
    """
    修改列名: 首字母大写
    """
    df = pd.read_csv(tax_file, sep='\t')
 
    # 改列名
    column_mapping = {
        'kingdom': 'Kingdom',
        'phylum': 'Phylum', 
        'subphylum': 'Subphylum',
        'class': 'Class',
        'order': 'Order',
        'suborder': 'Suborder',
        'family': 'Family',
        'subfamily': 'Subfamily',
        'genus': 'Genus'
    }
    
    df.rename(columns=column_mapping, inplace=True)

    df.to_csv(output_tax_file, sep='\t', index=False)

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


def main():
    """
    完整分析流程
    """
    # votus 路径
    votus = Path("/home/shijiabin/2025_2ME/03_virus_identification/votus.fna")

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
    tax_virify = virify_dir / votus.stem / "06-taxonomy" / f"{votus.stem}_renamed_filt0bp.fasta_prodigal_annotation_taxonomy.tsv"

    #run_virify(virify_path, votus, virify_dir, viphog_version, meta_version, db_dir)

    #### step2 优化分类注释表
    # 输入
    votus
    tax_virify
    # 输出
    final_tax = Path("votus_taxonomy.tsv")

    modify_column_names(tax_virify, final_tax)
    optimize_annotation_table(votus, final_tax, final_tax)

if __name__ == "__main__":
    main()
