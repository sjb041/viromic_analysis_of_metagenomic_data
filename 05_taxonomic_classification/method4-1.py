#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# Date: 2026-01-20
#
# Description:
#       Phage taxonomy annotation pipeline:
#       votu → MMseqs2(easy-taxonomy)
#       1. votu 输入 MMseqs2 (mmseqs easy-taxonomy)
#       2. 优化分类注释表
#   参考文献: 2021-Metagenomic compendium of 189,680 DNA viruses from the human gut microbiome
#   参考代码: https://github.com/snayfach/MGV/issues/1
# 
# Dependencies:
#       MMseqs2 Version: bd01c2229f027d8d8e61947f44d11ef1a7669212
######################################

import os
import subprocess
import pandas as pd
import logging
from pathlib import Path
from Bio import SeqIO

###################################### step1 MMseqs2 分类注释 
def run_mmseqs_easy_taxonomy(input_virus, db, output_prefix, threads):
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    tmp = os.path.join(os.path.dirname(output_prefix), "tmp")
    
    # 不加双引号会导致目录解析不正确
    cmd = [
        "mmseqs", "easy-taxonomy",
        input_virus,
        db,
        output_prefix,
        tmp,
        "--blacklist", "",
        "--tax-lineage", "0",
        "--lca-ranks", "kingdom,phylum,class,order,family,genus,species",
        "--threads", str(threads)
    ]

    subprocess.run(cmd, check=True)

###################################### step2 优化分类注释表
# 格式化 MMseqs2 分类注释表
def format_annotation_table(tax_file, output_tax_file): 
    # 读取 tax_file,提取第一列和第九列
    df = pd.read_csv(tax_file, sep="\t", header=None)
    df = df[[0, 8]]
    df.columns = ["vOTU", "taxonomy"]

    # 按 ; 拆分第九列, 按顺序装入 Kingdom,Phylum,Class,Order,Family,Genus,Species 列
    df[["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]] = df["taxonomy"].str.split(";", expand=True)

    # 把各列中以 "uc_" 开头的值修改为以 "Unclassified" 开头
    df = df.replace("uc_", "Unclassified", regex=True)
    df = df.drop(columns=["taxonomy"])
    df.to_csv(output_tax_file, sep="\t", index=False)
    return df

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

# =============================
# Main Workflow
# =============================

def main():
    """
    完整分析流程
    """
    # votus 路径
    votus = Path("/home/shijiabin/2025_2ME/03_virus_identification/votus.fna")

    # 参考数据库
    db = "/home/shijiabin/db/mmseqs_db/ictv_nr_db/ictv_nr_db"

    # cpus
    threads = 20

    ##### step1 MMseqs2 分类注释
    # 输入
    votus
    output_prefix = "01_mmseqs_taxonomy/votus"
    # 输出
    tax_assign = Path("01_mmseqs_taxonomy/votus_lca.tsv")

    run_mmseqs_easy_taxonomy(votus, db, output_prefix, threads)

    ##### step2 优化分类注释表
    # 输入
    tax_assign = Path("01_mmseqs_taxonomy/votus_lca.tsv")
    # 输出
    output_tax_file = Path("votus_taxonomy.tsv")

    format_annotation_table(tax_assign, output_tax_file)
    optimize_annotation_table(votus, output_tax_file, output_tax_file)

if __name__ == "__main__":
    main()
