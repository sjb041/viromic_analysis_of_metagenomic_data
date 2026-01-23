#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# Date: 2026-01-10
#
# Description:
#       Phage taxonomy annotation pipeline:
#       votus --> PhaBOX2-phagcn (phabox_db_v2_1)
#       1. votus 输入 phagcn (--task phagcn, --len 0)
#       2. 优化分类注释表
#   参考文献: 2024-The phageome of patients with ulcerative colitis treated with donor fecal microbiota reveals markers associated with disease remission
# 
# Dependencies:
#       PhaBOX2 --version 2.1.13
######################################

import os
import pandas as pd
from pathlib import Path
import subprocess
from Bio import SeqIO

###################################### step1 PhaBOX2-phagcn 分类注释
def run_phagcn(virus, db, min_len, output_dir, threads):
    os.makedirs(output_dir, exist_ok=True)
    cmd1 = f"micromamba run -n phabox2 phabox2 --task phagcn --contigs {virus} --dbdir {db} --outpth {output_dir} --len {min_len} --threads {threads}"
    subprocess.run(cmd1, shell=True, check=True)
###################################### step2 提取各分类等级, 并处理空值
def process_lineage(phagcn_prediction, output_tax_file, lineage_col="Lineage"):
    """
    从 Lineage 列解析出各分类等级，并处理空值。
    
    - 提取等级: kingdom, phylum, class, order, family, subfamily, genus
    - 空值处理逻辑: 从 genus 开始依次往左，
        如果该等级为空:
            - 上一级也为空 -> "unknown"
            - 上一级不为空 -> "Unclassified{上一级值}"

    修改列名:
    kingdom -> Kingdom
    phylum -> Phylum
    class -> Class
    order -> Order
    family -> Family
    subfamily -> Subfamily
    genus -> Genus
    """

    ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus']
    
    df = pd.read_csv(phagcn_prediction, sep='\t')
    # -------- 解析 Lineage 字符串为 dict --------
    def parse_lineage(lineage_str):
        items = str(lineage_str).split(";")
        rank_dict = {}
        for item in items:
            if ":" in item:
                k, v = item.split(":", 1)
                rank_dict[k.strip()] = v.strip()
        return rank_dict
    
    parsed = df[lineage_col].apply(parse_lineage)
    
    # -------- 提取目标等级 --------
    for r in ranks:
        df[r] = parsed.apply(lambda d: d.get(r, None))
    
    # -------- 空值回填处理 --------
    for i in range(len(ranks)-1, -1, -1):  # 从右往左
        rank = ranks[i]
        for idx in df.index:
            if pd.isna(df.at[idx, rank]) or df.at[idx, rank] == "":
                if i == 0:  # kingdom 特殊情况
                    df.at[idx, rank] = "unknown"
                else:
                    parent_rank = ranks[i-1]
                    parent_val = df.at[idx, parent_rank]
                    if pd.isna(parent_val) or parent_val == "":
                        df.at[idx, rank] = "unknown"
                    else:
                        df.at[idx, rank] = f"Unclassified{parent_val}"
    
    cols_to_keep = ["Accession", "Prokaryotic virus (Bacteriophages and Archaeal virus)", "Lineage", "kingdom", "phylum", "class", "order", "family", "subfamily", "genus"]
    df_out = df[cols_to_keep]

    # 改列名
    column_mapping = {
        'kingdom': 'Kingdom',
        'phylum': 'Phylum', 
        'class': 'Class',
        'order': 'Order',
        'family': 'Family',
        'subfamily': 'Subfamily',
        'genus': 'Genus'
    }
    
    df_out.rename(columns=column_mapping, inplace=True)

    df_out.to_csv(output_tax_file, sep='\t', index=False)

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

    # vcontact2使用的数据库版本
    db = Path("/home/shijiabin/db/phabox2_db/phabox_db_v2_1")

    # cpus
    threads = 10

    # 最小序列长度, 0 表示不限制
    min_len = 0    


    #### step1 PhaBOX2-phagcn 分类注释
    # 输入
    votus
    # 输出
    outdir_phagcn = Path("01_phagcn")
    tax_phagcn = outdir_phagcn / "final_prediction" / "phagcn_prediction.tsv"

    # 仅运行分类注释
    run_phagcn(votus, db, min_len, outdir_phagcn, threads)

    #### step2 提取各分类等级, 并处理空值
    # 输入
    tax_phagcn
    # 输出
    tax_formatted = outdir_phagcn / "final_prediction" / "taxonomy_formatted.tsv"

    # 提取各分类等级
    process_lineage(tax_phagcn, tax_formatted, "Lineage")

    #### step3 优化分类注释表
    # 输入
    votus
    tax_formatted
    # 输出
    final_tax = Path("votus_taxonomy.tsv")

    optimize_annotation_table(votus, tax_formatted, final_tax)


if __name__ == "__main__":
    main()
