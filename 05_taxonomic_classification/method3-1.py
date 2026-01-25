#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# Date: 2026-01-17
#
# Description:
#       Phage taxonomy annotation pipeline:
#       votus → Prodigal → DIAMOND (Virus-Host DB, release229) → taxonomy assignmen
#       1. 预测 ORF
#       2. 蛋白去冗余
#       3. 非冗余蛋白比对到 Virus-Host DB (使用带有数据库名称前缀的蛋白数据库)
#       4. 分配分类学注释 (virusTaxAnno.sh 内置)
#       5. 优化注释表
#   参考文献: (思路) 2025-The Chinese gut virus catalogue reveals gut virome diversity and disease-related viral signatures (与原文不同: 不跑 geNomad)
#             (代码) 2024-A multi-kingdom collection of 33,804 reference genomes for the human vaginal microbiome (基于比对结果的分类注释方法的具体代码)
# Dependencies:
#       prodigal --version 2.6.3
#       MMseqs2 Version: bd01c2229f027d8d8e61947f44d11ef1a7669212 
#       diamond --version 2.0.6
#       virusTaxAnno https://github.com/RChGO/VMGC/blob/main/Pipelines/virusTaxAnno.sh

##### 注意:
#    1. 根据比对结果分配注释的代码 https://github.com/RChGO/VMGC/blob/main/Pipelines/virusTaxAnno.sh
#    2. 真的需要蛋白质去冗余吗? (目前去冗余)

import os
import subprocess
import pandas as pd
import logging
from pathlib import Path
from Bio import SeqIO
import shutil

###################################### step1 预测 ORF
def run_prodigal(input_virus, output_protein):
    os.makedirs(os.path.dirname(output_protein), exist_ok=True)
    cmd = f"prodigal -i {input_virus} -a {output_protein} -p meta"
    subprocess.run(cmd, shell=True, check=True)

###################################### step2 蛋白去冗余
def run_mmseqs_linclust_protein(protein_fasta, out_dir, prefix_name, tmp_dir):
    """
    Run MMseqs2 easy-linclust for protein clustering and extract non-redundant proteins.

    Parameters
    ----------
    protein_fasta : Path
        Input protein FASTA file.
    out_dir : Path
        Output directory for MMseqs2 results.
    prefix_name : str
        Prefix name for MMseqs2 outputs.
    tmp_dir : Path
        Temporary directory for MMseqs2.
    """

    protein_fasta = Path(protein_fasta)
    out_dir = Path(out_dir)
    tmp_dir = Path(tmp_dir)

    prefix = out_dir / prefix_name
    rep_seq = out_dir / f"{prefix_name}_rep_seq.fasta"
    nonredundant_protein = out_dir / "nonredundant_protein.faa"

    # 创建输出目录
    out_dir.mkdir(parents=True, exist_ok=True)

    # MMseqs2 clustering
    cmd = ["mmseqs", "easy-linclust",
            str(protein_fasta),
            str(prefix),
            str(tmp_dir),
            "--min-seq-id", "0.9", 
            "--cov-mode", "1", 
            "-c", "0.8", 
            "--kmer-per-seq", "80"]

    subprocess.run(cmd, check=True)

    # 重命名代表序列为 nonredundant_protein.faa
    if not rep_seq.exists():
        raise FileNotFoundError(f"Expected MMseqs output not found: {rep_seq}")

    rep_seq.rename(nonredundant_protein)

    # 清理临时目录
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)

    return nonredundant_protein

###################################### step3 分配分类信息
def run_virusTaxAnno(virusTaxAnno, prot, outf, db, db_tax, db_tax_ref, threads):
    '''
    virusTaxAnno:   脚本的位置
    prot:           输入蛋白质序列文件(.faa)
    outf:           输出文件前缀
    db:             要比对的蛋白数据库(.faa, 脚本内会临时建库)
    db_tax:         蛋白质序列数据库分类注释表
    db_tax_ref:     蛋白质序列数据库拥有完整等级的分类注释表
    threads:        线程数
    '''
    # 创建输出目录
    os.makedirs(os.path.dirname(outf), exist_ok=True)

    # 使用 virusTaxAnno 脚本分配 family 注释
    cmd = f"{virusTaxAnno} {prot} {outf} {db} {db_tax} {threads}"
    subprocess.run(cmd, shell=True, check=True)

    # 读取分配结果
    df = pd.read_csv(f"{outf}.tax_family", sep=r"\s+", header=None)
    # 取第一列和最后一列
    df_sub = df.iloc[:, [0, -1]]
    # 重命名列名
    df_sub.columns = ["vOTU", "family"]

    # 根据 family 匹配 db_tax_ref 中对应的其他等级
    # 读取数据库的分类注释表
    df_tax_ref = pd.read_csv(db_tax_ref, sep="\t", header=0)
    # 对 family 去重，保留第一个
    df_family_map = (
        df_tax_ref
        .loc[:, ["family", "phylum", "class", "order"]]
        .drop_duplicates(subset="family", keep="first")
    )

    # 合并分配结果表和数据库分类注释表
    df_final = df_sub.merge(
        df_family_map,
        left_on="family",
        right_on="family",
        how="left"
    )

    # 整理列顺序 & 删除冗余列
    df_final = df_final[
        ["vOTU", "phylum", "class", "order", "family"]
    ]

    df_final.to_csv(f"{outf}.tax_all", sep='\t', index=False)

###################################### step4 优化分类注释表
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

    提取第一列和 Phylum~Family 之间的列(要找到这两个列名然后提取他们两以及之间的列)
    由高等级开始处理(第二列), 如果第二列为空,那么填充为 Unknown;
    然后处理下一列,如果为空或 Unknown,则判断上一列,若上一列为 Unknown,则填充为 Unknown;
    若上一列不为 Unknown 且不以 Unclassified 开头,则填充为 Unclassified + "上一列的值";
    若上一列以 Unclassified 开头,则直接复制上一列的值;若本列不为空或 Unknown,则保留本列的值. 
    再下一列的处理同第三列
    """

    #### step1 提取原分类注释表的第一列 + Phylum~Family 列
    # 读取分类注释表
    tax_df = pd.read_csv(tax_file, sep="\t")

    # 找到 Phylum 和 Family 的列索引
    cols = tax_df.columns.tolist()
    start_col = cols.index("Phylum")
    end_col = cols.index("Family")

    # 提取第一列 + Phylum~Species, 得到数据框用于合并
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
    # 添加 Genus 和 Species 列，值为空字符串
    merged_df["Genus"] = ""
    merged_df["Species"] = ""
    
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

    # 数据库
    virushostdb_dir = Path("/home/shijiabin/db/virushostdb_release229")
    virushostdb_path = virushostdb_dir / "virushostdb_formated.faa"       # 使用带有数据库名称前缀的蛋白数据库
    virushostdb_tax_path = virushostdb_dir / "taxonomy.tsv"               # 数据库蛋白id与family的映射
    virushostdb_tax_ref_path = virushostdb_dir / "taxonomy_ref.tsv"       # 数据库蛋白id与分类信息的映射,用于补充完整等级

    threads = 20

    virusTaxAnno = Path("/home/shijiabin/2025_2ME/scripts/VMGCPipelines/virusTaxAnno_modified.sh")

    #### step1 预测 ORF
    # 输入
    votus
    # 输出
    protein = Path("01_prodigal/protein.faa")

    run_prodigal(votus, protein)
    
    #### step2 蛋白去冗余
    # 输入
    protein
    # 输出
    dir_mmseqs = Path("02_mmseqs2")
    nonredundant_protein = dir_mmseqs / "nonredundant_protein.faa"

    run_mmseqs_linclust_protein(protein, dir_mmseqs, "protein", "tmp_dir")

    #### step3 分配分类信息
    # 输入
    nonredundant_protein    # 输入蛋白质序列文件(.faa)
    virushostdb_path        # 蛋白质序列数据库(.faa)
    virushostdb_tax_path    # 蛋白质序列数据库分类信息 
    # 输出
    outf = "03_assign/votus"	    # 输出文件前缀
    tax_assign = f"{outf}.tax_all"

    run_virusTaxAnno(virusTaxAnno, nonredundant_protein, outf, virushostdb_path, virushostdb_tax_path, virushostdb_tax_ref_path, threads)
    
    #### step4 优化分类注释表
    # 输入
    tax_assign
    votus
    # 输出
    final_tax = "votu_tax.tsv"

    modify_column_names(tax_assign, final_tax)
    optimize_annotation_table(votus, final_tax, final_tax)

if __name__ == "__main__":
    main()
