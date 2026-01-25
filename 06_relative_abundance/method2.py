#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date: 2026-01-26

Description:
    计算相对丰度,公式2:
        1. 使用 bowtie2 比对每一个样本的 reads 到 votus
        2. 从原 sam 文件中提取出 45 万个不同的 read 名,然后根据这些 reads 名筛选原 sam 文件
        3. 基于筛选后的 sam 文件重新计算每个 vOTU 的 read 数量
        4. 每个样本中每个 vOTU 的 相对丰度 定义为其 read 数量除以 450000
    参考文献: 2025-A metagenome-wide study of the gut virome in chronic kidney disease  

Dependencies:
    bowtie2 --version 2.3.5.1
    samtools --version 1.10

A:
{'F1_1A': 613338, 'F1_2A': 479386, 'F1_3A': 643799, 'F2_1A': 1334157, 'F2_2A': 662320, 'F2_3A': 882047, 'FG_1A': 3132991, 'FG_2A': 986908, 'FG_3A': 909400, 'L1_1A': 1347670, 'L1_2A': 825013, 'L1_3A': 1447137, 'L2_1A': 727399, 'L2_2A': 850688, 'L2_3A': 1415388, 'H_LX': 482631, 'H_O': 627098}
Recommended subsample size = 479386

E:
{'F1_1A': 704363, 'F1_2A': 558437, 'F1_3A': 816471, 'F2_1A': 1486359, 'F2_2A': 784834, 'F2_3A': 1012486, 'FG_1A': 3260976, 'FG_2A': 1115247, 'FG_3A': 1063891, 'L1_1A': 1572231, 'L1_2A': 1013624, 'L1_3A': 1752509, 'L2_1A': 1003769, 'L2_2A': 948206, 'L2_3A': 1642161, 'H_LX': 467012, 'H_O': 598872}
Recommended subsample size = 467012
"""

import os
import subprocess
import random
from collections import Counter
import pandas as pd
import argparse

# =============================
# Utility Functions
# =============================

def get_mapped_qnames(samfile):
    """从SAM文件中提取所有primary mapped read名"""
    qnames = set()
    with open(samfile, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            if rname == "*" or (flag & 0x4) or (flag & 0x100) or (flag & 0x800):
                continue
            qnames.add(qname)
    return qnames

def subsample_qnames(qnames, k, seed=42):
    """随机抽取k个read名"""
    random.seed(seed)
    if len(qnames) <= k:
        return qnames
    return set(random.sample(list(qnames), k))

def filter_sam_by_qnames(samfile, qnames, outfile):
    """根据QNAME集合筛选SAM文件(保留头部)"""
    with open(samfile, "r") as fin, open(outfile, "w") as fout:
        headers = []
        for line in fin:
            if line.startswith("@"):
                headers.append(line)
                continue
            fields = line.strip().split("\t")
            if not fields or len(fields) < 3:
                continue
            qname = fields[0]
            if qname in qnames:
                # 写入头部（仅在第一次遇到匹配时）
                if headers:
                    fout.writelines(headers)
                    headers = []  # 确保只写一次
                fout.write(line)

def reservoir_sample(iterator, k, seed=42):
    """Reservoir sampling: 从大数据流中随机抽取 k 条"""
    r = random.Random(seed)
    result = []
    for i, item in enumerate(iterator, start=1):
        if i <= k:
            result.append(item)
        else:
            j = r.randint(1, i)
            if j <= k:
                result[j-1] = item
    return result

def valid_reads(samfile):
    """迭代器: 从无头SAM文件读取 primary mapped read 的 reference_name"""
    with open(samfile, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue
            flag = int(fields[1])
            rname = fields[2]
            if rname == "*" or rname == "":
                continue
            if flag & 0x4 or flag & 0x100 or flag & 0x800:
                continue
            yield rname

def count_votu_reads(samfile, subsample=None, seed=42):
    counts = Counter()
    if subsample:
        sampled = reservoir_sample(valid_reads(samfile), subsample, seed=seed)
        counts.update(sampled)
        total_reads = subsample
    else:
        total_reads = 0
        for rname in valid_reads(samfile):
            counts[rname] += 1
            total_reads += 1
    return counts, total_reads

def count_mapped_reads(samfile):
    """统计单个SAM文件中唯一 mapped reads 的数量（去重QNAME）"""
    mapped_qnames = set()
    with open(samfile, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            if rname == "*" or (flag & 0x4) or (flag & 0x100) or (flag & 0x800):
                continue
            mapped_qnames.add(qname)
    return len(mapped_qnames)

def summarize_mapped_reads(samfiles):
    """打印每个样本的 mapped read 数，并找出最小值"""
    results = {}
    for sam in samfiles:
        sample_name = os.path.splitext(os.path.basename(sam))[0]
        mapped = count_mapped_reads(sam)
        results[sample_name] = mapped
        print(f"{sample_name:20s}  {mapped:,} mapped reads")
    min_mapped = min(results.values()) if results else 0
    print(f"\n>> Minimum mapped reads across samples: {min_mapped:,}\n")
    return results, min_mapped

# =============================
# Core Analysis Functions
# =============================

def run_bowtie2_mapping(votus, reads1, reads2, output_dir, outputfile_prefix, threads=16):
    """
    使用 Bowtie2 将 reads 比对到 vOTUs，并输出 SAM 文件
    """
    os.makedirs(output_dir, exist_ok=True)
    index = os.path.join(output_dir, "index")
    sam = os.path.join(output_dir, f"{outputfile_prefix}.sam")

    # 检查索引
    index_files = [f"{index}.{ext}.bt2" for ext in [1, 2, 3, 4]] + \
                  [f"{index}.rev.{ext}.bt2" for ext in [1, 2]]
    if not all(os.path.exists(f) for f in index_files):
        subprocess.run(f"bowtie2-build {votus} {index}", shell=True, check=True)

    # Bowtie2 比对
    cmd = (
        f"bowtie2 -x {index} -1 {reads1} -2 {reads2} "
        f"--end-to-end --fast --no-unal "
        f"-p {threads} -S {sam}"
    )
    subprocess.run(cmd, shell=True, check=True)
    return sam

def process_sam_to_bam(sam_file, output_dir, outputfile_prefix, threads=16):
    """
    将SAM文件转换为BAM格式, 排序, 建立索引并生成统计信息
    
    Args:
        sam_file (str): 输入SAM文件路径
        output_dir (str): 输出目录路径
        outputfile_prefix (str): 输出文件前缀
    
    Output files:
        - {output_dir}/{outputfile_prefix}.sorted.bam: 排序后的比对结果
        - {output_dir}/{outputfile_prefix}.idxstats.tsv: 比对统计信息, 格式: <refname> <ref_len> <mapped> <unmapped>
        - {output_dir}/{outputfile_prefix}.flagstat: 比对统计信息(基于reads数量)
    """
    
    # 构建输出路径
    bam = os.path.join(output_dir, f"{outputfile_prefix}.sorted.bam")
    idxstats = os.path.join(output_dir, f"{outputfile_prefix}.idxstats")
    flagstat = os.path.join(output_dir, f"{outputfile_prefix}.flagstat")
    
    # 将SAM转换为BAM格式并排序
    cmd1 = f"samtools view -bS {sam_file} | samtools sort -@ {threads} -o {bam}"
    subprocess.run(cmd1, shell=True, check=True)
    
    # 建立索引并生成idxstats统计信息
    cmd2 = f"samtools index {bam} && samtools idxstats {bam} > {idxstats}"
    subprocess.run(cmd2, shell=True, check=True)
    
    # 生成flagstat统计信息(基于reads数量)
    cmd3 = f"samtools flagstat {bam} > {flagstat}"
    subprocess.run(cmd3, shell=True, check=True)

def votu_abundance(samfile, output, sample, subsample=None, seed=42):
    """
    计算每个 vOTU 的 reads 计数与相对丰度(read_count / subsample 或 / total_reads)
    """
    counts, total_reads = count_votu_reads(samfile, subsample=subsample, seed=seed)

    with open(output, "w") as out:
        colname = f"rel_abun_{sample}"
        out.write("vOTU\tReadCount\t{}\n".format(colname))
        for votu, cnt in counts.items():
            ra = cnt / total_reads if total_reads > 0 else 0
            out.write(f"{votu}\t{cnt}\t{ra:.6f}\n")

    return output

# 合并 votus 在所有样本下的相对丰度表
def merge_relative_abundance(relative_abundance_files, output_file):
    """
    合并多个样本中vOTU的相对丰度表,要合并的相对丰度表格式为 <vOTU> <*> <*> <rel_abun_sample>
    
    Args:
        relative_abundance_files (list): 相对丰度文件路径列表
        output_file (str): 输出合并后的文件路径
    
    输出文件格式：
        vOTU rel_abun_sample1 rel_abun_sample2 rel_abun_sample3 .....
    """

    # 收集所有文件中的vOTU
    all_votus = set()
    all_dataframes = []

    # 读取所有文件并收集vOTU
    for file_path in relative_abundance_files:
        df = pd.read_csv(file_path, sep='\t')
        all_dataframes.append(df)
        all_votus.update(df['vOTU'].tolist())

    # 创建完整的vOTU列表
    all_votus = sorted(list(all_votus))
    # 创建结果DataFrame
    merged_df = pd.DataFrame({'vOTU': all_votus})

    # 依次处理每个文件的数据
    for df in all_dataframes:
        rel_cols = [c for c in df.columns if c.startswith("rel_abun_")]
        for col in rel_cols:
            df_indexed = df.set_index('vOTU')[col].reindex(all_votus, fill_value=0)
            merged_df[col] = df_indexed.values

    merged_df.to_csv(output_file, sep='\t', index=False)
    return output_file

# =============================
# Main Workflow
# =============================

def main():
    """完整分析流程"""

    # 样本名
    SAMPLES = [
        "F1_1A", "F1_2A", "F1_3A",
        "F2_1A", "F2_2A", "F2_3A",
        "FG_1A", "FG_2A", "FG_3A",
        "L1_1A", "L1_2A", "L1_3A",
        "L2_1A", "L2_2A", "L2_3A",
        "H_LX", "H_O"
    ]

    # votus 路径
    votus = "/home/shijiabin/2025_2ME/03_virus_identification/votus.fna"

    # reads 目录,该目录下的 reads 文件名格式为 sample_R1.fastq.gz, sample_R2.fastq.gz
    reads_dir = "/home/shijiabin/2025_2ME/01_raw_sequence_preprocessing/03_cleandata"

    threads = 20

    # 构建输出目录
    output_dir_map = "01_mapping"
    os.makedirs(output_dir_map, exist_ok=True)
    output_dir_abun = "03_abun"
    os.makedirs(output_dir_abun, exist_ok=True)

    rel_files = []
    sam_files = []
    subsampled_sam_files = []
    for sample in SAMPLES:

        # 步骤1. 使用 bowtie2 , 比对每一个样本的 reads 到 votus
        reads1 = os.path.join(reads_dir, f"{sample}_R1.fastq.gz")
        reads2 = os.path.join(reads_dir, f"{sample}_R2.fastq.gz")
        #sam = run_bowtie2_mapping(votus, reads1, reads2, output_dir_map, sample, threads)
        sam = os.path.join(output_dir_map, f"{sample}.sam")
        sam_files.append(sam)

        # 步骤2. 统计比对信息
        #process_sam_to_bam(sam, output_dir_map, sample, threads)

        # 步骤3. 子抽样 sam 文件
        qnames = get_mapped_qnames(sam)
        selected = subsample_qnames(qnames, 450000)
        output_dir_submapping = "02_submapping"
        os.makedirs(output_dir_submapping, exist_ok=True)
        subsampled_sam = os.path.join(output_dir_submapping, f"{sample}.subsampled.sam")
        filter_sam_by_qnames(sam, selected, subsampled_sam)
        subsampled_sam_files.append(subsampled_sam)

        # 步骤3. 计算 votus 在该 reads 样本中的相对丰度
        rel_file = os.path.join(output_dir_abun, f"{sample}_abun.tsv")
        votu_abundance(sam, rel_file, sample, subsample=450000)
        rel_files.append(rel_file)

    # 扫一遍所有样本，统计各自的 mapped read 数量
    stats, min_reads = summarize_mapped_reads(sam_files)
    print(stats)
    print(f"Recommended subsample size = {min_reads}")

    # 合并该样本的 votus 在所有样本下的相对丰度表
    merged_file = merge_relative_abundance(rel_files, "votu_abun.tsv")
if __name__ == "__main__":
    main()
