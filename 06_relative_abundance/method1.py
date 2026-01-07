#!/usr/bin/env python3
import subprocess
import os
import pandas as pd

######################################
# Date: 2026-01-07
#
# Description:
#   计算相对丰度,公式1:
#       每对样本随机采样 1500 对 reads,
#       使用 bowtie2 比对每一对样本的 reads 到 votus,
#       基因组长度归一化: count_per_bp = mapped_reads / genome_length,
#       相对丰度计算: rel_abun = count_per_bp / sum(all_count_per_bp)
#   参考文献: The Chinese gut virus catalogue reveals gut virome diversity and disease-related viral signatures   
# 
# Dependencies:
#   BBMap/reformat.sh --version 39.61
#   bowtie2 --version 2.3.5.1
#   samtools --version 1.10
######################################

###################################### step1 采样 1500 万条 reads
def sampling_reads(bbmap_reformat, in1, in2, out1, out2, num_reads, seed, threads):
    """
    使用 BBMap/reformat.sh 随机成对采样 reads
    """
    # 采样
    cmd = (
        f"{bbmap_reformat} in={in1} in2={in2} out={out1} out2={out2} "
        f"samplereadstarget={num_reads} "
        f"sampleseed={seed} threads={threads}"
    )
    subprocess.run(cmd, shell=True, check=True)
######################################

###################################### step2 bowtie2 比对
def run_bowtie2_mapping(votus, reads1, reads2, output_dir, outputfile_prefix, threads):
    """
    使用 Bowtie2 将 reads 比对到 vOTUs, 并输出 SAM 文件
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
def process_sam_to_bam(sam_file, output_dir, outputfile_prefix, threads):
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
######################################

###################################### step3 计算相对丰度
def parse_idxstats(idxstats_file, sample_name, abun_file):
    """
    解析 samtools idxstats 输出文件, 计算相对丰度
    输出两个表格: 
        相对丰度表, 2 列, vOTU    rel_abun_{sample_name}
        计算过程表, 6 列, vOTU    length    mapped    unmapped    count_per_bp    rel_abun
    """
    data = pd.read_csv(idxstats_file, sep='\t', header=None, 
                       names=['vOTU', 'length', 'mapped', 'unmapped'])
    # 过滤掉星号行
    data = data[data['vOTU'] != '*']

    # 基因组长度归一化: count_per_bp = mapped_reads / genome_length
    data['count_per_bp'] = data['mapped'] / data['length']
    
    # 计算总 count_per_bp 的和，用于计算相对丰度
    total_count_per_bp = data['count_per_bp'].sum()
    
    # 相对丰度计算: rel_abun = count_per_bp / sum(all_count_per_bp)
    if total_count_per_bp > 0:
        data['rel_abun'] = data['count_per_bp'] / total_count_per_bp
    else:
        data['rel_abun'] = 0  # 如果总和为0，所有相对丰度设为0
    
    # 构建输出列名
    rel_abun_col = f"rel_abun_{sample_name}"
    
    # 选择输出列
    result = data[['vOTU', 'rel_abun']].copy()
    result.columns = ['vOTU', rel_abun_col]
    
    # 输出到文件
    result.to_csv(abun_file, sep='\t', index=False)

    # 输出计算过程
    data.to_csv(f"{abun_file}.tmp", sep='\t', index=False)
def merge_relative_abundance(rel_abun_file_list, output_file):
    """
    合并多个样本中vOTU的相对丰度表,要合并的相对丰度表格式为 <vOTU> <*> <*> <rel_abun_sample>
    
    Args:
        rel_abun_file_list (list): 相对丰度文件路径列表
        output_file (str): 输出合并后的文件路径
    
    输出文件格式：
        vOTU rel_abun_sample1 rel_abun_sample2 rel_abun_sample3 .....
    """

    # 收集所有文件中的vOTU
    all_votus = set()
    all_dataframes = []

    # 读取所有文件并收集vOTU
    for file_path in rel_abun_file_list:
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

    # === 合并后统一检查列和是否为1 ===
    rel_cols = [c for c in merged_df.columns if c.startswith("rel_abun_")]

    for col in rel_cols:
        col_sum = merged_df[col].sum()

        # 明显错误：列和 > 1
        if col_sum > 1.0 + 1e-2:
            raise ValueError(
                f"合并后列 {col} 的相对丰度之和为 {col_sum:.6f} (>1)"
            )

        # 合理但需要注意的情况
        if col_sum < 1.0 - 1e-2:
            raise ValueError(
                f"合并后列 {col} 的相对丰度之和为 {col_sum:.6f} (<1)"
            )

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

    # BBMap/reformat.sh
    bbmap_reformat = "/home/shijiabin/opt/bbmap/reformat.sh"

    threads = 20

    # 逐个计算 votus 在各个样本下的相对丰度
    abun_files = [] 
    for sample in SAMPLES:

        #### step1 采样 1500 万对 reads 
        # 输入
        in1 = os.path.join(reads_dir, f"{sample}_R1.fastq.gz")
        in2 = os.path.join(reads_dir, f"{sample}_R2.fastq.gz")
        # 输出
        dir_sampling = "00_sampling_reads"
        os.makedirs(dir_sampling, exist_ok=True)
        out1 = f"{dir_sampling}/{sample}_R1.fastq.gz"
        out2 = f"{dir_sampling}/{sample}_R2.fastq.gz"
        # 采样
        sampling_reads(bbmap_reformat, in1, in2, out1, out2, num_reads=15000000, seed=42, threads=threads)

        #### step2 使用 bowtie2 , 比对每一个样本的 reads 到 votus
        # 输入
        reads1 = out1
        reads2 = out2
        # 输出
        dir_mapping = "01_mapping"
        os.makedirs(dir_mapping, exist_ok=True)
        # 比对
        sam = run_bowtie2_mapping(votus, reads1, reads2, dir_mapping, sample, threads)
        sam = os.path.join(dir_mapping, f"{sample}.sam")
        # 统计比对信息
        process_sam_to_bam(sam, dir_mapping, sample, threads)

        #### step3 计算 votus 在该 reads 样本中的相对丰度
        # 输入
        idxstats_file = f"{dir_mapping}/{sample}.idxstats"
        # 输出
        dir_abun = "02_abun"
        os.makedirs(dir_abun, exist_ok=True)
        abun_file = f"{dir_abun}/{sample}_abun.tsv"
        # 计算相对丰度
        parse_idxstats(idxstats_file, sample, abun_file)
        # 在该样本的下的相对丰度,添加进相对丰度文件列表
        abun_files.append(abun_file)

    # 合并 votus 在各个样本下的相对丰度    
    merge_relative_abundance(abun_files, "votu_abun.tsv")

if __name__ == "__main__":
    main()
