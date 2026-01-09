#!/usr/bin/env python3
import subprocess
import os
import pandas as pd

######################################
# Date: 2026-01-09
#
# Description:
#   计算 raw_counts, 使用公式3的 idxstats_file:
#       raw_counts 计算: raw_counts = mapped_reads
######################################

###################################### step3 计算 raw_counts
def parse_idxstats(idxstats_file, sample_name, abun_file):
    """
    解析 samtools idxstats 输出文件, 计算 raw_counts
    输出两个表格: 
        raw_counts 表, 2 列, vOTU  raw_counts_{sample_name}
    """
    data = pd.read_csv(idxstats_file, sep='\t', header=None, 
                       names=['vOTU', 'length', 'mapped', 'unmapped'])
    # 过滤掉星号行
    data = data[data['vOTU'] != '*']

    # 选择输出列
    result = data[['vOTU', 'mapped']].copy()
    result.columns = ['vOTU', f"raw_counts_{sample_name}"]
    
    # 输出到文件
    result.to_csv(abun_file, sep='\t', index=False)

def merge_raw_counts(raw_counts_file_list, output_file):
    """
    合并多个样本中vOTU的raw_counts表,要合并的相对丰度表格式为 <vOTU> <*> <*> <raw_counts_sample>
    
    Args:
        raw_counts_file_list (list): raw_counts 文件路径列表
        output_file (str): 输出合并后的文件路径
    
    输出文件格式：
        vOTU raw_counts_sample1 raw_counts_sample2 raw_counts_sample3 .....
    """

    # 收集所有文件中的vOTU
    all_votus = set()
    all_dataframes = []

    # 读取所有文件并收集vOTU
    for file_path in raw_counts_file_list:
        df = pd.read_csv(file_path, sep='\t')
        all_dataframes.append(df)
        all_votus.update(df['vOTU'].tolist())

    # 创建完整的vOTU列表
    all_votus = sorted(list(all_votus))
    # 创建结果DataFrame
    merged_df = pd.DataFrame({'vOTU': all_votus})

    # 依次处理每个文件的数据
    for df in all_dataframes:
        rel_cols = [c for c in df.columns if c.startswith("raw_counts_")]
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

    # idxstats 路径
    dir_mapping = "/home/shijiabin/2025_2ME/06_relative_abundance/method3/01_mapping"

    # 逐个计算 votus 在各个样本下的 raw_counts
    file_list = [] 
    for sample in SAMPLES:

        #### step1 计算 votus 在该 reads 样本中的 raw_counts
        # 输入
        idxstats_file = f"{dir_mapping}/{sample}.idxstats"
        # 输出
        dir = "01_raw_counts"
        os.makedirs(dir, exist_ok=True)
        raw_counts_file = f"{dir}/{sample}_raw_counts.tsv"
        # 计算 raw_counts
        parse_idxstats(idxstats_file, sample, raw_counts_file)
        # 在该样本的下的 raw_counts,添加进相对丰度文件列表
        file_list.append(raw_counts_file)

    # 合并 votus 在各个样本下的 raw_counts    
    merge_raw_counts(file_list, f"{dir}/votu_raw_counts.tsv")

if __name__ == "__main__":
    main()
