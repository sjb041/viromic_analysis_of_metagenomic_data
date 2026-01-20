#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# Date: 2026-01-17
#
# Description:
#       Phage taxonomy annotation pipeline:
#       votus → Prodigal → DIAMOND (Virus-Host DB, release229) → taxonomy assignmen
#         ↓
#       genomad
#       1. genomad 分类注释
#       2. 预测 ORF
#       3. 蛋白去冗余
#       4. 非冗余蛋白比对到 Virus-Host DB
#       5. 基于diamond结果分配科分类和完整分类
#       6. 合并基于diamond的分类信息和基于genomad的分类信息
#   参考文献: 2025-The Chinese gut virus catalogue reveals gut virome diversity and disease-related viral signatures
#
# Dependencies:
#       geNomad --version 1.11.2
#       prodigal --version 2.6.3
#       MMseqs2 Version: bd01c2229f027d8d8e61947f44d11ef1a7669212 
#       diamond --version 2.0.6

##### 待改进部分：
# 1. 是否需要选择最佳匹配？但是原文中设置了diamond --max-target-seqs 10,所以不需要选择最佳匹配? (目前不选择最佳匹配)
#    事实上,如果选择最佳匹配,不必自定义一个函数,只需 --max-target-seqs 1 即可
#    建议参考代码 https://github.com/RChGO/VMGC/blob/main/Pipelines/virusTaxAnno.sh (此脚本未参考)
# 2. 真的需要蛋白质去冗余吗? (目前去冗余)
# 3. genomad使用 --lenient-taxonomy 分类到科级以下的分类单元 (已使用)
# 4. 未按 lineage 分离各等级, 未优化各等级的 Unclassified 注释

import os
import subprocess
import pandas as pd
import logging
from pathlib import Path
from Bio import SeqIO
import shutil

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


###################################### step3 蛋白去冗余
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

###################################### step4 非冗余蛋白比对到 Virus-Host DB
def run_diamond(protein_fasta, db, diamond_results):
    """
    protein_fasta: 输入的蛋白序列
    db: 数据库前缀, 不带 .dmnd
    diamond_results: 结果
    """
    # 创建输出目录
    os.makedirs(os.path.dirname(diamond_results), exist_ok=True)

    # 构建 .dmnd
    dmnd_file = Path(f"{db}.dmnd")
    if not dmnd_file.exists():
        makedb_cmd = ["diamond", "makedb",
                    "--in", str(db),
                    "--db", str(db)]
        subprocess.run(makedb_cmd, check=True)

    cmd = ["diamond", "blastp",
            "-q", str(protein_fasta),
            "-d", str(db),
            "-o", str(diamond_results),
            "--id", "30",
            "--query-cover", "50",
            "--min-score", "50",
            "--max-target-seqs", "10",
            "--outfmt", "6", "qseqid", "stitle", "pident", "length", "evalue", "bitscore"]

    subprocess.run(cmd, check=True)

###################################### step5 基于diamond结果分配科分类和完整分类
# 选择最佳匹配 (实际不使用)
def select_best_matches(input_file_path, output_path):
    """
    选出相同查询序列的最佳匹配, 即 bitscore(最后一列的值)最高的那一行
    参数:
        input_file_path: 输入TSV文件路径
        output_path: 输出TSV文件路径
    返回:
        bool: 处理是否成功
    """
    try:
        # 读取文件内容并按查询序列分组
        query_groups = {}
        with open(input_file_path, 'r', encoding='utf-8') as f_in:
            for line in f_in:
                columns = line.strip().split('\t')
                if not columns or len(columns) < 6:
                    continue
                
                query_id = columns[0]
                bitscore = float(columns[-1])  # 最后一列是bitscore
                
                # 如果是新的查询序列或当前行的bitscore更高,则更新
                if query_id not in query_groups or bitscore > query_groups[query_id]['bitscore']:
                    query_groups[query_id] = {
                        'bitscore': bitscore,
                        'line': line
                    }
        
        # 写入最佳匹配结果
        with open(output_path, 'w', encoding='utf-8') as f_out:
            for query_id in query_groups:
                f_out.write(query_groups[query_id]['line'])
        
        return True
    except Exception as e:
        print(f"处理文件时发生错误: {e}")
        return False

# 格式化 diamond 结果
def format_diamond_results(diamond_results_path, output_path):
    """
    验证diamond_results.tsv中每行第二列按"|"分隔后的元素个数是否一致,如果一致,则将第二列修改为倒数第二个元素
    去掉diamond_results.tsv第一列中的最后一个_及之后的部分,得到蛋白对应的病毒ID
    输出新的TSV文件

    
    参数:
        diamond_results_path: diamond_results.tsv文件路径
        output_path: 输出文件路径,如果为None则不输出
    返回:
        bool: 如果所有行的第二列元素个数一致,则返回True,否则返回False
    """
    # 读取文件内容
    with open(diamond_results_path, 'r') as file:
        lines = file.readlines()

    first_line = True       # 标记是否为文件的第一行
    expected_count = None   # 一行第二列分割后的元素数量
    is_valid = True         # 标记是否所有行都有效
    modified_lines = []     # 存储修改后的行

    for line_num, line in enumerate(lines, 1):
    # - enumerate(lines, 1) : 内置函数,用于将可迭代对象 lines 转换为索引-元素对
    # - 第一个参数 lines ：要遍历的列表(这里是文件的所有行)
    # - 第二个参数 1 ：指定索引的起始值为1(默认是0)
    # - line_num : 变量名,存储当前元素的索引(行号)
    # - line : 变量名,存储当前遍历到的 lines 中的元素(文件行内容)
        
        columns = line.strip().split('\t')  # 对每一行,先按制表符 \t 分割成多列
        
        # 处理第一列,去掉最后一个_及之后的部分
        if columns:
            first_column = columns[0]
            last_underscore_index = first_column.rfind('_')
            if last_underscore_index != -1:
                columns[0] = first_column[:last_underscore_index]
        
        # 处理第二列,按"|"分隔
        second_column = columns[1]
        elements = second_column.split('|')
        
        # 先计算第一行第二列按“|”分隔后的元素个数
        if first_line:
            expected_count = len(elements)
            first_line = False
            # 对于第一行,修改第二列为倒数第二个元素
            if expected_count > 1:
                columns[1] = elements[-2]  # 取倒数第二个元素
            modified_line = '\t'.join(columns) + '\n'
            modified_lines.append(modified_line)
        # 后续行,检查元素个数是否与第一行相同
        else:   
            if len(elements) != expected_count:
                print(f"第{line_num}行: 第二列按'|'分隔后的元素个数与第一行不同")
                is_valid = False
                modified_lines.append(line)  # 添加原行
            else:
                # 元素个数一致,修改第二列为倒数第二个元素
                if len(elements) > 1:
                    columns[1] = elements[-2]  # 取倒数第二个元素
                modified_line = '\t'.join(columns) + '\n'
                modified_lines.append(modified_line)
    
    # 输出新文件(如果所有行都有效)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    if is_valid:
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.writelines(modified_lines)
        except Exception as e:
            print(f"写入文件时发生错误: {e}")
            is_valid = False
    
    return is_valid

# 合并 diamond 结果和 refseq 分类信息
def merge_diamond_refseq(refseq_family_path, diamond_results_processed_path, output_path = 'merge_diamond_refseq.tsv'):
    """
    根据预处理后的diamond结果表的第二列匹配refseq id 对应的分类信息表的对应行,(匹配不到则报错)
    将refseq id 对应的分类信息表的第2、3列添加到预处理后的diamond结果表中。

    参数:
        refseq_family_path: refseq id 对应的分类信息表
        diamond_results_processed_path: 预处理后的 diamond 结果表
        output_path: 输出文件的路径
    """
    # 读取文件refseq_family_path,创建ID到family和分类信息的映射
    refseq_family_map = {}
    try:
        with open(refseq_family_path, 'r', encoding='utf-8') as f2:
            for line in f2:
                columns = line.strip().split('\t')
                # 文件2的第一列作为键,第二列和第三列作为值
                refseq_family_map[columns[0]] = (columns[1], columns[2])
    except Exception as e:
        print(f"读取文件{refseq_family_path}时发生错误: {e}")
        return False

    # 处理文件diamond_results_processed_path并写入结果
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    try:
        with open(diamond_results_processed_path, 'r', encoding='utf-8') as f1, \
                open(output_path, 'w', encoding='utf-8') as out:
            for line in f1:
                columns = line.strip().split('\t')

                # 获取文件1第二列的值用于匹配
                id_to_match = columns[1]
                # 查找匹配项
                if id_to_match in refseq_family_map:
                    family, lineage = refseq_family_map[id_to_match]
                    # 添加文件2的第2、3列
                    columns.extend([family, lineage])
                else:
                    # 如果没有匹配项,抛出错误
                    raise ValueError(f"在{refseq_family_path}中未找到与{diamond_results_processed_path}中ID '{id_to_match}' 匹配的refseq id")

                # 写入结果
                out.write('\t'.join(columns) + '\n')
    except Exception as e:
        print(f"处理文件{diamond_results_processed_path}或写入输出文件时发生错误: {e}")
        return False

    return True

# 分配科分类和完整分类
def assign_family_classification(merged_file_path, output_path):
    """
    对于同一个查询序列的所有匹配,如果有至少20%匹配的科水平的分类信息相同,
    那么该查询序列的科级分类就是这个匹配的科级分类。
    对于未满足20%频率条件的查询序列,科级分类标记为Unclassified,
    但完整分类等级设置为频率最高的科首次匹配的完整分类。
    输出一个新的表格,第一列是查询序列,第二列是科,第三列是科对应的所有分类等级(原表格的最后一列)

    参数:
        merged_file_path: 输入TSV文件路径
        output_path: 输出TSV文件路径
    返回:
        bool: 处理是否成功
    """
    try:
        # 读取文件内容并按查询序列分组
        query_groups = {}
        with open(merged_file_path, 'r', encoding='utf-8') as f_in:
            for line in f_in:
                columns = line.strip().split('\t')
                if not columns or len(columns) < 8:
                    continue
                
                query_id = columns[0]
                family = columns[6]
                lineage = columns[7]
                
                if query_id not in query_groups:
                    query_groups[query_id] = []
                
                query_groups[query_id].append({
                    'family': family,
                    'lineage': lineage
                })
        
        # 处理每个查询序列,确定科级分类
        results = []
        for query_id, matches in query_groups.items():
            total_matches = len(matches)
            family_count = {}
            
            # 统计每个科出现的次数
            for match in matches:
                family = match['family']
                if family not in family_count:
                    family_count[family] = 0
                family_count[family] += 1
            
            # 找到频率最高的科及其首次匹配的完整分类
            max_count = 0
            most_frequent_family = 'Unclassified'
            most_frequent_lineage = 'Viruses; unclassified viruses'
            
            for family, count in family_count.items():
                if count > max_count:
                    max_count = count
                    most_frequent_family = family
                    # 获取该科首次出现的完整分类
                    for match in matches:
                        if match['family'] == family:
                            most_frequent_lineage = match['lineage']
                            break
            
            # 检查是否有科出现频率超过20%
            assigned_family = 'Unclassified'
            assigned_lineage = most_frequent_lineage  # 默认使用频率最高的科的分类等级
            
            for family, count in family_count.items():
                frequency = count / total_matches
                if frequency >= 0.2:
                    assigned_family = family
                    # 获取该科对应的分类等级
                    for match in matches:
                        if match['family'] == family:
                            assigned_lineage = match['lineage']
                            break
                    break
            
            results.append({
                'query_id': query_id,
                'family': assigned_family,
                'lineage': assigned_lineage
            })
        
        # 写入结果
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f_out:
            # 写入表头
            f_out.write('QueryID\tfamily\tlineage\n')
            # 写入数据行
            for result in results:
                f_out.write(f"{result['query_id']}\t{result['family']}\t{result['lineage']}\n")
        
        return True
    except Exception as e:
        print(f"处理文件时发生错误: {e}")
        return False

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
    db_dir = Path("/home/shijiabin/db")
    genomad_db_path = db_dir / "genomad_db" / "genomad_db_v1.9"
    virushostdb_dir = db_dir / "virushostdb_release229"
    virushostdb_path = virushostdb_dir / "virushostdb.cds.faa"
    refseq_family_path = virushostdb_dir / "refseq2family.tsv"

    threads = 10

    #### step1 genomad 分类注释
    # 输入
    votus
    # 输出
    outdir_genomad = Path("01_genomad")
    genomad_tax = outdir_genomad / "votus_annotate/votus_taxonomy.tsv"

    run_genomad(votus, outdir_genomad, genomad_db_path, threads)

    #### step2 预测 ORF
    # 输入
    votus
    # 输出
    protein = Path("02_prodigal/protein.faa")

    run_prodigal(votus, protein)
    
    #### step3 蛋白去冗余
    # 输入
    protein
    # 输出
    dir_mmseqs = Path("03_mmseqs2")
    nonredundant_protein = dir_mmseqs / "nonredundant_protein.faa"

    run_mmseqs_linclust_protein(protein, dir_mmseqs, "protein", "tmp_dir")

    #### step4 非冗余蛋白比对到 Virus-Host DB

    # 非冗余蛋白比对到Virus-Host DB
    # input
    nonredundant_protein
    # output
    diamond_results = Path("04_diamond/diamond_results.tsv")

    run_diamond(nonredundant_protein, virushostdb_path, diamond_results)

    #### step5 基于diamond结果分配科分类和完整分类
    # input
    diamond_results
    # output
    diamond_format = Path("04_diamond/diamond_results_format.tsv")
    merged_file_path = Path("05_diamond_classification/merge.tsv")
    diamond_classification = Path("05_diamond_classification/classification.tsv")

    # 格式化 diamond 结果
    format_diamond_results(diamond_results, diamond_format)
    # 合并 diamond 结果和 refseq 分类信息
    merge_diamond_refseq(refseq_family_path, diamond_format, merged_file_path)
    # 分配科分类和完整分类
    assign_family_classification(merged_file_path, diamond_classification)

    #### step6 合并基于diamond的分类信息和基于genomad的分类信息
    # input:
    diamond_classification
    genomad_tax
    # output:
    merged_path = Path("06_merged_classification/merge_diamond_genomad.tsv")

    # 读取两个TSV文件（带表头）
    diamond_classification = pd.read_csv(diamond_classification, sep='\t', header=0)
    genomad_tax = pd.read_csv(genomad_tax, sep='\t', header=0)

    # 删除 genomad_tax 的 n_genes_with_taxonomy 列，agreement 列，taxid 列
    genomad_tax = genomad_tax.drop(columns=['n_genes_with_taxonomy', 'agreement', 'taxid'])

    # 重命名第一列的列名为 vOTU, 以便直接合并
    diamond_classification = diamond_classification.rename(columns={diamond_classification.columns[0]: 'vOTU'})
    genomad_tax = genomad_tax.rename(columns={genomad_tax.columns[0]: 'vOTU'})

    # 使用 outer 连接合并两个表格
    merged_classification = pd.merge(diamond_classification, genomad_tax, on='vOTU', how='outer', suffixes=('_diamond', '_genomad'))

    # 保存合并结果（保留表头）
    os.makedirs(os.path.dirname(merged_path), exist_ok=True)
    merged_classification.to_csv(merged_path, sep='\t', index=False)

if __name__ == "__main__":
    main()
