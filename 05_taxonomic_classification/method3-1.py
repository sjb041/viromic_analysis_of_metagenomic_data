#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# Date: 2026-01-15
#
# Description:
#       Phage taxonomy annotation pipeline:
#       votus → Prodigal → DIAMOND (Virus-Host DB, release229_2025-06-03) → taxonomy assignmen
#         ↓
#       genomad
#       1. 
#       2. 
#   参考文献: 2025-The Chinese gut virus catalogue reveals gut virome diversity and disease-related viral signatures
#
# Dependencies:
#       prodigal --version 2.6.3
#       vcontact2 --version 0.11.3

##### 待改进部分：
# 1. 是否需要选择最佳匹配？但是原文中设置了diamond --max-target-seqs 10,所以不需要选择最佳匹配?
#    事实上,如果选择最佳匹配,不必自定义一个函数,只需 --max-target-seqs 1 即可
#    建议参考代码 https://github.com/RChGO/VMGC/blob/main/Pipelines/virusTaxAnno.sh
# 2. 真的需要蛋白质去冗余吗?
# 3. genomad使用 --lenient-taxonomy 分类到科级以下的分类单元

import os
import subprocess
import pandas as pd
import logging
from pathlib import Path
from Bio import SeqIO

###################################### step0 配置数据库
def setup_virushostdb(virushostdb_dir, virushostdb_path, release = "release229", release_date = "2025-06-03"):
    """
    Download and configure Virus-Host DB (KEGG) and build DIAMOND database.

    Parameters
    ----------
    virushostdb_dir : str
        Directory to store Virus-Host DB
    virushostdb_path : str
        Path to DIAMOND database prefix
    release : str
        Virus-Host DB release version
    release_date : str
        Release date for logging
    """

    flag_file = os.path.join(virushostdb_dir, "virushostdb_configured")

    # 1. Skip if already configured
    if os.path.isfile(flag_file):
        print(
            f"Virus-Host DB 数据库已配置, 位于 {virushostdb_path}, 跳过下载和处理步骤\n"
        )
        return

    print(
        f"正在配置 Virus-Host DB 数据库, 当前使用的版本是 {release} ({release_date})...\n"
    )

    os.makedirs(virushostdb_dir, exist_ok=True)

    def run(cmd):
        print(f"[CMD] {cmd}")
        subprocess.run(cmd, shell=True, check=True)

    os.chdir(virushostdb_dir)

    # 2. Download database files
    base_url = f"https://www.genome.jp/ftp/db/virushostdb/old/{release}"

    run(f'wget "{base_url}/virushostdb.tsv"')
    run(f'wget "{base_url}/virus_genome_type.tsv"')
    run(f'wget "{base_url}/virushostdb.cds.faa.gz"')

    # 3. Build DIAMOND database
    run("gunzip -f virushostdb.cds.faa.gz")
    run(f"diamond makedb --in virushostdb.cds.faa --db {virushostdb_path}")

    # 4. Process taxonomy / host mapping tables
    run(
        "cut -f 1,2,3,4 virushostdb.tsv | sort -r | uniq > virushostdb_nonredundant.tsv"
    )

    run(
        "awk -F'\t' -v OFS='\t' "
        "'{print $1, $4, $18, $15, $5, $7, $8, $10, $11, $12}' "
        "virus_genome_type.tsv > tmp"
    )

    run("head -n1 virushostdb_nonredundant.tsv > header1")
    run("tail -n +2 virushostdb_nonredundant.tsv | sort -k1,1 > sorted1")
    run("head -n1 tmp > header2")
    run("tail -n +2 tmp | sort -k1,1 > sorted2")

    run("cut -f2- header2 | paste header1 - > final_header")
    run("join -t $'\\t' -1 1 -2 1 -a 1 sorted1 sorted2 > joined_data")
    run("cat final_header joined_data > final_tax.tsv")

    run(
        "awk 'BEGIN {FS=OFS=\"\\t\"} NR>1 {"
        "split($4, refseqs, \",\"); "
        "for(i in refseqs) {"
        "gsub(/^[ \\t]+|[ \\t]+$/, \"\", refseqs[i]); "
        "print refseqs[i], $11, $3"
        "}}' final_tax.tsv > refseq_family.tsv"
    )

    run(
        "awk 'BEGIN {FS=OFS=\"\\t\"} {if($2 == \"\") $2=\"Unclassified\"} 1' "
        "refseq_family.tsv | sort | uniq > tmp && mv tmp refseq_family.tsv"
    )

    run(
        "rm -f header1 sorted1 header2 sorted2 final_header joined_data"
    )

    # 5. Mark as configured
    open(flag_file, "w").close()

    print("\nVirus-Host DB 配置完毕\n")

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
# 选择最佳匹配 ()
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
    votus = Path("/home/shijiabin/2025_2ME/03_virus_identification/methodA/votus.fna")

    # 数据库
    db_dir = Path("/home/shijiabin/db")
    genomad_db_path = db_dir / "genomad_db"
    virushostdb_dir = db_dir / "virushostdb_release229_2025-06-03"
    virushostdb_path = virushostdb_dir / "virushostdb"
    refseq_family_path = virushostdb_dir / "refseq_family.tsv"

    #### step1 genomad 分类注释
    # 输入
    votus
    # 输出
    outdir_genomad = Path("01_genomad")
    genomad_tax = outdir_genomad / "votus_annotate/votus_taxonomy.tsv"

    #run_genomad(votus, outdir_genomad, genomad_db_path, threads)

    #### step2 预测 ORF
    # 输入
    votus
    # 输出
    protein = Path("02_prodigal/protein.faa")

    #run_prodigal(votus, protein)
    
    #### step3 蛋白去冗余
    # 输入
    protein
    # 输出
    dir_mmseqs = Path("03_mmseqs2")
    prefix = dir_mmseqs / "protein" 
    rep_seq = dir_mmseqs / f"{prefix.name}_rep_seq.fasta"
    nonredundant_protein = dir_mmseqs / "nonredundant_protein.faa"

    os.makedirs(dir_mmseqs, exist_ok=True)
    
    cmd1 = ["mmseqs", "easy-linclust", 
           str(protein), 
           str(prefix), 
           "tmp_mmseqs", 
           "--min-seq-id", "0.9", 
           "--cov-mode", "1", 
           "-c", "0.8", 
           "--kmer-per-seq", "80"]
    subprocess.run(cmd1, check=True)

    cmd2 = ["mv", str(rep_seq), str(nonredundant_protein)]
    subprocess.run(cmd2, check=True)
    cmd3 = ["rm", "-r", "tmp_mmseqs"]
    subprocess.run(cmd3, check=True)

if __name__ == "__main__":
    main()
