# 手动版的virify
# 使用v3版数据库和v4版分类信息(2022版)
# 处理病毒识别方法A和E的结果, all_uvigs

# 参考文献 Massive expansion of human gut bacteriophage diversity

# 运行: snakemake -s pipeline_uvigs.smk --cores 5

import os
import argparse
import csv
from pathlib import Path
import pandas as pd

##################################### 定义全局变量
# 样本名 
SAMPLES = [
    "F1_1A", "F1_2A", "F1_3A",
    "F2_1A", "F2_2A", "F2_3A",
    "FG_1A", "FG_2A", "FG_3A",
    "L1_1A", "L1_2A", "L1_3A",
    "L2_1A", "L2_2A", "L2_3A",
    "H_LX", "H_O"
]

# 类型（A和E）
TYPES = ["A", "E"]

# uvigs目录,该目录下有两个目录A和E,分别是病毒识别方法A和E的结果
virus_dir = "/home/shijiabin/2025_2ME/05_taxonomic_classification/UVIGs"

# hmm数据库
hmmdb = "/home/shijiabin/db/virify_db/vpHMM_database_v3/vpHMM_database_v3.hmm"

# 数据库分类信息,版本 additional_data_vpHMMs_v{1,2,3,4}.tsv
additional_data = "/home/shijiabin/db/virify_db/models/additional_data_vpHMMs_v4.tsv"

# hmmscan中,同一蛋白的多个显著匹配,计算分类注释的比例时,是否需要去重?
deduplicate = False  # True

##################################### 定义流程规则
# 最终要生成的文件
rule all:
    input:
        expand("uvigs_{type}/{sample}_taxonomy.tsv", type=TYPES, sample=SAMPLES)

# 预测ORF
rule prodigal:
    input:
        virus = lambda wildcards: os.path.join(virus_dir, f"{wildcards.type}", f"{wildcards.sample}_uvigs.fna")
    output:
        protein = "uvigs_{type}/01_prodigal/{sample}_protein.faa",
        count_orf = "uvigs_{type}/01_prodigal/{sample}_count_orf.txt"
    log:
        "uvigs_{type}/log/prodigal/{sample}.log"
    shell:
        """
        prodigal -i {input.virus} -a {output.protein} -p meta > {log} 2>&1
        
        # 统计orf数量
        grep ">" {output.protein} | cut -d" " -f1 | sed 's/>//' | sed 's/_[0-9]*$//' | sort | uniq -c | awk '{{print $2 "\\t" $1}}' > {output.count_orf}
        """

# hmmscan
rule hmmscan:
    input:
        protein = "uvigs_{type}/01_prodigal/{sample}_protein.faa"
    output:
        tblout = "uvigs_{type}/02_hmmscan/{sample}.tbl"
    log:
        "uvigs_{type}/log/hmmscan/{sample}.log"
    shell:
        "hmmscan --cpu 2 -E 0.001 --domE 0.1 --tblout {output.tblout} {hmmdb} {input.protein} > {log} 2>&1"

# 格式化hmmscan输出
rule format:
    input:
        count_orf = "uvigs_{type}/01_prodigal/{sample}_count_orf.txt",
        tblout = "uvigs_{type}/02_hmmscan/{sample}.tbl"
    output:
        tblout_format = "uvigs_{type}/02_hmmscan/{sample}_format.tbl"
    run:
        # 取表格前四列,再添加 hit_count, orf_count, ratio, tax 列
        # 填充 target accession, query accession
        # 按 query accession 分组, 对每一组进行如下处理:
        # 1. 计算该组 query name 去重后的数量, 填充到 hit_count
        # 2. 根据 ORF 计数文件, 填充 orf_count
        # 3. 计算 ratio = hit_count/orf_count
        # 4. 根据 target accession 和数据库的分类信息表, 填充 tax

        # 定义表头
        _table_headers = [
            "target name",
            "target accession",
            "query name",
            "query accession",
            "hit_count",
            "orf_count",
            "ratio",
            "tax"
        ]

        input_table = input.tblout
        count_orf = input.count_orf
        output_table = output.tblout_format

        # 读取ORF计数文件，建立query accession到ORF数量的映射
        orf_counts = {}
        with open(count_orf, "r") as count_file:
            for line in count_file:
                parts = line.strip().split()
                if len(parts) == 2:
                    query_accession, orf_count = parts
                    orf_counts[query_accession] = int(orf_count)

        # 读取数据库的分类信息表
        taxa_df = pd.read_csv(additional_data, sep='\t')
        # 创建映射字典，Number列对应target accession，Associated列对应分类单元
        taxa_mapping = dict(zip(taxa_df['Number'], taxa_df['Associated']))

        # 读取输入表格并处理数据
        processed_queries = set()  # 存储已处理的query name,用于去重(如果需要的话)
        all_data = []
        with open(input_table, mode="r") as dt_reader:
            for line in dt_reader:
                if line.startswith("#"):  # 跳过以 # 开头的行
                    continue
                
                cols = line.split()
            
                data = cols[:4]  # 取前四列
                
                # 填充 target accession 列（第2列）
                target_name = data[0]
                target_accession = int(target_name[6:-4])  # 去掉 "ViPhOG" 前缀和 ".faa" 后缀
                data[1] = target_accession

                # 填充 query accession 列（第4列）
                query_name = data[2]
                query_accession = query_name.split("_")[0]  # 取下划线之前的部分
                data[3] = query_accession

                # 根据deduplicate布尔值决定是否去重
                if deduplicate:
                    # 只保留query name列中相同值的第一行
                    if query_name not in processed_queries:
                        processed_queries.add(query_name)
                        all_data.append(data)   
                else:
                    all_data.append(data)   # 保留所有行
                
        # 转换为DataFrame以便处理
        df = pd.DataFrame(all_data, columns=_table_headers[:4])

        # 按query accession分组处理
        result_data = []
        for query_accession, group in df.groupby('query accession'):
            # 计算hit_count
            hit_count = len(group['query name'].unique())
            
            # 获取orf_count
            if query_accession in orf_counts:
                orf_count = orf_counts[query_accession]
            else:
                raise ValueError(f"找不到'{query_accession}'对应的ORF数量")
            
            # 计算ratio
            ratio = hit_count / orf_count
            
            # 对于该组中的每一行数据
            for _, row in group.iterrows():
                # 根据target accession获取分类信息
                target_acc = row['target accession']
                tax = taxa_mapping.get(target_acc, 'unknown') if isinstance(target_acc, int) else 'unknown'
                    
                # 构建完整数据行
                full_row = [
                    row['target name'],
                    row['target accession'],
                    row['query name'],
                    row['query accession'],
                    hit_count,
                    orf_count,
                    ratio,
                    tax
                ]
                result_data.append(full_row)

        # 写入结果到输出文件
        with open(output_table, "w", newline="") as out_table:
            tsv_writer = csv.writer(out_table, delimiter="\t",
                                    quoting=csv.QUOTE_MINIMAL)
            tsv_writer.writerow(_table_headers)
            tsv_writer.writerows(result_data)

# 分配
rule assign:
    input:
        tblout_format = "uvigs_{type}/02_hmmscan/{sample}_format.tbl"
    output:
        assign_tax = "uvigs_{type}/03_assign/{sample}.tsv"
    run:
        # 按 query accession 分组, 保留 ratio >= 0.2 的组
        # 对每一组计算是否有相同的分类单元占比达到0.6,分配该组的分类

        input_table = input.tblout_format
        output_table = output.assign_tax

        df = pd.read_csv(input_table, sep='\t')

        # 筛选比值大于等于0.2的行
        df = df[df['ratio'] >= 0.2]

        # 按query accession分组，并计算每个分类单元的占比
        result = []
        for query_accession, group in df.groupby('query accession'):
            taxa_counts = group['tax'].value_counts()
            total_count = len(group)
            
            # 计算每个分类单元的占比
            taxa_ratios = taxa_counts / total_count
            
            # 查找占比超过60%的分类单元
            dominant_taxa = taxa_ratios[taxa_ratios >= 0.6]
            
            if len(dominant_taxa) > 0:
                # 如果有占比超过60%的分类单元，则使用占比最高的那个
                assigned_taxa = dominant_taxa.index[0]
            else:
                # 否则标记为unknown
                assigned_taxa = 'unknown'
            
            result.append({
                'query accession': query_accession,
                'assigned taxa': assigned_taxa
            })

        # 转换为DataFrame并保存结果，使用制表符作为分隔符
        result_df = pd.DataFrame(result)
        result_df.to_csv(output_table, sep='\t', index=False)

# 提取结果文件
rule postprocess:
    input:
        assign_tax = "uvigs_{type}/03_assign/{sample}.tsv"
    output:
        final_tax = "uvigs_{type}/{sample}_taxonomy.tsv"
    shell:
        "cp {input.assign_tax} {output.final_tax}"