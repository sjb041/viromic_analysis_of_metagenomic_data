# PhaBOX2-phagcn(类似phagcn3) + ICTV 2024
# 能分类到属级,分类结果现在包括分类信息的完整路径,注意：并非所有病毒都有完整的分类路径，所有分类路径都严格符合 ICTV 的分类标准
# 处理病毒识别方法A和E的结果, all_uvigs and votus

# 2024-The phageome of patients with ulcerative colitis treated with donor fecal microbiota reveals markers associated with disease remission

# 运行: snakemake -s run.smk --cores 5 --use-conda

##### 待改进部分：
# 1. 处理最终分类结果表,提取出科级分类另做一列

import os
import pandas as pd

### 定义函数,处理分类结果表
def extract_family(lineage):
    """
    从lineage字符串中提取family信息
    
    参数:
    lineage (str): 包含分类信息的字符串
    
    返回:
    str: family名称或'unknown'
    """
    # 处理空值和特殊值
    if pd.isna(lineage) or lineage == 'no hits to database':
        return 'unknown'
    
    # 按分号分割lineage字符串
    elements = lineage.split(';')
    
    # 查找以"family:"开头的元素
    for element in elements:
        if element.startswith('family:'):
            # 返回冒号后的部分（family名称）
            return element.split(':')[1]
    
    # 如果没有找到family信息，返回unknown
    return 'unknown'

def add_family_column(df):
    """
    在PhaGCN预测结果数据框中,
    Lineage列的内容按";"分隔,取以"family:"开头的元素,作为该行新的列Family的值,
    如果没有以"family:"开头的元素,则为"unknown",
    Family放在第3列位置

    参数:
    df (pandas.DataFrame): PhaGCN预测结果的数据框
    
    返回:
    pandas.DataFrame: 处理后的数据框
    """
    
    # 创建数据框的副本以避免修改原始数据
    df_copy = df.copy()
    
    # 创建Family列
    df_copy['Family'] = df_copy['Lineage'].apply(extract_family)
    
    # 重新排列列的顺序，将Family列放在第3列（索引为2）
    cols = df_copy.columns.tolist()
    
    # 找到Family列的当前位置
    family_col_index = cols.index('Family')
    
    # 将Family列移到第3列位置（索引为2）
    new_cols = cols[:2] + [cols[family_col_index]] + cols[2:family_col_index] + cols[family_col_index+1:]
    df_copy = df_copy[new_cols]
    
    return df_copy

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

# 数据库
db_dir = "/home/shijiabin/db/phabox2_db/phabox_db_v2"

threads = 16    # 线程数
min_len = 0    # 最小序列长度,0表示不限制

##################################### 定义流程规则
rule all:
    input:
        expand("votus_{type}/{sample}_taxonomy.tsv", type=TYPES, sample=SAMPLES),
        expand("uvigs_{type}/{sample}_taxonomy.tsv", type=TYPES, sample=SAMPLES)

# ------------------- votus 分类 ------------------- # 
# 聚类
rule cluster:
    input:
        virus = lambda wildcards: os.path.join(virus_dir, f"{wildcards.type}", f"{wildcards.sample}_uvigs.fna")
    output:
        dir = directory("votus_{type}/01_cluster/{sample}"),
        votus = "votus_{type}/01_cluster/{sample}/votus.fa"
    log:
        "votus_{type}/log/cluster/{sample}.log"
    shell:
        "cluster.sh -i {input.virus} -t {input.virus} -a 95 -o {output.dir} > {log} 2>&1"

# PhaBOX2-phagcn分类注释
rule PhaBOX2_phagcn_votus:
    input:
        votus = "votus_{type}/01_cluster/{sample}/votus.fa"
    output:
        dir = directory("votus_{type}/PhaBOX2_phagcn/{sample}"),
        tax = "votus_{type}/PhaBOX2_phagcn/{sample}/final_prediction/phagcn_prediction.tsv"
    conda:
        "phabox2"
    log:
        "votus_{type}/log/phabox2/{sample}.log"
    shell:
        # 仅运行分类注释
        "phabox2 --task phagcn --contigs {input.votus} --dbdir {db_dir} --outpth {output.dir} --len {min_len} --threads {threads} > {log} 2>&1"

# 添加Family列
rule postprocess_votus:
    input:
        tax = "votus_{type}/PhaBOX2_phagcn/{sample}/final_prediction/phagcn_prediction.tsv"
    output:
        final_tax = "votus_{type}/{sample}_taxonomy.tsv"
    run:
        df = pd.read_csv(input.tax, sep='\t')
        result_df = add_family_column(df)
        result_df.to_csv(output.final_tax, sep='\t', index=False)

# ------------------- uvigs 分类 ------------------- #
# PhaBOX2-phagcn分类注释
rule PhaBOX2_phagcn_uvigs:
    input:
        virus = lambda wildcards: os.path.join(virus_dir, f"{wildcards.type}", f"{wildcards.sample}_uvigs.fna")
    output:
        dir = directory("uvigs_{type}/PhaBOX2_phagcn/{sample}"),
        tax = "uvigs_{type}/PhaBOX2_phagcn/{sample}/final_prediction/phagcn_prediction.tsv"
    conda:
        "phabox2"
    log:
        "uvigs_{type}/log/phabox2/{sample}.log"
    shell:
        # 仅运行分类注释
        "phabox2 --task phagcn --contigs {input.virus} --dbdir {db_dir} --outpth {output.dir} --len {min_len} --threads {threads} > {log} 2>&1"

# 添加Family列
rule postprocess_uvigs:
    input:
        tax = "uvigs_{type}/PhaBOX2_phagcn/{sample}/final_prediction/phagcn_prediction.tsv"
    output:
        final_tax = "uvigs_{type}/{sample}_taxonomy.tsv"
    run:
        df = pd.read_csv(input.tax, sep='\t')
        result_df = add_family_column(df)
        result_df.to_csv(output.final_tax, sep='\t', index=False)
