# PhaBOX2-phagcn(类似phagcn3) + ICTV 2024
# 能分类到属级,分类结果现在包括分类信息的完整路径,注意：并非所有病毒都有完整的分类路径，所有分类路径都严格符合 ICTV 的分类标准
# 处理病毒识别方法A和E的结果, all_uvigs and votus

# 2024-The phageome of patients with ulcerative colitis treated with donor fecal microbiota reveals markers associated with disease remission

# 运行: snakemake -s run.smk --cores 5 --use-conda

import os
import pandas as pd

### 定义函数,处理分类结果表
# 提取各分类等级,并处理空值
def process_lineage(df, lineage_col="Lineage"):
    """
    从 Lineage 列解析出各分类等级，并处理空值。
    
    - 提取等级: kingdom, phylum, class, order, family, subfamily, genus
    - 空值处理逻辑: 从 genus 开始依次往左，
        如果该等级为空:
            - 上一级也为空 -> "unknown"
            - 上一级不为空 -> "Unclassified{上一级值}"
    """
    ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus']
    
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

    return df_out

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

# 添加各等级列
rule postprocess_votus:
    input:
        tax = "votus_{type}/PhaBOX2_phagcn/{sample}/final_prediction/phagcn_prediction.tsv"
    output:
        final_tax = "votus_{type}/{sample}_taxonomy.tsv"
    run:
        df = pd.read_csv(input.tax, sep='\t')
        result_df = process_lineage(df)
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

# 添加各等级列
rule postprocess_uvigs:
    input:
        tax = "uvigs_{type}/PhaBOX2_phagcn/{sample}/final_prediction/phagcn_prediction.tsv"
    output:
        final_tax = "uvigs_{type}/{sample}_taxonomy.tsv"
    run:
        df = pd.read_csv(input.tax, sep='\t')
        result_df = process_lineage(df)
        result_df.to_csv(output.final_tax, sep='\t', index=False)
