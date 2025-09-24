# Prodigal + vConTACT2 + 内置的RefSeq病毒蛋白参考数据库
# 处理病毒识别方法A和E的结果, votus

# 参考文献 2024-The gut ileal mucosal virome is disturbed in patients with Crohn’sdiseaseand exacerbates intestinal inflammation in mice

# 运行: snakemake -s run.smk --cores 1 --use-conda

import os
import pandas as pd
import logging
from pathlib import Path
from contextlib import redirect_stdout

##################################### 定义函数
# 设置日志记录器
def setup_logger(log_file):
    """设置日志记录器"""
    # 清除现有的handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    # 配置新的日志
    logging.basicConfig(
        filename=str(log_file),
        level=logging.INFO,
        format='%(message)s',
        force=True  # 覆盖之前的配置
    )
    return logging.getLogger()

# 根据vcontact2的输入g2g文件和输出表格,分配输入序列的分类等级
def assign_taxonomy(overview_file, g2g_file, overview_filtered, final_tax, log_file):
    """
    基于vConTACT2聚类结果为目标序列分配病毒分类信息
    
    该函数通过分析vConTACT2生成的基因组概览文件,为目标序列(来自FASTA文件)分配
    病毒分类信息。主要基于preVC科级别和VC属级别聚类结果,通过统计
    参考序列的分类信息频数来推断目标序列的分类。
    
    参数:
    overview_file (str): vConTACT2生成的genome_by_genome_overview.csv文件路径
    g2g_file (str): vConTACT2的输入g2g文件,用于获取目标序列的序列名
    output_overview_filtered (str): 输出筛选后的完整结果文件路径TSV格式
    output_tax (str): 输出目标序列分类结果文件路径TSV格式,仅前7列
    log_file (str): 日志文件路径,保存警告信息
    
    处理流程:
    1. 从FASTA文件提取目标序列ID列表
    2. 筛选overview文件,只保留同时包含目标序列和参考序列的preVC组
    3. 为每个preVC组的目标序列分配分类信息:
    - Family及以上级别:基于该组参考序列中出现频数最高的分类
    - Genus级别:基于VC小组参考序列中出现频数最高的分类
    """

    # 设置日志记录器
    logger = setup_logger(log_file)
    
    # 提取FASTA序列名
    seqid_list = pd.read_csv(g2g_file)['contig_id'].unique().tolist()
    # 读取overview文件
    overview = pd.read_csv(overview_file)

    ###### 筛选overview,只保留那些既有目标序列又有参考序列的preVC组
    # 按照Excel的筛选步骤
    # 1. 筛选那些目标序列preVC列不为空的行
    df1 = overview[(overview['Genome'].isin(seqid_list)) & 
                (overview['preVC'].notna()) & 
                (overview['preVC'] != '')]

    # 2. 使用1的preVC列筛选overview
    preVC_list = list(df1['preVC'].unique())     # preVC列去重,转化为列表
    df2 = overview[overview['preVC'].isin(preVC_list)]

    # 3. 只保留那些既有目标序列又有参考序列的preVC组
    # 找出满足条件的preVC组
    preVC_with_both = []
    for preVC in preVC_list:
        group = df2[df2['preVC'] == preVC]
        has_target = (group['Genome'].isin(seqid_list)).any()  # 组内是否有目标序列
        has_reference = (~group['Genome'].isin(seqid_list)).any()  # 组内是否有参考序列  
        # 只有当组内同时包含目标序列和参考序列时才保留该组
        if has_target and has_reference:
            preVC_with_both.append(preVC)
    # 筛选preVC组
    filtered_df = df2[df2['preVC'].isin(preVC_with_both)]

    ###### 为每个目标序列分配Family和Genus信息
    # Family分配：
    # 对每个preVC组，统计参考序列中各Family的出现频数
    # 排除'Unassigned'值
    # 选择出现频数最高的Family作为该组目标序列的Family
    # 如果有多个Family频数相同则报错
    # Genus分配：
    # 在preVC组内，根据VC进一步分组
    # 排除'Unassigned'值
    # 选择出现频数最高的Genus作为该VC小组目标序列的Genus
    # 如果有多个Genus频数相同则报错
    df = filtered_df.copy() 
    # 按preVC组处理
    for prevc in df['preVC'].unique():
        group_mask = (df['preVC'] == prevc)
        group_data = df[group_mask]
        
        # 获取参考序列用于Family分类
        reference_mask = ~group_data['Genome'].isin(seqid_list)
        reference_family_counts = group_data[reference_mask]['Family'].value_counts()
        
        # 移除'Unassigned'值
        valid_families = reference_family_counts[reference_family_counts.index != 'Unassigned']
        
        if len(valid_families) > 0:
            # 检查是否有多个最高频数
            max_count = valid_families.max()
            top_families = valid_families[valid_families == max_count]
            
            if len(top_families) > 1:
                warning_msg = f"警告: preVC {prevc} 中存在多个相同频数的Family: {list(top_families.index)}，跳过该组"
                print(warning_msg)
                logger.warning(warning_msg)
                continue
            
            predicted_family = top_families.index[0]
            
            # 为该组所有目标序列填充Family及更高级分类
            target_mask = group_mask & (df['Genome'].isin(seqid_list))
            
            # 填充Family及更高级分类信息
            reference_taxonomy = group_data[reference_mask].iloc[0]  # 取第一个参考序列的分类信息作为模板
            
            df.loc[target_mask, 'Family'] = predicted_family
            df.loc[target_mask, 'Order'] = reference_taxonomy['Order']
            df.loc[target_mask, 'Class'] = reference_taxonomy['Class']
            df.loc[target_mask, 'Phylum'] = reference_taxonomy['Phylum']
            df.loc[target_mask, 'Kingdom'] = reference_taxonomy['Kingdom']
            
            # 按VC小组处理Genus
            for vc in group_data['VC'].unique():
                if pd.notna(vc):
                    vc_mask = group_mask & (df['VC'] == vc)
                    vc_group_data = df[vc_mask]
                    
                    # 获取参考序列用于Genus分类
                    vc_reference_mask = ~vc_group_data['Genome'].isin(seqid_list)
                    reference_genus_counts = vc_group_data[vc_reference_mask]['Genus'].value_counts()
                    
                    # 移除'Unassigned'值
                    valid_genera = reference_genus_counts[reference_genus_counts.index != 'Unassigned']
                    
                    if len(valid_genera) > 0:
                        # 检查是否有多个最高频数
                        max_count = valid_genera.max()
                        top_genera = valid_genera[valid_genera == max_count]

                        if len(top_genera) > 1:
                            warning_msg = f"警告: VC {vc} 中存在多个相同频数的Genus: {list(top_genera.index)}，跳过该VC"
                            print(warning_msg)
                            logger.warning(warning_msg)
                            continue
                        
                        predicted_genus = top_genera.index[0]
                        
                        # 为该VC小组的所有目标序列填充Genus
                        vc_target_mask = vc_mask & (df['Genome'].isin(seqid_list))
                        df.loc[vc_target_mask, 'Genus'] = predicted_genus

    # 只保留目标序列的分类结果，并且只保留前7列
    target_classified = df[df['Genome'].isin(seqid_list)].iloc[:, :7]

    # 保存结果
    filtered_df.to_csv(overview_filtered, sep=',', index=False)
    target_classified.to_csv(final_tax, sep=',', index=False)

        


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

# vcontact2使用的数据库版本
db = "ProkaryoticViralRefSeq211-Merged"

##################################### 定义流程规则
# 最终要生成的文件
rule all:
    input:
        expand("votus_{type}/{sample}_taxonomy.csv", type=TYPES, sample=SAMPLES),
        expand("uvigs_{type}/{sample}_taxonomy.csv", type=TYPES, sample=SAMPLES)

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

# 预测蛋白
rule prodigal_votus:
    input:
        votus = "votus_{type}/01_cluster/{sample}/votus.fa"
    output:
        protein = "votus_{type}/02_prodigal/{sample}_protein.faa"
    log:
        "votus_{type}/log/prodigal/{sample}.log"
    shell:
        "prodigal -i {input.votus} -a {output.protein} -p meta > {log} 2>&1"

# 生成viral clusters (VCs)
rule vcontact2_votus:
    input:
        protein = "votus_{type}/02_prodigal/{sample}_protein.faa"
    output:
        dir = directory("votus_{type}/03_vcontact2/{sample}"),
        g2g = "votus_{type}/03_vcontact2/{sample}/g2g.csv",
        overview = "votus_{type}/03_vcontact2/{sample}/genome_by_genome_overview.csv"
    conda:
        "vc2"
    log:
        "votus_{type}/log/vcontact2/{sample}.log"
    shell:
        """
        echo "================= 开始执行 vcontact2_gene2genome 步骤 =================" > {log}
        vcontact2_gene2genome -s 'Prodigal-FAA' -p {input.protein} -o {output.g2g} >> {log} 2>&1
        echo "================= 开始执行 vcontact2 聚类步骤 =================" >> {log}
        vcontact2  --rel-mode 'Diamond' --db '{db}' --pcs-mode MCL --vcs-mode ClusterONE --raw-proteins {input.protein} --proteins-fp {output.g2g} --output-dir {output.dir} >> {log} 2>&1
        """

# 分配分类信息
rule assign_votus:
    input:
        g2g = "votus_{type}/03_vcontact2/{sample}/g2g.csv",
        overview = "votus_{type}/03_vcontact2/{sample}/genome_by_genome_overview.csv"
    output:
        overview_filtered = "votus_{type}/03_vcontact2/{sample}/overview_filtered.csv",
        final_tax = "votus_{type}/{sample}_taxonomy.csv"
    log:
        "votus_{type}/log/assign/{sample}.log"
    run:
        assign_taxonomy(input.overview, input.g2g, output.overview_filtered, output.final_tax, log)

# ------------------- uvigs 分类 ------------------- #
# 预测蛋白
rule prodigal_uvigs:
    input:
        virus = lambda wildcards: os.path.join(virus_dir, f"{wildcards.type}", f"{wildcards.sample}_uvigs.fna")
    output:
        protein = "uvigs_{type}/02_prodigal/{sample}_protein.faa"
    log:
        "uvigs_{type}/log/prodigal/{sample}.log"
    shell:
        "prodigal -i {input.virus} -a {output.protein} -p meta > {log} 2>&1"

# 生成viral clusters (VCs)
rule vcontact2_uvigs:
    input:
        protein = "uvigs_{type}/02_prodigal/{sample}_protein.faa"
    output:
        dir = directory("uvigs_{type}/03_vcontact2/{sample}"),
        g2g = "uvigs_{type}/03_vcontact2/{sample}/g2g.csv",
        overview = "uvigs_{type}/03_vcontact2/{sample}/genome_by_genome_overview.csv"
    conda:
        "vc2"
    log:
        "uvigs_{type}/log/vcontact2/{sample}.log"
    shell:
        """
        echo "================= 开始执行 vcontact2_gene2genome 步骤 =================" > {log}
        vcontact2_gene2genome -s 'Prodigal-FAA' -p {input.protein} -o {output.g2g} >> {log} 2>&1
        echo "================= 开始执行 vcontact2 聚类步骤 =================" >> {log}
        vcontact2  --rel-mode 'Diamond' --db '{db}' --pcs-mode MCL --vcs-mode ClusterONE --raw-proteins {input.protein} --proteins-fp {output.g2g} --output-dir {output.dir} >> {log} 2>&1
        """

# 分配分类信息
rule assign_uvigs:
    input:
        g2g = "uvigs_{type}/03_vcontact2/{sample}/g2g.csv",
        overview = "uvigs_{type}/03_vcontact2/{sample}/genome_by_genome_overview.csv"
    output:
        overview_filtered = "uvigs_{type}/03_vcontact2/{sample}/overview_filtered.csv",
        final_tax = "uvigs_{type}/{sample}_taxonomy.csv"
    log:
        "uvigs_{type}/log/assign/{sample}.log"
    run:
        assign_taxonomy(input.overview, input.g2g, output.overview_filtered, output.final_tax, log)

# 汇总要填充的分类信息
onsuccess:
    types = TYPES
    samples = SAMPLES

    for i in ["votus", "uvigs"]:
        for t in types:
            log_dir = Path(f"{i}_{t}")
            log_path = log_dir / "output.log"

            # 确保目录存在
            log_dir.mkdir(parents=True, exist_ok=True)

            # 覆盖写新日志
            with open(log_path, "w") as out_f:
                for s in samples:
                    sample_log = log_dir / "log" / "assign" / f"{s}.log"
                    out_f.write(f"---------------- {s} ----------------\n")
                    if sample_log.exists():
                        with open(sample_log) as in_f:
                            out_f.write(in_f.read())
                        out_f.write("\n")
                    else:
                        out_f.write(f"[WARNING] {sample_log} not found\n")

            # 打印日志内容到 stdout
            with open(log_path) as f:
                print(f.read())
