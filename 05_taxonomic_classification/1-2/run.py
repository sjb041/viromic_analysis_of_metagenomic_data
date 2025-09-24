# Prodigal + vConTACT2 + geNomad
# 处理病毒识别方法A和E的结果, all_uvigs and votus
# 参考文献 Metagenomic analysis reveals gut phage diversity across three mammalian models

# 运行: python run.py

# 不知道为什么,写的snakemake运行时一直卡在 Building DAG of jobs...,所以放弃,使用python将流程串起来

import os
import pandas as pd
import subprocess
from contextlib import redirect_stdout

# 分配分类信息的方法
def assign_classification(genomad_path, vcontact_path, merge_out, final_out):
    """
    输入：
        genomad_path（votus_taxonomy.tsv路径）
        vcontact_path（genome_by_genome_overview.csv路径）
    输出：
        merge_out（合并后文件路径）
        将votus_taxonomy的分类谱系按等级拆分，然后合并到genome_by_genome_overview中
    
        final_out（填充后文件路径）
        首先按VC编号分组，
        对每个VC组，进行如下处理：
        计算组内成员的Genus中非缺失值出现频数， 频数最高的成员，将它的各分类等级列的值赋值给其他Genus为缺失值的成员
        如果成员的Genus中没有非缺失值，则计算组内成员的Family中非缺失值出现频数....
    """
    # 定义分类等级列（从高到低）
    taxonomy_levels = ['Viruses', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

    # 读取输入文件
    df_genomad = pd.read_csv(genomad_path, sep='\t')
    df_vcontact = pd.read_csv(vcontact_path)

    # 将lineage列按分号进行拆分，分配拆分后的分类等级到各列
    split_result = df_genomad['lineage'].str.split(';', expand=True)
    df_genomad[taxonomy_levels] = split_result

    # 提取出 Genome, preVC, VC Status, VC 列
    df_vcontact = df_vcontact[['Genome', 'preVC', 'VC Status', 'VC']]

    # 将df_genomad的各分类等级添加到df_vcontact中
    df_vcontact = pd.merge(
        df_vcontact,
        df_genomad[['seq_name'] + taxonomy_levels],
        left_on='Genome',
        right_on='seq_name',
        how='left'
    )

    # 将空白串替换为缺失值
    df_vcontact = df_vcontact.replace('', pd.NA)

    # 保存合并结果
    df_vcontact.to_csv(merge_out, sep='\t', index=False)

    # 按VC分组
    groups = df_vcontact.groupby('VC')

    # 对每组成员的分类等级进行缺失值填充
    for vc, group in groups:
        #print(f"VC: {vc}")

        # 获取该组在原数据框中的索引
        group_indices = group.index
        
        # 按照分类等级从低到高进行处理（Genus -> Family -> Order ...）
        # 从Genus开始向上循环（不包括Species）
        genus_index = taxonomy_levels.index('Genus')  # 获取Genus在分类等级列表中的索引
        for level_idx in range(genus_index, -1, -1):  # 从Genus开始向上
            current_level = taxonomy_levels[level_idx]
            target_levels = taxonomy_levels[:level_idx+1]  # 包含当前等级及其以上所有等级
            
            # 获取该组中当前分类等级非缺失值
            non_na_values = group[current_level].dropna()

            if len(non_na_values) == 0:
                #print(f"{current_level} 等级没有分类信息，跳过当前等级的填充，向上一等级寻找...")
                continue
            
            if len(non_na_values) > 0:
                # 计算非缺失值出现频数
                value_counts = non_na_values.value_counts()
                
                # 获取频数最高的值
                most_frequent_value = value_counts.index[0]

                # 找到具有最高频数值的成员（取第一个即可）
                reference_member = group[group[current_level] == most_frequent_value].iloc[0]
                
                # 获取需要填充的成员（当前等级为缺失值的成员）
                members_to_fill = group[group[current_level].isna()]
                
                if len(members_to_fill) > 0:
                    # 打印需要填充的成员信息
                    print(f"{vc}")
                    print(f"{current_level} 等级有分类信息 {most_frequent_value}")
                    print(f"需要填充的成员为：")
                    print(members_to_fill[['Genome']].to_string(index=False, header=False))
                    print(f"\n")

                    # 对需要填充的成员进行填充
                    for idx in members_to_fill.index:
                        for level in target_levels:
                            # 只填充缺失值
                            if pd.isna(df_vcontact.loc[idx, level]):
                                df_vcontact.loc[idx, level] = reference_member[level]      
                #else:
                    #print(f"所有成员当前等级都有分类信息，不需要填充，{vc} 处理完毕")

                # 一旦找到可以处理的等级，就跳出循环
                break


    # 保存填充结果
    df_vcontact.to_csv(final_out, sep='\t', index=False)

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
METHODS = ["A", "E"]

# uvigs目录,该目录下有两个目录A和E,分别是病毒识别方法A和E的结果
virus_dir = "/home/shijiabin/2025_2ME/05_taxonomic_classification/UVIGs"

# 参考数据库的路径
genomad_db = "/home/shijiabin/db/genomad_db"

# 聚类并提取votu
def run_cluster(input_virus, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    cmd = f"cluster.sh -i {input_virus} -t {input_virus} -a 95 -o {output_dir}"
    subprocess.run(cmd, shell=True, check=True)

# genomad分类注释 
def run_genomad(input_virus, output_dir, db):
    os.makedirs(output_dir, exist_ok=True)
    cmd = f"conda run -n genomad genomad annotate {input_virus} {output_dir} {db} --lenient-taxonomy --threads 20"
    subprocess.run(cmd, shell=True, check=True)

# 预测蛋白
def run_prodigal(input_virus, output_protein):
    os.makedirs(os.path.dirname(output_protein), exist_ok=True)
    cmd = f"prodigal -i {input_virus} -a {output_protein} -p meta"
    subprocess.run(cmd, shell=True, check=True)

# 生成viral clusters (VCs)
def run_vcontact2(input_protein, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    g2g_file = f"{output_dir}/g2g.csv"
    cmd1 = f"conda run -n vc2 vcontact2_gene2genome -s 'Prodigal-FAA' -p {input_protein} -o {g2g_file}"
    cmd2 = f"conda run -n vc2 vcontact2 --rel-mode 'Diamond' --db 'None' --pcs-mode MCL --vcs-mode ClusterONE --raw-proteins {input_protein} --proteins-fp {g2g_file} --output-dir {output_dir}"
    subprocess.run(cmd1, shell=True, check=True)
    subprocess.run(cmd2, shell=True, check=True)

# 分类分配
def run_assign_classification(input_genomad_tax, input_overview, output_dir, sample, method, virus_type):
    os.makedirs(output_dir, exist_ok=True)
    merge_taxonomy = os.path.join(output_dir, f"{sample}_merge_taxonomy.tsv")
    final_taxonomy = os.path.join(output_dir, f"{sample}_final_taxonomy.tsv")
    
    log_file = f"{virus_type}_{method}/log/assign_classification/{sample}.log"
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    
    with open(log_file, "w") as f, redirect_stdout(f):
        assign_classification(input_genomad_tax, input_overview, merge_taxonomy, final_taxonomy)

# 对votus执行整个流程
for method in METHODS:
    for sample in SAMPLES:

        ### 聚类并提取votu
        input_virus = os.path.join(virus_dir, f"{method}", f"{sample}_uvigs.fna")
        output_dir_cluster = f"votus_{method}/01_cluster/{sample}"
        run_cluster(input_virus, output_dir_cluster)
        
        ### genomad分类注释
        input_votus = os.path.join(output_dir_cluster, "votus.fa")
        output_dir_genomad = f"votus_{method}/02-1_genomad/{sample}"
        db = genomad_db
        run_genomad(input_votus, output_dir_genomad, db)

        ### 预测蛋白
        output_protein = f"votus_{method}/02_prodigal/{sample}_protein.faa"
        run_prodigal(input_votus, output_protein)

        ### 生成viral clusters (VCs)
        output_dir_vc = f"votus_{method}/03_vcontact2/{sample}"
        run_vcontact2(output_protein, output_dir_vc)

        ### 分类分配
        input_genomad_tax = os.path.join(output_dir_genomad, "votus_annotate/votus_taxonomy.tsv")
        input_overview = os.path.join(output_dir_vc, "genome_by_genome_overview.csv")
        output_dir_assign = f"votus_{method}/04_assign_classification"
        run_assign_classification(input_genomad_tax, input_overview, output_dir_assign, sample, method, "votus")

        ### 提取最终结果
        final_taxonomy = merge_taxonomy = os.path.join(output_dir_assign, f"{sample}_merge_taxonomy.tsv")
        final_tax = f"votus_{method}/{sample}_taxonomy.tsv"
        if os.path.exists(final_taxonomy):
            subprocess.run(f"cp {final_taxonomy} {final_tax}", shell=True, check=True)

    # 合并日志,汇总有哪些是填充的
    log_dir = f"votus_{method}/log/assign_classification"
    merged_log_file = f"votus_{method}/output.log"
    os.makedirs(os.path.dirname(merged_log_file), exist_ok=True)  # 确保输出目录存在
    log_files = [os.path.join(log_dir, f"{sample}.log") for sample in SAMPLES]
    
    with open(merged_log_file, "w") as outfile:
        for i, log_file in enumerate(log_files):
            sample_name = SAMPLES[i]  # 获取对应的样本名
            if os.path.exists(log_file):
                outfile.write(f"=== Sample: {sample_name} ===\n")
                with open(log_file, "r") as infile:
                    outfile.write(infile.read())
                outfile.write("\n" + "="*50 + "\n\n")  # 添加更清晰的分隔符

# 对uvigs执行整个流程
for method in METHODS:
    for sample in SAMPLES:

        ### genomad分类注释
        input_virus = os.path.join(virus_dir, f"{method}", f"{sample}_uvigs.fna")
        output_dir_genomad = f"uvigs_{method}/02-1_genomad/{sample}"
        db = genomad_db
        run_genomad(input_virus, output_dir_genomad, db)

        ### 预测蛋白
        output_protein = f"uvigs_{method}/02_prodigal/{sample}_protein.faa"
        run_prodigal(input_virus, output_protein)

        ### 生成viral clusters (VCs)
        output_dir_vc = f"uvigs_{method}/03_vcontact2/{sample}"
        run_vcontact2(output_protein, output_dir_vc)

        ### 分类分配
        input_genomad_tax = os.path.join(output_dir_genomad, f"{sample}_uvigs_annotate/{sample}_uvigs_taxonomy.tsv")
        input_overview = os.path.join(output_dir_vc, "genome_by_genome_overview.csv")
        output_dir_assign = f"uvigs_{method}/04_assign_classification"
        run_assign_classification(input_genomad_tax, input_overview, output_dir_assign, sample, method, "uvigs")

        ### 提取最终结果
        final_taxonomy = merge_taxonomy = os.path.join(output_dir_assign, f"{sample}_merge_taxonomy.tsv")
        final_tax = f"uvigs_{method}/{sample}_taxonomy.tsv"
        if os.path.exists(final_taxonomy):
            subprocess.run(f"cp {final_taxonomy} {final_tax}", shell=True, check=True)

    # 合并日志,汇总有哪些是填充的
    log_dir = f"uvigs_{method}/log/assign_classification"
    merged_log_file = f"uvigs_{method}/output.log"
    os.makedirs(os.path.dirname(merged_log_file), exist_ok=True)  # 确保输出目录存在
    log_files = [os.path.join(log_dir, f"{sample}.log") for sample in SAMPLES]
    
    with open(merged_log_file, "w") as outfile:
        for i, log_file in enumerate(log_files):
            sample_name = SAMPLES[i]  # 获取对应的样本名
            if os.path.exists(log_file):
                outfile.write(f"=== Sample: {sample_name} ===\n")
                with open(log_file, "r") as infile:
                    outfile.write(infile.read())
                outfile.write("\n" + "="*50 + "\n\n")  # 添加更清晰的分隔符