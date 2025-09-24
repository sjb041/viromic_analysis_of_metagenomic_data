# Prodigal + vConTACT2 + geNomad
# 处理病毒识别方法A和E的结果, all_uvigs

# 参考文献 Metagenomic analysis reveals gut phage diversity across three mammalian models

# 运行: snakemake -s pipeline_uvigs.smk --cores 5 --use-conda

import os
import pandas as pd
from contextlib import redirect_stdout

##################################### 定义函数
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

# 参考数据库的路径
db_dir = "/home/shijiabin/db"
genomad_db_path = os.path.join(db_dir, "genomad_db")

##################################### 定义流程规则
# 最终要生成的文件
rule all:
    input:
        expand("uvigs_{type}/{sample}_taxonomy.tsv", type=TYPES, sample=SAMPLES)

# genomad分类注释
rule genomad:
    input:
        virus = lambda wildcards: os.path.join(virus_dir, f"{wildcards.type}", f"{wildcards.sample}_uvigs.fna")
    output:
        dir = directory("uvigs_{type}/02-1_genomad/{sample}"),
        genomad_tax = "uvigs_{type}/02-1_genomad/{sample}/{sample}_uvigs_annotate/{sample}_uvigs_taxonomy.tsv"
    conda:
        "genomad"
    log:
        "uvigs_{type}/log/genomad/{sample}.log"
    shell:
        # 仅运行分类注释,允许分类到科级以下的分类单元
        "genomad annotate {input.virus} {output.dir} {genomad_db_path} --lenient-taxonomy > {log} 2>&1"

# 预测蛋白
rule prodigal:
    input:
        virus = lambda wildcards: os.path.join(virus_dir, f"{wildcards.type}", f"{wildcards.sample}_uvigs.fna")
    output:
        protein = "uvigs_{type}/02_prodigal/{sample}_protein.faa"
    log:
        "uvigs_{type}/log/prodigal/{sample}.log"
    shell:
        "prodigal -i {input.virus} -a {output.protein} -p meta > {log} 2>&1"

# 生成viral clusters (VCs)
rule vcontact2:
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
        vcontact2  --rel-mode 'Diamond' --db 'None' --pcs-mode MCL --vcs-mode ClusterONE --raw-proteins {input.protein} --proteins-fp {output.g2g} --output-dir {output.dir} >> {log} 2>&1
        """
    
rule assign_classification:
    input:
        genomad_tax = "uvigs_{type}/02-1_genomad/{sample}/{sample}_uvigs_annotate/{sample}_uvigs_taxonomy.tsv",
        overview = "uvigs_{type}/03_vcontact2/{sample}/genome_by_genome_overview.csv"
    output:
        merge_taxonomy = "uvigs_{type}/04_assign_classification/{sample}_merge_taxonomy.tsv",
        final_taxonomy = "uvigs_{type}/04_assign_classification/{sample}_final_taxonomy.tsv"
    log:
        "uvigs_{type}/log/assign_classification/{sample}.log"
    run:
        with open(log[0], "w") as f, redirect_stdout(f):
            assign_classification(input.genomad_tax, input.overview, output.merge_taxonomy, output.final_taxonomy)

# 提取结果文件
rule postprocess:
    input:
        final_taxonomy = "uvigs_{type}/04_assign_classification/{sample}_final_taxonomy.tsv"
    output:
        final_tax = "uvigs_{type}/{sample}_taxonomy.tsv"
    shell:
        "cp {input.final_taxonomy} {output.final_tax}"