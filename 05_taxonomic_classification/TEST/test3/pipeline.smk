# snakemake -s pipeline.smk --cores 5 --use-conda
##### 待改进部分：
# ① 原文献中聚类votu时使用75%AF，此处使用更常用的85%，是否需要照抄原文？
# ② 目前不选择最佳匹配

import os
import pandas as pd
##################################### 定义函数
### 处理 virushostdb 比对结果
# 选择比对的最佳匹配
def select_best_matches_df(input_df):
    """
    选出相同查询序列的最佳匹配，即bitscore最高的那一行
    
    参数:
        input_df: 输入数据框，包含比对结果
        
    返回:
        pd.DataFrame: 包含最佳匹配行的数据框
    """
    try:
        # 创建数据框的副本以避免修改原始数据
        df = input_df.copy()
        
        # 获取bitscore列的索引（最后一列）
        bitscore_col = df.columns[-1]  # 最后一列是bitscore
        
        idx = df.groupby(df.columns[0])[bitscore_col].idxmax()
        best_matches = df.loc[idx].reset_index(drop=True)

        return best_matches
    
    except Exception as e:
        print(f"选择最佳匹配时发生错误: {e}")
        raise


# 从比对结果中的第二列提取出refseqid并替换为第二列的值
def format_diamond_results_df(input_df):
    """
    验证输入数据框中每行第二列按"|"分隔后的元素个数是否一致,如果一致,则将第二列修改为倒数第二个元素
    去掉数据框第一列中的最后一个_及之后的部分,得到蛋白对应的病毒ID
    返回处理后的数据框

    参数:
        input_df: 输入的数据框，包含diamond比对结果
        
    返回:
        pd.DataFrame: 处理后的数据框，如果所有行的第二列元素个数一致则返回处理后的数据框，否则返回原始数据框
    """
    try:
        # 创建数据框的副本以避免修改原始数据
        df = input_df.copy()
        
        # 处理第一列，去掉最后一个_及之后的部分
        df.iloc[:, 0] = df.iloc[:, 0].apply(lambda x: x[:x.rfind('_')] if '_' in str(x) else x)
        
        # 处理第二列，按"|"分隔并检查元素个数一致性
        second_column = df.iloc[:, 1].astype(str)
        split_elements = second_column.str.split('|')
        
        # 获取第一行第二列分割后的元素个数作为基准
        expected_count = len(split_elements.iloc[0])
        
        # 检查所有行的第二列元素个数是否一致
        element_counts = split_elements.apply(len)
        is_valid = (element_counts == expected_count).all()
        
        if not is_valid:
            print("警告: 数据中存在第二列按'|'分隔后的元素个数不一致的行")
            # 打印不一致的行号
            inconsistent_rows = element_counts[element_counts != expected_count].index.tolist()
            print(f"不一致的行号: {inconsistent_rows[:10]}{'...' if len(inconsistent_rows) > 10 else ''}")
            return df  # 返回未修改的数据框
        
        # 如果所有行元素个数一致，修改第二列为倒数第二个元素
        if expected_count > 1:
            df.iloc[:, 1] = split_elements.apply(lambda x: x[-2])
        else:
            print("警告: 第二列按'|'分隔后的元素个数小于等于1，无法提取倒数第二个元素")
        
        return df
    
    except Exception as e:
        print(f"格式化比对结果时发生错误: {e}")
        return input_df


# 合并diamond结果和refseq分类信息
def merge_diamond_refseq_df(refseq_df, diamond_results_df):
    """
    根据预处理后的diamond结果表的第二列匹配refseq id 对应的分类信息表的对应行,（匹配不到则报错）
    将refseq id 对应的分类信息表的第2、3列添加到预处理后的diamond结果表中。

    参数:
        refseq_df: refseq id 对应的分类信息数据框（包含至少3列）
        diamond_results_df: 预处理后的diamond结果数据框
        
    返回:
        pd.DataFrame: 合并后的数据框，包含原始diamond结果和对应的分类信息
    """
    try:
        # 创建refseq ID到family和分类谱系的映射字典
        refseq_family_map = {}
        for index, row in refseq_df.iterrows():
            if len(row) >= 3:
                refseq_family_map[row.iloc[0]] = (row.iloc[1], row.iloc[2])
        
        # 创建结果数据框的副本
        result_df = diamond_results_df.copy()
        
        # 为结果数据框添加新的列用于存储family和lineage信息
        family_list = []
        lineage_list = []
        
        # 遍历diamond结果数据框的每一行
        for index, row in result_df.iterrows():
            id_to_match = row.iloc[1]  # 第二列是refseq ID
            
            # 查找匹配项
            if id_to_match in refseq_family_map:
                family, lineage = refseq_family_map[id_to_match]
                family_list.append(family)
                lineage_list.append(lineage)
            else:
                # 如果没有匹配项，抛出错误
                raise ValueError(f"未找到与ID '{id_to_match}' 匹配的refseq id")
        
        # 将family和lineage信息添加到结果数据框
        result_df['family'] = family_list
        result_df['lineage'] = lineage_list

        # 仅保留前两列和倒数两列
        result_df = result_df.iloc[:, [0, 1, -2, -1]]

        # 添加第一列和第二列的标题行分别为 votu_id, refseq
        result_df.columns = ['votu_id', 'refseq', 'family', 'lineage']

        return result_df
        
    except Exception as e:
        print(f"合并数据框时发生错误: {e}")
        return None

### 处理 crasslike_db 比对结果
def process_crasslike_df(input_df):
    """
    处理crasslike数据框，提取所需列并添加分类信息
    
    参数:
        input_df: 输入的数据框，包含crasslike比对结果
        
    返回:
        pd.DataFrame: 处理后的数据框，包含votu_id, refseq, family, lineage列
    """
    try:
        # 创建数据框的副本以避免修改原始数据
        df = input_df.copy()
        
        # 仅保留1,2列并添加列名
        df = df.iloc[:, [0, 1]]
        df.columns = ['votu_id', 'refseq']
        
        # 修改第一列，去掉最后一个_及之后的部分，得到病毒ID
        df['votu_id'] = df['votu_id'].apply(lambda x: x[:x.rfind('_')] if '_' in str(x) else x)
        
        # 修改第二列，取第一个空格前的部分
        df['refseq'] = df['refseq'].str.split(' ').str[0]
        
        # 添加 family 列，值为 unclassified Crassvirales
        df['family'] = 'unclassified Crassvirales'
        
        # 添加 lineage 列
        df['lineage'] = 'Duplodnaviria; Heunggongvirae; Uroviricota; Caudoviricetes; Crassvirales; unclassified Crassvirales; ;'
        
        return df
    
    except Exception as e:
        print(f"处理crasslike数据框时发生错误: {e}")
        raise

### 处理 fromBenler_db 比对结果
def process_frombenler_df(input_df, fromBenler_formated):
    """
    处理fromBenler数据框，提取所需列并添加分类信息
    
    参数:
        input_df: 输入的数据框，包含fromBenler比对结果
        
    返回:
        pd.DataFrame: 处理后的数据框，包含votu_id, refseq, family, lineage列
    """
    # 创建数据框的副本以避免修改原始数据
    df = input_df.copy()
    
    # 仅保留1,2列并添加列名
    df = df.iloc[:, [0, 1]]
    df.columns = ['votu_id', 'refseq']
    
    # 修改第一列，去掉最后一个_及之后的部分，得到病毒ID
    df['votu_id'] = df['votu_id'].apply(lambda x: x[:x.rfind('_')] if '_' in str(x) else x)
    
    # 修改第二列，使用正则表达式提取"fromBenler_"与第二个下划线之间的部分
    df['refseq'] = df['refseq'].str.extract(r'fromBenler_([^_]+)_')
    
    # 根据refseq从fromBenler_formated中匹配，添加各分类等级
    # 合并
    df = pd.merge(df, fromBenler_formated, left_on='refseq', right_on='representative', how='left')
    # 检查
    if df['representative'].isna().any():
        unmatched_count = df['representative'].isna().sum()
        raise ValueError(f"发现 {unmatched_count} 个未匹配的 refseq 值")
    
    # 添加lineage列，
    # 如果Phylum列的值为 Uroviricota，则lineage值为 Duplodnaviria;Heunggongvirae;Uroviricota
    # 如果Phylum列的值为 Phixviricota，则lineage值为 Monodnaviria;Sangervirae;Phixviricota
    # 如果Phylum列的值为 Hofneiviricota，则lineage值为 Monodnaviria;Loebvirae;Hofneiviricota
    df['lineage'] = df['Phylum'].replace({
        'Uroviricota': 'Duplodnaviria; Heunggongvirae; Uroviricota',
        'Phixviricota': 'Monodnaviria; Sangervirae; Phixviricota',
        'Hofneiviricota': 'Monodnaviria; Loebvirae; Hofneiviricota'
    })
    
    # 将 Class Order Family Genus 列的值依次加入到 lineage 列中，（先添加“; ”再加入，如果为缺失值，则直接将“; ”加入）
    for col in ['Class', 'Order', 'Family', 'Genus']:
        df['lineage'] += '; ' + df[col].fillna('').astype(str)
    
    # 仅保留 votu_id refseq Family lineage 列
    df = df[['votu_id', 'refseq', 'Family', 'lineage']]

    # 将 Family 列中的缺失值替换为 Unclassified
    df['Family'] = df['Family'].fillna('Unclassified')

    # 修改 Family 为 family
    df = df.rename(columns={'Family': 'family'})
    
    return df

### 分配科分类和完整分类
# 分配科分类和完整分类
def assign_family_classification(merged_file_path, output_path):
    """
    对于同一个查询序列的所有匹配，如果有至少25%匹配的科水平的分类信息相同，
    那么该查询序列的科级分类就是这个匹配的科级分类。
    对于未满足25%频率条件的查询序列，科级分类标记为Unclassified，
    但完整分类等级设置为频率最高的科首次匹配的完整分类。
    输出一个新的表格，第一列是查询序列，第二列是科，第三列是科对应的所有分类等级（原表格的最后一列）

    参数:
        merged_file_path: 输入TSV文件路径
        output_path: 输出TSV文件路径
    返回:
        bool: 处理是否成功
    """
    try:
        # 使用pandas读取文件
        df = pd.read_csv(merged_file_path, sep='\t')
        
        # 按votu_id分组
        query_groups = {}
        for index, row in df.iterrows():
            votu_id = row['votu_id']
            family = row['family']
            lineage = row['lineage']
            
            if votu_id not in query_groups:
                query_groups[votu_id] = []
            
            query_groups[votu_id].append({
                'family': family,
                'lineage': lineage
            })
        
        # 处理每个查询序列，确定科级分类
        results = []
        for votu_id, matches in query_groups.items():
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
            
            # 检查是否有科出现频率超过25%
            assigned_family = 'Unclassified'
            assigned_lineage = most_frequent_lineage  # 默认使用频率最高的科的分类等级
            
            for family, count in family_count.items():
                frequency = count / total_matches
                if frequency >= 0.25:
                    assigned_family = family
                    # 获取该科对应的分类等级
                    for match in matches:
                        if match['family'] == family:
                            assigned_lineage = match['lineage']
                            break
                    break
            
            results.append({
                'votu_id': votu_id,
                'family': assigned_family,
                'lineage': assigned_lineage
            })
        
        # 创建结果DataFrame并保存
        result_df = pd.DataFrame(results)
        result_df.to_csv(output_path, sep='\t', index=False)
        
        return True
    except Exception as e:
        print(f"处理文件时发生错误: {e}")
        return False
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

# uvigs的父目录，要求文件名的格式为sample_uvigs.fa
virus_dir = "/home/shijiabin/2025_2ME/03_virus_identification/03_pipeline2-2/07_all_sample_results"

# 参考数据库的路径
db_dir = "/home/shijiabin/db/three_viral_databases"
# 蛋白序列
proteins_db_path = f"{db_dir}/proteins.faa"
# 分类信息
virushostdb_formated_path = f"{db_dir}/virushostdb_formated.tsv"
fromBenler_formated_path = f"{db_dir}/fromBenler_formated.csv"

##################################### 定义流程规则
onstart:
    shell("echo '\n开始处理样本 {SAMPLES}'")

# 最终要生成的文件
rule all:
    input:
        expand("04_classification/{sample}_final.tsv", sample=SAMPLES)

# 聚类
rule cluster:
    input:
        virus = lambda wildcards: os.path.join(virus_dir, f"{wildcards.sample}_uvigs.fa")
    output:
        dir = directory("01_cluster/{sample}"),
        votus = "01_cluster/{sample}/votus.fa"
    log:
        "log/cluster/{sample}.log"
    shell:
        "cluster.sh -i {input.virus} -t {input.virus} -a 95 -o {output.dir} > {log} 2>&1"

# votus预测蛋白
rule prodigal:
    input:
        votus = "01_cluster/{sample}/votus.fa"
    output:
        protein = "02_prodigal/{sample}_protein.faa"
    log:
        "log/prodigal/{sample}.log"
    shell:
        "prodigal -i {input.votus} -a {output.protein} -p meta > {log} 2>&1"

# 蛋白比对到three viral databases
rule diamond:
    input:
        protein = "02_prodigal/{sample}_protein.faa"
    output:
        diamond_results = "03_diamond/{sample}.tsv"
    log:
        "log/diamond/{sample}.log"
    shell:
        """
        diamond blastp \
            -q {input.protein} \
            -d {proteins_db_path} \
            --id 30 \
            --query-cover 50 \
            --min-score 50 \
            --max-target-seqs 10 \
            --outfmt 6 qseqid stitle pident length evalue bitscore \
            -o {output.diamond_results} > {log} 2>&1
        """

# 基于diamond结果分配科分类和完整分类
rule assign_classification:
    input:
        diamond_results = "03_diamond/{sample}.tsv"
    output:
        merge_results = "04_classification/{sample}_merge.tsv",
        final_results = "04_classification/{sample}_final.tsv"
    run:
        # 读入 diamond_results
        df = pd.read_csv(input.diamond_results, sep="\t", header=None)

        ### 处理 virushostdb 比对结果
        # 筛选出第二列以“virushostdb”开头的
        virushostdb_df = df[df.iloc[:, 1].str.startswith("virushostdb")]

        # 从比对结果中的第二列提取出refseqid并替换为第二列的值
        virushostdb_df_1 = format_diamond_results_df(virushostdb_df)

        # 合并diamond结果和refseq分类信息
        refseq_df = pd.read_csv(virushostdb_formated_path, sep="\t", header=None)
        virushostdb_df_2 = merge_diamond_refseq_df(refseq_df, virushostdb_df_1)

        ### 处理 crasslike_db 比对结果
        crasslike_df = df[df.iloc[:, 1].str.startswith("crasslike")]
        crasslike_df_1 = process_crasslike_df(crasslike_df)
        
        ### 处理 fromBenler_db 比对结果
        fromBenler_df = df[df.iloc[:, 1].str.startswith("fromBenler")]
        fromBenler_formated = pd.read_csv(fromBenler_formated_path)
        fromBenler_df_1 = process_frombenler_df(fromBenler_df, fromBenler_formated)

        ### 合并三个结果
        merged_df = pd.concat([virushostdb_df_2, crasslike_df_1, fromBenler_df_1], ignore_index=True)
        # 按 votu_id 重新排序
        merged_df = merged_df.sort_values(by='votu_id')

        # 保存结果
        merged_df.to_csv(output.merge_results, sep="\t", index=False)

        # 分配科级分类和完整分类
        assign_family_classification(output.merge_results, output.final_results)
