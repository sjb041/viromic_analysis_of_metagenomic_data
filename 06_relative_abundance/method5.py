#!/usr/bin/env python3

# 计算相对丰度,公式5:
# 使用 Prodigal 预测 UVIGs 的 ORFs，输出核苷酸序列
# 使用 CD-HIT(v4.6.6)在95%核苷酸同一性和85%的比对区域上进行聚类(cd-hit -c 0.95 -n 5 -aS 0.85 -M 6000)
# 将原始测序reads比对到聚类质心序列
# 计算RPK和CPM
# 计算UViGs相对丰度
# 计算物种级 vOTU相对丰度

# 参考文献 2024-Fecal microbiota transplantation alters gut phage communities in a clinical trial for obesity


import os
import pandas as pd
import re
import sys
from Bio import SeqIO
import subprocess

# 读取基因长度
def get_gene_lengths(fasta_file):
    gene_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_lengths[record.id] = len(record.seq)
    return gene_lengths

# 解析 idxstats 文件
def parse_idxstats(idxstats_file):
    """解析 samtools idxstats 输出文件"""
    data = pd.read_csv(idxstats_file, sep='\t', header=None, 
                       names=['gene_id', 'length', 'mapped', 'unmapped'])
    # 过滤掉星号行
    data = data[data['gene_id'] != '*']
    return data

# 解析 flagstat 文件
def parse_flagstat(flagstat_file):
    """解析 samtools flagstat 输出文件"""
    stats = {}
    with open(flagstat_file, 'r') as f:
        for i, line in enumerate(f):
            parts = line.strip().split()       # 将每行内容按空格分割成多个部分
            if i == 0 and 'in total' in line:  # 第一行总是总读段数
                stats['total'] = int(parts[0])
            elif i == 4 and 'mapped' in line:  # 第五行总是比对读段数
                stats['mapped'] = int(parts[0])
    return stats

def extract_uvig_from_gene_id(gene_id):
    """
    从基因 ID 提取 UViG ID
    """
    try:
        # 去掉最后一个下划线及之后的内容
        match = re.match(r'(.+)_\d+$', gene_id)
        if match:
            return match.group(1)
        return gene_id  # 如果无法匹配，返回原始ID
    except Exception as e:
        print(f"警告: 提取UViG ID时出错，基因ID: {gene_id}, 错误: {e}")
        return "unknown"

def extract_sample_from_uvig_id(uvig_id):
    """
    从 UViG ID 提取样本名
    """
    try:
        # 将UViG ID按下划线分割
        parts = uvig_id.split('_')
        if len(parts) >= 2:
            # 取倒数第二个下划线之前的内容
            return parts[0] + '_' + parts[1]
        return "unknown"
    except Exception as e:
        print(f"警告: 提取样本名时出错，UViG ID: {uvig_id}, 错误: {e}")
        return "unknown"



# 输出目录
os.makedirs('04_abundance', exist_ok=True)

# --------------------------------------------- all viruses orf centroids --------------------------------------------- #
#### 预测 orf
# 输入
virus = "/home/shijiabin/2025_2ME/03_virus_identification/methodA/all_viruses.fna"
# 输出
dir_centroids = "01_orfs_centroids"
os.makedirs(dir_centroids, exist_ok=True)
orf = f"{dir_centroids}/all_viruses.orf.fna"

cmd_orf = (f"prodigal -i {virus} -d {orf} -p meta")
#subprocess.run(cmd_orf, shell=True, check=True)

#### 聚类获取 centroids
# 输入
orf
# 输出
centroids = f"{dir_centroids}/all_viruses.orf.centroids.fna"

cmd_centroids = (f"cd-hit-est -i {orf} -o {centroids} -c 0.95 -n 10 -aS 0.85 -d 0")
#subprocess.run(cmd_centroids, shell=True, check=True)

# --------------------------------------------- reads mapping centroids --------------------------------------------- #
# 样本名
SAMPLES = [
    "F1_1A", "F1_2A", "F1_3A",
    "F2_1A", "F2_2A", "F2_3A",
    "FG_1A", "FG_2A", "FG_3A",
    "L1_1A", "L1_2A", "L1_3A",
    "L2_1A", "L2_2A", "L2_3A",
    "H_LX", "H_O"
]

# 输入的 reads 的目录
# reads 目录,该目录下的 reads 文件名格式为 sample_R1.fastq.gz, sample_R2.fastq.gz
reads_dir = "/home/shijiabin/2025_2ME/01_raw_sequence_preprocessing/03_cleandata"

# 输出目录
dir_mapping = "03_mapping"
os.makedirs(dir_mapping, exist_ok=True)

# 建索引
cmd_bwa_index = f"bwa-mem2 index {centroids}"
subprocess.run(cmd_bwa_index, shell=True, check=True)

for sample in SAMPLES:
    # 输入
    r1 = os.path.join(reads_dir, f"{sample}_R1.fastq.gz")
    r2 = os.path.join(reads_dir, f"{sample}_R2.fastq.gz")
    centroids
    # 输出
    bam = f"{dir_mapping}/{sample}.bam"
    
    cmd_bwa = (
        f"bwa-mem2 mem -t 20 {centroids} {r1} {r2} | "
        f"samtools view -bS - | "
        f"samtools sort -@ 6 -o {bam} -"
    )
    subprocess.run(cmd_bwa, shell=True, check=True)

# --------------------------------------------- calculate_gene_cpm --------------------------------------------- #
# 输入
centroids

# 输出
dir_gene_cpm = "04_gene_cpm"
os.makedirs(dir_gene_cpm, exist_ok=True)
cpm_table = f"{dir_gene_cpm}/gene_cpm.tsv"
alignment_stats = f"{dir_gene_cpm}/alignment_stats.tsv"

# 获取聚类代表序列的长度
gene_lengths = get_gene_lengths(centroids)

# 计算 RPK (Reads Per Kilobase)
all_data = {}
for sample in SAMPLES:
    idxstats_file = f'03_mapping/{sample}.idxstats'
    if os.path.exists(idxstats_file):
        data = parse_idxstats(idxstats_file)
        data['rpk'] = data['mapped'] / (data['length'] / 1000)
        all_data[sample] = data

# 创建基因丰度矩阵
gene_abundance = pd.DataFrame()
for sample in SAMPLES:
    if sample in all_data:
        # 提取 gene_id 和 rpk 列
        sample_data = all_data[sample][['gene_id', 'rpk']]
        sample_data.columns = ['gene_id', sample]
        
        if gene_abundance.empty:
            gene_abundance = sample_data
        else:
            gene_abundance = pd.merge(gene_abundance, sample_data, on='gene_id', how='outer')

# 填充缺失值为0
gene_abundance = gene_abundance.fillna(0)

# 设置 gene_id 为索引
gene_abundance = gene_abundance.set_index('gene_id')
 
# 计算 CPM (Copies Per Million)
for sample in SAMPLES:
    if sample in gene_abundance.columns:
        # 计算样本总RPK
        total_rpk = gene_abundance[sample].sum()
        # 标准化为CPM
        gene_abundance[f'{sample}_cpm'] = (gene_abundance[sample] / total_rpk) * 1000000

gene_abundance.to_csv(cpm_table, sep='\t')       

# 输出每个样本的比对率统计
alignment_stats = pd.DataFrame(index=SAMPLES, columns=['total_reads', 'mapped_reads', 'mapping_rate'])
for sample in SAMPLES:
    flagstat_file = f'03_mapping/{sample}.flagstat'
    if os.path.exists(flagstat_file):
        stats = parse_flagstat(flagstat_file)
        total = stats.get('total', 0)
        mapped = stats.get('mapped', 0)
        rate = mapped / total if total > 0 else 0
        alignment_stats.loc[sample] = [total, mapped, rate]

alignment_stats.to_csv(alignment_stats, sep='\t')

# --------------------------------------------- calculate_uvig_abundance --------------------------------------------- #
# 输入
cpm_table
# 输出
dir_uvig_abun = "05_uvig_abun"
os.makedirs(dir_uvig_abun, exist_ok=True)
uvig_abun_table = f"{dir_uvig_abun}/uvig_abun.tsv"

# 将基因 ID 映射到 UViGs
gene_to_uvig = {}
for i, gene_id in enumerate(gene_abundance.index):
    uvig = extract_uvig_from_gene_id(gene_id)
    gene_to_uvig[gene_id] = uvig
    # 打印前几个示例以供检查
    if i < 5:
        sample = extract_sample_from_uvig_id(uvig)
        print(f"样例 {i+1}: {gene_id} -> UViG: {uvig}, 样本: {sample}")

# 添加uvig列
gene_abundance['uvig'] = gene_abundance.index.map(lambda x: gene_to_uvig.get(x, 'unknown'))

# 计算每个UViG的基因丰度中位数
uvig_abundance = pd.DataFrame()

for sample in SAMPLES:
    sample_cpm = f"{sample}_cpm"
    if sample_cpm in gene_abundance.columns:
        # 按UViG分组并计算中位数CPM
        uvig_grouped = gene_abundance.groupby('uvig')[sample_cpm].median()
        
        # 排除unknown组
        if 'unknown' in uvig_grouped.index:
            uvig_grouped = uvig_grouped.drop('unknown')
            
        if uvig_abundance.empty:
            uvig_abundance = pd.DataFrame(uvig_grouped)
            uvig_abundance.columns = [sample]
        else:
            uvig_abundance[sample] = uvig_grouped

# 保存UViG丰度结果
uvig_abundance.to_csv(uvig_abun_table, sep='\t')

# --------------------------------------------- calculate_votu_abundance --------------------------------------------- #
# 输入
votu_clusters_file = "/home/shijiabin/2025_2ME/03_virus_identification/methodA/06_cluster/clusters.tsv"
uvig_abun_table
# 输出
votu_abun_table = "votu_abun.tsv"

# 读取 vOTU 聚类信息，使用 clusters.tsv
# 第一列是各组中最长的 UViG（即 vOTU 代表序列）
# 第二列是该组所有 UViG 成员，以逗号分隔

votu_clusters = {}  # 创建一个空字典来存储vOTU聚类信息
# 读取文件
with open(votu_clusters_file, 'r') as f:
    for line_num, line in enumerate(f, 1):  # line_num: 行号从1开始，line: 接收当前行的内容(字符串)
        parts = line.strip().split('\t')
        if len(parts) < 2:
            print(f"错误: 在{votu_clusters_file}文件第{line_num}行发现格式错误，该行没有两列: {line.strip()}")
            sys.exit(1)
        
        votu_id = parts[0]  # vOTU代表序列作为vOTU ID
        members = parts[1].split(',')  # 所有成员UViGs
        votu_clusters[votu_id] = members

# 检查哪些成员 UViGs 不在 UViG 丰度数据中
missing_uvigs = []
for votu_id, members in votu_clusters.items():
    for uvig in members:
        if uvig not in uvig_abundance.index:
            missing_uvigs.append(uvig)
if missing_uvigs:
    print(f"注意：有{len(missing_uvigs)}个UViG不在丰度数据中") 

# 计算vOTU丰度
valid_uvigs = set(uvig_abundance.index)
votu_abundance = {}

for votu_id, members in votu_clusters.items():
    votu_abundance[votu_id] = {}
    for sample in SAMPLES:
        abundance_sum = sum(uvig_abundance.loc[uvig, sample] 
                            for uvig in members
                            if uvig in valid_uvigs)
        votu_abundance[votu_id][sample] = abundance_sum

# for uvig in members: 遍历当前vOTU包含的所有UViG成员
# if uvig in uvig_abundance.index: 检查该UViG是否存在于丰度数据中
# uvig_abundance.loc[uvig, sample]: 获取特定UViG在特定样本中的丰度值
# loc是pandas中的定位器，用于通过标签获取数据
# uvig是行索引(UViG ID)
# sample是列名(样本名称)
# 返回的是该UViG在该样本中的丰度值(浮点数)
    
# 转换为DataFrame
votu_df = pd.DataFrame.from_dict(votu_abundance, orient='index')

# 重置索引，将vOTU ID从索引转换为普通列
votu_df = votu_df.reset_index()
# 重命名索引列名为'vOTU'
votu_df = votu_df.rename(columns={'index': 'vOTU'})

# 为所有丰度列添加rel_abun_前缀（除了vOTU列）
abundance_columns = [col for col in votu_df.columns if col != 'vOTU']
new_column_names = {'vOTU': 'vOTU'}
new_column_names.update({col: f'rel_abun_{col}' for col in abundance_columns})
votu_df = votu_df.rename(columns=new_column_names)

votu_df.to_csv("votu_cpm.tsv", sep='\t', index=False)

# 计算每个样本的总和
sample_sums = votu_df.set_index('vOTU').sum(axis=0)

# 归一化到百分比
votu_df_percent = votu_df.set_index('vOTU').div(sample_sums, axis=1)
votu_df_percent = votu_df_percent.reset_index()

# 保存结果
votu_df_percent.to_csv(votu_abun_table, sep='\t', index=False)
