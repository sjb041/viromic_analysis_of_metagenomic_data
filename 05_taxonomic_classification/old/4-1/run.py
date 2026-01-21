# MMseqs2构建蛋白质数据库进行分类分配
# 处理病毒识别方法A和E的结果, all_uvigs and votus

# 来源: 本想复现 Metagenomic compendium of 189,680 DNA viruses from the human gut microbiome 的病毒分类注释方法,
# 但是作者推荐使用MMseqs2, https://github.com/snayfach/MGV/issues/1

# 运行: python run.py

import os
import subprocess

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
db = "/home/shijiabin/db/mmseqs_db/ictv_nr_db/ictv_nr_db"

# 聚类并提取votu
def run_cluster(input_virus, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    cmd = f"cluster.sh -i {input_virus} -t {input_virus} -a 95 -o {output_dir}"
    subprocess.run(cmd, shell=True, check=True)

# MMseqs2分类注释 
def run_mmseqs_easy_taxonomy(input_virus, db, output_dir, prefix):
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = f"{output_dir}/{prefix}"
    tmp = f"{output_dir}/tmp"
    
    # 不加双引号会导致目录解析不正确
    cmd = f'mmseqs easy-taxonomy "{input_virus}" "{db}" "{output_prefix}" "{tmp}" --blacklist "" --tax-lineage 1 --lca-ranks family --threads 20'
    subprocess.run(cmd, shell=True, check=True)

# 对votus执行整个流程
for method in METHODS:
    for sample in SAMPLES:

        ### 聚类并提取votu
        input_virus = os.path.join(virus_dir, f"{method}", f"{sample}_uvigs.fna")
        output_dir_cluster = f"votus_{method}/01_cluster/{sample}"
        run_cluster(input_virus, output_dir_cluster)
        
        ### MMseqs2分类注释 
        input_votus = os.path.join(output_dir_cluster, "votus.fa")
        output_dir_mmseqs_easy_taxonomy = f"votus_{method}/02_mmseqs_easy_taxonomy"
        output_prefix = sample
        run_mmseqs_easy_taxonomy(input_votus, db, output_dir_mmseqs_easy_taxonomy, output_prefix)


# 对uvigs执行整个流程
for method in METHODS:
    for sample in SAMPLES:

        ### MMseqs2分类注释 
        input_virus = os.path.join(virus_dir, f"{method}", f"{sample}_uvigs.fna")
        output_dir_mmseqs_easy_taxonomy = f"uvigs_{method}/02_mmseqs_easy_taxonomy"
        output_prefix = sample
        run_mmseqs_easy_taxonomy(input_virus, db, output_dir_mmseqs_easy_taxonomy, output_prefix)
