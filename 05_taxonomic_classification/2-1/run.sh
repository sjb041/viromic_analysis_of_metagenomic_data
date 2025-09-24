# 自动版的virify
# 使用v3版数据库和v4版分类信息(2022版)
# 处理病毒识别方法A和E的结果, all_uvigs and votus

# 参考文献 VIRify: An integrated detection, annotation and taxonomic classification pipeline using virus-specific protein profile hidden Markov models

# 运行: bash run.sh

# virify 的位置
virify="/home/shijiabin/opt/emg-viral-pipeline-3.0.2/main.nf"

# 参考数据库的路径
db_dir="/home/shijiabin/db/virify_db"

# hmm 数据库版本 vpHMM_database_v{1,2,3}.hmm
viphog_version="v3"

# 数据库分类信息版本 additional_data_vpHMMs_v{1,2,3,4}.tsv
meta_version="v4"

# 样本集
samples=(
    "F1_1A" "F1_2A" "F1_3A"
    "F2_1A" "F2_2A" "F2_3A"
    "FG_1A" "FG_2A" "FG_3A"
    "L1_1A" "L1_2A" "L1_3A"
    "L2_1A" "L2_2A" "L2_3A"
    "H_LX" "H_O"
)

methods=("A", "E")

# 创建输出目录
for method in methods
do
    mkdir -p all_uvigs_$method votus_$method
done

# 处理A和E两种类型,分别是病毒识别方法A和E的结果
for method in methods
do
    ########################## virify注释votus ##########################
    # 聚类并提取votu,然后创建Samplesheet
    echo "id,assembly" > votus_${method}/samplesheet.csv
    for sample in "${samples[@]}"
    do
        # 输入文件
        uvigs="$HOME/2025_2ME/05_taxonomic_classification/UVIGs/${method}/${sample}_uvigs.fna"
        
        # 输出目录
        cluster_dir="$HOME/2025_2ME/05_taxonomic_classification/2-1/votus_${method}/01_cluster/${sample}"
        mkdir -p ${cluster_dir}

        # 执行聚类
        cluster.sh -i ${uvigs} -t ${uvigs} -a 95 -o ${cluster_dir}

        # 创建Samplesheet
        echo "${sample},${cluster_dir}/votus.fa" >> votus_${method}/samplesheet.csv
    done

    # 输出目录
    virify_dir="votus_${method}/02_virify_${viphog_version}${meta_version}_ncbi2022"
    mkdir -p ${virify_dir}    
    
    # 运行virify
    nextflow run ${virify} \
        --samplesheet  votus_${i}/samplesheet.csv \
        --output  ${virify_dir} \
        --viphog_version ${viphog_version} \
        --meta_version ${meta_version} \
        --databases ${db_dir} \
        --onlyannotate true \
        --length 0 \
        --publish_all
    
    # 提取分类表
    for sample in "${samples[@]}"
    do
        cp ${virify_dir}/${sample}/06-taxonomy/*_taxonomy.tsv ${virify_dir}/${sample}_taxonomy.tsv
    done
    
    ########################## virify注释all_uvigs ##########################
    # 创建Samplesheet,样本名与序列文件的映射,用于virify输入
    echo "id,assembly" > uvigs_${method}/samplesheet.csv
    for sample in "${samples[@]}"
    do
        # uvigs的路径
        uvigs="$HOME/2025_2ME/05_taxonomic_classification/UVIGs/${method}/${sample}_uvigs.fna"
        # 写入Samplesheet
        echo "${sample},${uvigs}" >> uvigs_${method}/samplesheet.csv
    done
    
    # 输出目录
    virify_dir="uvigs_${method}/virify_${viphog_version}${meta_version}_ncbi2022"
    mkdir -p ${virify_dir}
        
    # 运行virify
    nextflow run ${virify} \
        --samplesheet  uvigs_${i}/samplesheet.csv \
        --output  ${virify_dir} \
        --viphog_version ${viphog_version} \
        --meta_version ${meta_version} \
        --databases ${db_dir} \
        --onlyannotate true \
        --length 0 \
        --publish_all
        
    # 提取分类表
    for sample in "${samples[@]}"
    do
        cp ${virify_dir}/${sample}/06-taxonomy/*_taxonomy.tsv ${virify_dir}/${sample}_taxonomy.tsv
    done
done
