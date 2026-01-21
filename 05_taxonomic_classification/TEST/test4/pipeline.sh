#!/usr/bin/env bash

# 样本集
#SAMPLE=("F1_1A" "F1_2A" "F1_3A" \
#        "F2_1A" "F2_2A" "F2_3A" \
#        "FG_1A" "FG_2A" "FG_3A" \
#        "L1_1A" "L1_2A" "L1_3A" \
#        "L2_1A" "L2_2A" "L2_3A" \
SAMPLE=("H_LX" "H_O")

# uvigs的父目录，要求文件名的格式为sample_uvigs.fa
virus_dir="/home/shijiabin/2025_2ME/03_virus_identification/03_pipeline2-2/07_all_sample_results"

# virify 的位置
virify="/home/shijiabin/opt/emg-viral-pipeline-3.0.2/main.nf"

# 参考数据库的路径
db_dir="/home/shijiabin/db/virify_db"

# hmm 数据库版本 vpHMM_database_v{1,2,3}.hmm
viphog_version="v3"

# 数据库分类信息版本 additional_data_vpHMMs_v{1,2,3,4}.tsv
meta_version="v4"

# 日志
cluster_log="log/cluster" && mkdir -p $cluster_log
rename_log="log/rename" && mkdir -p $rename_log
virify_log="log/virify" && mkdir -p $virify_log


for sample in "${SAMPLE[@]}"
do
    echo -e "\n------------------- ${sample} 开始 -------------------"

    # ------------------------ cluster
    echo " 聚类..."
    # 输入
    virus="${virus_dir}/${sample}_uvigs.fa"

    # 输出:
    # 目录
    cluster_dir="01_cluster/${sample}" && mkdir -p ${cluster_dir}
    # 聚类后得到的votu
    votus="${cluster_dir}/votus.fa"

    # 执行聚类
    cluster.sh -i ${virus} -t ${virus} -a 95 -o ${cluster_dir} > ${cluster_log}/${sample}.log 2>&1

    # ------------------------ 序列重命名
    echo " 序列重命名..."
    # 输入：votus
    
    # 输出：
    # 目录
    rename_dir="02_rename" && mkdir -p ${rename_dir}
    # 映射文件
    mapping="${rename_dir}/${sample}_mapping.tsv"
    # 提取出的原名
    original_name="${rename_dir}/${sample}_original_name.txt"
    # 重命名后的votus
    renamed_fasta="${rename_dir}/${sample}_renamed.fa"

    # 执行重命名
    rename.py ${votus} -m ${mapping} -o ${renamed_fasta} > ${rename_log}/${sample}.log 2>&1
        
    # 取映射文件第1列
    cut -f 1 ${mapping} > ${original_name}

    # ------------------------ virify 分类注释
    echo " virify 分类注释..."
    # 输入：renamed_fasta

    # 输出：
    # 目录
    virify_dir="03_virify/${viphog_version}${meta_version}_ncbi2022/${sample}" && mkdir -p ${virify_dir}
    # 分类表
    taxonomy="${virify_dir}/renamed_taxonomy.tsv"

    # 执行分类注释
    nextflow run ${virify} \
        --fasta  ${renamed_fasta} \
        --output  ${virify_dir} \
        --viphog_version ${viphog_version} \
        --meta_version ${meta_version} \
        --databases ${db_dir} \
        --onlyannotate true \
        --length 0 \
        --publish_all
        # > ${virify_log}/${sample}.log 2>&1
    
    # 提取分类表
    mv ${virify_dir}/${sample}_renamed/06-taxonomy/*_taxonomy.tsv ${taxonomy}

    # ------------------------ 还原 taxonomy 中的序列名
    echo " 还原 taxonomy 中的序列名..."
    # 输入：taxonomy, original_name

    # 输出:
    # 恢复原名的分类表
    final_taxonomy="${virify_dir}/${sample}_final_taxonomy.tsv"

    # 按第一列排序分类表
    head -n 1 ${taxonomy} > ${virify_dir}/tmp
    tail -n +2 ${taxonomy} | sort -t$'\t' -k1,1V >> ${virify_dir}/tmp

    # 使用 original_name 替换 taxonomy 的第一列
    cut -f 2- ${virify_dir}/tmp | paste ${original_name} - > ${final_taxonomy}

    # 删除临时文件
    rm ${virify_dir}/tmp

    echo " ${sample} 完成"
done
