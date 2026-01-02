# 执行此脚本的方法
# bash 01_quality_assessment_rawdata.sh

# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"

# 创建文件夹存放fastqc的输出
mkdir -p 01_quality_assessment_rawdata

# 质量评估
(
    conda activate fastqc

    # 查找所有文件，并将文件名存储在变量 files 中
    files=$(find ./00_rawdata -maxdepth 2 -type f \( -name '*_L1_1.fq.gz' -o -name '*_L1_2.fq.gz' \))

    # 打印看看找到的文件
    echo "$files"

    # 使用 GNU Parallel 并行处理文件
    echo "$files" | parallel -j 30 'fastqc {} -o 01_quality_assessment_rawdata -f fastq -t 2'
)


# 汇总评估报告
(
    conda activate multiqc
    multiqc -d 01_quality_assessment_rawdata/ -o 01_quality_assessment_rawdata/
)


# 求(mean ± SD），file的第五列是每个样本的 Total Sequences
file="./01_quality_assessment_rawdata/multiqc_data/multiqc_fastqc.txt"
awk -F'\t' 'NR>1 {sum+=$5; sumsq+=$5*$5; count++} END {mean=sum/count; sd=sqrt(sumsq/count - (mean*mean)); printf "均值 ± 标准差: %.2f ± %.2f\n", mean, sd}' "$file"

# 检查测序双端序列标签是否唯一
zcat ./00_rawdata/F1.1A/F1.1A_L1_1.fq.gz|head
zcat ./00_rawdata/F1.1A/F1.1A_L1_2.fq.gz|head
