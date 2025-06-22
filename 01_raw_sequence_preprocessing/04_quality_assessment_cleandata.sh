mkdir -p 04_quality_assessment_cleandata/fastqc_paired
mkdir -p 04_quality_assessment_cleandata/fastqc_single

# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"

# 质量评估
(
    conda activate fastqc

    # 查找要评估的文件，并将文件名存储在变量中
    files=$(find ./03_cleandata -type f \( -name '*R1.fastq.gz' -o -name '*R2.fastq.gz' \))
    files_s=$(find ./03_cleandata -type f \( -name '*s1.fastq.gz' -o -name '*s2.fastq.gz' \))

    # 打印看看找到的文件
    echo "$files"
    # 使用 GNU Parallel 并行评估
    echo "$files" | parallel -j 30 'fastqc {} -o 04_quality_assessment_cleandata/fastqc_paired -f fastq -t 2'

    # 打印看看找到的文件
    echo "$files_s"
    # 使用 GNU Parallel 并行评估
    echo "$files_s" | parallel -j 30 'fastqc {} -o 04_quality_assessment_cleandata/fastqc_single -f fastq -t 2'
)

# 汇总评估报告
(
    conda activate multiqc
    multiqc -d 04_quality_assessment_cleandata/fastqc_paired -o 04_quality_assessment_cleandata/fastqc_paired/
    multiqc -d 04_quality_assessment_cleandata/fastqc_single -o 04_quality_assessment_cleandata/fastqc_single/
)

# 求(mean ± SD），file的第五列是每个样本的 Total Sequences
file="./04_quality_assessment_cleandata/fastqc_paired/multiqc_data/multiqc_fastqc.txt"
awk -F'\t' 'NR>1 {sum+=$5; sumsq+=$5*$5; count++} END {mean=sum/count; sd=sqrt(sumsq/count - (mean*mean)); printf "均值 ± 标准差: %.2f ± %.2f\n", mean, sd}' "$file"
