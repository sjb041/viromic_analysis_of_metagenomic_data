mkdir -p 04_quality_assessment_cleandata/fastqc_paired
mkdir -p 04_quality_assessment_cleandata/fastqc_single

# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"

# 样本名
sample="F1.1A F1.2A F1.3A F2.1A F2.2A F2.3A FG.1A FG.2A FG.3A L1.1A L1.2A L1.3A L2.1A L2.2A L2.3A H.LX H.O"

# 获取评估报告
for i in $sample
do
    mv 02_quality_control_rawdata/${i}/fastqc/${i}_L1_1_kneaddata_paired* 04_quality_assessment_cleandata/fastqc_paired/ &
    mv 02_quality_control_rawdata/${i}/fastqc/${i}_L1_1_kneaddata_unmatched* 04_quality_assessment_cleandata/fastqc_single/ &
done
wait

# 汇总评估报告
conda activate multiqc

multiqc -d 04_quality_assessment_cleandata/fastqc_paired/ -o 04_quality_assessment_cleandata/fastqc_paired/
multiqc -d 04_quality_assessment_cleandata/fastqc_single/ -o 04_quality_assessment_cleandata/fastqc_single/

# 求(mean ± SD），file的第五列是每个样本的 Total Sequences
file="./04_quality_assessment_cleandata/fastqc_paired/multiqc_data/multiqc_fastqc.txt"
awk -F'\t' 'NR>1 {sum+=$5; sumsq+=$5*$5; count++} END {mean=sum/count; sd=sqrt(sumsq/count - (mean*mean)); printf "均值 ± 标准差: %.2f ± %.2f\n", mean, sd}' "$file"
