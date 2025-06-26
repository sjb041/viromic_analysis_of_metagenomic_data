# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate quast

mkdir -p 03_quality_assessment
mkdir -p 03_quality_assessment_addS

sample=("F1_1A" "F1_2A" "F1_3A" "F2_1A" "F2_2A" "F2_3A" "FG_1A" "FG_2A" "FG_3A" "L1_1A" "L1_2A" "L1_3A" "L2_1A" "L2_2A" "L2_3A" "H_LX" "H_O")

# 评估组装质量（不加不成对的）
echo ${sample[@]} | tr ' ' '\n' | parallel --load 80% '
    metaquast.py 02_contigs/{}_contigs.fasta \
        --max-ref-number 0 \
        -t 8 \
        -o "03_quality_assessment_test/{}" \
        --contig-thresholds 0,500,1000,1500,2000,2200,3000,3500,4000,5000,10000,25000,50000
'

# 合并评估报告
cut -f1 03_quality_assessment/F1_1A/report.tsv > 03_quality_assessment/report.tsv
for i in "${sample[@]}"
do
    paste 03_quality_assessment/report.tsv <(cut -f2 03_quality_assessment/${i}/report.tsv) > tmp.tsv
    mv tmp.tsv 03_quality_assessment/report.tsv
done


# 评估组装质量（加不成对的）
for i in "${sample[@]}"
do
    r1="../01_raw_sequence_preprocessing/03_cleandata/${i}_R1.fastq.gz"
    r2="../01_raw_sequence_preprocessing/03_cleandata/${i}_R2.fastq.gz"
    rs="../01_raw_sequence_preprocessing/03_cleandata/${i}_S.fastq.gz"
    out="03_quality_assessment_addS/${i}"

    #metaquast.py 02_contigs_addS/${i}_contigs.fasta --pe1 ${r1} --pe2 ${r2} --single ${rs} --max-ref-number 0 -t 30 -o ${out} --contig-thresholds 0,500,1000,1500,2000,2200,3000,3500,4000,5000,10000,25000,50000
    metaquast.py 02_contigs_addS/${i}_contigs.fasta --max-ref-number 0 -t 30 -o ${out} \
        --contig-thresholds 0,500,1000,1500,2000,2200,3000,3500,4000,5000,10000,25000,50000
done

# 合并评估报告
cut -f1 ${out}/report.tsv > 03_quality_assessment_addS/report_addS.tsv
for i in "${sample[@]}"
do
    paste 03_quality_assessment_addS/report_addS.tsv <(cut -f2 03_quality_assessment_addS/${i}/report.tsv) > tmp.tsv
    mv tmp.tsv 03_quality_assessment_addS/report_addS.tsv
done

########################## 代码解释 ###########################

#cut -f1 quast_out_7/report.tsv > test.tsv

#for i in 7 8
#do
#    paste test.tsv <(cut -f2 quast_out_${i}/report.tsv) > tmp.tsv
#    mv tmp.tsv test.tsv
#done

# 提取file2的第二列（test8数据），然后与file1合并
#paste quast_out_7/report.tsv <(cut -f2 quast_out_8/report.tsv) > merged.tsv

# 进程替换（Process Substitution） < 
# 整个 <(command) 结构会创建一个临时文件描述符，将 command 的输出作为文件传递给 paste。
# 明确指定了 paste 的两个输入源。
# 最后使用 > 将 paste 的结果重定向写入 merged.tsv。(注意：重定向 > 会首先清空 merged.tsv ,所以代码错误)
#paste merged.tsv <(cut -f1 quast_out_8/report.tsv) > merged.tsv

# 如果不使用 < ，可以改用 - ,表示从标准输入读取数据。如果省略 -，paste 会只读取 merged.tsv，而不会处理 cut 命令的输出。
#cut -f1 quast_out_8/report.tsv | paste merged.tsv -






