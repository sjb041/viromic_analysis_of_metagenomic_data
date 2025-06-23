# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate DVF

contigspath="../02_assembly/02_contigs"
contigspath_addS="../02_assembly/02_contigs_addS"
out="01_dvf"
out_addS="01_dvf_addS"
sample=("F1_1A" "F1_2A" "F1_3A" "F2_1A" "F2_2A" "F2_3A" "FG_1A" "FG_2A" "FG_3A" "L1_1A" "L1_2A" "L1_3A" "L2_1A" "L2_2A" "L2_3A" "H_LX" "H_O")

#################### 运行DVF预测
# 注意使用双引号包裹整个命令字符串，确保变量在传递给parallel之前被解析。
echo ${sample[@]} | tr ' ' '\n' | parallel --load 80% "dvf.py -i ${contigspath}/{}_contigs.fasta -o ${out} -l 1500 -c 5"
echo ${sample[@]} | tr ' ' '\n' | parallel --load 80% "dvf.py -i ${contigspath_addS}/{}_contigs.fasta -o ${out_addS} -l 1500 -c 5"


#################### 根据阈值提取contigs
for i in ${sample[@]}
do
    # 筛选出符合的 contigID
    awk '$3 > 0.9 && $4 < 0.01 {print $1}' ${out}/${i}_contigs.fasta_gt1500bp_dvfpred.txt > ${out}/${i}.txt
 	awk '$3 > 0.9 && $4 < 0.01 {print $1}' ${out_addS}/${i}_contigs.fasta_gt1500bp_dvfpred.txt > ${out_addS}/${i}.txt
    
    # 提取contigs
    seqkit grep -f ${out}/${i}.txt ${contigspath}/${i}_contigs.fasta > ${out}/${i}_vir.fasta
    seqkit grep -f ${out_addS}/${i}.txt ${contigspath_addS}/${i}_contigs.fasta > ${out_addS}/${i}_vir.fasta

    # 错误检查
    expected=$(wc -l < ${out}/${i}.txt)
    real=$(grep -c '^>' ${out}/${i}_vir.fasta)
    if [ "$expected" -ne "$real" ]; then
        echo "检查 ${out}/${i}"
        exit 1
    fi

    expected_addS=$(wc -l < ${out_addS}/${i}.txt)
    real_addS=$(grep -c '^>' ${out_addS}/${i}_vir.fasta)
    if [ "$expected_addS" -ne "$real_addS" ]; then
        echo "检查 ${out_addS}/${i}"
        exit 1
    fi
done
