##################################### 旧版的预测，每个样本单独预测，并且序列ID未格式化

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
#echo ${sample[@]} | tr ' ' '\n' | parallel --load 80% "dvf.py -i ${contigspath_addS}/{}_contigs.fasta -o ${out_addS} -l 1500 -c 5"

#################### 格式化序列ID
# 替换 name
for i in ${sample[@]}
do
    cat ${i}_contigs.fasta_gt1500bp_dvfpred.txt | cut -f 1 | awk 'NR==1{print $0; next} {split($1,a,"_"); $1=a[1]"_"a[2]}1' > tmp1.txt
    cat ${i}_contigs.fasta_gt1500bp_dvfpred.txt | cut -f 2,3,4 > tmp2.txt
    paste tmp1.txt tmp2.txt > ${i}_gt1500bp.txt
    sed "s/NODE_/${i}_contig_/g" ${i}_gt1500bp.txt > ${i}_gt1500bp_rename.txt && rm ${i}_gt1500bp.txt
done
rm tmp*.txt

# 合并所有样本
head -n 1 ${i}_gt1500bp_rename.txt > all_gt1500bp.txt
for i in ${sample[@]}
do
    tail -n +2 ${i}_gt1500bp_rename.txt >> all_gt1500bp.txt
done

# 计算合并前的行数
wc -l *_gt1500bp_rename.txt

# 计算合并后的行数
wc -l all_gt1500bp.txt

# 备份原始的 dvfpred 文件
mkdir -p backup_gt1500bp_dvfpred
mv ./*gt1500bp_dvfpred.txt ./backup_gt1500bp_dvfpred/

# 备份替换 name 之后的 dvfpred 文件
mkdir -p backup_rename
mv ./*gt1500bp_rename.txt ./backup_rename/
