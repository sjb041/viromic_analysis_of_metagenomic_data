# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate DVF

contigspath="../02_assembly/02_contigs"
contigspath_addS="../02_assembly/02_contigs_addS"
out="01_dvf"
out_addS="01_dvf_addS"

sample=("F1_1A" "F2_1A" "FG_1A" "H_LX" "L1_1A" "L2_1A")

# 注意使用双引号包裹整个命令字符串，确保变量在传递给parallel之前被解析。
echo ${sample[@]} | tr ' ' '\n' | parallel --load 80% "dvf.py -i ${contigspath}/{}_contigs.fasta -o ${out} -l 1500 -c 5"

echo ${sample[@]} | tr ' ' '\n' | parallel --load 80% "dvf.py -i ${contigspath_addS}/{}_contigs.fasta -o ${out_addS} -l 1500 -c 5"
