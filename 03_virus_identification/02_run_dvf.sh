# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate DVF
dvf.py -i all_contigs.fna -o 02_run_dvf -l 4000 -c 13

# 根据分数阈值和长度阈值筛选
cd 02_run_dvf
awk 'NR==1 || ($2 >= 4000 &&  $3 > 0.9 && $4 < 0.05)' all_contigs.fna_gt4000bp_dvfpred.txt > gt4000bp_s0.9_p0.05.txt
# 提取序列
tail -n +2 gt4000bp_s0.9_p0.05.txt | cut -f 1 | seqtk subseq ../all_contigs.fna - > dvf.fna && grep -c '^>' dvf.fna

cd ..
