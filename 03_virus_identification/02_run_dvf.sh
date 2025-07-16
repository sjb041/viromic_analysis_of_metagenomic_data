# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate DVF

mkdir -p 02_run_dvf

dvf.py -i all_contigs.fna -o 02_run_dvf -l 1500 -c 15
mv all_contigs.fna_gt1500bp_dvfpred.txt all_gt1500bp.txt
