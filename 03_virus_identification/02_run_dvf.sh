# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate DVF
dvf.py -i all_contigs.fna -o 02_run_dvf -l 4000 -c 13
