# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate vs2

mkdir -p 01_run_vs2 && cd 01_run_vs2

time virsorter run -i ../00_rename_contigs/all_contigs.fna -w vs2_keep-l4k-s0.5 --keep-original-seq --include-groups dsDNAphage,ssDNA --min-length 4000 --min-score 0.5 -j 30 all

# real    987m31.932s
# user    17218m29.797s
# sys     502m17.363s

time virsorter run -i ../00_rename_contigs/all_contigs.fna -w vs2_l4k-s0.5 --include-groups dsDNAphage,ssDNA --min-length 4000 --min-score 0.5 -j 30 all
