# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"

# 获取所有样本的contigs
cat ../02_assembly/02_contigs/rename_contigs/*.fna > all_contigs.fna

mkdir -p 01_run_vs2 && cd 01_run_vs2

conda activate vs2
time virsorter run -i ../all_contigs.fna -w vs2_l4k-s0.5 --include-groups dsDNAphage,ssDNA --min-length 4000 --min-score 0.5 -j 30 all

# real    957m36.394s
# user    16872m10.353s
# sys     498m56.884s

# 如果只是想使用--keep-original-seq，其它参数不变，可以根据得分表从all_contigs.fna中获取那些||full序列，而不用重新跑一次
time virsorter run -i ../all_contigs.fna -w vs2_keep-l4k-s0.5 --keep-original-seq --include-groups dsDNAphage,ssDNA --min-length 4000 --min-score 0.5 -j 30 all

# real    987m31.932s
# user    17218m29.797s
# sys     502m17.363s

cd ..
