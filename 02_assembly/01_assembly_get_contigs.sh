# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate spades

mkdir -p 01_assembly_results
mkdir -p 02_contigs
mkdir -p 01_assembly_results_addS
mkdir -p 02_contigs_addS

sample="F1_1A F1_2A F1_3A F2_1A F2_2A F2_3A FG_1A FG_2A FG_3A L1_1A L1_2A L1_3A L2_1A L2_2A L2_3A H_LX H_O"

# 组装
#for i in "${sample[@]}"
for i in $sample
do
    in1="../01_raw_sequence_preprocessing/03_cleandata/${i}_R1.fastq.gz"
    in2="../01_raw_sequence_preprocessing/03_cleandata/${i}_R2.fastq.gz"
    ins="../01_raw_sequence_preprocessing/03_cleandata/${i}_S.fastq.gz"
    out="01_assembly_results/${i}"
    out_addS="01_assembly_results_addS/${i}"

    spades.py --meta -1 $in1  -2 $in2 -o $out -k 21,33,55,77 -m 100 -t 20
    spades.py --meta -1 $in1  -2 $in2 -s $ins -k 21,33,55,77 -m 100 -t 20 -o $out_addS 

    # 提取 contigs
    mv ${out}/contigs.fasta 02_contigs/${i}_contigs.fasta
    mv ${out}/scaffolds.fasta 02_contigs/${i}_scaffolds.fasta
    mv ${out_addS}/contigs.fasta 02_contigs_addS/${i}_contigs.fasta
    mv ${out_addS}/scaffolds.fasta 02_contigs_addS/${i}_scaffolds.fasta

    #rm -r $out
    #rm -r $out_addS
done
