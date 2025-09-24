# 统一重命名预测出的所有病毒,后续的所有分类注释方案均使用重命名后的序列

# 样本集
samples=("F1_1A" "F1_2A" "F1_3A" \
         "F2_1A" "F2_2A" "F2_3A" \
         "FG_1A" "FG_2A" "FG_3A" \
         "L1_1A" "L1_2A" "L1_3A" \
         "L2_1A" "L2_2A" "L2_3A" \
         "H_LX" "H_O")

# 3-2-A 得到的病毒
mkdir -p A

# 2-2-E 得到的病毒
mkdir -p E

dir="$HOME/2025_2ME/03_virus_identification"

# 序列重命名
for sample in "${samples[@]}"
do
    rename.py $dir/03_pipeline3-2/05_all_sample_results/${sample}_uvigs.fa -o A/${sample}_uvigs.fna -m A/${sample}_mapping.tsv
    
    rename.py $dir/03_pipeline2-2/07_all_sample_results/${sample}_uvigs.fa -o E/${sample}_uvigs.fna -m E/${sample}_mapping.tsv
done
