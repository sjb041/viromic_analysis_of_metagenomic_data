
sample="F1_1A F1_2A F1_3A F2_1A F2_2A F2_3A FG_1A FG_2A FG_3A L1_1A L1_2A L1_3A L2_1A L2_2A L2_3A H_LX H_O"

for i in $sample
do
    mkdir $i
    # 提取评分表
    awk -v prefix="$i" 'NR == 1 || $1 ~ "^" prefix' all-viral-score.tsv > $i/final-viral-score.tsv
    # 提取序列
    seqkit grep -r -n -p "^$i" all-viral-combined.fa > $i/final-viral-combined.fa && seqkit stats $i/final-viral-combined.fa
    
    # 判断评分表行数与序列数
    # 获取 TSV 中有效记录数（减去表头）
    tsv_lines=$(wc -l < "$i/final-viral-score.tsv")
    tsv_records=$((tsv_lines - 1))
    # 获取 FASTA 中序列数量
    fa_records=$(seqkit stat "$i/final-viral-combined.fa" | awk 'NR==2 {print $4}')
    # 判断数量是否一致
    if [ "$tsv_records" -ne "$fa_records" ]; then
        echo "Error: Number of records mismatch in $i" >&2
        echo "TSV records: $tsv_records, FASTA sequences: $fa_records" >&2
        exit 1
    fi
done
