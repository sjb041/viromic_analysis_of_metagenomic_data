mkdir vs2_l10k-highonly && cd vs2_l10k-highonly

# 计算原始序列中每条序列的长度, 筛选出 ≥10kb 的ID, 并加上 ||
seqkit fx2tab -nl ../../all_contigs.fna | awk '$2 >= 10000' | cut -f1 | sed 's/$/||/' > contigs_filter_length.txt

# 根据阈值筛选vs2评分表的行
awk -F'\t' 'NR==1 || ($4!="nan" && $4>=0.9) || ($4!="nan" && $4>=0.7 && $7>=1)' ../vs2_l4k-s0.5/final-viral-score.tsv > final-viral-score-filtered.tsv

# 提取出筛选后的vs2评分表的ID，并去掉 || 之后的部分
cut -f1 final-viral-score-filtered.tsv | sed '1d' | sed 's/||.*/||/' > vs2_ids.txt

# 找出交集，正常情况下，intersection.txt 不会有重复的行，因为 comm 命令在处理已排序的文件时，每个相同的ID只会输出一次。
comm -12 <(sort vs2_ids.txt) <(sort contigs_filter_length.txt) > intersection.txt

# 提取评分表中对应行，-F：把 intersection.txt 中的内容当作固定字符串（而不是正则表达式）来匹配。-f intersection.txt：从文件中读取要匹配的字符串列表。
head -n 1 ../vs2_l4k-s0.5/final-viral-score.tsv > final-viral-score.tsv && grep -f intersection.txt ../vs2_l4k-s0.5/final-viral-score.tsv >> final-viral-score.tsv
head -n 1 ../vs2_l4k-s0.5/final-viral-boundary.tsv > final-viral-boundary.tsv && grep -f intersection.txt ../vs2_l4k-s0.5/final-viral-boundary.tsv >> final-viral-boundary.tsv

# 提取序列
cut -f1 final-viral-score.tsv | sed '1d' | seqtk subseq ../vs2_l4k-s0.5/final-viral-combined.fa - > final-viral-combined.fa && grep -c '^>' final-viral-combined.fa

rm contigs_filter_length.txt final-viral-score-filtered.tsv vs2_ids.txt intersection.txt
