mkdir -p 00_rename_contigs && cd 00_rename_contigs

# 将组装后的 contig 复制过来
cp ../../02_assembly/02_contigs/*_contigs* ./

# 序列重命名
for sample in F1_1A F1_2A F1_3A F2_1A F2_2A F2_3A FG_1A FG_2A FG_3A L1_1A L1_2A L1_3A L2_1A L2_2A L2_3A H_LX H_O
do
 
    # 使用 seqkit 重命名序列
    seqkit replace -p ".*" -r "${sample}_contig_{nr}" "${sample}_contigs.fasta" -o "${sample}_contigs.fna"
    echo "${sample} done"

done

rm *.fasta

cat *.fna > all_contigs.fna
