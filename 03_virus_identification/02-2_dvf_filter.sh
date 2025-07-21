# 可根据分数阈值和长度阈值筛选

mkdir 02_dvf_l10k-s0.9-p0.05 && cd 02_dvf_l10k-s0.9-p0.05

awk 'NR==1 || ($3 > 0.9 && $4 < 0.05 && $2 >= 10000 )' ../02_run_dvf/all_gt1500bp.txt > all_gt10000bp_s0.9_p0.05.txt

# 提取序列
tail -n +2 all_gt10000bp_s0.9_p0.05.txt | cut -f 1 | seqtk subseq ../all_contigs.fna - > dvf.fna && grep -c '^>' dvf.fna

