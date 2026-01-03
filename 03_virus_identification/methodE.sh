#########################################
# VirSorter2 + DeepVirFinder + CheckV/geNomad
#########################################

# 合并所有样本的 contigs, 共 4747613 条 contig
cat /home/shijiabin/2025_2ME/02_assembly/02_contigs/rename_contigs/*_contigs.fna > all_contigs.fna && grep -c '^>' all_contigs.fna

# 计算 >=4kb 的 contig 数量, 共 112796 条 contig
seqkit seq -m 4000 all_contigs.fna -n | wc -l

######################################### step1 vs2 识别病毒,筛选高保守阈值结果
mkdir -p 01_vs2 && cd 01_vs2
#conda run -n vs2 virsorter run -i ../all_contigs.fna -w vs2_l4k-s0.5 --include-groups dsDNAphage,ssDNA --min-length 4000 --min-score 0.5 -j 30 all

# 如果只是想使用--keep-original-seq，其它参数不变，可以根据得分表从all_contigs.fna中获取那些||full序列，而不用重新跑一次
# --keep-original-seq：保留原始序列而不是修剪后的序列；默认情况下，已识别的病毒序列两端的未翻译区域会被修剪；环状序列会被修改以去除两端的重叠部分，并调整因基因分裂成两端的情况
# --keep-original-seq不影响评分【已验证】，建议不使用--keep-original-seq，如果想获取未被修剪的full序列，可以根据得分表从all_contigs.fna中获取那些full序列。

# 从 final-viral-score.tsv 中取出 high-confidence-only
# 根据阈值筛选vs2评分表的行, 共 4118 行
awk -F'\t' 'NR==1 || ($4!="nan" && $4>=0.9) || ($4!="nan" && $4>=0.7 && $7>=1)' vs2_l4k-s0.5/final-viral-score.tsv > viral_score_vs2.tsv && tail -n +2 viral_score_vs2.tsv | wc -l

# 根据评分表中的 seqname 提取对应的 high-confidence-only 序列, 共 4118 条 virus
cut -f1 viral_score_vs2.tsv | sed '1d' | seqkit grep -n -f - vs2_l4k-s0.5/final-viral-combined.fa > viral_vs2.fna && grep -c '^>' viral_vs2.fna

cd ..

######################################### step1 dvf 识别病毒,筛选高保守阈值结果
mkdir -p 01_dvf && cd 01_dvf
#conda run -n DVF dvf.py -i ../all_contigs.fna -o ./ -l 4000 -c 15

# 根据分数阈值和长度阈值筛选, 共 6004 行
awk 'NR==1 || ($2 >= 4000 &&  $3 > 0.9 && $4 < 0.05)' all_contigs.fna_gt4000bp_dvfpred.txt > viral_score_dvf.txt && tail -n +2 viral_score_dvf.txt | wc -l
# 从 all_contigs.fna 中提取 dvf 识别的病毒, 共 6004 条
tail -n +2 viral_score_dvf.txt | cut -f 1 | seqtk subseq ../all_contigs.fna - > viral_dvf.fna && seqkit seq -m 1 viral_dvf.fna -n | wc -l

cd ..

######################################### step2 vs2+dvf 去重
mkdir -p 02_dereplicate && cd 02_dereplicate

# 获取 vs2 序列 ID，并去掉 || 及之后的部分，然后去重, 共 4115 个 id, 为何要去重:（F1_3A_contig_403||1_partial；F1_3A_contig_403||2_partial；FG_3A_contig_3||1_partial） 
cut -f1 ../01_vs2/viral_score_vs2.tsv | sed '1d' | sed 's/||.*//' | sort | uniq > vs2_ids.txt && cat vs2_ids.txt | wc -l
# 获取 dvf 序列 ID, 共 6004 个 id
cut -f1 ../01_dvf/viral_score_dvf.txt | sed '1d' > dvf_ids.txt && cat dvf_ids.txt | wc -l

# 取交集，然后匹配 vs2 中的序列
# 交集共有 1175 个 id
grep -Fxf vs2_ids.txt dvf_ids.txt > intersection.txt && wc -l intersection.txt
# 能够匹配到 vs2 中的 1176 条序列, 因为: (FG_3A_contig_3||1_partial)
sed 's/^/^/; s/$/\\|\\|/' intersection.txt | seqkit grep -r -f - ../01_vs2/viral_vs2.fna > intersection.fna && grep -c '^>' intersection.fna

# 在 vs2 中不在 dvf 中的
# 2940
grep -vFxf dvf_ids.txt vs2_ids.txt > vs2only.txt && wc -l vs2only.txt
# 2942（F1_3A_contig_403||1_partial；F1_3A_contig_403||2_partial）
sed 's/^/^/; s/$/\\|\\|/' vs2only.txt | seqkit grep -r -f - ../01_vs2/viral_vs2.fna > vs2only.fna && grep -c '^>' vs2only.fna

# 在 dvf 中不在 vs2 中的, 4829
grep -vFxf vs2_ids.txt dvf_ids.txt > dvfonly.txt && wc -l dvfonly.txt
seqkit grep -n -f dvfonly.txt ../01_dvf/viral_dvf.fna > dvfonly.fna && grep -c '^>' dvfonly.fna

# 合并三个集合，8947
cat intersection.fna vs2only.fna dvfonly.fna > viral_dereplicated.fna && grep -c '^>' viral_dereplicated.fna

# 验证id无重复
# 首先取出 viral_dereplicated.fna 的id, 并去掉 || 及之后的部分，去重之后看行数与 8947 对比，预期少3行就说明无重复。
seqkit seq -m 1 viral_dereplicated.fna -n | sed 's/||.*//' | sort | uniq > viral_dereplicated_ids_uniq.txt && wc -l viral_dereplicated_ids_uniq.txt

cd ..

######################################### step3 genomad 处理 viral_dereplicated.fna
dir_genomad="03_genomad"

conda run -n genomad genomad end-to-end 02_dereplicate/viral_dereplicated.fna \
    $dir_genomad \
    $HOME/db/genomad_db \
    --enable-score-calibration \
    --min-score 0.7 --threads 30

# 总病毒数 4200
grep -c '^>' $dir_genomad/viral_dereplicated_summary/viral_dereplicated_virus.fna
# 找到的前噬菌体数 589
grep -c '^>' $dir_genomad/viral_dereplicated_find_proviruses/viral_dereplicated_provirus.fna
# 最终纳入的前噬菌体数 550 （序列名含有provirus）
grep -c '^>.*provirus' $dir_genomad/viral_dereplicated_summary/viral_dereplicated_virus.fna

# 质粒数 582
grep -c '^>' $dir_genomad/viral_dereplicated_summary/viral_dereplicated_plasmid.fna


######################################### step4 checkv 处理 viral_dereplicated_virus.fna, 修剪全部病毒
dir_checkv="04_checkv"

conda run -n checkv checkv end_to_end \
    $dir_genomad/viral_dereplicated_summary/viral_dereplicated_virus.fna \
    $dir_checkv \
    -t 30

# 拼接前噬菌体和非前噬菌体文件, 最终得到 4201 条病毒
cat $dir_checkv/proviruses.fna $dir_checkv/viruses.fna > $dir_checkv/viruses_checkv.fna && grep -c '^>' $dir_checkv/viruses_checkv.fna

# 评估的前噬菌体数 366
grep -w 'Yes' $dir_checkv/quality_summary.tsv | wc -l

# 评估的普通噬菌体数 3834
grep -w 'No'  $dir_checkv/quality_summary.tsv | wc -l

# 修剪后得到的前噬菌体数 367
grep -c '^>' $dir_checkv/proviruses.fna

######################################### (可选) 从 viruses_checkv.fna 中选取 high-quality 和 medium_quality
### 首先使用 all_viruses.fna 跑一遍 checkv 评估, 以下使用这次评估的 quality_summary.tsv

#### 从 all_viruses.fna 中选取 high-quality
# 用 Excel 筛选要保留的部分,然后手动复制第一列contig_id, 289
#quality_summary.tsv ——→ high_quality.xlsx ——→ high_ids.txt

# 提取 high, 289
#seqkit grep -f high_ids.txt all_viruses.fna > all_viruses_hq.fna && grep -c '^>' all_viruses_hq.fna

#### 从 viruses_checkv.fna 中选取 medium_quality
# 用 Excel 筛选要保留的部分,然后手动复制第一列contig_id, 376
#quality_summary.tsv ——→ medium_quality.xlsx ——→ medium_ids.txt

# 提取 mediums, 376
#seqtk subseq all_viruses.fna medium_ids.txt > all_viruses_mq.fna && grep -c '^>' all_viruses_mq.fna

######################################### step5 统一重命名预测出的所有病毒 viruses_checkv.fna, 后续的所有分类注释方案均使用重命名后的序列
### 分别提取出 viruses_checkv.fna 中各样本的序列, 然后重命名, 最后再合并所有样本到一个文件中
mkdir 05_separate_samples && cd 05_separate_samples
sample="F1_1A F1_2A F1_3A F2_1A F2_2A F2_3A FG_1A FG_2A FG_3A L1_1A L1_2A L1_3A L2_1A L2_2A L2_3A H_LX H_O"

# 按样本拆分评分表和序列文件
for i in $sample
do
    # 提取序列
    seqkit grep -r -p "^$i" ../$dir_checkv/viruses_checkv.fna > ${i}_viruses.fna

    # 重命名
    rename.py -i ${i}_viruses.fna -o ${i}_viruses_rename.fna -p ${i} -s virus -m ${i}_mapping.tsv
done

# 合并, 然后查看序列数是否是 4201
cat *_rename.fna > ../all_viruses.fna && grep -c '^>' ../all_viruses.fna

cd ..
