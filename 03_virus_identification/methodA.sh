#########################################
# VirSorter2 + CheckV + HVR
# 筛选对象是 checkv 评估并修剪出的所有病毒
# 去除假阳性 HVR > 0.386 & host genes > 4
#########################################

# 合并所有样本的 contigs, 共 4747613 条 contig
cat /home/shijiabin/2025_2ME/02_assembly/02_contigs/rename_contigs/*_contigs.fna > all_contigs.fna && grep -c '^>' all_contigs.fna

# 计算 >=4kb 的 contig 数量, 共 112796 条 contig
seqkit seq -m 4000 all_contigs.fna -n | wc -l

######################################### step1 vs2 识别病毒, 筛选宽松阈值结果, 使用 --keep-original-seq 参数
mkdir -p 01_vs2 && cd 01_vs2
#conda run -n vs2 virsorter run -i ../all_contigs.fna -w vs2_keep-l4k-s0.5 --include-groups dsDNAphage,ssDNA --min-length 4000 --min-score 0.5 -j 30 all

# 评分表 7647
mv vs2_keep-l4k-s0.5/final-viral-score.tsv viral_score_vs2.tsv && tail -n +2 viral_score_vs2.tsv | wc -l

# 总病毒数 7647
mv vs2_keep-l4k-s0.5/final-viral-combined.fa viral_vs2.fna && grep -c '^>' viral_vs2.fna

# 前噬菌体数 338
grep -c "^>.*_partial$" viral_vs2.fna

cd ..

######################################### step2 checkv 评估修剪
dir_checkv="02_checkv"
#conda run -n checkv checkv end_to_end 01_vs2/viral_vs2.fna $dir_checkv -t 30
#real    14m58.986s
#user    309m24.572s
#sys     4m39.187s

# 合并 checkv 的前噬菌体和非前噬菌体序列文件 7649
cat $dir_checkv/proviruses.fna $dir_checkv/viruses.fna > $dir_checkv/viruses_ckv.fna && grep -c '^>' $dir_checkv/viruses_ckv.fna

# 评估的前噬菌体数 895
grep -w 'Yes' $dir_checkv/quality_summary.tsv | wc -l

# 评估的普通噬菌体数 6752
grep -w 'No'  $dir_checkv/quality_summary.tsv | wc -l

# 修剪后得到的前噬菌体数 897
grep -c '^>' $dir_checkv/proviruses.fna

######################################### step3 筛选病毒, 从 step2 中的序列筛选, 根据 step1+step2 的评分表
# 序列满足以下任意条件的保留
# 1）checkv预测至少含有一个病毒基因。
# 2）checkv预测不含宿主基因。
# 3）vs2评分 ≥ 0.95。
# 4）vs2标记基因 ≥ 3。
mkdir 03_filter1 && cd 03_filter1

# 满足条件1或2, 5085
awk -F'\t' 'NR==1 || ($6>=1 || $7==0)' ../$dir_checkv/quality_summary.tsv > score_ckv.tsv && tail -n +2 ckv_score.tsv | wc -l

# 满足条件3或4, 2838
awk -F'\t' 'NR==1 || (($4!="nan" && $4>=0.95) || $7>=3)' ../01_vs2/viral_score_vs2.tsv > score_vs2.tsv && tail -n +2 score_vs2.tsv | wc -l

# 去重, 5214（checkv评估的是vs2的输出序列，所以check评分表的ID与vs2评分表的ID相同，可直接去重）
{ tail -n +2 score_ckv.tsv; tail -n +2 score_vs2.tsv; } | cut -f1 | sort | uniq > ids.txt && wc -l ids.txt
# 花括号 { ... ; }，shell中的“命令组”，可以把多条命令组合在一起，作为一个整体输出。

# 从病毒2中提取序列得到病毒3, 5216
# 因为 checkv 修剪后的序列 ID 与评估时的 ID 不同，所以需要正则表达式匹配
# 将 viruses3.tsv 改成正则表达式,使用 -r 正则提取
sed 's/|/\\|/g; s/^/^/' ids.txt | seqkit grep -r -f - ../$dir_checkv/viruses_ckv.fna > viruses_filtered1.fna && grep -c '^>' viruses_filtered1.fna

# 提取出病毒3的checkv评分, 5214
# 获取标题行
head -n 1 ../$dir_checkv/quality_summary.tsv > quality_summary_filtered1.tsv
# 获取评分行
grep -F -w -f ids.txt ../$dir_checkv/quality_summary.tsv >> quality_summary_filtered1.tsv && tail -n +2 quality_summary_filtered1.tsv | wc -l
# -F：固定字符串匹配（不用正则，速度快）-w：整词匹配（防止部分ID误匹配）-f：指定ID列表文件

cd ..

######################################### step4 根据 HVR 值去除假阳性, HVR > 0.386 & host genes > 4
mkdir 04_filter2 && cd 04_filter2

# 计算 HVR & host genes 并筛选出假阳性的 ids，679
cut -f 1,6,7 ../03_filter1/quality_summary_filtered1.tsv | awk 'BEGIN{OFS="\t"} NR==1{print $0, "HVR"} NR>1{if($2!=0){print $0, $3/$2}else{print $0, "NA"}}' | awk 'NR==1{print; next} ($3>4 && ($4=="NA" || $4>0.386))' > false_ids.txt && tail -n +2 false_ids.txt | wc -l

# 移除假阳性序列, 4535
# 因为 checkv 修剪后的序列 ID 与评估时的 ID 不同，viruses3中的序列是checkv修剪后的，所以需要正则表达式匹配
tail -n +2 false_ids.txt | cut -f 1 | sed 's/|/\\|/g; s/^/^/' | seqkit grep -r -v -f - ../03_filter1/viruses_filtered1.fna > viruses_filtered2.fna && grep -c '^>' viruses_filtered2.fna

# 提取出病毒4的checkv评分
# 因为是反选，所以不需要单独获取标题行, 4535
cut -f1 false_ids.txt | sed '1d' | grep -v -F -w -f - ../03_filter1/quality_summary_filtered1.tsv > quality_summary_filtered2.tsv && tail -n +2 quality_summary_filtered2.tsv | wc -l
# -F：固定字符串匹配（不用正则，速度快）-w：整词匹配（防止部分ID误匹配）-f：指定ID列表文件 -v：反选

cd ..

######################################### step5 统一重命名预测出的所有病毒 viruses_filtered2.fna, 后续的所有分类注释方案均使用重命名后的序列
### 分别提取出 viruses_filtered2.fna 中各样本的序列, 然后重命名, 最后再合并所有样本到一个文件中
mkdir 05_rename && cd 05_rename
sample="F1_1A F1_2A F1_3A F2_1A F2_2A F2_3A FG_1A FG_2A FG_3A L1_1A L1_2A L1_3A L2_1A L2_2A L2_3A H_LX H_O"

# 按样本拆分评分表和序列文件
for i in $sample
do
    # 提取序列
    seqkit grep -r -p "^$i" ../04_filter2/viruses_filtered2.fna > ${i}_viruses.fna

    # 重命名
    rename.py -i ${i}_viruses.fna -o ${i}_viruses_rename.fna -p ${i} -s virus -m ${i}_mapping.tsv
done

# 合并, 然后查看序列数是否是 4535
cat *_rename.fna > ../all_viruses.fna && grep -c '^>' ../all_viruses.fna

cd ..

######################################### step6 聚类得到 votu
dir_cluster="06_cluster"

cluster.sh -i all_viruses.fna -t all_viruses.fna -a 95 -o $dir_cluster

# 2821 条 votu 代表序列
cp $dir_cluster/votus.fna votus.fna && grep -c '^>' votus.fna
