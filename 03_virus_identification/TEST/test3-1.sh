# 测试 test3-1
# SOP修：VirSorter2+CheckV+HVR
# 筛选对象是vs2识别出的所有病毒

#------------------------- ① all_contigs → vs2宽松条件 → 评分表 vs2_1，病毒1：【前噬菌体338，总噬菌体7647】 -------------------------------------------#
mkdir test3-1 && cd test3-1

# vs2预测
conda activate vs2
time virsorter run -i ../all_contigs.fna -w 01_vs2_keep-l4k-s0.5 --keep-original-seq --include-groups dsDNAphage,ssDNA --min-length 4000 --min-score 0.5 -j 30 all

# 总病毒数 7647
grep -c '^>' 01_vs2_keep-l4k-s0.5/final-viral-combined.fa

# 前噬菌体数 338
grep -c "^>.*_partial$" 01_vs2_keep-l4k-s0.5/final-viral-combined.fa



#--------------------------- ② 病毒1 → checkv → 评分表cv_1：【评估结果：前噬菌体895，非前噬菌体6752】--------------------------------------------------#
conda activate checkv

time checkv end_to_end 01_vs2_keep-l4k-s0.5/final-viral-combined.fa 02_cv_vs2_keep-l4k-s0.5 -t 30
#real    14m58.986s
#user    309m24.572s
#sys     4m39.187s

# 评估的前噬菌体数 895
grep -w 'Yes' 02_cv_vs2_keep-l4k-s0.5/quality_summary.tsv | wc -l

# 评估的普通噬菌体数 6752
grep -w 'No'  02_cv_vs2_keep-l4k-s0.5/quality_summary.tsv | wc -l

# 修剪后得到的前噬菌体数 897 (不需要关注，因为不使用修剪)
grep -c '^>' 02_cv_vs2_keep-l4k-s0.5/proviruses.fna 



#--------------------------- ③ 病毒1 → 筛选（评分表vs2_1 + cv_1）→ 病毒2，评分表cv_2：【7647 → 5214】--------------------------------------------------#
mkdir 03_filter1 && cd 03_filter1

# 满足条件1或2, 5085
awk -F'\t' 'NR==1 || ($6>=1 || $7==0)' ../02_cv_vs2_keep-l4k-s0.5/quality_summary.tsv > viruses2_cv.tsv

# 满足条件3或4, 2838
awk -F'\t' 'NR==1 || (($4!="nan" && $4>=0.95) || $7>=3)' ../01_vs2_keep-l4k-s0.5/final-viral-score.tsv > viruses2_vs2.tsv

# 去重, 5214（checkv评估的是vs2的输出序列，所以check评分表的ID与vs2评分表的ID相同，可直接去重）
{ tail -n +2 viruses2_cv.tsv; tail -n +2 viruses2_vs2.tsv; } | cut -f1 | sort | uniq > viruses2.tsv
# 花括号 { ... ; }，shell中的“命令组”，可以把多条命令组合在一起，作为一个整体输出。

# 从病毒1中提取序列得到病毒2, 5214
seqtk subseq ../01_vs2_keep-l4k-s0.5/final-viral-combined.fa viruses2.tsv > viruses2.fna && grep -c '^>' viruses2.fna

# 提取出病毒2的checkv评分
# 获取标题行
head -n 1 ../02_cv_vs2_keep-l4k-s0.5/quality_summary.tsv > viruses2_quality_summary.tsv
# 获取评分行
grep -F -w -f viruses2.tsv ../02_cv_vs2_keep-l4k-s0.5/quality_summary.tsv >> viruses2_quality_summary.tsv
# -F：固定字符串匹配（不用正则，速度快）-w：整词匹配（防止部分ID误匹配）-f：指定ID列表文件

cd ..



#--------------------------- ④ 病毒2 → 去除假阳性（评分表cv_2， HVR）→ 病毒3，评分表cv_3：【5214 → 4584】----------------------------------------------#

#--------------------------------- 选取画图用的数据：88 complete putative viral contigs，53 complete vOTUs
mkdir 04_filter2 && cd 04_filter2

# Firstly, putative viral contigs predicted to be complete were selected.
# 从病毒2中选出那些完整的
awk -F'\t' 'NR==1 || $8=="Complete"' ../03_filter1/viruses2_quality_summary.tsv > complete.tsv
cut -f1 complete.tsv | sed '1d' | seqtk subseq ../03_filter1/viruses2.fna - > complete.fna && grep -c '^>' complete.fna
# 88

# These contigs were dereplicated into vOTUs (viral Operational Taxonomic Units; 95% Average Nucleotide Identity over 85% Alignment Fraction)
conda activate checkv
# 创建数据库
makeblastdb -in complete.fna -dbtype nucl -out blast_db/complete

# 比对
blastn -query complete.fna -db blast_db/complete -outfmt '6 std qlen slen' -max_target_seqs 10000 -out blast_out.tsv -num_threads 30

# 计算 ANI (anicalc.py已在$HOME/opt/中)
anicalc.py -i blast_out.tsv -o ani.tsv

# 使用推荐的 MIUVIG 参数（95% ANI + 85% AF）进行类似 CD-HIT 的聚类
aniclust.py --fna complete.fna --ani ani.tsv --out clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0


#---------------------------------- 画图select_HVR.Rmd
---
title: "select_HVR"
output: html_document
---

```{r}
library(ggplot2)
library(seqinr)
library(dplyr)

getwd()

################# 输入
# 聚类表
clusters_path <- "clusters.tsv"    
# complete序列
cpv_path <- "complete.fna"
# checkv评分表
quality_summary_path <- "../03_filter1/viruses2_quality_summary.tsv"


################# 输出
# votus
vOTUs_path <- "vOTUs.fa"
# 看图选阈值

# 定义函数，取fasta序列的ID
ids <- function(fasta) {
  # names 是取整个描述行，这里截取第一个空格前的部分（ID）
  return (sub(" .*", "", names(fasta)))
}
```

```{r}
# ------------------------------- 提取 vOTUs.fa ----------------------------
# 读入vOTUs的表clusters.tsv
clusters <- read.delim(file = clusters_path, header = F, stringsAsFactors = F, check.names = F)

# 读入complete.fa
cpv <- read.fasta(file = cpv_path,
                    as.string = TRUE,
                    forceDNAtolower = FALSE,
                    whole.header = TRUE)
                    
# 从complete.fa中提取出vOTUs 
vOTUs_fasta <- cpv[names(cpv) %in% clusters[,1]]   # 不能使用 clusters[1] 因为它是一个数据框

# 判断是否能全部提取出clusters里的序列
if (dim(clusters)[1] != length(vOTUs_fasta)){
    stop("不正确的筛选！")
}
message("Number complete vOTUs:",length(vOTUs_fasta))

# 输出vOTUs 
write.fasta(sequences = vOTUs_fasta, names = names(vOTUs_fasta), file.out = vOTUs_path)
```

```{r}
# --------------------------------计算 HVR ----------------------------------------
# 读入 checkv 的评分表
df1 <- read.delim(file = quality_summary_path, header = T, stringsAsFactors = F, check.names = F)

# 筛选评分表，提取出 vOTUs 的评分
df1 <- dplyr::filter(df1, contig_id %in% names(vOTUs_fasta))

# 判断是否每一个vOTUs都找到了自己的评分表
if (dim(df1)[1] != length(vOTUs_fasta)){
    stop("存在某个vOTU的评分表未找到！")
}


# 计算 HVR 值
df1 <- mutate(df1,
    HVR = case_when(
      host_genes != 0 & viral_genes == 0 ~ Inf,  # 正无穷大
      host_genes == 0 & viral_genes != 0 ~ NA,   # 缺失值
      host_genes == 0 & viral_genes == 0 ~ NaN,  # 无效的运算
      TRUE ~ host_genes / viral_genes            # 否则，HVR 为 host_genes / viral_genes
    )
)

# 提取出有用的列
df2 <- df1[, c("contig_id", "provirus", "host_genes", "HVR")]

# 找到 HVR 的正常最大值
max_hvr <- max(df2$HVR, na.rm = TRUE, finite = TRUE)
# 将 Inf 设置为最大 HVR
df2$HVR <- ifelse(is.infinite(df2$HVR), max_hvr, df2$HVR)

# 找到 HVR 的正常最小值
min_hvr <- min(df2$HVR, na.rm = TRUE, finite = TRUE)
# 将 NA 和 NaN 设置为最小 HVR
df2$HVR <- ifelse(is.na(df2$HVR), min_hvr, df2$HVR)
    
# 将 df2 分为前噬菌体和非前噬菌体两组
pro_vOTUs <- dplyr::filter(df2, provirus == "Yes")
non_pro_vOTUs <- dplyr::filter(df2, provirus == "No")
# 确保分组后的数量
if ((dim(pro_vOTUs)[1] + dim(non_pro_vOTUs)[1]) != length(vOTUs_fasta)){
    stop("pro + non_pro != df2")
}
```

```{r}
# 创建绘图 a)
p <- ggplot(df2, aes(x = provirus, y = log10(HVR))) +
  #geom_boxplot(aes(fill = provirus), outlier.shape = NA, alpha = 0.5) + # 绘制箱形图 
  geom_jitter(width = 0.2, size = 2, alpha = 0.25) + # 散点图，增加点的随机抖动，width = 0.2控制抖动的宽度，size = 2设置点的大小，alpha = 0.25设置点的透明度
  # 不显式设置 limits，让 ggplot 自动调整
  # 假设 df2$HVR 的值为 c(1, 10, 100, NA, 1000)，经过 log10() 转换后，得到的值为 c(0, 1, 2, NA, 3)。计算范围后，得到 c(0, 3)，然后 pretty(c(0, 3)) 可能返回 c(0, 1, 2, 3)。
  scale_y_continuous(breaks = pretty(range(log10(df2$HVR), na.rm = TRUE))) +   
  labs(x = "Provirus", y = expression(paste("Log"[10], "(HVR)"))) +
  theme_bw() +   # 设置主题，使用黑白主题
  theme(         # 自定义主题元素
      panel.grid.major = element_blank(),   # 不需要网格线
      panel.grid.minor = element_blank(),   # 不需要网格线
      panel.border = element_blank(),       # 移除面板边界
      axis.text = element_text(size = 12),  # 调整坐标轴刻度文本的大小
      axis.title = element_text(size = 14), # 调整坐标轴标题的大小
      axis.line = element_line(color = "black", linewidth = 1), # 设置坐标轴线条
      text = element_text(size = 12)
  )

# 使用 wilcox.test 进行检验
# 比较 log10(HVR) 在不同 provirus 组之间的差异。log10(HVR) ~ provirus 表示以 provirus 作为分组变量。
wilcox_result <- wilcox.test(log10(HVR) ~ provirus, data = df2)
wilcox_p_value <- format.pval(wilcox_result$p.value, digits = 3)   # 格式化 Wilcoxon 检验的 p 值，保留三位有效数字。

# 在图中标记 p 值
p <- p + 
  # 在图中添加文本注释。
  annotate("text", x = 1.5, y = max(log10(df2$HVR), na.rm = TRUE) - 0.2,    # 文本的放置位置
           label = paste0("Wilcoxon, p <", wilcox_p_value), size = 7)

# 展示图形
print(p)

```

```{r}
# 绘制对数密度图 b)
#p <- ggplot(non_pro_vOTUs, aes(x = log10(HVR))) +
p <- ggplot(non_pro_vOTUs, aes(x = (HVR))) +
  geom_density(color = "black", size = 0.8) +            # 密度图层
  geom_rug(sides = "b", color = "black", size = 0.8) +   # rug 图层，显示数据的分布情况
  labs(title = "Logarithmic Density Plot of HVR", x = "Log10 of HVR", y = "Density") +
  theme_minimal()   # 使用最小主题
p + #scale_x_continuous(breaks = seq(-2, 2, by = 0.5)) +
  #scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  theme(panel.grid.major = element_blank(), # 不需要网格线
        panel.grid.minor = element_blank(), # 不需要网格线
        axis.text = element_text(size = 12),  # 调整坐标轴刻度文本的大小
        axis.title = element_text(size = 14), # 调整坐标轴标题的大小
        axis.line = element_line(color = "black", linewidth = 1)  # 设置坐标轴线条
       )
#non_pro_vOTUs$HVR
#log1p(non_pro_vOTUs$HVR)

```

```{r}
# 绘制对数密度图 c)
#non_pro_vOTUs$host_genes
#log10(non_pro_vOTUs$host_genes)
# 创建一个新的列，用于表示 HVR 是否小于或等于 0.333 (log10 = -0.477)  0.231 (log10 = -0.637)  
non_pro_vOTUs$HVR_group <- ifelse(non_pro_vOTUs$HVR <= 0.5, "Yes", "No")

# 绘制对数密度图
pc <- ggplot(non_pro_vOTUs, aes(x = log10(host_genes), fill = HVR <= 0.5)) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = c("#E69F00", "#D5B4AD")) +  # 设置填充颜色
      labs(x = expression(Log[10]~"of the number of host genes"), y = "Density") +
      theme_minimal()

pc + 
  #scale_x_continuous(limits = c(0, 2), breaks = seq(0, 2, by = 0.5)) +  # 设置x轴范围并增加刻度
  #scale_y_continuous(limits = c(0, NULL), breaks = seq(0, 1, by = 0.25)) +  # 设置y轴下限为0，上限自动确定
  theme(panel.grid.major = element_blank(), # 不需要网格线
        panel.grid.minor = element_blank(), # 不需要网格线
        axis.text = element_text(size = 12),  # 调整坐标轴刻度文本的大小
        axis.title = element_text(size = 14), # 调整坐标轴标题的大小
        axis.line = element_line(color = "black", linewidth = 1)  # 设置坐标轴线条
       )
```


#--------------------------------- 根据阈值去除假阳性 HVR > 0.386 & host genes > 5
# 计算 HVR & host genes 并筛选出假阳性的 ids，630
cut -f 1,6,7 ../03_filter1/viruses2_quality_summary.tsv | awk 'BEGIN{OFS="\t"} NR==1{print $0, "HVR"} NR>1{if($2!=0){print $0, $3/$2}else{print $0, "NA"}}' | awk 'NR==1{print; next} ($3>5 && ($4=="NA" || $4>0.386))' > false_ids.txt

# 移除假阳性序列, 4584
tail -n +2 false_ids.txt | cut -f 1 | seqkit grep -n -v -f - ../03_filter1/viruses2.fna > viruses3.fna && grep -c '^>' viruses3.fna

# 提取出病毒3的checkv评分
# 因为是反选，所以不需要单独获取标题行
cut -f1 false_ids.txt | sed '1d' | grep -v -F -w -f - ../03_filter1/viruses2_quality_summary.tsv > viruses3_quality_summary.tsv
# -F：固定字符串匹配（不用正则，速度快）-w：整词匹配（防止部分ID误匹配）-f：指定ID列表文件 -v：反选

cd ..
