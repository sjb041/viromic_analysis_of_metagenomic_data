################################## 使用03_quast_out_test/test9的contigs.fasta测试各个病毒识别工具或管道 ################################## 

################################## contigs_test9.fasta的来源：
# 提取 F1.1A 子集，100万条reads
#zcat F1.1A_L1_1.fq.gz | seqtk sample -s100 - 1000000 | gzip > R1.fastq.gz
#zcat F1.1A_L1_2.fq.gz | seqtk sample -s100 - 1000000 | gzip > R2.fastq.gz
# test9 测试
# ①使用kneaddata 0.7.4质控 + 去宿主，参数：
#- 接头：TruSeq3-PE
#- 滑动窗口：SLIDINGWINDOW:4:20
#- 保留最小长度：MINLEN:50
#--bowtie2-options '--end-to-end --sensitive --reorder'
#spades组装，再评估


################################## vs1
conda activate vs1

# 没有限制 contig 长度的参数，在输入前根据需要自行过滤掉过短的 contig
wrapper_phage_contigs_sorter_iPlant.pl -f contigs_test9.fasta --db 1 --wdir vs1 --ncpu 16 --data-dir $HOME/db/vs1_db/virsorter-data


################################## vs2
conda activate vs2

# 增加 hmmsearch 的线程数，默认是2
virsorter config --set HMMSEARCH_THREADS=4

# 设置数据库的位置
virsorter config --init-source --db-dir $HOME/db/vs2_db

# 使用所有分类器，最小长度3000，最小得分0.5
virsorter run -w vs2 -i contigs_test9.fasta --include-groups "dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae" --min-length 3000 --min-score 0.5 --label l3000s0.5 all
