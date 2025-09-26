#!/bin/bash

##################################################
# Prodigal + DIAMOND + [Virus-Host DB(release207,2021-10-07) + crAss-like噬菌体蛋白 + Benler研究的蛋白目录]
# 处理任意病毒序列文件

# 参考文献: A multi-kingdom collection of 33,804 reference genomes for the human vaginal microbiome

##################################################


# 检查输入参数数量，必须为3个参数
if [ "$#" -ne "5" ]; then
    echo -e "\nusage: sh $0 <path/proteins.faa> <path/output/prefix> <path/db/db.faa> <path/db/db.tsv> <threads>\n" 
    exit 2
fi


# 定义输入输出文件和参数
prot=$1                    # 输入蛋白质序列文件(.faa)
outf=$2					   # 输出文件前缀
db=$3                      # 蛋白质序列数据库
db_tax=$4                  # 蛋白质序列数据库分类信息    
threads=$5                 # 使用的线程数

# 添加依赖脚本的父目录到PATH
export PATH=/home/shijiabin/2025_2ME/bin/VMGCPipelines:$PATH

#### 基于蛋白序列的科水平注释

# 使用diamond进行blastp比对，比对蛋白质序列到蛋白质数据库
diamond blastp --threads $threads --max-target-seqs 10 --db $db --query $prot --outfmt 6 --out $outf.pep.bt --quiet

# 整理diamond比对结果中的TargetID
# 如果以"fromBenler_"开头,则去掉这个前缀然后再去掉第一个""及之后的部分; 如果以"virushostdb"开头,则去掉这个前缀,然后再去掉去掉第一个"|"及之后的部分; 如果以"crasslike_"开头,则去掉这个前缀
awk 'BEGIN{OFS="\t"} {
    if ($2 ~ /^fromBenler_/) {
        sub(/^fromBenler_/, "", $2);
        gsub(/_.*/, "", $2);
    } else if ($2 ~ /^virushostdb_/) {
        sub(/^virushostdb_/, "", $2);
        gsub(/\|.*$/, "", $2);
    } else if ($2 ~ /^crasslike_/) {
        sub(/^crasslike_/, "", $2);
    }
    print $0
}' $outf.pep.bt > $outf.pep.bt1

# 获取蛋白质序列长度信息
get_pep.length.pl $prot $outf.pep.len

# 过滤blast结果，设置查询覆盖度50%，身份相似度30%，tops 40，score 50等参数
filter_blast -i $outf.pep.bt1 -o $outf.pep.bt2 --qfile $outf.pep.len --qper 50 --identity 30 --tops 40 --score 50 

# 提取蛋白质序列ID并统计每个基因的蛋白质数量
grep '^>' $prot | sed 's/^>//' |perl -pne 's/_(\d+) #.*//' |sort |uniq -c |awk '{print $2"\t"$1}' > $outf.ngenes

# 使用vctg_stat.v2.pl脚本统计病毒contig信息
vctg_stat.v2.pl $outf.pep.bt2 $db_tax $outf.ngenes $outf.pep.bt2.f

# 过滤掉分类为NA的结果
awk '$2!="NA"' $outf.pep.bt2.f > $outf.pep.bt2.f2

# 按照特定规则排序并过滤结果，生成科水平分类文件
sort -k1,1 -k3,3nr -k6,6nr $outf.pep.bt2.f2 | perl -ne 'chomp;@s=split /\s+/;next if exists $h{$s[0]}; $pct=$s[2]/$s[4]; next unless ($pct>=0.2 or $s[2]>=10) and $s[5]>=30; printf "$s[0]\t$s[2]/$s[4]\t%.2f\t$s[1]\n",$s[5];$h{$s[0]}=1;' > $outf.tax_family
