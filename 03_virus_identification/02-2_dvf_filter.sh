#!/usr/bin/env bash

######################################## 根据阈值提取contigs ########################################
# chmod +x 02-2_dvf_filter.sh
# ./02-2_dvf_filter.sh -h

# 或者 bash 02-2_dvf_filter.sh

################## 参数检查
show_help() {
    echo "用法: $0 -f <序列文件> -t <dvf预测文件> -s <score> -p <pvalue> -o <输出目录>"
    echo ""
    echo "参数说明："
    echo "  -f    序列文件"
    echo "  -t    dvf预测文件"
    echo "  -s    score"
    echo "  -p    pvalue"
    echo "  -o    输出目录"
    echo "  -h    显示帮助信息"
    echo ""
    echo "示例:"
    echo "  $0  -f seq.fasta -t pred.txt -s 0.9 -p 0.05 -o ."
}

# 初始化变量
f_file=""
t_file=""
score=""
pvalue=""
out=""

# 解析参数
while getopts "f:t:s:p:o:h" opt; do
    case $opt in
        f) f_file=$OPTARG ;;
        t) t_file=$OPTARG ;;
        s) score=$OPTARG ;;
        p) pvalue=$OPTARG ;;
        o) out=$OPTARG ;;
        h)
            show_help
            exit 0
            ;;
        *)
            show_help
            exit 1
            ;;
    esac
done

# 检查必需参数
if [[ -z $f_file || -z $t_file || -z $score || -z $pvalue || -z $out ]]; then
    echo "错误：缺少必要参数！"
    show_help
    exit 1
fi


################## 筛选出符合的contigID
t_file_base=$(basename "$t_file")   # 移除父路径
t_file_base="${t_file_base%.*}"     # 移除扩展名

t_file_filter="${out}/${t_file_base}_s${score}_p${pvalue}.txt"
awk '$3 > '"$score"' && $4 < '"$pvalue"' {print $1}' $t_file > $t_file_filter

################## 提取contigs
f_file_base=$(basename "$f_file")
f_file_base="${f_file_base%.*}"

f_file_filter="${out}/${f_file_base}_s${score}_p${pvalue}.fna"
seqkit grep -f $t_file_filter $f_file > $f_file_filter

# 错误检查
expected=$(wc -l < $t_file_filter)
real=$(grep -c '^>' $f_file_filter)
if [ "$expected" -ne "$real" ]; then
    echo "${t_file_filter}中某些序列未提取"
    exit 1
fi
