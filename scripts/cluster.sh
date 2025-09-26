#!/usr/bin/env bash

# 参数解析
usage() {
    echo "Usage: $0 -i <query_sequence> -t <target_sequence> -a <min_ani_value> -o <output_directory> " >&2
    exit 1
}

# 解析命令行参数
while getopts ":i:t:a:o:" opt; do
    case $opt in
        i) seq="$OPTARG" ;;
        t) target="$OPTARG" ;;
        a) min_ani="$OPTARG" ;;
        o) out="$OPTARG" ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# 检查必要参数
if [ -z "$seq" ] || [ -z "$target" ] || [ -z "$min_ani" ] || [ -z "$out" ] ; then
    echo "Error: Missing required parameters!" >&2
    echo "All of -i, -t, -a and -o must be specified" >&2
    usage
fi

# 建库
makeblastdb -in ${target} -dbtype nucl -out ${out}/db

# 比对
blastn -query ${seq} -db ${out}/db -outfmt '6 std qlen slen' -max_target_seqs 10000 -out ${out}/blast.tsv -num_threads 20

# 计算 ANI
anicalc.py -i ${out}/blast.tsv -o ${out}/ani.tsv

# 聚类
aniclust.py --fna ${seq} --ani ${out}/ani.tsv --out ${out}/clusters.tsv --min_ani ${min_ani} --min_tcov 85 --min_qcov 0

# 根据 cluster.tsv 获取 votu 代表序列
if [[ "$(realpath "$seq")" == "$(realpath "$target")" ]]; then
  # 提取序列
  cut -f1 ${out}/clusters.tsv | seqkit grep -f - ${seq} > ${out}/votus.fna

  # 获取 TSV 中行数
  tsv_lines=$(wc -l < "${out}/clusters.tsv")
  # 获取 FASTA 中序列数量
  fa_records=$(grep -c "^>" "${out}/votus.fna")
  
  # 判断数量是否一致
  if [ "$tsv_lines" -eq "$fa_records" ]; then
    echo "All sequences were successfully extracted."
  else
    echo "Error: Number of sequences in TSV and FASTA do not match!" >&2
    exit 1
  fi
fi

