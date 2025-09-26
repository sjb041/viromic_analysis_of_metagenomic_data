#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

_table_headers = [
    "target name",
    "target accession",
    "query name",
    "query accession",
    "full sequence E-value",
    "full sequence score",
    "full sequence bias",
    "orf"
]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Format hmmer domain hits table.")
    parser.add_argument("-i", dest="input_table",
                        help="hmmer domain hits table")
    parser.add_argument("-o", "--output_table",
                        dest="output_table", help="Output table name")
    parser.add_argument("-c", "--count_orf",
                        dest="count_orf", help="ORF count file")
    args = parser.parse_args()

    domain_table = Path(args.input_table)
    if not domain_table.is_file():
        raise Exception(
            "Input domain hits table missing. Path: " + args.input_table)

    # 读取ORF计数文件，建立query accession到ORF数量的映射
    orf_counts = {}
    with open(args.count_orf, "r") as count_file:
        for line in count_file:
            parts = line.strip().split()
            if len(parts) == 2:
                query_accession, orf_count = parts
                orf_counts[query_accession] = orf_count

    # 存储已处理的query name
    processed_queries = set()

    with open(args.output_table, "w", newline="") as out_table:
        tsv_writer = csv.writer(out_table, delimiter="\t",
                                quoting=csv.QUOTE_MINIMAL)
        tsv_writer.writerow(_table_headers)
        with open(domain_table, mode="r") as dt_reader:
            for line in dt_reader:
                if line.startswith("#"):
                    continue
                cols = line.split()
                data = cols[:7]
                
                # 填充 target accession 列（第2列）
                target_name = data[0]
                data[1] = target_name[6:-4]  # 去掉 "ViPhOG" 前缀和 ".faa" 后缀

                # 填充 query accession 列（第4列）
                query_name = data[2]
                query_accession = query_name.split("_")[0]  # 取下划线之前的部分
                data[3] = query_accession

                # 添加orf数量列
                if query_accession in orf_counts:
                    data.append(orf_counts[query_accession])
                else:
                    data.append("0")  # 如果找不到对应的ORF数量，默认为0

                # 只保留query name列中相同值的第一行
                if query_name not in processed_queries:
                    processed_queries.add(query_name)
                    tsv_writer.writerow(data)