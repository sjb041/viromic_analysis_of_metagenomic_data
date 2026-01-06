#!/usr/bin/env bash

######################################
# Date: 2026-01-06
# Description:
#     从 vOTU 层级的相对丰度表和分类注释表计算 family 层级的相对丰度，
#     并绘制 family 层级相对丰度的热图。
# Dependencies:
#     family_abun.py
######################################

method_viral_list=(A E)                                     # 病毒识别方法
method_abun_list=(method1 method2 method3 method4 method5)  # 相对丰度计算方法
method_tax_list=(1-1 1-2 1-3 2-1 3-1 3-2 4-1)               # 分类注释方法

for method_viral in "${method_viral_list[@]}"; do
    for method_abun in "${method_abun_list[@]}"; do
        for method_tax in "${method_tax_list[@]}"; do
            # 输入
            abun_file="/home/shijiabin/2025_2ME/06_relative_abundance/$method_viral/$method_abun/votus_relative_abundance.tsv"
            tax_file="/home/shijiabin/2025_2ME/05_taxonomic_classification/votus_taxonomy_${method_viral}_${method_tax}.tsv"
            
            # 输出文件放在输入文件父目录
            abun_dir=$(dirname "$abun_file")
            family_abun_file="$abun_dir/family_relative_abundance.${method_tax}.tsv"
            png_file="$abun_dir/heatmap.${method_tax}.png"
            
            # 调用处理脚本（示例）
            python family_abun.py --abun $abun_file --tax $tax_file --family_abun $family_abun_file --png $png_file

            # 打印处理信息
            echo "Processing viral=$method_viral, abun=$method_abun, tax=$method_tax"
            echo "Output file: $family_abun_file"
        done
    done
done
