#!/usr/bin/env python3

"""
Date: 2026-01-06
Description:
    从 vOTU 层级的相对丰度表和分类注释表计算 family 层级的相对丰度，
    并绘制 family 层级相对丰度的热图。
"""

import os
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# 汇总科水平的相对丰度
def calculate_family_abundance(votu_abundance_file, taxonomy_file, family_file):
    """
    按 family 层级汇总相对丰度
    1. 将vOTU相对丰度表与分类注释表进行合并
    2. 没有Family分类信息的记录填 'unknown'
    3. 按Family分组并对相对丰度进行汇总
    """
    votu_df = pd.read_csv(votu_abundance_file, sep='\t')
    taxonomy_df = pd.read_csv(taxonomy_file, sep='\t')

    # 合并
    merged = pd.merge(votu_df, taxonomy_df, left_on=votu_df.columns[0], right_on=taxonomy_df.columns[0], how="left")

    # 统一列名
    if 'family' in merged.columns and 'Family' not in merged.columns:
        merged = merged.rename(columns={'family': 'Family'})

    # 将缺失 Family 填充为 'Unknown'
    merged['Family'] = merged['Family'].fillna('Unknown')

    # 将 unclassified / Unclassified / 'unknown' 等统一改为 'Unknown'
    merged['Family'] = merged['Family'].astype(str).str.strip()  # 去掉前后空格
    merged['Family'] = merged['Family'].replace(
        to_replace=['unclassified', 'Unclassified', 'unknown', "'unknown'"],
        value='Unknown'
    )

    #merged.to_csv(merged_file, sep='\t', index=False)

    # 按 Family 汇总相对丰度
    rel_cols = [c for c in merged.columns if c.startswith("rel_abun_")]
    family_df = merged.groupby("Family")[rel_cols].sum().reset_index()
    family_df.to_csv(family_file, sep='\t', index=False)

    return family_df

# 作图
def plot_family_heatmap(family_df, output_png):
    """
    绘制 family 层级的相对丰度热图
    """
    sns.set_style("whitegrid")
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    # 排除 Family 为 'Unknown' 或者以 'uc_' 开头的
    family_df = family_df[
        (family_df['Family'] != 'Unknown') &
        (~family_df['Family'].str.startswith('uc_'))
    ]
    heatmap_data = family_df.set_index('Family')
    samples = [c.replace('rel_abun_', '') for c in heatmap_data.columns]
    heatmap_data.columns = samples

    plt.figure(figsize=(20, 12))
    sns.heatmap(heatmap_data, cmap='YlOrRd', cbar_kws={'label': 'Relative Abundance'})
    plt.title('Family Level Relative Abundance Heatmap', fontsize=14, fontweight='bold')
    plt.xlabel('Samples')
    plt.ylabel('Family')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')


def main():
    parser = argparse.ArgumentParser(description="计算并绘制 family 层级的相对丰度热图")
    parser.add_argument("--abun", help="vOTU 相对丰度文件路径")
    parser.add_argument("--tax", help="vOTU 分类信息文件路径")
    parser.add_argument("--family_abun", help="family 层级相对丰度结果输出路径")
    parser.add_argument("--png", help="热图输出文件路径")

    args = parser.parse_args()

    # 确保输出目录存在
    for path in [args.family_abun, args.png]:
        out_dir = os.path.dirname(os.path.abspath(path))
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

    # 汇总科水平的相对丰度
    family_df = calculate_family_abundance(args.abun, args.tax, args.family_abun)

    # 作图
    plot_family_heatmap(family_df, args.png)

if __name__ == "__main__":
    main()
