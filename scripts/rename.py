#!/usr/bin/env python3

import argparse
import os

def rename_fasta_headers(input_fasta, output_fasta, prefix, postfix, output_mapping):
    """
    重命名FASTA文件的描述行并生成映射文件
    
    参数:
    input_fasta: 输入FASTA文件路径 (必须设置)
    output_fasta: 输出FASTA文件路径 (必须设置)
    prefix: 新标题前缀，默认为输入文件名 (不包含路径和扩展名)
    postfix: 新标题后缀，默认为'seq'
    output_mapping: 输出映射文件路径, 默认为output_fasta文件名加上"_mapping.tsv"
    """
    # 初始化计数器
    count = 0
    
    
    # 打开输入文件和输出文件
    with open(input_fasta, 'r') as fin, \
            open(output_fasta, 'w') as fout_fasta, \
            open(output_mapping, 'w') as fout_map:
        
        # 写入映射文件表头
        fout_map.write("original_name\trenamed\n")
        
        for line in fin:
            # 如果是描述行（以'>'开头）
            if line.startswith('>'):
                count += 1
                original_header = line.strip()[1:]  # 移除'>'和换行符
                new_header = f"{prefix}_{postfix}{count}"
                
                # 写入新的FASTA文件
                fout_fasta.write(f">{new_header}\n")
                
                # 写入映射文件
                fout_map.write(f"{original_header}\t{new_header}\n")
            else:
                # 如果是序列行，直接写入
                fout_fasta.write(line)


def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='重命名FASTA文件描述行并生成映射文件')
    parser.add_argument('-i', '--input', required=True,
						help='输入FASTA文件路径')
    parser.add_argument('-o', '--output', required=True,
                       help='输出FASTA文件路径')
    parser.add_argument('-p', '--prefix', default=None,
                       help='新标题前缀，默认为输入文件名')
    parser.add_argument('-s', '--postfix', default='seq',
                       help='新标题后缀，默认为seq')
    parser.add_argument('-m', '--mapping', default=None,
                       help='映射文件路径，默认为output_fasta文件名加上"_mapping.tsv"')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        print(f"错误: 输入文件 '{args.input}' 不存在")
        return
    
    # 设置默认前缀为输入文件名(不含路径和扩展名)
    if args.prefix is None:
        args.prefix = os.path.splitext(os.path.basename(args.input))[0]
    
    # 设置默认映射文件名为输出文件名加上"_mapping.tsv"
    if args.mapping is None:
        args.mapping = os.path.splitext(args.output)[0] + "_mapping.tsv"
    
    # 执行重命名操作
    rename_fasta_headers(args.input, args.output, args.prefix, args.postfix, args.mapping)

if __name__ == "__main__":
    main()