用自己的数据，复现参考文献的流程, 相当于手动版的virify

参考文献：[*Massive expansion of human gut bacteriophage diversity - ScienceDirect*](https://www.sciencedirect.com/science/article/pii/S0092867421000726)

Viral taxonomic assignment of contigs was performed using a custom database of phylogenetically informative profile HMMs (ViPhOG v1, available here: ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/viral-pipeline/hmmer_databases), where each model is specific to one viral taxon. We used ‘hmmscan’ from HMMER v3.1b2 (Eddy, 1998) to query each protein sequence against the ViPhOG database, setting a full-sequence E-value reporting threshold of 10−3 and a per-domain independent E-value threshold of 0.1. Resulting hits were analyzed to predict the most likely and specific taxon for the whole contig based on the following criteria: (i) a minimum of 20% of genes with hits against the ViPhOG database, or at least two genes if the contig had less than 10 total genes; and (ii) among those with hits against the ViPhOG database, a minimum of 60% assigned to the same viral taxon.

#### 流程步骤

- 步骤1: 预测ORF (`prodigal`, `-p meta`) 文献提示**输入uvigs** **,**实际输入**uvigs and votus**

```shell
# 预测orf
prodigal -i {input.renamed_fasta} -a {output.protein} -p meta 

# 统计每条序列的orf数
grep ">" {output.protein} | cut -d" " -f1 | sed 's/>//' | sed 's/_[0-9]*$//' | sort | uniq -c | awk '{print $2 "\\t" $1}' > {output.count_orf}
```

- 步骤2: 使用hmmscan将蛋白比对到hmm数据库 (`-E 0.001 --domE 0.1`)

```shell
# hmmdb 是 vpHMM_database_v3.hmm
# We used ‘hmmscan’ from HMMER v3.1b2 (Eddy, 1998) to query each protein sequence against the ViPhOG database, setting a full-sequence E-value reporting threshold of 10−3 and a per-domain independent E-value threshold of 0.1.
hmmscan -E 0.001 --domE 0.1 --tblout {output.tblout} {hmmdb} {input.protein}
```

- 步骤3: 手动版_根据比对结果注释分类信息 (使用`additional_data_vpHMMs_v4.tsv`)

```
# Resulting hits were analyzed to predict the most likely and specific taxon for the whole contig based on the following criteria: (i) a minimum of 20% of genes with hits against the ViPhOG database, or at least two genes if the contig had less than 10 total genes; and (ii) among those with hits against the ViPhOG database, a minimum of 60% assigned to the same viral taxon.
```

① 使用去重后的比对结果, 计算是否满足0.2

② 若满足0.2, 使用不去重的比对结果注释每一个匹配, 再计算是否满足0.6

#### 运行此流程

`snakemake -s pipeline.smk --cores 5 --use-conda`

#### 输出示例

结果很差,估计是没能复现