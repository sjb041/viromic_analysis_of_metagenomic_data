#!/usr/bin/env bash

######################################
# Date: 2026-01-09
# Description:
#     Cluster genomes to Genus-clusters based on AAI
#     https://github.com/snayfach/MGV/tree/master/aai_cluster
# Dependencies:
#     prodigal --version 2.6.3
#     diamond --version 2.0.6
#     python --version 2.7.15 --> numpy biopython
#     mcl --version 14-137
######################################

#### step1 预测 ORF
dir="02_genus_clusters" && mkdir -p $dir && cd $dir

votus="/home/shijiabin/2025_2ME/03_virus_identification/votus.fna"

prodigal -i $votus -a votu_orfs.faa -p meta

#### step2 DIAMOND 进行成对蛋白质序列比对
mkdir diamond

# 建库
diamond makedb --in votu_orfs.faa --db diamond/db --threads 10

# 执行 all-vs-all 比对（很快）
diamond blastp -q votu_orfs.faa -d diamond/db -e 1e-5 --max-target-seqs 99999  --query-cover 50 --subject-cover 50 --outfmt 6 -o votu_orfs_diamond_results.tsv
rm -r diamond

#### step3 计算每对vOTU之间共享基因的百分比和平均氨基酸相似性（AAI）
cat > amino_acid_identity.py << 'EOF'
#!/usr/bin/env python

import subprocess as sp, os, numpy as np, sys, Bio.SeqIO, csv, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--in_faa', required=True, metavar='PATH')
parser.add_argument('--in_blast', required=True, metavar='PATH')
parser.add_argument('--out_tsv', required=True, metavar='PATH')
args = vars(parser.parse_args())

from collections import defaultdict
genome_size = defaultdict(int)
gene_to_genome = {}
for index, r in enumerate(Bio.SeqIO.parse(args['in_faa'], 'fasta')):
	genome_id = r.id.rsplit('_', 1)[0]
	genome_size[genome_id] += 1
	gene_to_genome[r.id] = genome_id
genome_alns = {}
for genome_id in genome_size:
	genome_alns[genome_id] = {}

print "parse"
for index, line in enumerate(open(args['in_blast'])):
	
	aln = line.split()
	gene = aln[0]
	query = gene_to_genome[aln[0]]
	target = gene_to_genome[aln[1]]
	score = float(aln[-1])
	
	if target not in genome_alns[query]:
		genome_alns[query][target] = {}
	if gene not in genome_alns[query][target]:
		genome_alns[query][target][gene] = aln
	elif score > float(genome_alns[query][target][gene][-1]):
		genome_alns[query][target][gene] = aln

print "compute"
rows = []
for query in genome_alns:
	query_genes = genome_size[query]
	for target in genome_alns[query]:
		target_genes = genome_size[target]
		alns = genome_alns[query][target].values()
		shared_genes = len(alns)
		qcov = 100.0*shared_genes/query_genes
		tcov = 100.0*shared_genes/target_genes
		aai = np.mean([float(_[2]) for _ in alns])
		row = [query, target, query_genes, target_genes, shared_genes, qcov, tcov, aai]
		rows.append(row)

print "write"
with open(args['out_tsv'], 'w') as out:
	header = ['qname', 'tname', 'qgenes', 'tgenes', 'sgenes', 'qcov', 'tcov', 'aai']
	out.write('\t'.join(header)+'\n')
	for row in rows:
		out.write('\t'.join([str(_) for _ in row])+'\n')
EOF

python amino_acid_identity.py --in_faa votu_orfs.faa --in_blast votu_orfs_diamond_results.tsv --out_tsv votu_aai.tsv
rm votu_orfs.faa votu_orfs_diamond_results.tsv amino_acid_identity.py

#### step4 建立votu间的连接网络
cat > filter_aai.py << 'EOF'
import csv, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--in_aai', required=True, metavar='PATH')
parser.add_argument('--min_num_shared', type=int, required=True, metavar='INT')
parser.add_argument('--min_percent_shared', type=float, required=True, metavar='FLOAT')
parser.add_argument('--min_aai', type=float, required=True, metavar='FLOAT')
parser.add_argument('--out_tsv', required=True, metavar='PATH')
args = vars(parser.parse_args())

rows = []
for i, r in enumerate(csv.DictReader(open(args['in_aai']), delimiter='\t')):
	qcov, tcov = float(r['qcov']), float(r['tcov'])
	cov = min([qcov, tcov])
	aai = float(r['aai'])
	shared = int(r['sgenes'])
	score = cov/100.0 * aai/100.0
	if not (
			(cov>=args['min_percent_shared'] or shared>=args['min_num_shared'])
			and aai >= args['min_aai']):
		continue
	row = [r['qname'], r['tname'], str(score)]
	rows.append(row)

with open(args['out_tsv'], 'w') as out:
	for row in rows:
		out.write('\t'.join([str(_) for _ in row])+'\n')
EOF

python filter_aai.py --in_aai votu_aai.tsv --min_percent_shared 20 --min_num_shared 0 --min_aai 50 --out_tsv genus_edges.tsv
rm votu_aai.tsv filter_aai.py

#### step5 基于vOTU之间的连接使用MCL进行聚类
mcl genus_edges.tsv -te 8 -I 2.0 --abc -o genus_clusters.txt
rm genus_edges.tsv
