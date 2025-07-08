## 使用 kneaddata 0.7.4 版本，miniconda3 安装
## 需要自行构建参考基因组 hg_38_p14，修改“db_path”到自己的数据库位置

## wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz
## gunzip GCA_000001405.29_GRCh38.p14_genomic.fna.gz
## conda activate kd0.7.4
## bowtie2-build GCA_000001405.29_GRCh38.p14_genomic.fna $HOME/db/kneaddata_db/hg_38_p14 --threads 16

## wget https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip
## unzip GRCm39.zip
############################################################################################

mkdir -p 02_quality_control_rawdata

# 动态获取Conda路径
CONDA_BASE=$(conda info --base)
# 加载Conda基础配置
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate kd0.7.4

# 样本名
sample="F1.1A F1.2A F1.3A F2.1A F2.2A F2.3A FG.1A FG.2A FG.3A L1.1A L1.2A L1.3A L2.1A L2.2A L2.3A H.LX H.O"

for i in $sample
do
	in1="00_rawdata/${i}/${i}_L1_1.fq.gz"
	in2="00_rawdata/${i}/${i}_L1_2.fq.gz"
	outdir="02_quality_control_rawdata/${i}"
	mkdir -p ${outdir}

    # H样本去除宿主 hg_38_p14，其他样本去除宿主 GRCm39
    if [ "$i" = "H.LX" ] || [ "$i" = "H.O" ]; then
        db_path="$HOME/db/kneaddata_db/hg_38_p14"
    else
        db_path="$HOME/db/kneaddata_db/GRCm39"
    fi

	kneaddata -i ${in1} -i ${in2} -o ${outdir} \
		-db ${db_path} \
		--trimmomatic $HOME/miniconda3/envs/kd0.7.4/share/trimmomatic \
		--trimmomatic-options 'ILLUMINACLIP:$HOME/miniconda3/envs/kd0.7.4/share/trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50' \
		--bowtie2-options '--end-to-end --sensitive --reorder' \
		--remove-intermediate-output \
		-t 30

	kneaddata_read_count_table --input ${outdir} --output ${outdir}/sum.tsv
done

# 合并sum.tsv
head -n 1 02_quality_control_rawdata/F1.1A/sum.tsv > 02_quality_control_rawdata/sum.tsv
for i in F1.1A F1.2A F1.3A F2.1A F2.2A F2.3A FG.1A FG.2A FG.3A L1.1A L1.2A L1.3A L2.1A L2.2A L2.3A
do
    # 修改第二行的第一列为样本名，然后添加该行到 sum 中
    sed -n '2s/^[^[:space:]]\+/'"$i"'/p' "02_quality_control_rawdata/${i}/sum.tsv" >> 02_quality_control_rawdata/sum.tsv
done

head -n 1 02_quality_control_rawdata/H.LX/sum.tsv >> 02_quality_control_rawdata/sum.tsv
for i in H.LX H.O
do
    # 修改第二行的第一列为样本名，然后添加该行到 sum 中
    sed -n '2s/^[^[:space:]]\+/'"$i"'/p' "02_quality_control_rawdata/${i}/sum.tsv" >> 02_quality_control_rawdata/sum.tsv
done
