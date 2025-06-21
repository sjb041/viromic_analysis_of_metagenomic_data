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

    # H样本去除宿主 hg_38_p14，其他样本去除宿主 mouse_C57BL_6NJ
    if [ "$i" = "H.LX" ] || [ "$i" = "H.O" ]; then
        db_path="$HOME/db/kneaddata_db/hg_38_p14"
    else
        db_path="$HOME/db/kneaddata_db/mouse_C57BL_6NJ"
    fi

	kneaddata -i ${in1} -i ${in2} -o ${outdir} \
		-db ${db_path} \
		--trimmomatic $HOME/miniconda3/envs/kd0.7.4/share/trimmomatic \
		--trimmomatic-options 'ILLUMINACLIP:$HOME/miniconda3/envs/kd0.7.4/share/trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50' \
		--bowtie2-options '--end-to-end --sensitive --reorder' \
		--remove-intermediate-output \
		--run-fastqc-end \
		-t 30

	kneaddata_read_count_table --input ${outdir} --output ${outdir}/sum.tsv
done
