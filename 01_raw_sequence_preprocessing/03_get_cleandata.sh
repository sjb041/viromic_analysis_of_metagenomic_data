mkdir -p 03_cleandata

# 样本名
sample="F1.1A F1.2A F1.3A F2.1A F2.2A F2.3A FG.1A FG.2A FG.3A L1.1A L1.2A L1.3A L2.1A L2.2A L2.3A H.LX H.O"

# 提取 cleandata ，文件重命名，压缩
for i in $sample
do
    samplei=$i
    sampleo=${samplei/./_}

    mv 02_quality_control_rawdata/${samplei}/${samplei}_L1_1_kneaddata_paired_1.fastq 03_cleandata/${sampleo}_R1.fastq 
    mv 02_quality_control_rawdata/${samplei}/${samplei}_L1_1_kneaddata_paired_2.fastq 03_cleandata/${sampleo}_R2.fastq 
    mv 02_quality_control_rawdata/${samplei}/${samplei}_L1_1_kneaddata_unmatched_1.fastq 03_cleandata/${sampleo}_s1.fastq 
    mv 02_quality_control_rawdata/${samplei}/${samplei}_L1_1_kneaddata_unmatched_2.fastq 03_cleandata/${sampleo}_s2.fastq 

    cat 03_cleandata/${sampleo}_s1.fastq 03_cleandata/${sampleo}_s2.fastq > 03_cleandata/${sampleo}_S.fastq 
done

pigz -p 8 03_cleandata/*.fastq
