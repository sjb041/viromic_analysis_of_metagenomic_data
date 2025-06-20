mkdir 03_cleandata

# 获取样本名
sample=$(ls 02_quality_control_rawdata)

# 提取 cleandata ，文件重命名，压缩
for i in $sample
do
    samplei=$i
    sampleo=${samplei/./_}

    mv 02_quality_control_rawdata/${samplei}/${samplei}_L1_1_kneaddata_paired_1.fastq 03_cleandata/${sampleo}_R1.fastq &
    mv 02_quality_control_rawdata/${samplei}/${samplei}_L1_1_kneaddata_paired_2.fastq 03_cleandata/${sampleo}_R2.fastq &
    mv 02_quality_control_rawdata/${samplei}/${samplei}_L1_1_kneaddata_unmatched_1.fastq 03_cleandata/${sampleo}_s1.fastq &
    mv 02_quality_control_rawdata/${samplei}/${samplei}_L1_1_kneaddata_unmatched_2.fastq 03_cleandata/${sampleo}_s2.fastq &
    wait

    cat 03_cleandata/${sampleo}_s1.fastq 03_cleandata/${sampleo}_s2.fastq > 03_cleandata/${sampleo}_S.fastq &
    wait

    pigz -p 8 03_cleandata/*.fastq &
    wait
done
