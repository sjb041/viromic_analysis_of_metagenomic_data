nextflow run /home/shijiabin/opt/emg-viral-pipeline-3.0.2/main.nf \
    --fasta 02_rename/H_LX_renamed.fa \
    --output 03_virify/v3v4_ncbi2022_hex_bex/H_LX \
    --viphog_version v3 \
    --meta_version v4 \
    --databases /home/shijiabin/db/virify_db \
    --onlyannotate true \
    --publish_all \
    --hmmextend \
    --blastextend \
    --length 0 \
    -profile local,docker












