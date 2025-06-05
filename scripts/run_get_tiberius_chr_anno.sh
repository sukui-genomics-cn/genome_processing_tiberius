#!/bin/bash


spec_dir="/home/share/huadjyin/home/s_liulin4/datasets/referenceGenome/T2T"
#anno_dir=${spec_dir}/"anno_tiberius"
out_dir="/home/share/huadjyin/home/yinpeng/refseq/Homo_sapiens_T2T"

fasta="${spec_dir}/*.fna.gz"
gtf="${spec_dir}/*.gtf.gz"
python ./get_tiberius_chr_anno.py \
    --fasta $fasta \
    --anno $gtf \
    --filter_inframestop \
    --filter_short 90 \
    --type tiberius \
    --transition \
    --save_chr \
    --output $out_dir 

