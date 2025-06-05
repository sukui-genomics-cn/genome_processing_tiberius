#!/bin/bash

data_dir="/home/share/huadjyin/home/yinpeng/refseq/Homo_sapiens_T2T"


for spec_dir in "$data_dir"/*/; do
    anno_dir=${spec_dir}"anno_tiberius"
    fasta_dir=${spec_dir}"fasta"
    spec=$(basename "$spec_dir")
    ls ${anno_dir}/*.pkl | while read -r file; do
        pkl=$file
        pkl_file=$(basename "$file")
        chr=$(echo "$pkl_file" | cut -d '_' -f 1,2)
        fasta=${fasta_dir}/${spec}_${chr}"_forward.txt"

        python ./get_tiberius_chunks.py \
            --fasta ${fasta} \
            --pkl ${pkl} \
            --chunk_size 9999 \
            --overlap 0 \
            --out_dir ${spec_dir} 
    done
done


for spec_dir in "$data_dir"/*/; do
    anno_dir=${spec_dir}"anno_tiberius"
    fasta_dir=${spec_dir}"fasta"
    spec=$(basename "$spec_dir")
    ls ${anno_dir}/*.pkl | while read -r file; do
        pkl=$file
        pkl_file=$(basename "$file")
        chr=$(echo "$pkl_file" | cut -d '_' -f 1,2)
        fasta=${fasta_dir}/${spec}_${chr}"_forward.txt"

        python ./get_tiberius_chunks.py \
            --fasta ${fasta} \
            --pkl ${pkl} \
            --chunk_size 50004 \
            --overlap 0 \
            --out_dir ${spec_dir}
    done
done



for spec_dir in "$data_dir"/*/; do
    anno_dir=${spec_dir}"anno_tiberius"
    fasta_dir=${spec_dir}"fasta"
    spec=$(basename "$spec_dir")
    ls ${anno_dir}/*.pkl | while read -r file; do
        pkl=$file
        pkl_file=$(basename "$file")
        chr=$(echo "$pkl_file" | cut -d '_' -f 1,2)
        fasta=${fasta_dir}/${spec}_${chr}"_forward.txt"

        python ./get_tiberius_chunks.py \
            --fasta ${fasta} \
            --pkl ${pkl} \
            --chunk_size 200016 \
            --overlap 0 \
            --out_dir ${spec_dir}
    done
done


for spec_dir in "$data_dir"/*/; do
    anno_dir=${spec_dir}"anno_tiberius"
    fasta_dir=${spec_dir}"fasta"
    spec=$(basename "$spec_dir")
    ls ${anno_dir}/*.pkl | while read -r file; do
        pkl=$file
        pkl_file=$(basename "$file")
        chr=$(echo "$pkl_file" | cut -d '_' -f 1,2)
        fasta=${fasta_dir}/${spec}_${chr}"_forward.txt"

        python ./get_tiberius_chunks.py \
            --fasta ${fasta} \
            --pkl ${pkl} \
            --chunk_size 800064 \
            --overlap 0 \
            --out_dir ${spec_dir}
    done
done



for spec_dir in "$data_dir"/*/; do
    anno_dir=${spec_dir}"anno_tiberius"
    fasta_dir=${spec_dir}"fasta"
    spec=$(basename "$spec_dir")
    ls ${anno_dir}/*.pkl | while read -r file; do
        pkl=$file
        pkl_file=$(basename "$file")
        chr=$(echo "$pkl_file" | cut -d '_' -f 1,2)
        fasta=${fasta_dir}/${spec}_${chr}"_forward.txt"

        python ./get_tiberius_chunks.py \
            --fasta ${fasta} \
            --pkl ${pkl} \
            --chunk_size 1600128 \
            --overlap 0 \
            --out_dir ${spec_dir}
    done
done



























#    python ./get_tiberius_chunks.py \
