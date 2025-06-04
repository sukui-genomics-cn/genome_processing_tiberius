#!/bin/bash

#data_dir="/home/share/huadjyin/home/yinpeng/refseq/model_species_tiberius"

source /home/HPCBase/tools/anaconda3/etc/profile.d/conda.sh
conda init bash


conda activate /home/share/huadjyin/home/s_liulin4/miniconda3/envs/Evo

data_dir="/home/share/huadjyin/home/yinpeng/refseq/30_verbertate_others"

for spec_dir in "$data_dir"/*/; do
#    if [[ "$spec_dir" == "$data_dir/Homo_sapiens/" ]]; then
#        continue
#    fi
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
