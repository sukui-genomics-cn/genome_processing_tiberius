#!/bin/bash

data_dir="/home/share/huadjyin/home/yinpeng/refseq/model_species_tiberius"
out_dir="/home/share/huadjyin/home/yinpeng/refseq/model_species_tiberius/Homo_sapiens/inter_gene_chunks_chrom"

for spec_dir in "$data_dir"/Homo_sapiens/; do
    anno_dir=${spec_dir}"anno_tiberius"
    fasta_dir=${spec_dir}"fasta"
    spec=$(basename "$spec_dir")
    ls ${anno_dir}/NC_000002.12_forward.pkl | while read -r file; do
        pkl=$file
        pkl_file=$(basename "$file")
        chr=$(echo "$pkl_file" | cut -d '_' -f 1,2)
        fasta=${fasta_dir}/${spec}_${chr}"_forward.txt"

        echo "------Cut chromosomes of ${spec}-------"
        python ./get_intergenic_gene_chunks.py \
            --fasta ${fasta} \
            --pkl ${pkl} \
            --chunk_size 6144 \
            --overlap 0 \
            --out_dir ${out_dir} 
    done
done
