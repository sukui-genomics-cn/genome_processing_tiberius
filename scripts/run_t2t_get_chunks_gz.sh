#!/bin/bash


data_dir=/home/nvme01/data/refseq_fish_30/
output_root=$data_dir

project_dir=/home/nvme01/sukui/03.project/genome_processing_tiberius/
cd $project_dir

for spec_dir in "$data_dir"*/; do
    anno_dir=${spec_dir}"anno_tiberius"
    fasta_dir=${spec_dir}"fasta"
    spec=$(basename "$spec_dir")
    out_dir=${output_root}${spec}
    echo "Processing species: $spec"
    echo "Annotation directory: $anno_dir"
    echo "Fasta directory: $fasta_dir"
    echo "Output directory: $out_dir"
    if [ ! -d "$out_dir" ]; then
        mkdir -p "$out_dir"
        echo "Created output directory: $out_dir"
    fi

    ls ${anno_dir}/*.pkl.gz | while read -r file; do
        pkl=$file
        pkl_file=$(basename "$file")
        chr=$(echo "$pkl_file" | cut -d '_' -f 1,2)
        fasta=${fasta_dir}/${spec}_${chr}"_forward.txt.gz"

        python ./src/get_tiberius_chunks_gz.py \
            --fasta ${fasta} \
            --pkl ${pkl} \
            --chunk_size 50004 \
            --overlap 1 \
            --out_dir ${out_dir} 
    done
done