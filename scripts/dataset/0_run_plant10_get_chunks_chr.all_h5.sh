#!/bin/bash

data_dir="/home/nvme01/data/dicots/"
output_root=$data_dir

root=/home/nvme01/sukui/03.project/genome_processing_tiberius/
cd $root

for spec_dir in "$data_dir"/*/; do
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
        # chr=$(echo "$pkl_file" | cut -d '_' -f 1,2)
        chr=${pkl_file%_*} # Extract everything before the last underscore
        fasta=${fasta_dir}/${spec}_${chr}"_forward.txt.gz"

        python ./src/get_tiberius_chr.all_h5.py \
            --fasta ${fasta} \
            --pkl ${pkl} \
            --out_dir ${spec_dir} 
    done
done





















#    python ./get_tiberius_chunks.py \
