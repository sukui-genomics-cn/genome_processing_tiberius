#!/bin/bash

root=/home/nvme01/sukui/03.project/genome_processing_tiberius
cd $root

python src/dnaseq_preprocess.v2.py \
    --input_dir /home/nvme01/T2T \
    --output_dir /home/nvme01/sukui/01.data/T2T/t2t_chr.all \
    --chunk_size 1000




















#    python ./get_tiberius_chunks.py \
