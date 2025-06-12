data_root=/home/nvme01/data/refseq_fish_30
filter_file=/home/nvme01/sukui/01.data/GeneStructure/refseq_fish.30_50k_filter.1M/filter.txt
save_path=/home/nvme01/sukui/01.data/GeneStructure/refseq_fish.30_50k_filter.1M
project_root=/home/nvme01/sukui/03.project/genome_processing_tiberius

cd $project_root

eccho "data_root: $data_root"
eccho "save_path: $save_path"
eccho "filter_file: $filter_file"

python src/split_dataset.py \
    -dp $data_root \
    -sp $save_path \
    -tr 0.97\
    -vr 0.01 \
    -ttr 0.02 \
    -sd "chunks_50004_1" \
    -fa "*.pkl" \
    -rc "true" \
    -ffp $filter_file
