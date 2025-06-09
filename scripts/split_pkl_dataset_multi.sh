data_root=/home/nvme01/data/refseq_fish_30
filter_file=/home/nvme01/sukui/01.data/GeneStructure/refseq_fish.30_50k/filter.txt
project_root=/home/nvme01/sukui/03.project/genome_processing_tiberius

cd $project_root
python src/split_dataset.py \
    -dp $data_root \
    -sp /home/nvme01/sukui/01.data/GeneStructure/refseq_fish.30_50k \
    -tr 0.95\
    -vr 0.01 \
    -ttr 0.01 \
    -fa "*.pkl" \
    -rc "true" \
    -ffp $filter_file