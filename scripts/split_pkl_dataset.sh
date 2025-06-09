data_root=/home/nvme01/data/refseq_fish_30
filter_file=/home/nvme01/sukui/01.data/GeneStructure/refseq_fish.30_50k/filter.txt
project_root=/home/nvme01/sukui/03.project/genome_processing_tiberius

python src/split_dataset.py \
    -dp $data_root \
    -sp /home/nvme01/sukui/01.data/GeneStructure/refseq_fish.30_50k \
    -ffp $filter_file


