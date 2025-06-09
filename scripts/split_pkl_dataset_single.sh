data_root=/home/nvme01/data/refseq_fish_30/Danio_rerio/chunks_50004_0/
project_root=/home/nvme01/sukui/03.project/genome_processing_tiberius

cd $project_root
python src/split_dataset.py \
    -dp $data_root \
    -sp /home/nvme01/sukui/01.data/GeneStructure/Danio.rerio_50k/chunks_50004_0/ \
    -tr 0.95\
    -vr 0.01 \
    -ttr 0.01 \
    -fa "*.pkl" \