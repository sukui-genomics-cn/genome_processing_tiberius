#!/bin/bash


source /home/HPCBase/tools/anaconda3/etc/profile.d/conda.sh
conda init bash

conda activate /home/share/huadjyin/home/s_liulin4/miniconda3/envs/Evo

spec_dir="/home/share/huadjyin/home/s_liulin4/datasets/refseq/vertebrate_others/down"
#anno_dir=${spec_dir}/"anno_tiberius"
#fish_list=("Acanthochromis_polyacanthus" "Carassius_auratus" "Chanos_chanos" "Gadus_morhua" "Mugil_cephalus" "Oncorhynchus_mykiss" "Salmo_salar" "Thunnus_albacares" "Xiphias_gladius" "Clupea_harengus")

#bird_list=("Aquila_chrysaetos" "Accipiter_gentilis" "Corvus_brachyrhynchos" "Haliaeetus_albicilla" "Anas_acuta" "Catharus_ustulatus" "Pelecanus_crispus" "Egretta_garzetta" "Struthio_camelus" "Chamaea_fasciata")

#bird_list=("Harpia_harpyja")

species_file="/home/share/huadjyin/home/yinpeng/refseq/30_verbertate_others/add.txt"

mapfile -t species_list < $species_file

echo "species_list=("
for species in "${species_list[@]}"; do
    echo "  \"$species\""
done
echo ")"


out_dir="/home/share/huadjyin/home/yinpeng/refseq/30_verbertate_others"

export script_PATH=/home/share/huadjyin/home/s_liulin4/data_preprocess/genome_processing_tiberius_lin
cd $script_PATH



for spec in "$spec_dir"/*/; do
    spec_name=$(basename "$spec")
    if [[ " ${species_list[@]} " =~ " ${spec_name} " ]]; then
        echo "$spec_name 匹配成功"
        fasta="${spec}*.fna.gz"
        gtf="${spec}*.gtf.gz"
        python ./get_tiberius_chr_anno.py \
        --fasta $fasta \
        --anno $gtf \
        --filter_inframestop \
        --filter_short 90 \
        --type tiberius \
        --transition \
        --save_chr \
        --output $out_dir 
    fi
done








#    python ./get_tiberius_chr_anno.py \
#    --fasta /home/share/huadjyin/home/s_liulin4/datasets/refseq/tiberius_chunks/Homo_sapiens/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
#    --anno /home/share/huadjyin/home/s_liulin4/datasets/refseq/tiberius_chunks/Homo_sapiens/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz \
#    --filter_inframestop \
#    --filter_short 90 \
 #   --type tiberius \
#    --transition \
#    --save_chr \
#    --output /home/share/huadjyin/home/s_liulin4/datasets/refseq/tiberius_chunks/Homo_sapiens/test
