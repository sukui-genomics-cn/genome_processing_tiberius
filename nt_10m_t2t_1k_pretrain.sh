#export DNA_LLM_PATH=/home/share/huadjyin/home/s_sukui/03_project/01_GeneLLM/DNA_LLM
#export T2T_PATH=/home/share/huadjyin/home/s_sukui/02_data/07_genomics_data/T2T/pretrain

export DNA_LLM_PATH=/home/STOmics/01_user/sukui/03_project/DNA_LLM
export T2T_PATH=/data/sukui_data/01_data/01_genomics_data/T2T/pretrain
export WANDB_WATCH="all"
cd $DNA_LLM_PATH

CUDA_VISIBLE_DEVICES=1 torchrun --nproc_per_node=1 --master_port=29539 pretrain.py \
    experiment=pretrain/t2t_nt_10M \
    dataset.dest_path=$T2T_PATH \
    dataset.dataset_name=8K \
    dataset.batch_size=2 \
    dataset.max_length=8192 \
    train.gradient_accumulation_steps=1 \
    train.logging_steps=100 \
    train.max_steps=1000 \
    paths.output_root=$DNA_LLM_PATH/outputs/debugs \
    wandb.mode=offline \
    wandb.job_type=debug

#CUDA_VISIBLE_DEVICES=1 torchrun --nproc_per_node=1 --master_port=29509 pretrain.py \
#    experiment=pretrain/t2t_nt_250M \
#    dataset.dest_path=$T2T_PATH \
#    dataset.dataset_name=1K \
#    dataset.batch_size=128 \
#    dataset.max_length=1024 \
#    train.gradient_accumulation_steps=2 \
#    train.logging_steps=1000 \
#    train.max_steps=-1 \
#    paths.output_root=$DNA_LLM_PATH/outputs/pretrain \
#    wandb.mode=online \
#    wandb.job_type=H100_T2T_10M