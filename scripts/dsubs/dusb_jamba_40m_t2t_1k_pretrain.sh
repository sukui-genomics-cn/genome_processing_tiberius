#!/bin/bash
#DSUB -n Jamba40MT2TPretrain
#DSUB -A root.project.P24Z10200N0983
#DSUB -R 'cpu=64;gpu=2;mem=100000'
#DSUB -N 2
#DSUB -eo /home/share/huadjyin/home/s_sukui/03_project/01_GeneLLM/DNA_LLM/logs/pretrain/Jamba40MT2TPretrain.%J.%I.err
#DSUB -oo /home/share/huadjyin/home/s_sukui/03_project/01_GeneLLM/DNA_LLM/logs/pretrain/Jamba40MT2TPretrain.%J.%I.out

## Set scripts
RANK_SCRIPT="./scripts/tasks/pretrain/jamba_40m_t2t_1k_pretrain.sh"

###Set Start Path
JOB_PATH="/home/share/huadjyin/home/s_sukui/03_project/01_GeneLLM/DNA_LLM"

## Set NNODES
NNODES=2
GPUS=2

## Create nodefile

JOB_ID=${BATCH_JOB_ID}
NODEFILE=${JOB_PATH}/outputs/tmp/${JOB_ID}.nodefile
# touch ${NODEFILE}
touch $NODEFILE
cat $CCS_ALLOC_FILE | grep ^cyclone001-agent | awk '{print $1,"slots="$2}' > ${JOB_PATH}/outputs/tmp/${JOB_ID}.nodefile
cat $CCS_ALLOC_FILE | grep ^cyclone001-agent | awk '{print $1}' > ${NODEFILE}
cat ${CCS_ALLOC_FILE} > :q!${JOB_PATH}/outputs/tmp/CCS_ALLOC_FILE

cd ${JOB_PATH};/usr/bin/bash ${RANK_SCRIPT} ${NNODES} ${GPUS} ${NODEFILE}
