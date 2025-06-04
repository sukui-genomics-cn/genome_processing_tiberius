#!/bin/bash
#DSUB -n Getchunks
#DSUB -A root.project.P24Z10200N0983
#DSUB -R 'cpu=64;gpu=0;mem=100000'
#DSUB -N 2
#DSUB -eo /home/share/huadjyin/home/s_liulin4/data_preprocess/genome_processing_tiberius_lin/logs/Cuttingseqs.%J.%I.err
#DSUB -oo /home/share/huadjyin/home/s_liulin4/data_preprocess/genome_processing_tiberius_lin/logs/Cuttingseqs.%J.%I.out

## Set scripts
RANK_SCRIPT="/home/share/huadjyin/home/s_liulin4/data_preprocess/genome_processing_tiberius_lin/run_get_chunks_all.sh"

###Set Start Path
JOB_PATH="/home/share/huadjyin/home/s_liulin4/data_preprocess/genome_processing_tiberius_lin"

## Set NNODES
NNODES=1
GPUS=0

## Create nodefile

JOB_ID=${BATCH_JOB_ID}
NODEFILE=${JOB_PATH}/outputs/${JOB_ID}.nodefile
# touch ${NODEFILE}
touch $NODEFILE
cat $CCS_ALLOC_FILE | grep ^cyclone001-agent | awk '{print $1,"slots="$2}' > ${JOB_PATH}/outputs/${JOB_ID}.nodefile
cat $CCS_ALLOC_FILE | grep ^cyclone001-agent | awk '{print $1}' > ${NODEFILE}
cat ${CCS_ALLOC_FILE} > :q!${JOB_PATH}/outputs/$CCS_ALLOC_FILE

cd ${JOB_PATH};/usr/bin/bash ${RANK_SCRIPT} ${NNODES} ${GPUS} ${NODEFILE}
