#!/bin/bash

CONSOL_FQ_DIR="/cluster/home/bhenderson/ALS_biomarkers/Sequencing/ALS_Group3_SeqCenter/fq_paths"

CONSOL_FQ_LIST="/cluster/home/bhenderson/ALS_biomarkers/Sequencing/ALS_Group3_SeqCenter/bwh_04.02.2024_ALS_biomarker_Group3_fq_gz_list.txt"

ls $CONSOL_FQ_DIR/*'.fastq.gz' > $CONSOL_FQ_LIST

TRIM_OUT_DIR="/cluster/home/bhenderson/ALS_biomarkers/Sequencing/ALS_Group3_SeqCenter/TRIMMED"

if [ ! -e "$TRIM_OUT_DIR" ]; then
 mkdir $TRIM_OUT_DIR
fi

while read FQ_IN; do
 BASE=$(echo $FQ_IN | awk -F "/" '{print$(NF)}' | awk -F ".fastq" '{print$1}')
 echo $BASE

 LOG_TRIM=$TRIM_OUT_DIR'/'$BASE'_trim_log.txt'

 sbatch -n 1 --time-min=01:00:00 --job-name $BASE --output $LOG_TRIM /cluster/home/bhenderson/BWH_Utility_Scripts/sRNA_pipeline_scripts/process_umi_new_homebrew_smallRNAseq_5pv3_3pv5.sh $FQ_IN $TRIM_OUT_DIR

done < $CONSOL_FQ_LIST
