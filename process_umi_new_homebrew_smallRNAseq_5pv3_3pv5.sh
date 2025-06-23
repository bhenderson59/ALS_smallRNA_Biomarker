#!/bin/bash
#this script take raw SE75 fastq.gz form qiagen lib prep and trims off adaptor and adds UMI to fastq info line

RAW_FQ_GZ=$1
OUT_DIR=$2

BASE=$(echo $RAW_FQ_GZ | awk -F "/" '{print $(NF)}' | awk -F ".fastq" '{print$1}') 

TOO_SHORT_FILE=$OUT_DIR"/"$BASE"_too_short.fastq.gz"
NO_3AD_FILE=$OUT_DIR"/"$BASE"_no_3AD.fastq.gz"
NO_5AD_FILE=$OUT_DIR"/"$BASE"_no_5AD.fastq.gz"

OUT_FINAL_FQ_GZ=$OUT_DIR"/"$BASE"_trim_read_and_umi.fastq.gz"

module load cluster/cutadapt/2.8
module load cluster/perl-modules

HOMEBREW_3P_AD="ACGGGCTAATATTTATCGGTGGAGCATCACGATCTCGTAT"

FIVE_PRIME_AD1_SEQ="^CAGTCG"
FIVE_PRIME_AD2_SEQ="^TGACTC"
FIVE_PRIME_AD3_SEQ="^GCTAGA"
FIVE_PRIME_AD4_SEQ="^ATCGAT"

cutadapt -a $HOMEBREW_3P_AD -e 0.3 --match-read-wildcards --untrimmed-output=$NO_3AD_FILE $RAW_FQ_GZ | cutadapt -e 0.5 --match-read-wildcards -m 15 -O 6 -n 1 -g $FIVE_PRIME_AD1_SEQ -g $FIVE_PRIME_AD2_SEQ -g $FIVE_PRIME_AD3_SEQ -g $FIVE_PRIME_AD4_SEQ --untrimmed-output=$NO_5AD_FILE --too-short-output=$TOO_SHORT_FILE - | /cluster/home/broberts/ben_miRNA_seq_analysis_example/process_umi_homebrew_3p_UMI_only.pl | gzip -c > $OUT_FINAL_FQ_GZ

