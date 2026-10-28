#!/bin/bash
# ---------------------------------------------------------------------------
# Paths: edit these before running. This file is a step-by-step command log,
# meant to be run block by block rather than executed as one script.
export PROJECT_DIR=/path/to/project      # working directory for this analysis
export SOFTWARE_DIR=/path/to/software    # dorado, usearch, etc.
# ---------------------------------------------------------------------------

batch=1
mkdir ${PROJECT_DIR}/00_SUP_$batch        # basecalled fastq files
mkdir ${PROJECT_DIR}/00_SUP_demux_$batch  # demultiplexed fastq files

#----------------------------------- 1) Basecalling and demultiplexing -----------------------------------
# GPU: this step was run on GPU nodes.
# Dorado is heavily optimised for NVIDIA A100 and H100 GPUs and delivers its best performance on them.
# Input: pod5 files in ${PROJECT_DIR}/00_Pod5_$batch/
dorado="${SOFTWARE_DIR}/dorado-1.1.1-linux-x64/bin/dorado"

time $dorado \
    basecaller ${SOFTWARE_DIR}/dorado-1.3.1-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0  \
     ${PROJECT_DIR}/00_Pod5_$batch/ --recursive  --emit-fastq --min-qscore 8 \
    > ${PROJECT_DIR}/00_SUP_$batch/SUP_raw_$batch.fastq

time  $dorado demux --kit-name SQK-NBD114-24 --emit-fastq \
    --output-dir  ${PROJECT_DIR}/00_SUP_demux_$batch/ \
    ${PROJECT_DIR}/00_SUP_$batch/SUP_raw_$batch.fastq

#----------------------------------- 2) Read length and quality filtering -----------------------------------
minLen=3500
maxLen=6000
Q=10
mkdir ${PROJECT_DIR}/00_SUP_demux_Q${Q}_${minLen}_${maxLen}_$batch
cd ${PROJECT_DIR}/00_SUP_demux_Q${Q}_${minLen}_${maxLen}_$batch

conda activate chopper
LOG_FILE="${PROJECT_DIR}/00_SUP_demux_Q${Q}_${minLen}_${maxLen}_$batch/chopper_summary.log"
> "$LOG_FILE"

bar_list="01 02 03 04 05 06 07 08 09 10 11 12 13"
for bn in $bar_list
do
  br=barcode$bn
  echo "=== Processing $br at $(date) ===" >> "$LOG_FILE"
  chopper -q ${Q} --minlength $minLen --maxlength $maxLen \
          -i ${PROJECT_DIR}/00_SUP_demux_$batch/$br.fastq \
      > ${PROJECT_DIR}/00_SUP_demux_Q${Q}_${minLen}_${maxLen}_$batch/${br}_QC.fastq \
      2>> "$LOG_FILE"

  echo "=== Completed $br at $(date) ===" >> "$LOG_FILE"
  echo "" >> "$LOG_FILE"
done
