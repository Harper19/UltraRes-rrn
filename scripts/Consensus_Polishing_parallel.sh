#!/bin/bash
# ---------------------------------------------------------------------------
# Paths: edit these before running (this file is a step-by-step command log,
# meant to be run block by block rather than executed as one script).
export PROJECT_DIR=/path/to/project      # working directory for this analysis
export SOFTWARE_DIR=/path/to/software    # dorado, usearch, etc.
# ---------------------------------------------------------------------------
# GPU: this step (racon + medaka polishing) was run on a GPU node (NVIDIA V100, 32 GB memory per GPU in this study);
# medaka inference uses the GPU.

cd ${PROJECT_DIR}/
batch=1
bar_list=""
minLen=3500
maxLen=6000  
Q=10
OUT_DIR=${PROJECT_DIR}/01_UMI_ref_paired_end_Q${Q}_${minLen}_${maxLen}_$batch

conda activate rrUMI_2
# Barcode to process. This step was originally run as a job array (one barcode per task):
bn=01
br=barcode${bn}
echo $br
batch=1

############################################################## 0. consens UMI ############################################################## 

minLen=3500
maxLen=6000  
Q=10
CUTADAPT=cutadapt
BWA=bwa
SAMTOOLS=samtools
SEQTK=seqtk
USEARCH=${SOFTWARE_DIR}/usearch
OUT_DIR=${PROJECT_DIR}/01_UMI_ref_paired_end_Q${Q}_${minLen}_${maxLen}_$batch
UMI_DIR=$OUT_DIR

  
### Bin reads based on UMIs ----------------------------------------------------
# conda activate seqtk
BINNING_DIR=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch
mkdir $BINNING_DIR
cd $BINNING_DIR

START_READ_CHECK=110
END_READ_CHECK=110
SEQTK=seqtk
GAWK=gawk

UMI_MATCH_ERROR=4 # maybe 4, lets test
UMI_MATCH_ERROR_single=$(($UMI_MATCH_ERROR/2))
UMI_MATCH_ERROR_SD=$(echo "3/2" | bc -l) # this is for Nanopore R10.3. Maybe need to change to 3/3 for R10.4. but the PNAS nexus paper did not change. 
RO_FRAC=0.3 #-O  Normalize read orientation fraction to 'O' if < 'O' reads are either +/- strand orientation.
MAX_BIN_SIZE=10000 # -N  Max number of reads with +/- orientation. [Default = 10000]
BIN_CLUSTER_RATIO=20 #-S  UMI bin size/UMI cluster size cutoff. [Default = 10]
# Checks whether UMI bin coverage is reasonable: a very small ratio means few reads in the UMI bin (possible noise or low coverage); a very large ratio means the bin takes up too many reads (possibly PCR amplification bias or duplicates)
# conda install conda-forge::parallel -y
GNUPARALLEL=parallel
THREADS=4


conda activate rrUMI_2
FW1=CAAGCAGAAGACGGCATACGAGAT # universal_fw
FW2=AGRRTTYGATYHTDGYTYAG # 16S
RV1=AATGATACGGCGACCACCGAGATC # universal_rv
RV2=CCRAMCTGTCTCACGACG # 23S

############################## 1. RACON polish ##############################
CUTADAPT=cutadapt
BWA=bwa
SAMTOOLS=samtools
SEQTK=seqtk
USEARCH=${SOFTWARE_DIR}/usearch11.0.667_i86linux32
GNUPARALLEL=parallel
RACON=racon
MINIMAP2=minimap2

minLen=3500
maxLen=6000
Q=10
OUT_DIR=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch
UMI_DIR=$OUT_DIR

conda activate seqtk
cd ${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch
# Wrapper
# Select a representative seed sequence from the raw reads of one UMI bin and polish it with several rounds of racon to generate a high-accuracy consensus sequence for downstream analysis.
seed_racon () {
  # Input
  local RB=$1 # $IN Directory containing UMI read bins in the format'umi*bins.fastq'.
  local OUT=$2 #$OUT  output directory
  local ROUNDS=$3 #$ROUNDS Number of racon polishing rounds.
  local PRESET=$4 #$PRESET. Minimap2 preset. 'map-ont' for Nanopore
  local RACON_ARG=$5 # $RACON_ARG Additional arguments for racon.

  # Name format
  local UMINO=${RB##*/}
  local UMINO=${UMINO%.*}
  local OUT=$OUT/$UMINO

  # Create dirs
  mkdir $OUT

  # Count binsize
  local BINSIZE=$(gawk 'NR%4==1 {N++} END {print N}' $RB)

  # Find seed read
  echo "usearch" $USEARCH
  ${SOFTWARE_DIR}/usearch11.0.667_i86linux32 -cluster_fast $RB -id 0.9 -strand both\
    -sizeout -centroids $OUT/${UMINO}_centroids.fa
  # ${SOFTWARE_DIR}/usearch11.0.667_i86linux32 -cluster_fast $RB -id 0.75 -strand both\
  #   -sizeout -centroids $OUT/${UMINO}_centroids.fa # old setting: very different reads could end up in the same bin
  ${SOFTWARE_DIR}/usearch11.0.667_i86linux32 -sortbysize $OUT/${UMINO}_centroids.fa -topn 1 \
    -fastaout $OUT/${UMINO}_sr.fa -relabel seed

  # Racon polishing
  echo "start Racon polishing"
  for i in `seq 1 $ROUNDS`; do
    minimap2 \
      -t 32 \
      -x $PRESET \
      $OUT/${UMINO}_sr.fa \
      $RB > $OUT/ovlp.paf
   # mapping all
    racon \
      -t 64 \
      -m 8 \
      -x -6 \
      -g -8 \
      -w 500 \
      $RACON_ARG \
      $RB \
      $OUT/ovlp.paf \
      $OUT/${UMINO}_sr.fa > $OUT/${UMINO}_tmp.fa

    mv $OUT/${UMINO}_tmp.fa $OUT/${UMINO}_sr.fa
  done

  # Rename and add binsize
  sed -i "s/^>.*/>$UMINO;ubs=$BINSIZE/" $OUT/${UMINO}_sr.fa
}
export -f seed_racon

ROUNDS=3
PRESET="map-ont"
RACON_ARG="--no-trimming"

  echo $OUT_DIR
  IN=$OUT_DIR/${br}_bin
  echo $IN
  OUT=$OUT_DIR/${br}_consensus
  echo $OUT
  mkdir -p $OUT

  echo "Assembly..."
  # Perform assembly in parallel
  find $IN -name 'umi*bins.fastq'  |\
    ( [[ -f "${SAMPLE}" ]] && grep -Ff $SAMPLE || cat ) |\
    $GNUPARALLEL -j 32 \
      --env seed_racon \
      --progress \
      "seed_racon {} $OUT $ROUNDS $PRESET $RACON_ARG $IN"


  #Collect seed-racon consensus sequences
  echo "Collect consensus sequences"
  OUT_NAME=${OUT##*/} # keep only the file name, not the full path
  echo $OUT_NAME
  OUT_NAME=${OUT_NAME%.*} # remove the file extension
  echo $OUT_NAME

  find $OUT/ \
    -mindepth 2 \
    -maxdepth 2 \
    -name "*_sr*.fa" \
    -exec cat {} \; \
    > $OUT/consensus_${OUT_NAME}.fa

echo $OUT/consensus_${OUT_NAME}.fa


############################## ############################## ##############################
############################## 2. Medaka polish ##############################
############################## ############################## ##############################
conda activate rrUMI_2
THREADS=32
MEDAKA_JOBS=16
# Format names
OUT_NAME=${OUT_DIR##*/}
echo "OUT_NAME:"$OUT_NAME
# Medaka jobs
MEDAKA_THREADS=$(( THREADS/MEDAKA_JOBS ))
CHUNK_SIZE=6000

medaka_align() {
  # Input
  local IN=$(cat)  # br_bin/job/umi_name.fastq
  local BINNING_DIR=$1 #br_bin
  local OUT_DIR=$2 #br_mapping

  # Name format
  local UMI_NAME=$(echo "$IN" | grep -o "umi.*bins")
  local UMI_BIN=$(find $BINNING_DIR -name ${UMI_NAME}.fastq)
  echo ${UMI_NAME}
  # Setup working directory
  mkdir $OUT_DIR/$UMI_NAME
  echo "$IN" > $OUT_DIR/$UMI_NAME/$UMI_NAME.fa

  # Map UMI reads to consensus
  mini_align  \
    -i $UMI_BIN \
    -r $OUT_DIR/$UMI_NAME/$UMI_NAME.fa \
    -m \
    -p $OUT_DIR/$UMI_NAME/${UMI_NAME} \
    -t 32
}
export -f medaka_align

  # Align the reads of one UMI bin (UMI_BIN) to its initial consensus (UMI_NAME.fa)
  # Output: a .bam file stored under $OUT_DIR/mapping/UMI_NAME/

# Calculate consensus probabilities--------------------------
echo "----------------------------------------------Calculate consensus probabilities--------------------------"
echo "-------------------------------------------------------------------------------------------------------------------------"

consensus_wrapper() {
  # Input
  local JOB_NR=$1
  local OUT_DIR=$2
  local MEDAKA_MODEL=$3
  local MEDAKA_THREADS=$4
  local CHUNK_SIZE=$5

  # Merge bam files

  ### Custom merge function required for > 1024 files
  bam_merge() {
    # Input
    local OUT_DIR=$1
	  local JOB_NR=$2
    awk -v OUT_BAM="$OUT_DIR/${JOB_NR}" '
      # Print HD line
      (NR==1 && $1 ~ /^@HD$/){print $0 > OUT_BAM ".header"}
      # Print reference lines
      ($1 ~ /^@SQ$/){print $0 > OUT_BAM ".header"}
      # Print reads lines
      ($1 ~ !/^@/){print $0 > OUT_BAM ".body"}
    '
  }
  export -f bam_merge

  ### View bam files in parallel and pipe into merge function
  cat |\
    parallel -j $MEDAKA_THREADS \
      --env bam_merge \
      samtools view -h {} |\
      bam_merge $OUT_DIR $JOB_NR
# ${PROJECT_DIR}/02_read_binning_Q10_3500_6000/barcode20_mapping/umi137170bins
  ### Build bam file
  cat $OUT_DIR/${JOB_NR}.header $OUT_DIR/${JOB_NR}.body  |\
    samtools view -b - > $OUT_DIR/${JOB_NR}.bam
  rm $OUT_DIR/${JOB_NR}.header $OUT_DIR/${JOB_NR}.body

  # Index bam
  samtools index \
    $OUT_DIR/${JOB_NR}.bam

  # Medaka consensus
  echo "------------- Medaka consensus"
  medaka inference \
    $OUT_DIR/${JOB_NR}.bam \
    $OUT_DIR/${JOB_NR}_consensus.hdf \
    --threads $MEDAKA_THREADS \
    --model $MEDAKA_MODEL \
	  --chunk_len $CHUNK_SIZE
}
export -f consensus_wrapper


# ${SOFTWARE_DIR}/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0: r1041_e82_400bps_sup_v5.2.0
MEDAKA_MODEL=r1041_e82_400bps_sup_v5.2.0
# 'r1041_e82_400bps_sup_variant_v5.0.0',
MEDEKA_ROUND=1

    BINNING_DIR=${OUT_DIR}/${br}_bin
    echo $BINNING_DIR
    CONSENSUS_FILE=${OUT_DIR}/${br}_consensus/consensus_${br}_consensus.fa
    ll $CONSENSUS_FILE
  for MR in `seq 1 $MEDEKA_ROUND`
  do
    mkdir $OUT_DIR/${br}_mapping_${MR}
    echo "Mapping UMI bins to consensus..."
    echo "Barcode: "$br", medeka polishing round:" $MR
    echo "Consensus file:"$CONSENSUS_FILE
    echo "Output_dir:"$OUT_DIR/${br}_mapping_${MR}

    # Individual mapping of UMI bins to consensus
    cat $CONSENSUS_FILE |\
      seqtk seq -l0 - |\
      ( [[ -f "${SAMPLE}" ]] && grep -A1 -Ff $SAMPLE | sed '/^--$/d' || cat ) |\
      parallel \
        --env medaka_align \
        --progress  \
        --recstart ">" \
        -N 1 \
        --pipe \
        "medaka_align \
           $BINNING_DIR \
           $OUT_DIR/${br}_mapping_${MR}
        "

  # output: bamw file under $OUT_DIR/mapping/UMI_NAME/
  echo "Medaka polishing..."
  mkdir $OUT_DIR/${br}_consensus_hdf_${MR}
  echo $br
  echo "INPUT: "
  echo "OUTPUT polished consensus file: "$OUT_DIR/${br}_consensus_hdf_${MR}/consensus_${br}_stitch.fa
  
  find $OUT_DIR/${br}_mapping_${MR}/ \
    -mindepth 2 \
    -maxdepth 2 \
    -type f \
    -name "umi*bins.bam" |\
  parallel -j $MEDAKA_JOBS \
    --env consensus_wrapper \
    --progress \
    -N1 \
    --roundrobin \
    --pipe \
    "consensus_wrapper \
       {#} \
       $OUT_DIR/${br}_consensus_hdf_${MR}  \
       $MEDAKA_MODEL \
       $MEDAKA_THREADS \
       $CHUNK_SIZE"

  # Stitch consensus sequences
  echo "--------- Medaka stitch"
  medaka sequence --threads 32 $OUT_DIR/${br}_consensus_hdf_${MR}/*_consensus.hdf $CONSENSUS_FILE  $OUT_DIR/${br}_consensus_hdf_${MR}/consensus_${br}_stitch.fa

  # Clean consensus header
    echo "--------- stitch head clean"
    sed -i -e "s/:.*//" -e "s/_segment.*//" $OUT_DIR/${br}_consensus_hdf_${MR}/consensus_${br}_stitch.fa

    # CONSENSUS_FILE=$OUT_DIR/${br}_consensus_hdf_${MR}/consensus_${br}_stitch.fa
  done




# compare raw read, consensus, medaka polished consensus.
for bn in $bar_list
do
  br=barcode$bn
  echo "-----------------------raw read count: " $br
  rq="${PROJECT_DIR}/00_SUP_demux_Q${Q}_${minLen}_${maxLen}_$batch/${br}_QC.fastq"
  grep "^@" $rq | wc -l
  
  echo "-----------------------consensus read count: " $br
  cq="${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/${br}_consensus/consensus_${br}_consensus.fa"
  grep "^>" $cq | wc -l
  
  echo "-----------------------medaka polished consensus read count: " $br
  mq="${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/${br}_consensus_hdf_1/consensus_${br}_stitch.fa"
  grep "^>" $mq | wc -l
done
