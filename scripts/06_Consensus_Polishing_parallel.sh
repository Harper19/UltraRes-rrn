#!/bin/bash
# ---------------------------------------------------------------------------
# Paths: edit these before running (this file is a step-by-step command log,
# meant to be run block by block rather than executed as one script).
export PROJECT_DIR=/path/to/project      # working directory for this analysis
export SOFTWARE_DIR=/path/to/software    # dorado, usearch, etc.
export GROND_DIR=/path/to/GROND          # GROND database checkout (07_* scripts)
export LOCAL_DATA=/path/to/local/data    # raw pod5 storage (03_* only)
# SLURM: replace YOUR_ACCOUNT and your_email@example.com in the #SBATCH lines.
# ---------------------------------------------------------------------------
# # parallel script
# salloc -q inferno -A YOUR_ACCOUNT -N4 --ntasks-per-node=8 --time=2:30:00 --mem=64G
# module load anaconda3/2022.05
# pace-quota

cd ${PROJECT_DIR}/
batch=1
bar_list="01 02 03 04 05 06 07 08 09 10 11 12 13"
minLen=3500
maxLen=6000  
Q=10
OUT_DIR=${PROJECT_DIR}/01_UMI_ref_paired_end_Q${Q}_${minLen}_${maxLen}_$batch

mkdir ${PROJECT_DIR}/Scripts/
salloc -A YOUR_ACCOUNT -N1 --mem-per-gpu=12G -qinferno -t0:15:00 --gres=gpu:H100:1
module load anaconda3/2022.05                   # Load module dependencies
conda activate rrUMI_2
# SLURM_ARRAY_TASK=1
{
script_polish_array="${PROJECT_DIR}/Scripts/6_Consensus_polish_array_$batch.sbatch"
echo -n "" > $script_polish_array # 清空文件
touch $script_polish_array #创建新文件 or 更新时间戳

# 输入文件内容
cat > $script_polish_array << \EOF
#!/bin/bash
#SBATCH -J02_racon_medaka_polish_GPU                            # Job name
#SBATCH --account=YOUR_ACCOUNT                 # charge account
#SBATCH -N1            # Number of nodes and cores per node required
#SBATCH --gres=gpu:V100:1                # GPU type (V100) and number of GPUs 
#SBATCH --mem-per-gpu=32G                        # Memory per core
#SBATCH -t3:00:00                              # Duration of the job (Ex: integrate:15 mins,D-HH:MM:SS for days, hours, minutes, and seconds requested)
#SBATCH -qinferno                               # QOS Name
#SBATCH -oReport-polish_drep%j.out                         # Combined output and error messages file
#SBATCH --mail-type=START,END,FAIL              # Mail preferencespac
#SBATCH --mail-user=your_email@example.com        # E-mail address for notifications
#SBATCH --array=1-13%4 # 最多同时运行4个barcode
# module load gcc
module load cuda
module load anaconda3/2022.05

# #---------------------------------------------------------------- Variables --------------------------------------------------------------- 
# echo "SLURM_ARRAY_TASK_ID:"$SLURM_ARRAY_TASK_ID
# br=barcode$(($SLURM_ARRAY_TASK_ID))
echo "SLURM_ARRAY_TASK_ID" $SLURM_ARRAY_TASK_ID
bn=$(printf "%02d" $(($SLURM_ARRAY_TASK_ID)))
br=barcode${bn}

echo $br
# br=barcode20
batch=1

############################################################## 0. consens UMI ############################################################## 

minLen=3500
maxLen=6000  
Q=10
CUTADAPT=cutadapt
BWA=bwa
SAMTOOLS=samtools
SEQTK=seqtk
USEARCH=${SOFTWARE_DIR}/usearch11.0.667_i86linux32
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
# 衡量 UMI bin 的覆盖度是否合理：如果 ratio 很小：说明这个 UMI bin 的 reads 很少 → 可能噪声或低覆盖；如果 ratio 很大：说明 UMI bin 过度占用 reads → 可能是 PCR 放大偏差或重复
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
# 从一个 UMI bin 的原始 reads 中自动挑选“代表性种子序列”，并通过多轮 racon 矫正生成高精度的 consensus（共识序列），用于后续分析。
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
  #   -sizeout -centroids $OUT/${UMINO}_centroids.fa # 旧的参数，会有差异很大的reads bin在一起
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
  OUT_NAME=${OUT##*/} # 只保留最后的文件名not full path
  echo $OUT_NAME
  OUT_NAME=${OUT_NAME%.*} # 去掉后缀
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
# conda create -n  rrUMI_2 python=3.12 minimap2 seqtk racon -y
# conda activate rrUMI_2
# conda install parallel cutadapt gawk samtools bwa -y
# pip install medaka==2.1.1
# pip install pyabpoa
# medaka --version

conda activate rrUMI_2

#
# # olishing of UMI consensus sequences (ONT data only). The racon-polished ONT consensus sequences were further polished individually by using all reads in each UMI bin and one round of Medaka (v.1.0.1) (https://github.com/nanoporetech/ medaka) with the commands Medaka mini_align -m and Medaka consensus --model r103_min_high_g360 --chunk_len 6,000 (ONT R10.3) or Medaka consensus --model r941_min_high_g360 --chunk_len 6,000 (ONT R9.4.1).
#     -c  File containing consensus sequences.CONSENSUS_FILE=$OPTARG;;
#     -m  Medaka model. MEDAKA_MODEL=$OPTARG;;
#     -l  Expected minimum chunk size. [Default = 6000]CHUNK_SIZE=$OPTARG;;
#     -d  Directory containing UMI read bins in the format
#         'umi*bins.fastq'. Recursive search. BINNING_DIR=$OPTARG;;
#     -o  Output directory.OUT_DIR=$OPTARG;;
#     -t  Number of threads to use.THREADS=$OPTARG;;
#     -n  Process n number of bins. If not defined all bins
#         are processed. SAMPLE=$OPTARG;;
#     -T  Number of Medaka jobs to run. [Default = 1].MEDAKA_JOBS=$OPTARG;;
#
#
### Medaka polishing assembly -------------------------------------------------
echo "----------------------------------------------1. Medaka polishing assembly -------------------------------------------------"
echo "-------------------------------------------------------------------------------------------------------------------------"
THREADS=32
MEDAKA_JOBS=16
# Format names
OUT_NAME=${OUT_DIR##*/}
echo "OUT_NAME:"$OUT_NAME
# Medaka jobs
MEDAKA_THREADS=$(( THREADS/MEDAKA_JOBS ))
CHUNK_SIZE=6000



# # Start medaka environment if relevant
# eval "$MEDAKA_ENV_START"


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

  # 把某个 UMI bin 的 reads (UMI_BIN) 比对到它的初始 consensus (UMI_NAME.fa)
  # 输出 .bam 文件，存在 $OUT_DIR/mapping/UMI_NAME/ 下

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



EOF
}
# rm Report*
cd ${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch
sbatch $script_polish_array  # 4028940 (3-13) 4029178 (1-2) 4030052 (2-13) 4030061(1) 4085943 (25 Gb) 4246381 (25Gb-failed) 4255405
squeue -u $USER
pace-quota

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














