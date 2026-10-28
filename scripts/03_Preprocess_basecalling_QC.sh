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

cd  ${LOCAL_DATA}/
# failed experiment was stored in sractch folder: ${PROJECT_DIR}/scratch_failed_run/

####################  install dorado in scratch  ################### 
# salloc -A YOUR_ACCOUNT -N1 --mem-per-gpu=24G -qinferno -t4:15:00 --gres=gpu:V100:1
# module load cuda
# mkdir ${SOFTWARE_DIR}
# cd ${SOFTWARE_DIR}
# wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.3.1-linux-x64.tar.gz
# tar -xf dorado-1.3.1-linux-x64.tar.gz
# dorado="${SOFTWARE_DIR}/dorado-1.3.1-linux-x64/bin/dorado"
# $dorado -v #1.3.1+7c84b01de
# mkdir ${SOFTWARE_DIR}/dorado-1.3.1-linux-x64/models
# cd ${SOFTWARE_DIR}/dorado-1.3.1-linux-x64/models
# $dorado download --model all


################### transfer data ################### 
batch=1
#------------------ batch 01 ------------------
# fist transfer Feb16 2026 5:31PM
# scp -r ${LOCAL_DATA}/52_PdNA_all_samples_02152026/no_sample_id/20260215_1711_MD-105372_FBF40675_925eb2a5/pod5/*.pod5 \
#   USER@CLUSTER_HOST:${PROJECT_DIR}/00_Pod5_$batch/
rsync -avP ${LOCAL_DATA}/52_PdNA_all_samples_02152026/no_sample_id/20260215_1711_MD-105372_FBF40675_925eb2a5/pod5/ \
  USER@CLUSTER_HOST:${PROJECT_DIR}/00_Pod5_$batch/

##########################################
#  first batch
##########################################

batch=1


cd  ${PROJECT_DIR}/
mkdir ${PROJECT_DIR}/00_SUP_$batch # all fastq files
mkdir ${PROJECT_DIR}/00_SUP_demux_$batch # all demultiplex files
mkdir ${PROJECT_DIR}/Scripts/
{
script_basecalling1="${PROJECT_DIR}/Scripts/1_2_Basecalling_demultiplexing_SUP_Q7.sbatch"
echo -n "" > $script_basecalling1 # 清空文件
touch $script_basecalling1 #创建新文件，更新时间戳
  # 输入文件内容
cat > $script_basecalling1 << \EOF
#!/bin/bash
#SBATCH -JbasecallingSUP_1                            # Job name
#SBATCH -A YOUR_ACCOUNT                            # Charge account
#SBATCH -N1 --gres=gpu:H200:2       # Number of nodes and GPUs required
#SBATCH --mem-per-gpu=142G                           # Memory per gpu
#SBATCH -t8:00:00                                        # Duration of the job (Ex: 15 mins)
#SBATCH -qinferno                                   # QOS name
#SBATCH -oReport-Superhigh-basecalling2%j.out                             # Combined output and error messages file
#SBATCH --mail-type=BEGIN,END,FAIL                  # Mail preferences
#SBATCH --mail-user=your_email@example.com            # e-mail address for notifications

module load cuda
# module load tensorflow-gpu/2.9.0 
# --gres=gpu:RTX_6000:3
# -N2  --gres=gpu:V100:2
# --gres=gpu:A100:2
#SBATCH -N2 --gres=gpu:H100:2        # Number of nodes and GPUs required
#Dorado is heavily-optimised for Nvidia A100 and H100 GPUs and will deliver maximal performance on systems with these GPUs.
#SBATCH -N2 --gres=gpu:V100:2        # Number of nodes and GPUs required

batch=1
# dorado="${SOFTWARE_DIR}/dorado-1.3.1-linux-x64/bin/dorado"
dorado="${SOFTWARE_DIR}/dorado-1.1.1-linux-x64/bin/dorado"

time $dorado \
    basecaller ${SOFTWARE_DIR}/dorado-1.3.1-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0  \
     ${PROJECT_DIR}/00_Pod5_$batch/ --recursive  --emit-fastq --min-qscore 8 \
    > ${PROJECT_DIR}/00_SUP_$batch/SUP_raw_$batch.fastq
   
time  $dorado demux --kit-name SQK-NBD114-24 --emit-fastq \
    --output-dir  ${PROJECT_DIR}/00_SUP_demux_$batch/ \
    ${PROJECT_DIR}/00_SUP_$batch/SUP_raw_$batch.fastq   

EOF
}
cd  ${PROJECT_DIR}/00_SUP_$batch
sbatch $script_basecalling1 # 4010185 (latest model 1.3.1 for basecalling) #4013663 (1.1.1 model for demux) #4014384 (1.1.1 for basecalling and demux) 4062695 (25Gb batch, model 1.1.1)
squeue -u $USER
pace-quota

#-------------------------------------------------- 1) rename demux file ----------------------------------------
# rename all demux files
{
  bar_list="01 02 03 04 05 06 07 08 09 10 11 12 13"
  for bn in $bar_list
  do
    echo $bn
    # mv ${PROJECT_DIR}/00_SUP_demux_1/no_sample/20260215_2211_0_FBF40675_925eb2a5/fastq_pass/barcode$bn/FBF40675_pass_barcode${bn}_925eb2a5_00000000_0.fastq \
    mv ${PROJECT_DIR}/00_SUP_demux_1/925eb2a5-3894-4aa8-b816-0c11fb7a0ae5_SQK-NBD114-24_barcode$bn.fastq \
    ${PROJECT_DIR}//00_SUP_demux_1/barcode${bn}.fastq
  done
}

#----------------------------------- 2) Quality check for 4TH batch SUP data -----------------------------------
salloc -q inferno -A YOUR_ACCOUNT -N4 --ntasks-per-node=8 --time=8:30:00 --mem=64G
module load anaconda3/2022.05
{
# quality check for raw reads
conda activate np
batch=1
# set name list (put duplicate together)
bar_list="01 02 03 04 05 06 07 08 09 10 11 12 13"
names=()
for bn in $bar_list; do
  names+=(${PROJECT_DIR}//00_SUP_demux_$batch/barcode${bn}.fastq)
done
echo ${names[@]}

# set output name set
ID=("PCE_0" "Methanol" "Glycerol_1" "Glycerol_2" "Glycerol_3" "Glycerol_4" "Glycerol_5" "Glycerol_6" "PCE_1" "PCE_2" "PCE_3" "PCE_4" "Blank")
echo ${ID[@]}
output=()
for id in "${ID[@]}"; do
  output+=("${id}")
done
echo "${output[@]}"

# set colors
colors=(
  "#8C33FF" 
  "#33FFF5"
  "#FF5733" "#FF5733" "#FF5733" "#FF5733" "#FF5733" "#FF5733"
  "#3357FF" "#3357FF" "#3357FF" "#3357FF"
  "#D3D3D3"
  # "#33FFF5" #"#33FFF5" "#33FFF5"
  # "#F5FF33" #"#F5FF33" "#F5FF33"
  #"#FF8C33" "#FF8C33" "#FF8C33"
  # "#33FF8C" "#33FF8C"
  # "#FF3333" "#FF3333"
  # "#3333FF" "#3333FF"
)
echo "${colors[@]}"

# run the analysis
# check the number of elements
echo "${#colors[@]}"
echo "${#ID[@]}"
echo "${#output[@]}"
echo "${#colors[@]}"

mkdir ${PROJECT_DIR}//00_SUP_NanoComp_$batch
NanoComp -t 24 --fastq ${names[@]} \
  --outdir  ${PROJECT_DIR}//00_SUP_NanoComp_$batch \
  -n ${output[@]} \
  -c ${colors[@]}
  
}


#----------------------------------- 3) Quality control and merge to all data -----------------------------------
# read quality and length
{
minLen=3500 
maxLen=6000  
Q=10
mkdir ${PROJECT_DIR}//00_SUP_demux_Q${Q}_${minLen}_${maxLen}_$batch
cd ${PROJECT_DIR}//00_SUP_demux_Q${Q}_${minLen}_${maxLen}_$batch

conda activate chopper
LOG_FILE="${PROJECT_DIR}//00_SUP_demux_Q${Q}_${minLen}_${maxLen}_$batch/chopper_summary.log"
> "$LOG_FILE"

for bn in $bar_list
do
  br=barcode$bn
  echo "=== Processing $br at $(date) ===" >> "$LOG_FILE"
  chopper -q ${Q} --minlength $minLen --maxlength $maxLen \
          -i ${PROJECT_DIR}//00_SUP_demux_$batch/$br.fastq \
      > ${PROJECT_DIR}//00_SUP_demux_Q${Q}_${minLen}_${maxLen}_$batch/${br}_QC.fastq \
      2>> "$LOG_FILE"
      
  echo "=== Completed $br at $(date) ===" >> "$LOG_FILE"
  echo "" >> "$LOG_FILE"
done

}
#----------------------------------- 4) Qquality check of filtered dataset-----------------------------------
# visualize the quality
conda activate np
{
minLen=3500 
maxLen=6000
Q=10
names=()
for bn in $bar_list; do
  br=barcode$bn
  echo $br
  names+=(${PROJECT_DIR}//00_SUP_demux_Q${Q}_${minLen}_${maxLen}_$batch/${br}_QC.fastq)
done
echo ${names[@]}

# set output name set
# set output name set
ID=("PCE_0" "Methanol" "Glycerol_1" "Glycerol_2" "Glycerol_3" "Glycerol_4" "Glycerol_5" "Glycerol_6" "PCE_1" "PCE_2" "PCE_3" "PCE_4" "Blank")
echo ${ID[@]}
output=()
for id in "${ID[@]}"; do
  output+=("${id}")
done
echo "${output[@]}"

# set colors
colors=(
  "#8C33FF" 
  "#33FFF5"
  "#FF5733" "#FF5733" "#FF5733" "#FF5733" "#FF5733" "#FF5733"
  "#3357FF" "#3357FF" "#3357FF" "#3357FF"
  "#D3D3D3"
  # "#33FFF5" #"#33FFF5" "#33FFF5"
  # "#F5FF33" #"#F5FF33" "#F5FF33"
  #"#FF8C33" "#FF8C33" "#FF8C33"
  # "#33FF8C" "#33FF8C"
  # "#FF3333" "#FF3333"
  # "#3333FF" "#3333FF"
)
echo "${colors[@]}"

# run the analysis
# check the number of elements
echo "${#colors[@]}"
echo "${#ID[@]}"
echo "${#output[@]}"
echo "${#colors[@]}"


mkdir ${PROJECT_DIR}//00_SUP_QC_Q${Q}_${minLen}_${maxLen}_NanoComp_$batch
NanoComp --fastq ${names[@]} \
  --outdir  ${PROJECT_DIR}//00_SUP_QC_Q${Q}_${minLen}_${maxLen}_NanoComp_$batch \
  -n ${output[@]} \
  -c ${colors[@]} \
  -t 24
}



