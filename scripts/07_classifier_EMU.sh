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
salloc -q inferno -A YOUR_ACCOUNT -N4 --ntasks-per-node=8 --time=12:30:00 --mem=64G
module load anaconda3/2022.05

# pip install osfclient
# export EMU_DATABASE_DIR="${SOFTWARE_DIR}/emu/emu_database"
# cd ${EMU_DATABASE_DIR}
# osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar
# tar -xvf emu.tar
# 
# conda create --name emu python=3.7 
# conda activate emu
# 
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda install emu
# conda install taxonkit
# # generate seq2taxid.map --by onDemand R - 09_RDP_classifier_generation
# META_TABLE="${GROND_DIR}/01_Github_file/GROND/Grond_GTDB_207/sub_Meta_region_genome/META_nospace.tsv"
# head $META_TABLE
# db_genome_16S23S=${PROJECT_DIR}/database_genome_16S23S


# Install RDP and training
# 1. generate files on R
# 2. extract fasta into the rdp
#!/bin/bash
# conda create -n rdp -c bioconda rdp_classifier

############################# ----------------------- ############################# 
############################# RDP classifier - 16S23S ############################# 
############################# ----------------------- ############################# 
# 设置路径
RAW_FASTA="${PROJECT_DIR}/database_genome_16S23S/db_genome_16S23S_all.fa"  # 你的原始序列文件
META_MAP="${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_16S23S_Classifier/rdp_headers.map" # 由R生成
head $META_MAP #16S_23S_6       Root;Bacteria;Pseudomonadota;Gammaproteobacteria;Xanthomonadales;Xanthomonadaceae;Xanthomonas;Xanthomonas_oryzae
head $RAW_FASTA #>16S_23S_1421 rrn_product_1421 AGAGTTTGATCCTGGCTCAGGATGAATGCTGGCGGCGTGCCTAACACATGCAAGTCGAGCGGTAACAGGTCTAGTTTACT
# 
# # 3. 使用 seqkit 根据 map 修改 FASTA Header
# # 将 Header 修改为: >ID <Tab> Root;Domain;...
# conda activate seqkit
# echo "Updating FASTA headers for RDP..."
# 
# awk -F'\t' 'NR==FNR{a[$1]=$2; next} {if(/^>/){split($1, id_part, " "); id=substr(id_part[1],2); if(id in a){print ">"id"\t"a[id]}else{print $0}}else{print $0}}' ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_16S23S_Classifier/rdp_headers.map $RAW_FASTA > ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_16S23S_Classifier/rdp_train.fasta
# head ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_16S23S_Classifier/rdp_train.fasta
# 
# # 4. 训练 RDP 模型 (创建自定义分类器)
# # 注意：RDP 训练需要特定的目录结构或使用 train 命令
# echo "Training RDP Classifier..."
# conda activate rdp
# mkdir ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_16S23S_Classifier/my_custom_model
# export _JAVA_OPTIONS="-Xmx100g"
# rdp_classifier train -o ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_16S23S_Classifier/my_custom_model/ \
#                     -s ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_16S23S_Classifier/rdp_train.fasta \
#                     -t ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_16S23S_Classifier/taxid.txt
# 
# cat > ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_16S23S_Classifier/my_custom_model/my_train.properties << 'EOF'
# classifierVersion=2.14
# bergeyTree=bergeyTrainingTree.xml
# probabilityList=genus_wordConditionalProbList.txt
# wordPrior=logWordPrior.txt
# probabilityIndex=wordConditionalProbIndexArr.txt
# EOF

# 5. classifer
conda activate rdp
minLen=3500
maxLen=6000  
Q=10
batch=1
bar_list="01 02 03 04 05 06 07 08 09 10 11 12"
region=16S23S
export _JAVA_OPTIONS="-Xmx100g"
for bn in $bar_list
do
  seq_dir=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1
  echo $bn
  rdp_classifier classify \
      -t ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_16S23S_Classifier/my_custom_model/my_train.properties \
     -o ${PROJECT_DIR}/04_rdp/rdp_classification_results_ccs_${region}_$bn.txt \
      -q   $seq_dir/consensus_barcode${bn}_stitch_${region}.fasta 
done
less ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/rdp_classification_results_16S23S_$bn.txt 


# compared to rrn results 
# attention: cannot directly use raw ccs reads, they have the universal primers. need to extract rrn region or chopper the universal primers
############################# ----------------------- ############################# 
############################# RDP classifier - ccs rrn.   ############################# 
############################# ----------------------- ############################# 
# 设置路径
RAW_FASTA="${PROJECT_DIR}/database_genome_rrn/db_genome_rrn_all.fa"  # 你的原始序列文件
META_MAP="${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_rrn_Classifier/rdp_headers.map"
head $RAW_FASTA # >rrn_product_1421 seq GCA_001029695__CP011211.1:filter(unmasked) start 846363 length 6424 AGAGTTTGATCCTGGCTCAGGATGAATGCTGGCGGCGTGCCTAACACATGCAAGTCGAGCGGTAACAGGTCTAGTTTACT
head $META_MAP #head $META_MAP rrn_product_13  Root;Bacteria;Fusobacteriota;Fusobacteriia;Fusobacteriales;Fusobacteriaceae;Fusobacterium;Fusobacterium_nucleatum

# 3. 使用 seqkit 根据 map 修改 FASTA Header
# 将 Header 修改为: >ID <Tab> Root;Domain;...
conda activate seqkit
echo "Updating FASTA headers for RDP..."
awk -F'\t' 'NR==FNR{a[$1]=$2; next} {if(/^>/){split($1, id_part, " "); id=substr(id_part[1],2); if(id in a){print ">"id"\t"a[id]}else{print $0}}else{print $0}}' $META_MAP $RAW_FASTA > ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_rrn_Classifier/rdp_train.fasta
head ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_rrn_Classifier/rdp_train.fasta


# 4. 训练 RDP 模型 (创建自定义分类器)
# 注意：RDP 训练需要特定的目录结构或使用 train 命令
echo "Training RDP Classifier..."
conda activate rdp
mkdir ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_rrn_Classifier/my_custom_model
export _JAVA_OPTIONS="-Xmx100g"
rdp_classifier train -o ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_rrn_Classifier/my_custom_model/ \
                    -s ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_rrn_Classifier/rdp_train.fasta \
                    -t ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_rrn_Classifier/taxid.txt

cat > ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_rrn_Classifier/my_custom_model/my_train.properties << 'EOF'
classifierVersion=2.14
bergeyTree=bergeyTrainingTree.xml
probabilityList=genus_wordConditionalProbList.txt
wordPrior=logWordPrior.txt
probabilityIndex=wordConditionalProbIndexArr.txt
EOF
ll -h ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_rrn_Classifier/my_custom_model/
# 5. classifer
conda activate rdp
# bn=
minLen=3500
maxLen=6000  
Q=10
batch=1
bar_list="01 02 03 04 05 06 07 08 09 10 11 12"
region=rrn
for bn in $bar_list
do
  seq_dir=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1
  echo $bn
  ll -h $seq_dir/consensus_barcode${bn}_stitch_${region}.fasta 
  rdp_classifier classify \
      -t ${PROJECT_DIR}/03_coverm_Q${Q}_${minLen}_${maxLen}_$batch/RDP_${region}_Classifier/my_custom_model/my_train.properties \
      -o ${PROJECT_DIR}/04_rdp/rdp_classification_results_ccs_${region}_$bn.txt \
      -q   $seq_dir/consensus_barcode${bn}_stitch_${region}.fasta 
done
less ${PROJECT_DIR}/04_rdp/rdp_classification_results_ccs_${region}_$bn.txt


############################# ----------------------- ############################# 
############################# RDP classifier - 16S/ITS/23S ########################
############################# ----------------------- ############################# 
# 设置路径

region_list="16S ITS 23S"
# region=16S
for region in $region_list
do
  echo $region
  RAW_FASTA="${PROJECT_DIR}/database_genome_${region}/db_genome_${region}_all.fa"  # 你的原始序列文件
  META_MAP="${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/rdp_headers.map" # 由R生成
  head $META_MAP 
  head $RAW_FASTA 
  
  # 3. 使用 seqkit 根据 map 修改 FASTA Header
  # 将 Header 修改为: >ID <Tab> Root;Domain;...
  conda activate seqkit
  echo "Updating FASTA headers for RDP..."
  
  awk -F'\t' 'NR==FNR{a[$1]=$2; next} {if(/^>/){split($1, id_part, " "); id=substr(id_part[1],2); if(id in a){print ">"id"\t"a[id]}else{print $0}}else{print $0}}' ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/rdp_headers.map $RAW_FASTA > ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/rdp_train.fasta
  head ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/rdp_train.fasta
  
  # 4. 训练 RDP 模型 (创建自定义分类器)
  # 注意：RDP 训练需要特定的目录结构或使用 train 命令
  echo "Training RDP Classifier..."
  conda activate rdp
  mkdir ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/my_custom_model
  export _JAVA_OPTIONS="-Xmx100g"
  rdp_classifier train -o ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/my_custom_model/ \
                      -s ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/rdp_train.fasta \
                      -t ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/taxid_${region}.txt
  
cat > ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/my_custom_model/my_train.properties << 'EOF'
classifierVersion=2.14
bergeyTree=bergeyTrainingTree.xml
probabilityList=genus_wordConditionalProbList.txt
wordPrior=logWordPrior.txt
probabilityIndex=wordConditionalProbIndexArr.txt
EOF
done

# 5. classifer
conda activate rdp
for region in $region_list
do
  echo $region
  echo "Running Classifier......"
  minLen=3500
  maxLen=6000  
  Q=10
  batch=1
  bar_list="01 02 03 04 05 06 07 08 09 10 11 12"
  region=${region}
  export _JAVA_OPTIONS="-Xmx100g"
  
  for bn in $bar_list
  do
    seq_dir=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1
    echo $bn
    rdp_classifier classify \
        -t ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/my_custom_model/my_train.properties \
        -o ${PROJECT_DIR}/04_rdp/rdp_classification_results_ccs_${region}_$bn.txt \
        -q   $seq_dir/consensus_barcode${bn}_stitch_${region}.fasta 
  done
done  
  


############################# ----------------------- ############################# 
############################# RDP classifier - V3V4   ############################# 
############################# ----------------------- ############################# 
salloc -q inferno -A YOUR_ACCOUNT -N4 --ntasks-per-node=8 --time=5:30:00 --mem=64G
module load anaconda3/2022.05

# 1. Extract V3V4 region from rrn 
# CCTACGGGNBGCWSCAG
# GACTACNVGGGTATCTAATCC
input_file_V3V4="${GROND_DIR}/02_Database/input_rrn_full_V3V4.txt"
id="V3V4"
FP_V3V4="CCTACGGGNBGCWSCAG"
RP_V3V4="GACTACNVGGGTATCTAATCC"
min_len="400"
max_len="500"
echo "$id $FP_V3V4 $RP_V3V4 $min_len $max_len" > $input_file_V3V4
head $input_file_V3V4


conda activate ipcress
mismatch=2
minLen=3500
maxLen=6000  
Q=10
batch=1
bar_list="01 02 03 04 05 06 07 08 09 10 11 12"
region=V3V4

for bn in $bar_list
do
  echo $bn
    raw_ccs=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/consensus_barcode${bn}_stitch.fa
    ccs_V3V4=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/consensus_barcode${bn}_stitch_${region}.fa
    OUT_DIR=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/
    echo $region
    echo "Mismatch:" $mismatch
    # input_primer_file=${GROND_DIR}/02_Database/input_rrn_full_${region}_multi.txt
    input_primer_file=${GROND_DIR}/02_Database/input_rrn_full_${region}.txt
    echo "IPCRess file: "  $input_primer_file
    ipcress $input_primer_file $raw_ccs \
      -m $mismatch -M 64 -S 106 -p F -P \
      > $ccs_V3V4   #fasta
    sed -i '/^ipcress/d' $ccs_V3V4
    sed -i '/^-- completed ipcress analysis/d' $ccs_V3V4
done
grep -c ">" $ccs_V3V4 #931
grep -c ">" consensus_barcode12_stitch_16S23S.fasta # 907
grep -c ">" consensus_barcode12_stitch_rrn.fasta # 931


# 2. Generate Classifier
region=V3V4
echo $region
mkdir ${PROJECT_DIR}/database_genome_${region}/
RAW_FASTA="${PROJECT_DIR}/database_genome_${region}/db_genome_${region}_all.fa"  # 你的原始序列文件
# cp ${GROND_DIR}/01_Github_file/GROND/Grond_GTDB_207/Regions_complete/extraction_V3V4.fasta $RAW_FASTA
META_MAP="${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/rdp_headers.map" # 由R生成
head $META_MAP 
head $RAW_FASTA 
  
# 3. 使用 seqkit 根据 map 修改 FASTA Header
  # 将 Header 修改为: >ID <Tab> Root;Domain;...
conda activate seqkit
echo "Updating FASTA headers for RDP..."
awk -F'\t' 'NR==FNR{a[$1]=$2; next} {if(/^>/){split($1, id_part, " "); id=substr(id_part[1],2); if(id in a){print ">"id"\t"a[id]}else{print $0}}else{print $0}}' ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/rdp_headers.map $RAW_FASTA > ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/rdp_train.fasta
head ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/rdp_train.fasta
  
# 4. 训练 RDP 模型 (创建自定义分类器)
# 注意：RDP 训练需要特定的目录结构或使用 train 命令
echo "Training RDP Classifier..."
conda activate rdp
mkdir ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/my_custom_model
export _JAVA_OPTIONS="-Xmx100g"
rdp_classifier train -o ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/my_custom_model/ \
                      -s ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/rdp_train.fasta \
                      -t ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/taxid_${region}.txt
  
cat > ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/my_custom_model/my_train.properties << 'EOF'
classifierVersion=2.14
bergeyTree=bergeyTrainingTree.xml
probabilityList=genus_wordConditionalProbList.txt
wordPrior=logWordPrior.txt
probabilityIndex=wordConditionalProbIndexArr.txt
EOF


# 5. Classify
region=V3V4
echo $region
echo "Running Classifier......"
minLen=3500
maxLen=6000  
Q=10
batch=1
bar_list="01 02 03 04 05 06 07 08 09 10 11 12"
export _JAVA_OPTIONS="-Xmx100g"
  
for bn in $bar_list
do
    seq_dir=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1
    echo $bn
    rdp_classifier classify \
        -t ${PROJECT_DIR}/03_coverm_Q10_3500_6000_1/RDP_${region}_Classifier/my_custom_model/my_train.properties \
        -o ${PROJECT_DIR}/04_rdp/rdp_classification_results_ccs_${region}_$bn.txt \
        -q   $seq_dir/consensus_barcode${bn}_stitch_${region}.fa

done



