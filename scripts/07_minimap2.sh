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
salloc -q inferno -A YOUR_ACCOUNT -N4 --ntasks-per-node=8 --time=4:30:00 --mem=64G
module load anaconda3/2022.05
pace-quota

# coverm genome cannot work, it will take the genus file as a whole, not single read 

############### ############### ###############   
  ############### minimap2 - each read level ############### 
############### ############### ############### 
batch=1
# bar_list="05 06 07 08 09 10 11 12 13"
bar_list="01 02 03 04 05 06 07 08 09 10 11 12"

minLen=3500
maxLen=6000  
Q=10

#------------- 16S_23S 
#-------- 0. prepare database (need to merge all fna files into one fa file for minimap2)
# 16S23S
genus_16S23S=${GROND_DIR}/01_Github_file/GROND/Grond_GTDB_207/Regions_complete/split_by_genus_16S_23S/
db_genus_16S23S=${PROJECT_DIR}/database_genus_16S23S
db_genome_16S23S=${PROJECT_DIR}/database_genome_16S23S
# # mkdir $db_genome_16S23S
# mkdir $db_genus_16S23S
# for f in "$genus_16S23S"/*.fasta; do
#     # 提取文件名（不带路径和后缀）
#     base=$(basename "$f" .fasta)
#     # 复制到新目录并改名为 .fna
#     cp "$f" "$db_genus_16S23S/${base}.fna"
# done
# cat $db_genus_16S23S/*.fna > $db_genome_16S23S/db_genome_16S23S_all.fa

#----- 1. prepare ITS omitting reads
conda activate ipcress
mismatch=2

for bn in $bar_list
do
  echo $bn
  raw_ccs=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/consensus_barcode${bn}_stitch.fa
  # ccs_16S_23S=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/consensus_barcode${bn}_stitch_16S_23S.fa
  OUT_DIR=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/
  # for region in "16S" "23S" "ITS" "rrn"
  # # for region in "rrn"
  # do
  #   echo $region
  #   echo "Mismatch:" $mismatch
  #   input_primer_file=${GROND_DIR}/02_Database/input_rrn_full_${region}_multi.txt
  #   echo "IPCRess file: "  $input_primer_file
  #   ipcress $input_primer_file $raw_ccs \
  #     -m $mismatch -M 64 -S 106 -p F -P \
  #     > $OUT_DIR/consensus_barcode${bn}_stitch_${region}.fasta   #fasta
  #   sed -i '/^ipcress/d' $OUT_DIR/consensus_barcode${bn}_stitch_${region}.fasta
  #   sed -i '/^-- completed ipcress analysis/d' $OUT_DIR/consensus_barcode${bn}_stitch_${region}.fasta
  # done

  for region in "16S" "23S" "ITS" "rrn" "V3V4"
  # for region in "V3V4"
  do
    echo $region
    awk 'BEGIN{
      OFS="\t";
      print "rrn_number","accession","start_position","length"
    }
    # header example:
    # >rrn_product_1 seq GCA_000017645__CP000781.1:filter(unmasked) start 1013944 length 4839
    /^>/{
      rrn="NA"; acc="NA"; st="NA"; ln="NA";
      if (match($0, /^>([^[:space:]]+)/, m))              rrn=m[1];
      if (match($0, /seq[[:space:]]+([^:]+):/, m))        acc=m[1];
      if (match($0, /start[[:space:]]+([0-9]+)/, m))      st=m[1];
      if (match($0, /length[[:space:]]+([0-9]+)/, m))     ln=m[1];
      print rrn, acc, st, ln;
    }' $OUT_DIR/consensus_barcode${bn}_stitch_${region}.fasta > $OUT_DIR/Summary_extraction_${region}.tsv
  done

# merge 16S_23S together-- order of 16S and 23S is according to the physical position
awk '
/^>/{
    # 提取 UMI 编号
    if (match($0, /umi([0-9]+)bins/, m)) {
        u_num = m[1]
    } else { next }

    # 提取 start 坐标
    pos = ""
    if (match($0, /start ([0-9]+)/, s)) {
        pos = s[1] + 0
    }

    # 判断类型
    type = ($0 ~ /16S/) ? "16S" : "23S"

    # 记录信息
    start_pos[u_num, type] = pos
    full_header[u_num, type] = $0
    found_umi[u_num] = 1
    if (u_num > max) max = u_num
    next
}

{
    # 累加序列
    gsub(/[ \t\r\n]/, "", $0)
    seqs[u_num, type] = seqs[u_num, type] $0
}

END {
    for (i = 0; i <= max; i++) {
        # 只有 16S 和 23S 都存在时才输出
        if (found_umi[i] && seqs[i, "16S"] != "" && seqs[i, "23S"] != "") {
            s16 = start_pos[i, "16S"]
            s23 = start_pos[i, "23S"]

            # 固定 16S 在前，23S 在后
            final_seq = seqs[i, "16S"] seqs[i, "23S"]

            printf ">16S_23S_merged_umi%d_order_16S_then_23S 16S_start:%s 23S_start:%s\n", i, s16, s23
            print final_seq
        }
    }
}' \
"$OUT_DIR/consensus_barcode${bn}_stitch_16S.fasta" \
"$OUT_DIR/consensus_barcode${bn}_stitch_23S.fasta" \
> "$OUT_DIR/consensus_barcode${bn}_stitch_16S23S.fasta"
done
# summarize 16S23S region
# >16S_23S_merged_umi27_order_16S_then_23S 16S_start:3200 23S_start:42  (example header)
bar_list="01 02 03 04 05 06 07 08 09 10 11 12"
minLen=3500
maxLen=6000  
Q=10
batch=1
for bn in $bar_list
do
  wkdir=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/
  cd "$wkdir"
  echo "$bn"

  awk '
  BEGIN {
      FS="\t"; OFS="\t";
      print "UMI_ID", "ubs", "rrn_header", "16S_header", "23S_header", "ITS_header", "16S23S_header", "V3V4_header";
  }

  # 跳过各个 summary 文件表头
  $1 ~ /_number$/ { next }

  {
      # 提取 umi 编号，例如 umi3499bins -> umi3499
      if (match($2, /umi[0-9]+/)) {
          uid = substr($2, RSTART, RLENGTH);
      } else {
          next
      }

      # 提取 ubs
      if (match($2, /ubs=[0-9]+/)) {
          split(substr($2, RSTART, RLENGTH), a, "=");
          ubs_val[uid] = a[2];
      }

      # 根据文件名写入对应列
      if (FILENAME ~ /Summary_extraction_rrn.tsv$/)      rrrn[uid]    = $1;
      else if (FILENAME ~ /Summary_extraction_16S.tsv$/) r16s[uid]    = $1;
      else if (FILENAME ~ /Summary_extraction_23S.tsv$/) r23s[uid]    = $1;
      else if (FILENAME ~ /Summary_extraction_ITS.tsv$/) rits[uid]    = $1;
      else if (FILENAME ~ /Summary_extraction_16S23S.tsv$/) r16s23s[uid] = $1;
      else if (FILENAME ~ /Summary_extraction_V3V4.tsv$/) rv3v4[uid]  = $1;

      uids[uid] = 1;
  }

  END {
      for (id in uids) {
          print id, \
                (id in ubs_val   ? ubs_val[id]   : "N/A"), \
                (id in rrrn      ? rrrn[id]      : "N/A"), \
                (id in r16s      ? r16s[id]      : "N/A"), \
                (id in r23s      ? r23s[id]      : "N/A"), \
                (id in rits      ? rits[id]      : "N/A"), \
                (id in r16s23s   ? r16s23s[id]   : "N/A"), \
                (id in rv3v4     ? rv3v4[id]     : "N/A");
      }
  }
  ' Summary_extraction_rrn.tsv \
    Summary_extraction_16S.tsv \
    Summary_extraction_23S.tsv \
    Summary_extraction_ITS.tsv \
    Summary_extraction_16S23S.tsv \
    Summary_extraction_V3V4.tsv > Summary_merged_regions.tsv
done



# merge all summary table together
for bn in $bar_list
do
  wkdir=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/
  cd $wkdir
  echo $bn
  awk '
  BEGIN {
      FS="\t"; OFS="\t";
      # 打印总表表头
      print "UMI_ID", "ubs", "16S_header", "23S_header", "ITS_header", "16S23S_header", "V3V4_header";
  }
  
  # 忽略所有原始文件的表头行
  $1 ~ /rrn_number/ { next }
  
  {
      # 1. 提取 UMI 数字标识 (例如 从 umi3499bins 提取 umi3499)
      match($2, /umi[0-9]+/);
      uid = substr($2, RSTART, RLENGTH);
  
      # 2. 提取 ubs 数值 (如果存在)
      if (match($2, /ubs=[0-9]+/)) {
          split(substr($2, RSTART, RLENGTH), a, "=");
          ubs_val[uid] = a[2];
      }
  
      # 3. 根据文件名存入对应的列 (ARGIND 是当前处理的文件序号)
      if (ARGIND == 1) r16s[uid] = $1;
      if (ARGIND == 2) r23s[uid] = $1;
      if (ARGIND == 3) rits[uid] = $1;
      if (ARGIND == 4) r16s23s[uid] = $1;
      if (ARGIND == 5) rv3v4[uid] = $1;
  
      # 记录所有出现过的 UMI ID
      uids[uid] = 1;
  }
  
  END {
      for (id in uids) {
          # 如果某个 UMI 在某个区域缺失，则用 "N/A" 填充
          print id, \
                (ubs_val[id] ? ubs_val[id] : "N/A"), \
                (r16s[id] ? r16s[id] : "N/A"), \
                (r23s[id] ? r23s[id] : "N/A"), \
                (rits[id] ? rits[id] : "N/A"), \
                (r16s23s[id] ? r16s23s[id] : "N/A"), \
                (rv3v4[id] ? rv3v4[id] : "N/A");
      }
  }
  ' Summary_extraction_16S.tsv \
    Summary_extraction_ITS.tsv \
    Summary_extraction_23S.tsv \
    Summary_extraction_16S23S.tsv \
    Summary_extraction_V3V4.tsv > Summary_merged_regions.tsv
done 
head Summary_merged_regions.tsv 


# check merged file lines
bn=12
OUT_DIR=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/

grep -c ">" $OUT_DIR/consensus_barcode${bn}_stitch_16S.fasta #933 multi_mm2:931 
grep -c ">" $OUT_DIR/consensus_barcode${bn}_stitch_ITS.fasta #1361 multi_mm2:807
grep -c ">" $OUT_DIR/consensus_barcode${bn}_stitch_23S.fasta #1335 multi_mm2:907
grep -c ">" $OUT_DIR/consensus_barcode${bn}_stitch_rrn.fasta #927 multi_mm2:931
grep -c ">" $OUT_DIR/consensus_barcode${bn}_stitch_16S23S.fasta #929 multi_mm2:907

head $OUT_DIR/consensus_barcode${bn}_stitch_16S.fasta
head $OUT_DIR/consensus_barcode${bn}_stitch_23S.fasta
head $OUT_DIR/consensus_barcode${bn}_stitch_rrn.fasta
head $OUT_DIR/consensus_barcode${bn}_stitch_16S23S.fasta
head $OUT_DIR/consensus_barcode${bn}_stitch_V3V4.fa

# summarize# extracted reads
csv_out=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/summary_read_count.tsv
echo "bn,stitch_16S,stitch_ITS,stitch_23S,stitch_rrn,stitch_16S23S" > "$csv_out"

for bn in $bar_list
do
    OUT_DIR=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/
    f16S="$OUT_DIR/consensus_barcode${bn}_stitch_16S.fasta"
    fITS="$OUT_DIR/consensus_barcode${bn}_stitch_ITS.fasta"
    f23S="$OUT_DIR/consensus_barcode${bn}_stitch_23S.fasta"
    frrn="$OUT_DIR/consensus_barcode${bn}_stitch_rrn.fasta"
    f16S23S="$OUT_DIR/consensus_barcode${bn}_stitch_16S23S.fasta"

    c16S=$( [ -f "$f16S" ] && grep -c ">" "$f16S" || echo "NA" )
    cITS=$( [ -f "$fITS" ] && grep -c ">" "$fITS" || echo "NA" )
    c23S=$( [ -f "$f23S" ] && grep -c ">" "$f23S" || echo "NA" )
    crrn=$( [ -f "$frrn" ] && grep -c ">" "$frrn" || echo "NA" )
    c16S23S=$( [ -f "$f16S23S" ] && grep -c ">" "$f16S23S" || echo "NA" )

    echo "${bn},${c16S},${cITS},${c23S},${crrn},${c16S23S}" >> "$csv_out"
done



# if merged sequence << 16S\23S, lets check
# 提取 16S 的所有 UMI 编号
grep -oP "umi[0-9]+bins" $OUT_DIR/consensus_barcode${bn}_stitch_16S.fasta | sort -u > 16S_umi.list
# 提取 23S 的所有 UMI 编号
grep -oP "umi[0-9]+bins" $OUT_DIR/consensus_barcode${bn}_stitch_23S.fasta | sort -u > 23S_umi.list
# 计算交集数量
comm -12 16S_umi.list 23S_umi.list | wc -l #929
comm -23 16S_umi.list 23S_umi.list | head -n 20  # 只出现在16S umi1154bins umi938bins
comm -23 23S_umi.list 16S_umi.list | head -n 20  # 只出现在23S umi1154bins umi938bins
cat 16S_umi.list 23S_umi.list | sort -u | wc -l # 所有特异性read #933

# 拽出原始的reads
grep -A 1 "umi1062bins" ${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/consensus_barcode${bn}_stitch.fa > umi1062_23Sread.fasta
grep -A 1 "umi1062bins" $OUT_DIR/consensus_barcode${bn}_stitch_16S.fasta

#------ 2. coverm - 16S23S
conda activate np23
identity=99
coverage=99
mkdir ${PROJECT_DIR}/03_coverm_Q${Q}_${minLen}_${maxLen}_$batch
for bn in $bar_list
do
  echo $bn
  OUT_DIR=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/
  coverm contig --single $OUT_DIR/consensus_barcode${bn}_stitch_16S23S.fasta \
    --reference $db_genome_16S23S/db_genome_16S23S_all.fa  \
    --bam-file-cache-directory ${PROJECT_DIR}/03_coverm_Q${Q}_${minLen}_${maxLen}_${batch} \
    -m mean length covered_bases \
    -o  ${PROJECT_DIR}/03_coverm_Q${Q}_${minLen}_${maxLen}_$batch/test_16S23S_Mutli_${bn}_$identity$coverage.tsv \
    -p minimap2-ont \
    --min-covered-fraction 0 \
    --min-read-percent-identity $identity \
    --min-read-aligned-percent $coverage \
    -t 32  
  # summarize taxonomic information
  # 设置变量
  COVERM_RESULT="${PROJECT_DIR}/03_coverm_Q${Q}_${minLen}_${maxLen}_$batch/test_16S23S_Mutli_${bn}_$identity$coverage.tsv"
  META_TABLE="${GROND_DIR}/01_Github_file/GROND/Grond_GTDB_207/sub_Meta_region_genome/META_nospace.tsv"
gawk -F'\t' '
    # 1. 读取 META 表，建立 16S23S 到物种信息的映射
    (NR==FNR) {
        if(NR==1) next;
        # 索引 2 为 16S_23S_product_xxx，存储对应的 Accession, Genus, Species
        # $1=GCA_id, $24=Genus, $25=Species
        meta[$2] = $1 "\t" $24 "\t" $25;
        next;
    }

    # 2. 处理 CoverM 结果表 (排除表头)
    ($1 ~ /^16S_23S_product_/) {
        # 如果该 rrn_product 在 META 表中存在且有覆盖 (第四列 > 0)
        if ($1 in meta) {
            # 以 Accession+Genus+Species 作为 Key 累加覆盖度
            # 这里累加第四列 Covered Bases ($4)
            tax_info = meta[$1];
            sum_coverage[tax_info] += $4;
            # 如果想累加 Mean 丰度，可以改为 += $2
            sum_mean[tax_info] += $2;
        }
    }

    # 3. 输出结果
    END {
        print "Accession\tGenus\tSpecies\tTotal_Covered_Bases\tTotal_Mean_Abundance";
        for (t in sum_coverage) {
            if (sum_coverage[t] > 0) {
                print t "\t" sum_coverage[t] "\t" sum_mean[t];
            }
        }
    }
' "$META_TABLE" "$COVERM_RESULT" > ${PROJECT_DIR}/03_coverm_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_genus_species_summary_16S23S.tsv
done
# 
# 样本编号 (Barcode),Mapped Reads,Total Reads,Mapping Rate (%)
# 01,85,940,9.04%
# 02,14,231,6.06%
# 03,49,"1,027",4.77%
# 04,55,981,5.61%
# 05,137,"2,083",6.58%
# 06,152,"1,476",10.30%
# 07,113,"1,150",9.83%
# 08,82,"1,155",7.10%
# 09,15,482,3.11%
# 10,45,"1,268",3.55%
# 11,32,815,3.93%
# 12,35,907,3.86%
# 总计 (Total),814,"12,515",6.50%

#################################### full length rrn ####################################
# rrn database
genus_rrn=${GROND_DIR}/01_Github_file/GROND/Grond_GTDB_207/Regions_complete/split_by_genus_rrn/
db_genus_rrn=${PROJECT_DIR}/database_genus_rrn
db_genome_rrn=${PROJECT_DIR}/database_genome_rrn
# mkdir $db_genome_rrn
# mkdir $db_genus_rrn
# for f in "$genus_rrn"/*.fasta; do
#     # 提取文件名（不带路径和后缀）
#     base=$(basename "$f" .fasta)
#     # 复制到新目录并改名为 .fna
#     cp "$f" "$db_genus_rrn/${base}.fna"
# done
# cat $db_genus_rrn/*.fna > $db_genome_rrn/db_genome_rrn_all.fa


#--------- 1. coverm 
conda activate np23
identity=95
coverage=95
# rrn minimap2
for bn in $bar_list
do
  echo $bn
  OUT_DIR=${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/
  coverm contig --single $OUT_DIR/consensus_barcode${bn}_stitch_rrn.fasta \
    --reference $db_genome_rrn/db_genome_rrn_all.fa  \
    --bam-file-cache-directory ${PROJECT_DIR}/03_coverm_Q${Q}_${minLen}_${maxLen}_$batch \
    -m mean length covered_bases \
    -o  ${PROJECT_DIR}/03_coverm_Q${Q}_${minLen}_${maxLen}_$batch/rrn_Mutli_${bn}_$identity$coverage.tsv \
    -p minimap2-ont \
    --min-covered-fraction 0 \
    --min-read-percent-identity $identity \
    --min-read-aligned-percent $coverage \
    -t 32  
  # summarize taxonomic information
  # 设置变量
  COVERM_RESULT="${PROJECT_DIR}/03_coverm_Q${Q}_${minLen}_${maxLen}_$batch/rrn_Mutli_${bn}_$identity$coverage.tsv"
  META_TABLE="${GROND_DIR}/01_Github_file/GROND/Grond_GTDB_207/sub_Meta_region_genome/META_nospace.tsv"
gawk -F'\t' '
    # 1. 读取 META 表，建立 rrn 到物种信息的映射
    (NR==FNR) {
        if(NR==1) next;
        # 索引 2 为 rrn_product_xxx，存储对应的 Accession, Genus, Species
        # $1=GCA_id, $24=Genus, $25=Species
        meta[$2] = $1 "\t" $24 "\t" $25;
        next;
    }

    # 2. 处理 CoverM 结果表 (排除表头)
    ($1 ~ /^rrn_product_/) {
        # 如果该 rrn_product 在 META 表中存在且有覆盖 (第四列 > 0)
        if ($1 in meta) {
            # 以 Accession+Genus+Species 作为 Key 累加覆盖度
            # 这里累加第四列 Covered Bases ($4)
            tax_info = meta[$1];
            sum_coverage[tax_info] += $4;
            # 如果想累加 Mean 丰度，可以改为 += $2
            sum_mean[tax_info] += $2;
        }
    }

    # 3. 输出结果
    END {
        print "Accession\tGenus\tSpecies\tTotal_Covered_Bases\tTotal_Mean_Abundance";
        for (t in sum_coverage) {
            if (sum_coverage[t] > 0) {
                print t "\t" sum_coverage[t] "\t" sum_mean[t];
            }
        }
    }
' "$META_TABLE" "$COVERM_RESULT" > ${PROJECT_DIR}/03_coverm_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_genus_species_summary_rrn_extracted_$identity$coverage.tsv

done

#################################### 16S/23S/ITS ####################################
# 1. Generate database for each region
region_list="16S ITS 23S"
for region in $region_list
do
  echo $region
  genus_region=${GROND_DIR}/01_Github_file/GROND/Grond_GTDB_207/Regions_complete/split_by_genus_$region/
  # db_genus_region=${PROJECT_DIR}/database_genus_$region
  db_genome_region=${PROJECT_DIR}/database_genome_$region
  mkdir $db_genome_region
  # mkdir $db_genus_region
  cat $genus_region/*.fasta > $db_genome_region/db_genome_${region}_all.fa
done





# less ${PROJECT_DIR}/03_coverm_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_genus_species_summary.tsv 
# 
# 
# 
# #--------- 2. minimap2
# # 比对
# minimap2 -ax map-ont \
#       $db_genome_rrn/db_genome_rrn_all.fa \
#       ${PROJECT_DIR}/02_read_binning_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_consensus_hdf_1/consensus_barcode${bn}_stitch.fa \
#       -t 16 --cs |\
#       samtools view -F 2308 - |\
#       cut -f1-9,12,21 \
#       > ${PROJECT_DIR}/03_minimap2_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_rrn.sam
# # 统计
# # 定义变量（确保这些变量与你上一步比对时的一致）
# SAM_FILE="${PROJECT_DIR}/03_minimap2_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_rrn.sam"
# DB_FILE="$db_genome_rrn/db_genome_rrn_all.fa"
# OUT_TAX="${PROJECT_DIR}/03_minimap2_Q${Q}_${minLen}_${maxLen}_$batch/barcode${bn}_rrn_tax.txt"
# gawk \
#   -v SAM="$SAM_FILE" \
#   -v DB="$DB_FILE" \
#   '
#     # 第一步：读取精简后的 SAM 文件
#     (FILENAME == SAM){
#       # 在你上一步 cut -f1-9,12,21 之后，参考 ID 依然在第 3 列
#       TARGETS[$3]++ 
#     }
#     # 第二步：读取数据库文件获取 Taxonomy 信息
#     (FILENAME == DB && $0 ~ /^>/){
#       # 提取 Header ID（去掉 >）
#       TARGET=substr($1,2)
#       if (TARGET in TARGETS){
#         # 打印 ID 和 对应的 Taxonomy 描述（$2 是 ID 后的第一个单词）
#         # 如果你想打印整行描述，可以将 $2 改为 substr($0, index($0,$2))
#         print TARGET, $2, "Count:" TARGETS[TARGET]
#       }
#     }
#   ' \
#   "$SAM_FILE" \
#   "$DB_FILE" > "$OUT_TAX"
#   
  
  
# 物种信息对照表：${GROND_DIR}/01_Github_file/GROND/Grond_GTDB_207/sub_Meta_region_genome/META_nospace.tsv
