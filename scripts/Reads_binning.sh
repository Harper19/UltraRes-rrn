#!/bin/bash
# ---------------------------------------------------------------------------
# Paths: edit these before running (this file is a step-by-step command log,
# meant to be run block by block rather than executed as one script).
export PROJECT_DIR=/path/to/project      # working directory for this analysis
export SOFTWARE_DIR=/path/to/software    # dorado, usearch, etc.
# ---------------------------------------------------------------------------
# reads binning using UMIs
minLen=3500 
maxLen=6000  
Q=10
batch=1
CUTADAPT=cutadapt
BWA=bwa
SAMTOOLS=samtools
SEQTK=seqtk
USEARCH=${SOFTWARE_DIR}/usearch
OUT_DIR=${PROJECT_DIR}/01_UMI_ref_paired_end_Q${Q}_${minLen}_${maxLen}_$batch
UMI_DIR=$OUT_DIR
bar_list=""

  
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
MAX_BIN_SIZE=10000000 # -N  Max number of reads with +/- orientation. [Default = 10000]
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


umi_binning() {
  # Input
  local UMIMAP=$1
  local OUT=$2

  # Binning
  gawk -v out="$OUT" '
    BEGIN {g=1; outsub=out"/"g; system("mkdir \047" outsub "\047");}
    NR==FNR {
      # Get read name
      sub(";.*", "", $1);
      # Associate read name and umi match
      bin[$2]=$1;
      # Assign umi to a folder group if it has none
      if (foldergrp[$1] == ""){
        j++;
        if (j <= 4000){
          foldergrp[$1]=g;
        } else {
          j = 0;
          g++;
          foldergrp[$1]=g;
          outsub=out"/"g;
          system("mkdir \047" outsub "\047");
        }
      }
      next;
    }
    FNR%4==1 {
      read=substr($1,2);
      bin_tmp=bin[read]
      if ( bin_tmp != "" ){
        binfile=out"/"foldergrp[bin_tmp]"/"bin_tmp"bins.fastq";
        print > binfile;
        getline; print > binfile;
        getline; print > binfile;
        getline; print > binfile;
      }
    }
  ' $UMIMAP -
}
export -f umi_binning


aggregate_bins() {
  # Input
  local IN=$1
  local OUTDIR=$2
  local OUTNAME=$3
  local JOB=$4

  # Determine output folder
  local BIN=$(( ($JOB - 1)/4000 ))
  
  mkdir -p $OUTDIR/$BIN

  # Aggregate data
  cat $IN > $OUTDIR/$BIN/$OUTNAME  
}
export -f aggregate_bins


{
for bn in $bar_list
do
  br=barcode$bn
  echo $br
  echo "Extract UMI region from files..."
  # Extract UMI region
  $GAWK -v BD="$BINNING_DIR" -v TL="$START_READ_CHECK" -v br="$br" '
    NR%4==1{
      print ">" substr($1,2) > BD"/"br"_reads_tf_umi1.fa";
    }
    NR%4==2{
      print substr($0, 1, TL) > BD"/"br"_reads_tf_umi1.fa";
    }
  ' $UMI_DIR/${br}_reads_tf_start.fq
  # 1387968
  $GAWK -v BD="$BINNING_DIR" -v TL="$END_READ_CHECK" -v br="$br" '
    NR%4==1{
       print ">" substr($1,2) > BD"/"br"_reads_tf_umi2.fa";
     }
     NR%4==2{
       print substr($0, length($0) - TL + 1, TL)  > BD"/"br"_reads_tf_umi2.fa";
     }
  ' $UMI_DIR/${br}_reads_tf_end.fq
  # 1387968

  echo "Devide UMIs from UMI_ref into 2 parts according to the length (18 bp)..."
  # Divide in barcode1 and barcode2 files
  cat $UMI_DIR/${br}_umi_ref.fa <($SEQTK seq -r $UMI_DIR/${br}_umi_ref.fa |\
    $GAWK 'NR%2==1{print $0 "_rc"; getline; print};') |\
    $GAWK -v BD="$BINNING_DIR" -v br="$br" 'NR%2==1{
         print $0 > BD"/"br"_umi_ref_b1.fa";
         print $0 > BD"/"br"_umi_ref_b2.fa";
       }
       NR%2==0{
         print substr($0, 1, 18) > BD"/"br"_umi_ref_b1.fa";
         print substr($0, 19, 18)  > BD"/"br"_umi_ref_b2.fa";
       }'
  # Map UMIs to UMI references
  echo "UMIs to UMI references ......................................"
  echo "Indexing (end 1+2)..."
  $BWA index $BINNING_DIR/${br}_reads_tf_umi1.fa
  $BWA index $BINNING_DIR/${br}_reads_tf_umi2.fa

  echo "Aln (end 1)..."
  $BWA aln $BINNING_DIR/${br}_reads_tf_umi1.fa $BINNING_DIR/${br}_umi_ref_b1.fa \
    -n 3 -N -t 16 > $BINNING_DIR/${br}_umi1_map.sai
  echo "Mapping (end 1)..."

  $BWA samse -n 10000000 $BINNING_DIR/${br}_reads_tf_umi1.fa $BINNING_DIR/${br}_umi1_map.sai\
    $BINNING_DIR/${br}_umi_ref_b1.fa | $SAMTOOLS view -F 20 - > $BINNING_DIR/${br}_umi1_map.sam
  echo "Aln (end 2)..."
  $BWA aln -t 16 $BINNING_DIR/${br}_reads_tf_umi2.fa $BINNING_DIR/${br}_umi_ref_b2.fa \
    -n 3 -N > $BINNING_DIR/${br}_umi2_map.sai
  echo "Mapping (end 2)..."
  $BWA samse -n 10000000 $BINNING_DIR/${br}_reads_tf_umi2.fa $BINNING_DIR/${br}_umi2_map.sai\
    $BINNING_DIR/${br}_umi_ref_b2.fa | $SAMTOOLS view -F 20 - > $BINNING_DIR/${br}_umi2_map.sam

  echo "UMI binning and filtering ......................................"
  $GAWK \
  -v BD="$BINNING_DIR" \
  -v UME_MATCH_ERROR="$UMI_MATCH_ERROR" \
  -v UME_MATCH_ERROR_SINGLE="$UMI_MATCH_ERROR_single" \
  -v UME_MATCH_ERROR_SD="$UMI_MATCH_ERROR_SD" \
  -v RO_FRAC="$RO_FRAC" \
  -v MAX_BIN_SIZE="$MAX_BIN_SIZE"  \
  -v BIN_CLUSTER_RATIO="$BIN_CLUSTER_RATIO" \
  -v br="$br" \
  '
## Process the first file (UMI1 matches)
  NR==1 {
    print "[" strftime("%T") "] ### Read-UMI match filtering ###" > "/dev/stderr";
    print "[" strftime("%T") "] Reading UMI1 match file..." > "/dev/stderr";
  }
  
  # Read UMI match file
  NR==FNR{
    # Extract data from optional fields
    for (i=12; i <= NF; i++){
      # Find NM field and remove prefix (primary hit err)
      if($i ~ /^NM:i:/){sub("NM:i:", "", $i); perr = $i};
      # Find secondary hit field, remove prefix and split hits
      if($i ~ /^XA:Z:/){sub("XA:Z:", "", $i); split($i, shits, ";")};
    }
    # Add primary mapping to hit list
    err1[$1][$3]=perr;
    # Add secondary mapping to hit list
    #Iterate over each hit
    for (i in shits){
      # Split hit in subdata (read, pos, cigar, err)  
      split(shits[i], tmp, ",");
      # Add hit if non-empty, not seen before and not target reverse strand
      if (tmp[1] != "" && !(tmp[1] in err1[$1]) && tmp[2] ~ "+"){
        err1[$1][tmp[1]] = tmp[4];
      }
    }
    next;
  }

## Process the second file (UMI2 matches) 
  FNR==1 {
   print "[" strftime("%T") "] Reading UMI2 match file..." > "/dev/stderr";
  }
  
  # Read UMI match file
  {
    # Extract data from optional fields
    for (i=12; i <= NF; i++){
      # Find NM field and remove prefix (primary hit err)
      if($i ~ /^NM:i:/){sub("NM:i:", "", $i); perr = $i};
      # Find secondary hit field and remove prefix
      if($i ~ /^XA:Z:/){sub("XA:Z:", "", $i); split($i, shits, ";")};
    }
    # Add primary mapping to hit list
    err2[$1][$3]=perr;
    # Add secondary mapping to hit list
    # Split list of hits 
    #Iterate over each hit
    for (i in shits){
      # Split hit in subdata (read, pos, cigar, err)
      split(shits[i], tmp, ",");
      # Add hit if non-empty, not seen before and not target reverse strand
      if (tmp[1] != "" && !(tmp[1] in err2[$1]) && tmp[2] ~ "+"){
        err2[$1][tmp[1]] = tmp[4];
      }
    }
  } END {
    print "[" strftime("%T") "] UMI match filtering..." > "/dev/stderr"; 
    
    
## Start filtering
  # 1. Filter by UMI match error
    # Filter reads based on UMI match error
    for (umi in err1){    
      for (read in err1[umi]){
        # Define vars
        e1 = err1[umi][read];
        e2 = err2[umi][read];
        # Filter reads not matching both UMIs
        if (e1 != "" && e2 != ""){
          # Filter based on mapping error 
          if (e1 + e2 <= UME_MATCH_ERROR && e1 <= UME_MATCH_ERROR_SINGLE && e2 <= UME_MATCH_ERROR_SINGLE){
            # Add read to bin list or replace bin assignment if error is lower
            if (!(read in match_err)){
              match_umi[read] = umi;
              match_err[read] = e1 + e2;
            } else if (match_err[read] > e1 + e2 ){
              match_umi[read] = umi;
              match_err[read] = e1 + e2;
            } 
          }
        }
      }
    }
    print "[" strftime("%T") "] Read orientation filtering..." > "/dev/stderr";
    
  # 2. Filter by read orientation
    # Count +/- strand reads
    for (s in match_umi){
      UM=match_umi[s]
      sub("_rc", "", UM)
      # Read orientation stats
      ROC=match(match_umi[s], /_rc/)
      if (ROC != 0){
        umi_ro_plus[UM]++
        roc[s]="+"
      } else {
        umi_ro_neg[UM]++
        roc[s]="-"
      }
      # Count reads per UMI bin
      umi_n_raw[UM]++;
    }
    
    # Calculate read orientation fraction
    for (u in umi_ro_plus){
      # Check read orientation fraction
      if (umi_ro_plus[u] >1 && umi_ro_neg[u] > 1){
#########  changed from 1 to 0
      # if (umi_ro_plus[u] >0 && umi_ro_neg[u] > 0){
        if (umi_ro_plus[u]/(umi_ro_neg[u]+umi_ro_plus[u]) < RO_FRAC ){
          rof_check[u]="rof_subset"
          rof_sub_neg_n[u] = umi_ro_plus[u]*(1/RO_FRAC-1)
          rof_sub_pos_n[u] = rof_sub_neg_n[u]
        } else if (umi_ro_neg[u]/(umi_ro_neg[u]+umi_ro_plus[u]) < RO_FRAC ){
          rof_check[u]="rof_subset"
          rof_sub_neg_n[u]=umi_ro_neg[u]*(1/RO_FRAC-1)
          rof_sub_pos_n[u]=rof_sub_neg_n[u]
        } else {
          rof_check[u]="rof_ok"
          rof_sub_neg_n[u]=MAX_BIN_SIZE
          rof_sub_pos_n[u]=MAX_BIN_SIZE
        }
      } else {
        rof_check[u]="rof_fail"
      }
    }
    
    # Subset reads
    for (s in match_umi){
      UMI_NAME=match_umi[s]
      sub("_rc", "", UMI_NAME)
      if(roc[s] == "+"){
        if(rof_sub_pos_n[UMI_NAME]-- > 0){
          ror_filt[s]=UMI_NAME
        }
      } else if (roc[s] == "-"){
        if(rof_sub_neg_n[UMI_NAME]-- > 0){
          ror_filt[s]=UMI_NAME
        }
      }
    }

    print "[" strftime("%T") "] UMI match error filtering..." > "/dev/stderr";
    
  # 3. Filter by the SD of the UMI match error
    # Calculate UME stats
    for (s in ror_filt){
      UM=ror_filt[s]
      # Count matching reads
      umi_n[UM]++;
      # UMI match error stats
      umi_me_sum[UM] += match_err[s]
      umi_me_sq[UM] += (match_err[s])^2
    }

    # Check UMI match error
    for (u in umi_n){
      UME_MEAN[u] = umi_me_sum[u]/umi_n[u]
      UME_SD[u] = sqrt((umi_me_sq[u]-umi_me_sum[u]^2/umi_n[u])/umi_n[u])
      if (UME_MEAN[u] > UME_MATCH_ERROR || UME_SD[u] > UME_MATCH_ERROR_SD){
        ume_check[u] = "ume_fail"
      } else {
        ume_check[u] = "ume_ok"
      }
    }

  # 4. Filter by UMI bin / cluster size ratio
    print "[" strftime("%T") "] UMI bin/cluster size ratio filtering..." > "/dev/stderr";
    for (u in umi_n){
      CLUSTER_SIZE=u
      sub(".*;", "", CLUSTER_SIZE)
      bcr[u]=umi_n_raw[u]/CLUSTER_SIZE
      if (bcr[u] > BIN_CLUSTER_RATIO){
        bcr_check[u] = "bcr_fail"
      } else {
        bcr_check[u] = "bcr_ok"
      }
    }

## Write binning statistics
    # Print filtering stats
    print "umi_name", "read_n_raw", "read_n_filt", "read_n_plus", "read_n_neg", \
      "read_max_plus", "read_max_neg", "read_orientation_ratio", "ror_filter", \
      "umi_match_error_mean", "umi_match_error_sd", "ume_filter", "bin_cluster_ratio", \
      "bcr_filter" > BD"/"br"_umi_binning_stats.txt"
    for (u in umi_n){
      print u, umi_n_raw[u], umi_n[u], umi_ro_plus[u], umi_ro_neg[u], \
        rof_sub_pos_n[u] + umi_ro_plus[u], rof_sub_neg_n[u] + umi_ro_neg[u], rof_check[u], \
        UME_MEAN[u], UME_SD[u], ume_check[u], bcr[u], bcr_check[u]\
        > BD"/"br"_umi_binning_stats.txt"
    }
	
    print "[" strftime("%T") "] Print UMI matches..." > "/dev/stderr"; 
    for (s in ror_filt){
      UMI_NAME=ror_filt[s]
      if( \
          ume_check[UMI_NAME] == "ume_ok" && \
          rof_check[UMI_NAME] == "rof_ok" && \
          bcr_check[UMI_NAME] == "bcr_ok" \
      ){print UMI_NAME, s, match_err[s]}
    }
    # Print to terminal
    print "[" strftime("%T") "] Done." > "/dev/stderr"; 
  }
' $BINNING_DIR/${br}_umi1_map.sam $BINNING_DIR/${br}_umi2_map.sam > $BINNING_DIR/${br}_umi_bin_map.txt

  echo "Extract binned reads..........."
  BINNING_DIR_BR=${BINNING_DIR}/${br}_bin
      cat ${PROJECT_DIR}/00_SUP_demux_Q${Q}_${minLen}_${maxLen}_$batch/${br}_QC.fastq | \
      parallel \
      -L4 \
      --block 300M \
      --pipe \
      "mkdir -p $BINNING_DIR_BR/job{#}; cat | umi_binning ${BINNING_DIR}/${br}_umi_bin_map.txt $BINNING_DIR_BR/job{#}"
      
      
  # prepare binning files
  echo $br
  BINNING_DIR_BR=${BINNING_DIR}/${br}_bin
  echo "create br_bin directory"
  
  find $BINNING_DIR_BR/*/*/ -name "*bins.fastq" -printf "%f\n" |\
    sort |\
    uniq |\
    $GNUPARALLEL \
      -j 10 \
        "aggregate_bins \"$BINNING_DIR_BR/*/*/{1}\" \"$BINNING_DIR_BR\" \"{1}\" {#}" 
  rm -r $BINNING_DIR_BR/job*
  
  # Trimming and Filtering based on sub UMIs
  echo "Filtering based on sub UMIs Output:"$UMI_DIR/${br}_sumi.fa # 18bp single end UMI
  cat \
    $UMI_DIR/${br}_reads_tf_start.fq \
    <($SEQTK seq -r $UMI_DIR/${br}_reads_tf_end.fq) |\
  cutadapt -j 32 \
    -e 0.1 \
    -O 11 \
    -m 18 \
    -M 18 \
    --discard-untrimmed \
    -g $FW1...$FW2 \
    -g $RV1...$RV2 \
    - 2> $UMI_DIR/${br}_sumi_trim.log |\
    gawk '
      NR%4==1{print ">" substr($1, 2)}
      NR%4==2{print $0}
    ' > $UMI_DIR/${br}_sumi.fa


  $USEARCH -threads 32 \
    -fastx_uniques $UMI_DIR/${br}_sumi.fa \
    -fastaout $UMI_DIR/${br}_sumi_u.fa \
    -relabel sumi \
    -strand plus \
    -sizeout \
    -minuniquesize 1

  PATTERN="[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}"
  grep -B1 -E "$PATTERN" $UMI_DIR/${br}_sumi_u.fa |\
    sed '/^--$/d' > $UMI_DIR/${br}_sumi_f.fa

  $SEQTK seq -r $UMI_DIR/${br}_sumi_f.fa > $UMI_DIR/${br}_sumi_frc.fa
  
  $GAWK \
    -v UF="$UMI_DIR/${br}_sumi_f.fa" \
    -v UR="$UMI_DIR/${br}_sumi_frc.fa" \
    -v UFR="$UMI_DIR/${br}_umi12u.fa" \
    '
    (FILENAME == UF && FNR%2==1){
      SIZE=$1
      sub(".*=|;", "", SIZE)
      getline
      UMI_FWD[$0]=SIZE
    }
    (FILENAME == UR && FNR%2==1){
      SIZE=$1
      gsub(".*=|;", "", SIZE)
      getline
      UMI_RV[$0]=SIZE
    }
    (FILENAME == UFR && FNR%2==1){
      NAME=$1
      getline
      # Extract UMI1 and UMI2 in both orientations
      s1 = substr($1, 1, 18);
      s2 = substr($1, 19, 36);
      if (s1 in UMI_FWD && s2 in UMI_RV){
        # print NAME UMI_FWD[s1] ";" UMI_RV[s2] ";" s1 ";" s2 ";\n" $1 
        print NAME UMI_FWD[s1] UMI_RV[s2] ";" s1 ";" s2 ";\n" $1 
  
      }
    }
   ' \
   $UMI_DIR/${br}_sumi_f.fa \
   $UMI_DIR/${br}_sumi_frc.fa \
   $UMI_DIR/${br}_umi12u.fa \
   > $UMI_DIR/${br}_umi12uf.fa    # a UMI pair is retained only if both umi1 and its reverse complement occur among the unique UMI pairs
done
}


# test
for bn in $bar_list
do
  br=barcode$bn
  echo $br
  awk '{print $1}' ${br}_umi1_map.sam | sort | uniq -c | sort -nr | head
  grep -v "^@" ${br}_umi1_map.sam | wc -l
  awk '{for(i=12;i<=NF;i++){if($i ~ /^NM:i:/){split($i,a,":"); print a[3]}}}' ${br}_umi1_map.sam | sort -n | uniq -c
done 



# cehck the mapping results (count error_value)
br=barcode12
for f in ${br}_umi1_map.sam ${br}_umi2_map.sam
do
  echo "=== $f ==="
  awk '{
    for (i=12; i<=NF; i++){
      if ($i ~ /^NM:i:/){
        split($i,a,":");
        print a[3]
      }
    }
  }' $f | sort -n | uniq -c
done
