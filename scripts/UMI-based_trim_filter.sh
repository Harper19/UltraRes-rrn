#!/bin/bash
# ---------------------------------------------------------------------------
# Paths: edit these before running (this file is a step-by-step command log,
# meant to be run block by block rather than executed as one script).
export PROJECT_DIR=/path/to/project      # working directory for this analysis
export SOFTWARE_DIR=/path/to/software    # dorado, usearch, etc.
# ---------------------------------------------------------------------------
#############################################################
################ 1. trim and filter raw data - all################ 
minLen=3500 
maxLen=6000  
Q=10
bar_list="" # Barcode used in Nanopore sequencing

#------------------------- 1) Prepare UMI ref (actual UMI sequence in the dataset) ------------------------- 
# extract UMI pairs from sequencing datasets --> generate a database for UMI for downstream analysis (concensus reads generation)
# prepare pair end reads
conda activate seqtk
# generate primers and their rc
{
TERMINAL_LENGTH=200
FW1=CAAGCAGAAGACGGCATACGAGAT # universal_fw 24bp
FW2=AGRRTTYGATYHTDGYTYAG # 16S
RV1=AATGATACGGCGACCACCGAGATC # universal_rv
RV2=CCRAMCTGTCTCACGACG # 23S
# generate reverse complement 
revcom() {
  echo $1 |\
  gawk '{print ">dummy\n" $0}' |\
  seqtk seq -r - |\
  gawk '!/^>/'  
}
FW1R=$(revcom "$FW1")
FW2R=$(revcom "$FW2")
RV1R=$(revcom "$RV1") 
RV2R=$(revcom "$RV2")
  
}
echo $FW1R $FW2R $RV1R $RV2R
echo $RV1 $RV1R
# extract the UMI from each reads and primer position check for downstream length set up
UMI_DIR=${PROJECT_DIR}/01_UMI_ref_paired_end_Q${Q}_${minLen}_${maxLen}_$batch
OUT_DIR=$UMI_DIR
mkdir $OUT_DIR
cd  $OUT_DIR
{
for bn in $bar_list
do
  br=barcode$bn
  echo $br
  echo $OUT_DIR
  INPUT_FASTQ="${PROJECT_DIR}/00_SUP_demux_Q${Q}_${minLen}_${maxLen}_$batch/${br}_QC.fastq"
  ll -h $INPUT_FASTQ
# extract 200 bp from each end and generate fasta - might not need
  echo "Extract 200 bp from each end to fasta files"
  gawk -v BD="$OUT_DIR" -v TL="$TERMINAL_LENGTH" -v br="$br" '
    NR%4==1{
       print ">" substr($1,2) > BD"/"br"_reads_t1.fa";
       print ">" substr($1,2) > BD"/"br"_reads_t2.fa";
     }
     NR%4==2{
       print substr($0, 1, TL) > BD"/"br"_reads_t1.fa";
       print substr($0, length($0) - TL + 1, TL)  > BD"/"br"_reads_t2.fa";
     }
' "$INPUT_FASTQ"
  seqtk seq -r $OUT_DIR/${br}_reads_t2.fa > $OUT_DIR/${br}_reads_t2_rc.fa


# primer position check (use 16S and 23S primer, include either of them will be included)
  echo "primer position check"
  cutadapt -j 24 \
      -O 10 \
      -e 0.2 \
      -a $FW2 \
      -a $RV2 \
      $OUT_DIR/${br}_reads_t1.fa \
      --discard-untrimmed \
      2>> $OUT_DIR/${br}_umi_endpos_log.txt | \
    awk '
      NR%2==0{
        # Position
        p=length($0)+0
        c[p]++
        if(pm+0 < p){pm = p}
        # SD and MEAN
        sum+=p
        sums+=p*p
        n++
      } END {
        for(j=0; j <= pm; j++){
          print j, c[j]+0
        }
        print "\nMean\tsd"
        print sum/n"\t"sqrt(sums/n -(sum/n)^2)
      }
    '
  # Position Mean: 42
  
  # Extract UMI terminal region to fastq
  echo "Extract 200 bp from each end to fq files"
  UMI_DIR=$OUT_DIR
  gawk -v UD="$UMI_DIR" -v br="$br" 'NR%4==1{
       print $0 > UD"/"br"_reads_tf_start.fq";
       print $0 > UD"/"br"_reads_tf_end.fq";
     }
     NR%4==2{
       print substr($0, 1, 200) > UD"/"br"_reads_tf_start.fq";
       print substr($0, length($0) - 199, 200)  > UD"/"br"_reads_tf_end.fq";
     }
     NR%4==3{
       print $0 > UD"/"br"_reads_tf_start.fq";
       print $0 > UD"/"br"_reads_tf_end.fq";
     }
     NR%4==0{
       print substr($0, 1, 200) > UD"/"br"_reads_tf_start.fq";
       print substr($0, length($0) - 199, 200)  > UD"/"br"_reads_tf_end.fq";
     }
  ' $INPUT_FASTQ

# extract UMI from both ends
  echo "Extract UMI from each end"
  cutadapt -j 16 -e 0.2 -O 11 -m 18 -M 18 \
    --discard-untrimmed \
    -g $FW1...$FW2 -g $RV1...$RV2 \
    -G $RV2R...$RV1R -G $FW2R...$FW1R \
    -o $UMI_DIR/${br}_umi1.fq -p $UMI_DIR/${br}_umi2.fq \
    $UMI_DIR/${br}_reads_tf_start.fq $UMI_DIR/${br}_reads_tf_end.fq \
    > $UMI_DIR/${br}_perfect_trim.log
  
  echo "DONE"
done

}

# check data summary (how many reads have both end UMI and summary of UMI)
for bn in $bar_list
do
  br=barcode$bn
  echo $br
  head $UMI_DIR/${br}_perfect_trim.log
done
less $UMI_DIR/${br}_perfect_trim.log
# check the perfect trim log ()
cat $UMI_DIR/*_perfect_trim.log >  $UMI_DIR/sum_perfect_trim.log
less $UMI_DIR/sum_perfect_trim.log

# summarize the UMI number
{
echo -e " File, total, min, max, median" > $OUT_DIR/PairEnd_summary.csv
echo -e " File, total, min, max, median" > $OUT_DIR/UMI_summary.csv

for bn in $bar_list
do
  br=barcode$bn
  echo "Processing" $br
  # summarize Pair end read number
  for end in {"start","end"}
  do
  echo "Processing pair-end file" $br $end
  Pairend_file=$OUT_DIR/${br}_reads_tf_${end}.fq
  [ -e "$Pairend_file" ] || continue  # skip if no match
    awk -v File="$(basename "$Pairend_file")" 'NR%4==2{
                 len=length($0);
                 count++;
                 lengths[len]++;
                 if(len<min || min==0) min=len;
                 if(len>max) max=len
               }
         END {
           total = count
           mid = (total+1)/2
           csum=0
           for(l in lengths){all_lens[l]=lengths[l]}
           n=asorti(all_lens, idx)
           for(i=1; i<=n; i++){
             l=idx[i]
             csum+=lengths[l]
             if(csum>=mid){median=l; break}
           }
           printf("%s,%d,%d,%d,%d\n", File, total, min, max, median)
         }' "$Pairend_file" >> $OUT_DIR/PairEnd_summary.csv
  done

  for end in {"1","2"}
  do 
    echo "Processing UMI file" $br $end
    UMI_fq=$OUT_DIR/${br}_umi${end}.fq
    [ -e "$UMI_fq" ] || continue  # skip if no match
    awk -v File="$(basename "$UMI_fq")" 'NR%4==2{
                 len=length($0); 
                 count++; 
                 lengths[len]++; 
                 if(len<min || min==0) min=len; 
                 if(len>max) max=len
               } 
         END {
           total = count
           mid = (total+1)/2
           csum=0
           for(l in lengths){all_lens[l]=lengths[l]}
           n=asorti(all_lens, idx)
           for(i=1; i<=n; i++){
             l=idx[i]
             csum+=lengths[l]
             if(csum>=mid){median=l; break}
           }
           printf("%s,%d,%d,%d,%d\n", File, total, min, max, median)
         }' "$UMI_fq" >> $OUT_DIR/UMI_summary.csv

  done 
done
  
}
less $OUT_DIR/PairEnd_summary.csv 
less $OUT_DIR/UMI_summary.csv


# merge two ends together 
{
for bn in $bar_list
do
  br=barcode$bn
  echo $br
  paste -d "" <( sed -n '1~4s/^@/>/p;2~4p' $OUT_DIR/${br}_umi1.fq ) \
            <( sed -n '1~4s/^@/>/p;2~4p' $OUT_DIR/${br}_umi2.fq ) |\
  cut -d " " -f1 > $OUT_DIR/${br}_umi12.fa
  LN=$(cat "$OUT_DIR/${br}_umi12.fa" | wc -l)
  RN=$(($LN/2))
  echo $br $RN
done
}

{
PATTERN="[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{6}\
[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}"
for bn in $bar_list
do
  br=barcode$bn
  echo $br
  grep -B1 -E "$PATTERN" $OUT_DIR/${br}_umi12.fa |\
    sed '/^--$/d' > $OUT_DIR/${br}_umi12f.fa
  LN=$(cat $OUT_DIR/${br}_umi12f.fa | wc -l)
  echo $(($LN/2))
done

}


########## Cluster UMI pairs. ########## 
cd $UMI_DIR
$USEARCH --version #usearch v11.0.667_i86linux32
{
for bn in $bar_list
do
  br=barcode$bn
  echo $br
  echo "--------------------------------- Usearch unique:" $br
  # Collapse identical sequences into unique sequences and count their occurrences
  $USEARCH  -threads 16 \
    -fastx_uniques \
    $OUT_DIR/${br}_umi12f.fa\
    -fastaout $UMI_DIR/${br}_umi12u.fa \
    -sizeout \
    -minuniquesize 1 \
    -relabel umi \
    -strand both
  echo "--------------------------------- Usearch cluster:" $br
  $USEARCH  -threads 16 \
    -cluster_fast $OUT_DIR/${br}_umi12u.fa \
    -id 0.9 \
    -centroids $OUT_DIR/${br}_umi12c.fa \
    -uc $UMI_DIR/${br}_umi12c.txt \
    -sizein \
    -sizeout \
    -strand both \
    -minsize 1
done
}


# summarize the UMI_unique and UMI_ref
{
for bn in $bar_list
do
  br=barcode$bn
  LN_u=$(cat "$OUT_DIR/${br}_umi12u.fa" | wc -l)
  RN_u=$(($LN_u/2))
  LN_c=$(cat "$OUT_DIR/${br}_umi12c.fa" | wc -l)
  RN_c=$(($LN_c/2))
  echo $br $RN_u $RN_c
done 
  
}

##########  Extract putative UMI pairs ########## 
CUTADAPT=cutadapt
echo $UMI_DIR
{
for bn in $bar_list
do
  br=barcode$bn
  echo "Extract putative UMI pairs:" $br
$CUTADAPT -j 16 -e 0.2 -O 11 -m 18 -l 18 \
  --discard-untrimmed \
  -g $FW1 -g $RV1 \
  -G $RV2R -G $FW2R \
  -o $UMI_DIR/${br}_umi1p.fq -p $UMI_DIR/${br}_umi2p.fq \
  $UMI_DIR/${br}_reads_tf_start.fq $UMI_DIR/${br}_reads_tf_end.fq \
  > $UMI_DIR/${br}_putative_trim.log

  paste -d "" <( sed -n '1~4s/^@/>/p;2~4p' $UMI_DIR/${br}_umi1p.fq ) \
            <( sed -n '1~4s/^@/>/p;2~4p' $UMI_DIR/${br}_umi2p.fq ) |\
  cut -d " " -f1 > $UMI_DIR/${br}_umi12p.fa
done 
}

# generate the UMI ref database
BWA=bwa
SAMTOOLS=samtools
GAWK=gawk 
{
for bn in $bar_list
do
  br=barcode$bn
  echo $br
  echo "Indexing..."
  $BWA index $UMI_DIR/${br}_umi12c.fa 
  echo "Aligning..."
  $BWA aln \
    $UMI_DIR/${br}_umi12c.fa \
    $UMI_DIR/${br}_umi12p.fa \
    -n 6 \
    -N > $UMI_DIR/${br}_umi12p_map.sai
  echo "Mapping..."
  $BWA samse \
    -n 10000000 \
    $UMI_DIR/${br}_umi12c.fa \
    $UMI_DIR/${br}_umi12p_map.sai \
    $UMI_DIR/${br}_umi12p.fa|\
    $SAMTOOLS view -F 4 - \
    > $UMI_DIR/${br}_umi12p_map.sam
    
  echo $br
  echo "retain only cluster size >2"
  $GAWK \
    -v UMS="$UMI_DIR/${br}_umi12p_map.sam" \
    -v UC="$UMI_DIR/${br}_umi12c.fa" \
    '
    (FILENAME == UMS){
        CLUSTER[$3]++
    }
    (FILENAME == UC && FNR%2==1){
      NAME=substr($1,2)
      if (NAME in CLUSTER){
        if (CLUSTER[NAME]+0 > 1){
          SIZE=CLUSTER[NAME]
          gsub(";.*", "", NAME)
          print ">" NAME ";size=" SIZE ";"
          getline; print
        }
      }
    }
    ' \
    $UMI_DIR/${br}_umi12p_map.sam \
    $UMI_DIR/${br}_umi12c.fa \
    > $UMI_DIR/${br}_umi12cf.fa

# remove potential chimeras
  echo "Remove potential chimeras"
  paste <(cat $UMI_DIR/${br}_umi12cf.fa | paste - - ) \
    <($GAWK '!/^>/{print}' $UMI_DIR/${br}_umi12cf.fa | rev | tr ATCG TAGC) |\
    $GAWK -v UD="$UMI_DIR" -v br="$br" 'NR==FNR {
        #Format columns
        split($1, a, /[>;]/);
        sub("size=", "", a[3]);
        # Extract UMI1 and UMI2 in both orientations
        s1 = substr($2, 1, 18);
        s2 = substr($2, 19, 36);
        s1rc= substr($3, 1, 18);
        s2rc= substr($3, 19, 36);
        # Register UMI1 size if larger than current highest or if empty
        if ((g1n[s1]+0) <= (a[3]+0) || g1n[s1] == ""){
          g1n[s1] = a[3];
          g1[s1] = a[2];
        }
        # Register UMI2 size if larger than current highest or if empty
        if ((g2n[s2]+0) <= (a[3]+0) || g2n[s2] == ""){
          g2n[s2] = a[3];
          g2[s2] = a[2];
        }
        # Register UMI1rc size if larger than current highest or if empty
        if ((g1n[s1rc]+0) <= (a[3]+0) || g1n[s1rc] == ""){
          g1n[s1rc] = a[3];
          g1[s1rc] = a[2];
        }
        # Register UMI2rc size if larger than current highest or if empty
        if ((g2n[s2rc]+0) <= (a[3]+0) || g2n[s2rc] == ""){
          g2n[s2rc] = a[3];
          g2[s2rc] = a[2];
        }
        # Register UMI1 and UMI matches for current UMI
        u[a[2]] = a[3];
        s1a[a[2]] = s1;
        s2a[a[2]] = s2;
        s1arc[a[2]] = s1rc;
        s2arc[a[2]] = s2rc;
      } END {
        for (i in u){
          keep="no";
          if (g1[s1a[i]] == i && g2[s2a[i]] == i && g1[s1arc[i]] == i && g2[s2arc[i]] == i && s1a[i] != s1arc[i]){
            keep="yes";
            print ">"i";"u[i]"\n"s1a[i]s2a[i] > UD"/"br"_umi_ref.fa";
          } else if (s1a[i] == s1arc[i]){
            keep="tandem"
            print ">"i";"u[i]"\n"s1a[i]s2a[i] > UD"/"br"_umi_ref.fa";
          }
          print i, n[i], s1a[i], s2a[i], keep, g1[s1a[i]]"/"g2[s2a[i]]"/"g1[s1arc[i]]"/"g2[s2arc[i]], u[i]
        }  
      }' > $UMI_DIR/${br}_umi_ref.txt
done
}

#  considering change to bwa-mem2 or minimap2 (not minimap2 cuz this part is short reads, i.e., UMI)
br=barcode12
grep "^umi1 " $UMI_DIR/${br}_umi_ref.txt | wc -l
grep "^umi " $UMI_DIR/${br}_umi_ref.txt | wc -l


# Remove potential chimeras
{
# summarize the putative UMI
# br_name="barcode22 barcode23 barcode24"
for bn in $bar_list
do
  br=barcode$bn
  echo $br
  LN_12p=$(less $UMI_DIR/${br}_umi12p.fa | wc -l)
  RN_12p=$(($LN_12p/2))

  LN_12cf=$(less $UMI_DIR/${br}_umi12cf.fa | wc -l)
  RN_12cf=$(($LN_12cf/2))
  # echo "cf_UMI" $RN_12cf

  LN_ref=$(less $UMI_DIR/${br}_umi_ref.fa | wc -l)
  RN_ref=$(($LN_ref/2))
  echo $br "purative_UMI cf_UMI ref_UMI:" $RN_12p $RN_12cf $RN_ref

  awk '{
    if ($4=="yes") yes++
    else if ($4=="tandem") tandem++
    else chimera++
  } END {
    print "yes:",yes, "tandem:",tandem, "chimera:",chimera
  }' ${br}_umi_ref.txt
done
}