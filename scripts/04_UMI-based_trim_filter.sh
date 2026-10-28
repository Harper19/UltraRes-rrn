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
# 04_trim based on UMI

salloc -q inferno -A YOUR_ACCOUNT -N1 --ntasks-per-node=12 --time=12:30:00 --mem=32G
module load anaconda3/2022.05
batch=1

#############################################################
################ 1. trim and filter raw data - all################ 
minLen=3500 
maxLen=6000  
Q=10
bar_list="01 02 03 04 05 06 07 08 09 10 11 12 13"

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
# check 
head $UMI_DIR/barcode12_perfect_trim.log

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
#   10Gb
# barcode01
# barcode01 77629
# barcode02
# barcode02 17703
# barcode03
# barcode03 56658
# barcode04
# barcode04 77108
# barcode05
# barcode05 135402
# barcode06
# barcode06 160421
# barcode07
# barcode07 95672
# barcode08
# barcode08 114365
# barcode09
# barcode09 31282
# barcode10
# barcode10 64024
# barcode11
# barcode11 58429
# barcode12
# barcode12 38595
# barcode13
# barcode13 44
}
# 25Gb
# barcode01
# barcode01 153611
# barcode02
# barcode02 35332
# barcode03
# barcode03 112256
# barcode04
# barcode04 152206
# barcode05
# barcode05 268204
# barcode06
# barcode06 315164
# barcode07
# barcode07 189374
# barcode08
# barcode08 225635
# barcode09
# barcode09 61098
# barcode10
# barcode10 126816
# barcode11
# barcode11 115139
# barcode12
# barcode12 76047
# barcode13
# barcode13 92
# Extract UMI pairs with correct patterns 
# Pattern: (NNNYRNNNYRNNNYRNNN NNNYRNNNYRNNNYRNNN)
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
# 10 Gb
# barcode01
# 75426
# barcode02
# 17269
# barcode03
# 55222
# barcode04
# 75004
# barcode05
# 131703
# barcode06
# 155905
# barcode07
# 93191
# barcode08
# 111262
# barcode09
# 30298
# barcode10
# 62403
# barcode11
# 56761
# barcode12
# 37461
# barcode13
# 42

# 20 Gb
# barcode01
# 148700
# barcode02
# 34365
# barcode03
# 108868
# barcode04
# 147544
# barcode05
# 260207
# barcode06
# 305283
# barcode07
# 183757
# barcode08
# 218707
# barcode09
# 59062
# barcode10
# 123174
# barcode11
# 111567
# barcode12
# 73662
# barcode13
# 85




# test 
br=barcode12
grep -v "^>" $OUT_DIR/${br}_umi12.fa | sort | uniq -c | sort -nr > $OUT_DIR/${br}_umi12_counts.txt # check the size of each umi pairs
less $OUT_DIR/${br}_umi12_counts.txt 
# head $OUT_DIR/${br}_umi12_counts.txt 
grep -v "^>" $OUT_DIR/${br}_umi12.fa | sort | uniq -c | awk '$1==1' | wc -l # check how many singlton: 710404 barcode24
awk 'NR % 4 == 2' $OUT_DIR/${br}_umi1.fq | sort | uniq -c | sort -nr > $OUT_DIR/${br}_umi1_counts.txt
less $OUT_DIR/${br}_umi1_counts.txt 
awk '{freq[$1]++} END {for (size in freq) print size, freq[size] | "sort -n"}' $OUT_DIR/${br}_umi1_counts.txt > $OUT_DIR/${br}_umi1_depth_freq.txt
less $OUT_DIR/${br}_umi1_depth_freq.txt
# 10 Gb
# 1 163976
# 2 90960
# 3 58037
# 4 31548
# 5 14478
# 6 6058
# 7 2228
# 8 790
# 9 253
# 10 113
# 11 36
# 12 14
# 13 5
# 14 1
# 15 3

# 20 Gb
# 4 49
# 5 33
# 6 16
# 7 20
# 8 17
# 9 18
# 10 12
# 11 12
# 12 15
# 13 15
# 14 19
# 15 27
# 16 16
# 17 27
# 18 21
# 19 19
# 20 32
# 21 27
# 22 38
# 23 30
# 24 46
# 25 47
# 26 42
# 27 46
# 28 40
# 29 46
# 30 50
# 31 64
# 32 70
# 33 81
# 34 71
# 35 66
# 36 59
# 37 75
# 38 57
# 39 70
# 40 49
# 41 71
# 42 63
# 43 56
# 44 54
# 45 51

########## Cluster UMI pairs. ########## 
# install usearch
# cd ${SOFTWARE_DIR}
# wget https://drive5.com/downloads/usearch11.0.667_i86linux32.gz
# gunzip usearch11.0.667_i86linux32.gz
# chmod a+rx ${SOFTWARE_DIR}/usearch11.0.667_i86linux32
cd $UMI_DIR
USEARCH=${SOFTWARE_DIR}/usearch11.0.667_i86linux32
$USEARCH --version #usearch v11.0.667_i86linux32
{
for bn in $bar_list
do
  br=barcode$bn
  echo $br
  echo "--------------------------------- Usearch unique:" $br
  # 把相同序列合并为一条唯一序列，并统计出现次数
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
# test 
br=barcode11
grep -c "^H" $OUT_DIR/${br}_umi12c.txt # 10Gb 1235 20Gb 3465
grep -c "^S" $OUT_DIR/${br}_umi12c.txt # 10Gb 2405 20 Gb 5538
# 如果 S ≈ 总数 且 H ≈ 0 → 阈值太高。
# 如果两者都很多，但数量接近 → 序列本身差异大。
# 如果 .uc 文件里所有 size 都是 1 → 忘记加 -sizein。
# 如果 .fa 文件行数一模一样 → 输入文件错了。
# check the distribution of UMI unique size.
grep "^>" $OUT_DIR/${br}_umi12u.fa | wc -l # 10 Gb 3628 20Gb 8919
grep "^>" $OUT_DIR/${br}_umi12c.fa | wc -l # 10 Gb 0.9: 2405 #0.8: ?  ## 20 Gb 0.9 -> 5538
head $OUT_DIR/${br}_umi12c.fa



cd $OUT_DIR/

echo ${br}
br=barcode12
grep -o "size=[0-9]*" ${br}_umi12u.fa | cut -d= -f2 | sort -n | uniq -c
grep -o "size=[0-9]*" ${br}_umi12c.fa | cut -d= -f2 | sort -n | uniq -c
head ${br}_umi12u.fa
head ${br}_umi12c.fa
# # id=0.9 10 Gb 
# >umi3499;size=90;
# TTGCATTCTGATCCATGTATGCGCGTTACGACGTCA
# >umi2652;size=89;
# AGCTAATATGTTATACTCCTATAGTTCATTACGGAT
# >umi3;size=86;
# AGGTGTCACGAACTACTGCACTGCTGCGAGGTACCC
# >umi1368;size=84;
# ACCCGCCATACGCCGAAGGAATACACCGCCACACCT
# >umi7;size=83;
# TGCTGTTGTAGTTTGTTACACTGAATCATAGCACGC
# id=0.8
# ?

## id 0.9 20 Gb 
# >umi1;size=176;
# TGCTGTTGTAGTTTGTTACACTGAATCATAGCACGC
# >umi2;size=176;
# AGGTGTCACGAACTACTGCACTGCTGCGAGGTACCC
# >umi3;size=171;
# TTGCATTCTGATCCATGTATGCGCGTTACGACGTCC
# >umi4;size=163;
# AGCTAATATGTTATACTTCTATAGTTCATTACGGAT
# >umi5;size=162;
# AGGTGTGGCGGTGTATTCGGTCGGCGTATGGCGGGT


# generate umi12c 代表多少条原始reads
# 因为用了size in 和size out 所以原始reads-》 unique -》 clsuter已经传递了数量


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
## 10 Gb 
# barcode01 6084 3881
# barcode02 1332 819
# barcode03 4647 2913
# barcode04 5865 3599
# barcode05 11647 7651
# barcode06 12130 7413
# barcode07 6641 3915
# barcode08 8161 4941
# barcode09 2962 1934
# barcode10 5550 3592
# barcode11 5404 3718
# barcode12 3628 2405
# barcode13 42 42

## 20 Gb
# barcode01 10136 5797
# barcode02 2215 1194
# barcode03 7719 4335
# barcode04 10068 5539
# barcode05 19690 11912
# barcode06 20853 11885
# barcode07 11502 6155
# barcode08 14021 7851
# barcode09 4586 2589
# barcode10 9254 5266
# barcode11 8919 5538
# barcode12 5852 3451
# barcode13 83 83

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
## 10 Gb 
# barcode01
# barcode01 purative_UMI cf_UMI ref_UMI: 101990 1622 1284
# yes: 1284 tandem:  chimera: 338
# barcode02
# barcode02 purative_UMI cf_UMI ref_UMI: 25522 365 303
# yes: 303 tandem:  chimera: 62
# barcode03
# barcode03 purative_UMI cf_UMI ref_UMI: 74020 1493 1274
# yes: 1274 tandem:  chimera: 219
# barcode04
# barcode04 purative_UMI cf_UMI ref_UMI: 103158 1636 1314
# yes: 1314 tandem:  chimera: 322
# barcode05
# barcode05 purative_UMI cf_UMI ref_UMI: 174012 3340 2763
# yes: 2763 tandem:  chimera: 577
# barcode06
# barcode06 purative_UMI cf_UMI ref_UMI: 208022 2958 2248
# yes: 2248 tandem:  chimera: 710
# barcode07
# barcode07 purative_UMI cf_UMI ref_UMI: 124582 1854 1526
# yes: 1526 tandem:  chimera: 328
# barcode08
# barcode08 purative_UMI cf_UMI ref_UMI: 150623 2154 1680
# yes: 1680 tandem:  chimera: 474
# barcode09
# barcode09 purative_UMI cf_UMI ref_UMI: 44714 816 652
# yes: 652 tandem:  chimera: 164
# barcode10
# barcode10 purative_UMI cf_UMI ref_UMI: 83857 1896 1631
# yes: 1631 tandem:  chimera: 265
# barcode11
# barcode11 purative_UMI cf_UMI ref_UMI: 78137 1483 1139
# yes: 1139 tandem:  chimera: 344
# barcode12
# barcode12 purative_UMI cf_UMI ref_UMI: 51798 1275 1085
# yes: 1085 tandem:  chimera: 190
# barcode13
# barcode13 purative_UMI cf_UMI ref_UMI: 117 0 0
# yes:  tandem:  chimera:
  
## 20 Gb
# barcode01
# barcode01 purative_UMI cf_UMI ref_UMI: 203083 2088 1196
# yes: 1196 tandem:  chimera: 892
# barcode02
# barcode02 purative_UMI cf_UMI ref_UMI: 48725 429 288
# yes: 288 tandem:  chimera: 141
# barcode03
# barcode03 purative_UMI cf_UMI ref_UMI: 148470 1803 1193
# yes: 1193 tandem:  chimera: 610
# barcode04
# barcode04 purative_UMI cf_UMI ref_UMI: 203196 2092 1241
# yes: 1241 tandem:  chimera: 851
# barcode05
# barcode05 purative_UMI cf_UMI ref_UMI: 349427 4187 2568
# yes: 2568 tandem:  chimera: 1619
# barcode06
# barcode06 purative_UMI cf_UMI ref_UMI: 414373 3976 2053
# yes: 2053 tandem:  chimera: 1923
# barcode07
# barcode07 purative_UMI cf_UMI ref_UMI: 249630 2399 1456
# yes: 1456 tandem:  chimera: 943
# barcode08
# barcode08 purative_UMI cf_UMI ref_UMI: 297682 2886 1536
# yes: 1536 tandem:  chimera: 1350
# barcode09
# barcode09 purative_UMI cf_UMI ref_UMI: 85244 989 603
# yes: 603 tandem:  chimera: 386
# barcode10
# barcode10 purative_UMI cf_UMI ref_UMI: 165568 2300 1518
# yes: 1518 tandem:  chimera: 782
# barcode11
# barcode11 purative_UMI cf_UMI ref_UMI: 152786 1877 1059
# yes: 1059 tandem:  chimera: 818
# barcode12
# barcode12 purative_UMI cf_UMI ref_UMI: 100789 1506 1076
# yes: 1076 tandem:  chimera: 430
# barcode13
# barcode13 purative_UMI cf_UMI ref_UMI: 190 1 1
# yes: 1 tandem:  chimera: 
  
