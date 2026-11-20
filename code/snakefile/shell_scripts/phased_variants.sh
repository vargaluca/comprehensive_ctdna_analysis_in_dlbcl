#!/usr/bin/env bash

# parameters and process specifications
if [ $# -eq 0 ]; then
        echo "No parameters specified."
        echo " "
        echo "phased variant capture tool"
        echo "github: https://github.com/vargaluca/comprehensive_ctdna_analysis_in_dlbcl"
        echo " "
        echo "usage:"
        echo "phased_variants.sh -o output_dir -i input_file -p call/filter"
        echo " "
        echo "options:"
        echo "-h, --help                show brief help"
        echo "-i, --input               specify input tsv file"
	echo "-b, --bam			specify corresponding input bam file"
        echo "-o, --outdir              specify output directory"
        echo "-p, --process             'call' or 'filter'"
exit 1
else
while test $# -gt 0; do
  case "$1" in
    -h|--help)
        echo " "
        echo "phased variant capture tool"
        echo "github: https://github.com/vargaluca/comprehensive_ctdna_analysis_in_dlbcl"
        echo " "
        echo "usage:"
        echo "phased_variants.sh -o output_dir -i input_file -p call/filter"
        echo " "
        echo "options:"
        echo "-h, --help                show brief help"
        echo "-i, --input               specify input tsv file"
        echo "-b, --bam                 specify corresponding input bam file"
        echo "-o, --outdir              specify output directory"
        echo "-p, --process             'call' or 'filter'"
      exit 0
      ;;
    -i|--input)
      shift
      if test $# -gt 0; then
        export INPUT=$1
      else
        echo "no input tsv specified"
        exit 1
      fi
      shift
      ;;
    -o|--outdir)
      shift
      if test $# -gt 0; then
        export OUTDIR=$1
      else
        echo "no output directory specified"
        exit 1
      fi
      shift
      ;;
    -b|--bam)
      shift
      if test $# -gt 0; then
        export BAM=$1
      else
        echo "no bam file specified"
        exit 1
      fi
      shift
      ;;
    -t|--ct-variants)
      shift
      if test $# -gt 0; then
        export CTDNA=$1
      else
        echo "no ctDNA variant list specified"
        exit 1
      fi
      shift
      ;;
    -p|--process)
      shift
      if test $# -gt 0; then
        export PROCESS=$1
        else
        echo "no process specified"
        exit 1
      fi
      shift
      ;;
    *)
      break
      ;;
  esac
done
fi

#configuration
wd=$(pwd)
config=$wd/config.txt

Rscript=$(grep "Rscript" $config | cut -f 2 -d ":")
samtools=$(grep "samtools" $config | cut -f 2 -d ":")

# functions
call() {
mkdir -p $OUTDIR

table=$(basename "$INPUT")
sample_id=$(echo $table | cut -f 1 -d ".")
outfile=$OUTDIR/$sample_id.unfiltered_pv.csv

if [ -f $INPUT ]; then
  echo "processing $sample_id"
  $Rscript $wd/R/phased.R $table $outfile $sample_id #running phased.R script

      if [ -f $outfile ]; then #check whether unfiltered file exists and isn't empty, return appropriate message

          if [ -s $outfile ]; then
            unfiltered=$(tail -n +2 $outfile | wc -l)
            echo "$outfile variants successfully created with $unfiltered variant pairs."
          fi

      else
        echo "$sample_id unfiltered variants file does not exist."
      fi

else
  echo "variant table not found"
fi
}

filter() {
table=$(basename "$INPUT")
sample_id=$(echo $table | cut -f 1 -d ".")
outfile=$OUTDIR/$sample_id.filtered_pv.csv
bamfile=$BAM

if [ -f $INPUT ]; then

    if [ -f $BAM ]; then
      echo "filtering cis variants from $table"
      touch $outfile
      echo "$(head -n 1 $table) reads_from_bam" >> $outfile

          cat $table | tail -n +2 | while read -r line #for every line in unfiltered phased file perform this whole thing
          do
            pattern=$(echo $line | cut -f 30 -d ",") #mutation pattern
            md_tag=$(echo $line | cut -f 30 -d ",") #md-tag
            chrom=$(echo $line | cut -f 4 -d ",") #chromosome
            start=$(($(echo $line | cut -f 5 -d ",") - 100)) #genomic coordinate of first variant, minus 100
            end=$(($(echo $line | cut -f 6 -d ",") + 100)) #genomic coordinate of second variant, plus 100
            region="$chrom:$start-$end" #region to search in bam file for cis variants

                if [[ $($samtools view -f 1 $bamfile $region | grep $pattern | cut -f 1-25 | grep $md_tag | grep -v "$pattern.*$pattern" | wc -l) > 0 ]]; then #if grep pattern returns at
                  reads=$($samtools view -f 1 $bamfile $region | grep $pattern | cut -f 1-25 | grep $md_tag | grep -v "$pattern.*$pattern" | wc -l)
                  echo "$line,$reads" >> $outfile
                fi
          done

              if [ -f $outfile ]; then #check whether filtered file exists and isn't empty, return appropriate message
                  if [ -s $outfile ]; then
                    unfiltered=$(tail -n +2 $table | wc -l)
                    filtered=$(tail -n +2 $outfile | wc -l)
                    out=$(($unfiltered - $filtered)) #number of variants thrown out
                    echo "$sample filtered variants successfully created with $filtered variant pairs."
                    echo "$out variant pairs thrown out."
                  else
                    echo "$sample filtered variants file does not exist."
              fi

    else
      echo "bam file not found"
    fi

else
  echo "unfiltered pv file not found"
fi
}

extract() {
sample=$(echo $INPUT | cut -f 1 -d ".")
bamfile=$($BAM)
liq_table=$($CTDNA)

if [ -f $bamfile ]; then
    if [ -f $INPUT ]; then
      $samtools view -H $bamfile > $OUTDIR/$sample.ctdna.sam

          cat $INPUT | tail -n +2 | while read -r line
          do
            pattern=$(echo $line | cut -f 30 -d " ")
            md_tag=$(echo $line | cut -f 31 -d " ")
            chrom=$(echo $line | cut -f 4 -d " ")
            start=$(($(echo $line | cut -f 5 -d " ") - 100))
            end=$(($(echo $line | cut -f 6 -d " ") + 100))
            region="$chrom:$start-$end"
            depth1=$(echo $line | cut -f 13 -d " ")
            depth2=$(echo $line | cut -f 14 -d " ")

              if [[ $(awk -v depth1="$depth1" -v depth2="$depth2" 'BEGIN {if (depth1 == "NA" && depth2 == "NA") {print "TRUE"} else {print "FALSE"}}') == "FALSE" ]]; then
                $samtools view -f 1 -F 0x400 $bamfile $region | grep $pattern | grep $md_tag >> $outfolder/$sample.ctdna.sam
              fi
          done

      $samtools view -h $outfolder/$sample.wt.sam | awk '!seen[$0]++' | samtools view -b - > $outfolder/$sample.wt.bam
      $samtools sort -o $outfolder/$sample.sorted.wt.bam $outfolder/$sample.wt.bam
      $samtools index $outfolder/$sample.sorted.wt.bam
      rm $outfolder/$sample.sam $outfolder/$sample.wt.bam

    else
      echo "input pv file not found"
else
  echo "input bam file not found"
fi
}

# starting defined process function

if [[ $PROCESS == "call" ]]; then
  call
  filter
  exit 0
elif [[ $PROCESS == "filter" ]]; then
  filter
  exit 0
elif [[ $PROCESS == "extract" ]]; then
  extract
  exit 0
elif [[ $PROCESS == "recall" ]]; then
  recall
  exit 0
else
  echo "$PROCESS not recognised"
  exit 1
fi

