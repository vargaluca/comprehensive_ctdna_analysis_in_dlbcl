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
    -o|--bam)
      shift
      if test $# -gt 0; then
        export BAM=$1
      else
        echo "no output bam specified"
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

# unfiltered variant calling

Rscript='/disk/work/shared/tools/R-4.0.2/bin/Rscript'
scripts='/disk/work/users/lv1/DLBCL_genclass_K21/code/R'

if [[ $PROCESS == "call" ]]; then

	echo "variant calling"

	mkdir -p $OUTDIR

	table=$(basename "$INPUT")
	sample_id=$(echo $table | cut -f 1 -d ".")
	outfile=$OUTDIR/$sample_id-unfiltered_pv.csv

        if [ ! -f $path/tables/phased_variants/$run/$sample_id.hard_unfiltered_phased.tsv ]
        then
        echo "processing $sample_id"
#        $Rscript $scripts/phased.R $table $outfile $sample_id #running phased.R script
        if [ -f $outfile ] #check whether unfiltered file exists and isn't empty, return appropriate message
        then
                if [ -s $outfile ]
                then
                        unfiltered=$(tail -n +2 $outfile | wc -l)
                        echo "$sample_id unfiltered variants successfully created with $unfiltered variant pairs."
                fi
        else
                echo "$sample_id unfiltered variants file does not exist."
        fi
        fi
exit 0
elif [[ $PROCESS == "filter" ]]; then
	echo "variant filtering yay"
exit 0
else
echo "you stupid fuck"
exit 0
fi
Rscript='/disk/work/shared/tools/R-4.0.2/bin/Rscript'
scripts='/disk/work/users/lv1/DLBCL_genclass_K21/code/R'

cd $path/maf_files/$run/subsets_local_liq_ranges/
mkdir -p $path/tables/phased_variants/$run

table=$path/maf_files/$run/subsets_local_liq_ranges/$maf
sample_id=$(echo $maf | cut -f 1 -d "_")
outfile=$path/tables/phased_variants/$run/$sample_id.hard_unfiltered_phased.tsv
sample=$(echo $maf | cut -f 1 -d "_")

	if [ ! -f $path/tables/phased_variants/$run/$sample_id.hard_unfiltered_phased.tsv ]
	then
	echo "processing $maf"

	$Rscript $scripts/phased.R $maffile $outfile $sample #running phased.R script
	if [ -f $outfile ] #check whether unfiltered file exists and isn't empty, return appropriate message
	then
		if [ -s $outfile ]
		then
			unfiltered=$(tail -n +2 $outfile | wc -l)
			echo "$sample_id unfiltered variants successfully created with $unfiltered variant pairs."
		fi
	else
		echo "$sample_id unfiltered variants file does not exist."
	fi
	fi


cd $path/tables/phased_variants/$run

module load samtools
for table in *hard_unfiltered*tsv #filtering cis variants
do
sample=$(echo $table | cut -f 1 -d ".")
if [ ! -f $path/tables/phased_variants/$run/$sample.hard_filtered_phased.tsv ]
then
echo "filtering cis variants from $table"
touch $path/tables/phased_variants/$run/$sample.hard_filtered_phased.tsv
echo "$(head -n 1 $table) reads_from_bam" >> $path/tables/phased_variants/$run/$sample.hard_filtered_phased.tsv

	cat $table | tail -n +2 | while read -r line #for every line in unfiltered phased file perform this whole thing
	do
		pattern=$(echo $line | cut -f 30 -d " ") #mutation pattern
		md_tag=$(echo $line | cut -f 30 -d " ") #md-tag
		chrom=$(echo $line | cut -f 4 -d " ") #chromosome
		start=$(($(echo $line | cut -f 5 -d " ") - 100)) #genomic coordinate of first variant, minus 100
		end=$(($(echo $line | cut -f 6 -d " ") + 100)) #genomic coordinate of second variant, plus 100
		region="$chrom:$start-$end" #region to search in bam file for cis variants
		bamfile=$(ls $path/recal_aligned_reads/$run/$sample*.recalibrated.bam) #bam file to be searched
		if [[ $(samtools view -f 1 $bamfile $region | grep $pattern | cut -f 1-25 | grep $md_tag | grep -v "$pattern.*$pattern" | wc -l) > 0 ]] #if grep pattern returns at least one line, write the variants to a new file (plus number of reads)
			then
			reads=$(samtools view -f 1 $bamfile $region | grep $pattern | cut -f 1-25 | grep $md_tag | grep -v "$pattern.*$pattern" | wc -l)
			echo "$line $reads" >> $path/tables/phased_variants/$run/$sample.hard_filtered_phased.tsv
		fi
	done
fi

if [ -f $path/tables/phased_variants/$run/$sample.hard_filtered_phased.tsv ] #check whether filtered file exists and isn't empty, return appropriate message
	then
	if [ -s $path/tables/phased_variants/$run/$sample.hard_filtered_phased.tsv ]
		then
		unfiltered=$(tail -n +2 $table | wc -l)
		filtered=$(tail -n +2 $path/tables/phased_variants/$run/$sample.hard_filtered_phased.tsv | wc -l)
		out=$(($unfiltered - $filtered)) #number of variants thrown out
		echo "$sample filtered variants successfully created with $filtered variant pairs."
		echo "$out variant pairs thrown out."
	fi
	else
        echo "$sample filtered variants file does not exist."
fi
done
module unload samtools
