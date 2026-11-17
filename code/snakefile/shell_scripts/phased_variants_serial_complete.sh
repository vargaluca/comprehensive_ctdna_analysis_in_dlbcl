#!/usr/bin/env bash

echo path to tissue phased files:
read tis_files
echo path to liquid bams:
read liq_bam
echo path to liquid mafs:
read liq_maf
echo output folder:
read outfolder

Rscript='/disk/work/shared/tools/R-4.0.2/bin/Rscript'
scripts='/disk/work/users/lv1/DLBCL_genclass_K21/code/R'

cd $tis_files
mkdir -p $outfolder

module load samtools

for file in *_filtered*
do

	sample=$(echo $file | cut -f 1 -d ".")
	echo $sample
	patient=$(echo $file | cut -f 2 -d "-")
	echo $patient
	bamfile=$(ls $liq_bam/*$patient*bam)
	echo $bamfile
	maffile=$(ls $liq_maf/*$patient*maf)
	echo $maffile
	liq_sample=$(echo $maffile | cut -f 12 -d "/" | cut -f 1 -d "_")
	echo $liq_sample

	if [ -f $bamfile ]
		then
#		if [ ! -f $outfolder/$liq_sample.recurring-phased.tsv ]
#		then
		echo "getting recurring phased variants from $liq_sample"
		headers=$(head -n 1 $file | cut -f 1-30)

		touch $outfolder/$liq_sample.recurring-phased.tsv
		echo "$headers reads_from_liquid_bam bam_depth" >> $outfolder/$liq_sample.recurring-phased.tsv

		cat $file | tail -n +2 | while read -r line
		do
			pattern=$(echo $line | cut -f 30 -d " ")
			md_tag=$(echo $line | cut -f 31 -d " ")
			chrom=$(echo $line | cut -f 4 -d " ")
			start=$(($(echo $line | cut -f 5 -d " ") - 100))
			end=$(($(echo $line | cut -f 6 -d " ") + 100))
			region="$chrom:$start-$end"
				if [[ $(samtools view -f 1 -F 0x400 $bamfile $region | grep $pattern | cut -f 1-25 | grep $md_tag | grep -v $pattern.*$pattern | wc -l) > 0 ]]
				then
#					echo "getting $region from $bamfile"
					line_start=$(echo $line | cut -f 2-12 -d " ")
					line_mid=$(echo $line | cut -f 21-26 -d " ")
					line_end=$(echo $line | cut -f 29-32 -d " ")

					gene=$(echo $line | cut -f 2 -d " ")
                                        start1=$(echo $line | cut -f 5 -d " ")
                                        start2=$(echo $line | cut -f 6 -d " ")
					ref1=$(echo $line | cut -f 7 -d " ")
                                        ref2=$(echo $line | cut -f 8 -d " ")
                                        alt1=$(echo $line | cut -f 9 -d " ")
                                        alt2=$(echo $line | cut -f 10 -d " ")

					variant_line1=$(awk -v gene="$gene" -v start1="$start1" -v ref1="$ref1" -v alt1="$alt1" -F \
					'\t' '{if ($2 == gene && $6 == start1 && $12 == ref1 && $13 == alt1) print $0}' $maffile)
					t_depth1=$(awk -v gene="$gene" -v start1="$start1" -v ref1="$ref1" -v alt1="$alt1" -F \
                                        '\t' '{if ($2 == gene && $6 == start1 && $12 == ref1 && $13 == alt1) print $40}' $maffile)
					t_ref1=$(awk -v gene="$gene" -v start1="$start1" -v ref1="$ref1" -v alt1="$alt1" -F \
                                        '\t' '{if ($2 == gene && $6 == start1 && $12 == ref1 && $13 == alt1) print $41}' $maffile)
					t_alt1=$(awk -v gene="$gene" -v start1="$start1" -v ref1="$ref1" -v alt1="$alt1" -F \
                                        '\t' '{if ($2 == gene && $6 == start1 && $12 == ref1 && $13 == alt1) print $42}' $maffile)
					filter1=$(awk -v gene="$gene" -v start1="$start1" -v ref1="$ref1" -v alt1="$alt1" -F \
                                        '\t' '{if ($2 == gene && $6 == start1 && $12 == ref1 && $13 == alt1) print $101}' $maffile)
#					echo $t_depth1 $t_ref1 $t_alt1 $filter1

					vaf1=$(awk -v gene="$gene" -v start1="$start1" -v ref1="$ref1" -v alt1="$alt1" -F \
                                        '\t' '{if ($2 == gene && $6 == start1 && $12 == ref1 && $13 == alt1) print ($42/$40)*100}' $maffile)
#					echo $vaf1

					variant_line2=$(awk -v gene="$gene" -v start2="$start2" -v ref2="$ref2" -v alt2="$alt2" -F \
                                        '\t' '{if ($2 == gene && $6 == start2 && $12 == ref2 && $13 == alt2) print $0}' $maffile)
                                        t_depth2=$(awk -v gene="$gene" -v start2="$start2" -v ref2="$ref2" -v alt2="$alt2" -F \
                                        '\t' '{if ($2 == gene && $6 == start2 && $12 == ref2 && $13 == alt2) print $40}' $maffile)
                                        t_ref2=$(awk -v gene="$gene" -v start2="$start2" -v ref2="$ref2" -v alt2="$alt2" -F \
                                        '\t' '{if ($2 == gene && $6 == start2 && $12 == ref2 && $13 == alt2) print $41}' $maffile)
                                        t_alt2=$(awk -v gene="$gene" -v start2="$start2" -v ref2="$ref2" -v alt2="$alt2" -F \
                                        '\t' '{if ($2 == gene && $6 == start2 && $12 == ref2 && $13 == alt2) print $42}' $maffile)
                                        filter2=$(awk -v gene="$gene" -v start2="$start2" -v ref2="$ref2" -v alt2="$alt2" -F \
                                        '\t' '{if ($2 == gene && $6 == start2 && $12 == ref2 && $13 == alt2) print $101}' $maffile)
					echo $t_depth2 $t_ref2 $t_alt2 $filter2

                                        vaf2=$(awk -v gene="$gene" -v start2="$start2" -v ref2="$ref2" -v alt2="$alt2" -F \
                                        '\t' '{if ($2 == gene && $6 == start2 && $12 == ref2 && $13 == alt2) print ($42/$40)*100}' $maffile)
#					echo $vaf2

					if [[ -z $t_depth1 ]]
					then
						if [[ -z $t_depth2 ]]
						then
						t_depth1="NA" t_depth2="NA" t_ref1="NA" t_ref2="NA" t_alt1="NA" t_alt2="NA" filter1="NA" filter2="NA" vaf1=0 vaf2=0
						elif [[ ! -z $t_depth2 ]]
						then
						t_depth1="NA" t_ref1="NA" t_alt1="NA" filter1="NA" vaf1=0
						fi

					elif [[ ! -z $t_depth1 ]]
					then
						if [[ -z $t_depth2 ]]
						then
						t_depth2="NA" t_ref2="NA" t_alt2="NA" filter2="NA" vaf2=0
						fi
					fi

					reads=$(samtools view -f 1 -F 0x400 $bamfile $region | grep $pattern | cut -f 1-25 | grep $md_tag | grep -v $pattern.*$pattern | wc -l)
					bam_depth=$(samtools view -f 1 -F 0x400 $bamfile $region | wc -l)
#					echo $reads
					full_line="${liq_sample} ${line_start} ${t_depth1} ${t_depth2} ${t_ref1} ${t_ref2} ${t_alt1} ${t_alt2} ${vaf1} ${vaf2} ${line_mid} ${filter1} ${filter2} ${line_end} ${reads} ${bam_depth}"
					fields=$(echo $full_line | grep -o ' ' | wc -l)

#					if [[ $(awk -v vaf1="$vaf1" -v vaf2="$vaf2" 'BEGIN {if (vaf1 <= 25 && vaf2 <= 25) print "TRUE"}') == "TRUE" ]]
#					then
#						if [[ $(awk -v vaf1="$vaf1" -v vaf2="$vaf2" 'BEGIN {if (vaf1 == 0 && vaf2 == 0) print "TRUE"}') == "TRUE" ]]
#						then
#							if [[ $(awk -v reads="$reads" -v bam_depth="$bam_depth" 'BEGIN {if (((reads/bam_depth)*100) <= 25) print "TRUE" }') == "TRUE" ]]
#							then
#							echo "none are there"
#							awk -v reads="$reads" -v bam_depth="$bam_depth" 'BEGIN {if (((reads/bam_depth)*100) <= 25) print (reads/bam_depth)*100 }'
#							echo "$vaf1 $vaf2"
							echo $full_line >> $outfolder/$liq_sample.recurring-phased.tsv
							echo $full_line
#							else
#							echo "does not fit the criteria"
#							fi
#						elif [[ $(awk -v vaf1="$vaf1" -v vaf2="$vaf2" 'BEGIN {if (vaf1 > 0 || vaf2 > 0) print "TRUE"}') == "TRUE" ]]
#						then
#						echo "at least one is there $vaf1 $vaf2"
#						echo $full_line
#                                               echo $full_line >> $outfolder/$liq_sample.recurring-phased.tsv
#						fi
#					elif [[ $(awk -v reads="$reads" -v bam_depth="$bam_depth" 'BEGIN {if (((reads/bam_depth)*100) <= 25) print "TRUE" }') == "TRUE" ]]
#					then
#					echo $full_line >> $outfolder/$liq_sample.recurring-phased.tsv
#					fi
				fi
		done
	fi
#	fi
#	if [ -f $outfolder/$liq_sample.recurring-phased.tsv ]
#		then
#		if [ -s $outfolder/$liq_sample.recurring-phased.tsv ]
#			then
#			tissue=$(tail -n +2 $file | wc -l)
#			recurring=$(tail -n +2 $outfolder/$liq_sample.recurring-phased.tsv | wc -l)
#			out=$(($tissue - $recurring))
#			echo "$recurring number of phased variants detected in $liq_sample"
#			echo "$out variant pairs not detected."
#		fi
#		else
#	       	echo "$sample filtered variants file does not exist."
#	fi
done
module unload samtools
