#!/usr/bin/env bash

#this script is for the extraction of ctdna fragments using phased variant data

echo path to liq phased files:
read tis_files
echo path to liquid bams:
read liq_bam
echo output folder:
read outfolder

cd $tis_files
mkdir -p $outfolder

module load samtools

for file in LIQ*tsv
do

        sample=$(echo $file | cut -f 1 -d ".")
        echo $sample
        patient=$(echo $file | cut -f 2 -d "-")
        echo $patient
        bamfile=$(ls $liq_bam/*$patient*bam)
        echo $bamfile
        liq_sample=$(echo $maffile | cut -f 12 -d "/" | cut -f 1 -d "_")
        echo $liq_sample

        if [ -f $bamfile ]
                then
		samtools view -H $bamfile > $outfolder/$sample.wt.sam
		cat $file | tail -n +2 | while read -r line
                do
                        pattern=$(echo $line | cut -f 30 -d " ")
                        md_tag=$(echo $line | cut -f 31 -d " ")
                        chrom=$(echo $line | cut -f 4 -d " ")
                        start=$(($(echo $line | cut -f 5 -d " ") - 100))
                        end=$(($(echo $line | cut -f 6 -d " ") + 100))
                        region="$chrom:$start-$end"
			depth1=$(echo $line | cut -f 13 -d " ")
			depth2=$(echo $line | cut -f 14 -d " ")

			if [[ $(awk -v depth1="$depth1" -v depth2="$depth2" 'BEGIN {if (depth1 == "NA" && depth2 == "NA") {print "TRUE"} else {print "FALSE"}}') == "FALSE" ]]
				then
				samtools view -f 1 -F 0x400 $bamfile $region | grep -v $pattern | grep -v $md_tag >> $outfolder/$sample.wt.sam
			fi
		done

		samtools view -h $outfolder/$sample.wt.sam | awk '!seen[$0]++' | samtools view -b - > $outfolder/$sample.wt.bam
		samtools sort -o $outfolder/$sample.sorted.wt.bam $outfolder/$sample.wt.bam
		samtools index $outfolder/$sample.sorted.wt.bam
		rm $outfolder/$sample.sam $outfolder/$sample.wt.bam
	fi
done
