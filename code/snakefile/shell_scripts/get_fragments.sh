#!/usr/bin/env bash

# this script is for getting specific fragment data after extracting ctdna from larger bam files

echo input folder:
read path
echo output folder:
read out

mkdir -p $out
cd $path

module load samtools
for bamfile in *bam
do
id=$(echo $bamfile | cut -f 1 -d ".")
touch $out/$id.fragments.tsv

#bcl6_regions=$(samtools view -f 1 -F 0x400 $bamfile "chr3:187459936-187463582" | cut -f 9 | awk '{print sqrt($0^2)}')
#cxcr4_regions=$(samtools view -f 1 -F 0x400 $bamfile "chr2:136873818-136875814" | cut -f 9 | awk '{print sqrt($0^2)}')
#bcl2_regions=$(samtools view -f 1 -F 0x400 $bamfile "chr18:60981723-60986853" | cut -f 9 | awk '{print sqrt($0^2)}')
#bcl7a_regions=$(samtools view -f 1 -F 0x400 $bamfile "chr12:122459690-122461321" | cut -f 9 | awk '{print sqrt($0^2)}')
#tmsb4x_regions=$(samtools view -f 1 -F 0x400 $bamfile "chrX:12993125-12994482" | cut -f 9 | awk '{print sqrt($0^2)}')
#pim1_regions=$(samtools view -f 1 -F 0x400 $bamfile "chr6:37137828-37140123" | cut -f 9 | awk '{print sqrt($0^2)}')
#socs1_regions=$(samtools view -f 1 -F 0x400 $bamfile "chr16:11348988-11350159" | cut -f 9 | awk '{print sqrt($0^2)}')
#ighm upstream: chr14:106322324-106330947

samtools view -h -q 30 -f 81 -F 0x400 -L /disk/work/diagnostics/DLBC/lv1/DLBCL_liquid_K21/bed_files/promoters_bed.bed $bamfile > $out/tmp_bam1.bam
echo "tmp1 created"
samtools view -h -q 30 -f 93 -F 0x400 -L /disk/work/diagnostics/DLBC/lv1/DLBCL_liquid_K21/bed_files/promoters_bed.bed $bamfile > $out/tmp_bam2.bam
echo "tmp2 created"
samtools view -h -q 30 -f 97 -F 0x400 -L /disk/work/diagnostics/DLBC/lv1/DLBCL_liquid_K21/bed_files/promoters_bed.bed $bamfile > $out/tmp_bam3.bam
echo "tmp3 created"
samtools view -h -q 30 -f 99 -F 0x400 -L /disk/work/diagnostics/DLBC/lv1/DLBCL_liquid_K21/bed_files/promoters_bed.bed $bamfile > $out/tmp_bam4.bam
echo "tmp4 created"
samtools view -h -q 30 -f 145 -F 0x400 -L /disk/work/diagnostics/DLBC/lv1/DLBCL_liquid_K21/bed_files/promoters_bed.bed $bamfile > $out/tmp_bam5.bam
echo "tmp5 created"
samtools view -h -q 30 -f 147 -F 0x400 -L /disk/work/diagnostics/DLBC/lv1/DLBCL_liquid_K21/bed_files/promoters_bed.bed $bamfile > $out/tmp_bam6.bam
echo "tmp6 created"
samtools view -h -q 30 -f 161 -F 0x400 -L /disk/work/diagnostics/DLBC/lv1/DLBCL_liquid_K21/bed_files/promoters_bed.bed $bamfile > $out/tmp_bam7.bam
echo "tmp7 created"
samtools view -h -q 30 -f 163 -F 0x400 -L /disk/work/diagnostics/DLBC/lv1/DLBCL_liquid_K21/bed_files/promoters_bed.bed $bamfile > $out/tmp_bam8.bam
echo "tmp8 created"

#samtools view -h -q 30 -f 81 -F 0x400 $bamfile > $out/tmp_bam1.bam
#echo "tmp1 created"
#samtools view -h -q 30 -f 93 -F 0x400 $bamfile > $out/tmp_bam2.bam
#echo "tmp2 created"
#samtools view -h -q 30 -f 97 -F 0x400 $bamfile > $out/tmp_bam3.bam
#echo "tmp3 created"
#samtools view -h -q 30 -f 99 -F 0x400 $bamfile > $out/tmp_bam4.bam
#echo "tmp4 created"
#samtools view -h -q 30 -f 145 -F 0x400 $bamfile > $out/tmp_bam5.bam
#echo "tmp5 created"
#samtools view -h -q 30 -f 147 -F 0x400 $bamfile > $out/tmp_bam6.bam
#echo "tmp6 created"
#samtools view -h -q 30 -f 161 -F 0x400 $bamfile > $out/tmp_bam7.bam
#echo "tmp7 created"
#samtools view -h -q 30 -f 163 -F 0x400 $bamfile > $out/tmp_bam8.bam
#echo "tmp8 created"

samtools merge $out/$id.filtered.bam $out/tmp_bam1.bam \
$out//tmp_bam2.bam $out/tmp_bam3.bam $out/tmp_bam4.bam \
$out/tmp_bam5.bam $out/tmp_bam6.bam $out/tmp_bam7.bam \
$out/tmp_bam8.bam

samtools view -h $out/$id.filtered.bam | awk '!seen[$0]++' | samtools view -b - > $out/$id.filtered_fixed.bam
samtools index $out/$id.filtered_fixed.bam

#mv $path/$id.bed_filtered.bam $path/$id.filtered_fixed.bam
#mv $path/$id.bed_filtered.bam.bai $path/$id.filtered_fixed.bam.bai

echo "bam files merged"

rm $out/tmp_bam*
rm $out/*filtered.bam

	samtools view -f 1 -F 0x400 $out/$id.filtered_fixed.bam "chr3:187459936-187463582" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
	do
	echo "$id	BCL6	$line" >> $out/$id.fragments.tsv
	done

	samtools view -f 1 -F 0x400 $out/$id.filtered_fixed.bam "chr2:136873818-136875814" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo "$id	CXCR4	$line" >> $out/$id.fragments.tsv
        done

	samtools view -f 1 -F 0x400 $out/$id.filtered_fixed.bam "chr18:60981723-60986853" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo "$id	BCL2	$line" >> $out/$id.fragments.tsv
        done

	samtools view -f 1 -F 0x400 $out/$id.filtered_fixed.bam "chr12:122459690-122461321" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo "$id	BCL7A	$line" >> $out/$id.fragments.tsv
        done

	samtools view -f 1 -F 0x400 $out/$id.filtered_fixed.bam "chrX:12993125-12994482" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo "$id	TMSB4X	$line" >> $out/$id.fragments.tsv
        done

	samtools view -f 1 -F 0x400 $out/$id.filtered_fixed.bam "chr6:37137828-37140123" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo "$id	PIM1	$line" >> $out/$id.fragments.tsv
        done

	samtools view -f 1 -F 0x400 $out/$id.filtered_fixed.bam "chr16:11348988-11350159" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo "$id	SOCS1	$line" >> $out/$id.fragments.tsv
        done

        samtools view -f 1 -F 0x400 $out/$id.filtered_fixed.bam "chr14:106322324-106330947" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo "$id	IGHM	$line" >> $out/$id.fragments.tsv
        done

        samtools view -f 1 -F 0x400 $out/$id.filtered_fixed.bam "chr8:128748093-128752594" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo "$id	MYC	$line" >> $out/$id.fragments.tsv
        done


#        samtools view -f 1 -F 0x400 $out/$id.filtered_fixed.bam "chr14:106322324-106330947" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
#        do
#        echo "$id	$line" >> $out/$id.fragments.tsv
#        done

echo "$id finished"

done

