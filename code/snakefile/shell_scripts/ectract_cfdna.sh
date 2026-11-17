#!/usr/bin/env bash

echo path to liquid bams:
read liq_bam
echo output folder:
read outfolder

cd $liq_bam
mkdir -p $outfolder

module load samtools
module load sambamba
for bamfile in *bam
do
id=$(echo $bamfile | cut -f 1 -d "_")
sambamba view -h $bamfile -L /disk/work/diagnostics/DLBC/lv1/DLBCL_liquid_K21/bed_files/promoters_bed.bed > $outfolder/regions_tmp.bam

echo "bed filtered bam complete"
samtools view -h -q 30 -f 81 -F 0x400 $outfolder/regions_tmp.bam > $outfolder/tmp_bam1.bam
echo "tmp1 created"
samtools view -h -q 30 -f 93 -F 0x400 $outfolder/regions_tmp.bam > $outfolder/tmp_bam2.bam
echo "tmp2 created"
samtools view -h -q 30 -f 97 -F 0x400 $outfolder/regions_tmp.bam > $outfolder/tmp_bam3.bam
echo "tmp3 created"
samtools view -h -q 30 -f 99 -F 0x400 $outfolder/regions_tmp.bam > $outfolder/tmp_bam4.bam
echo "tmp4 created"
samtools view -h -q 30 -f 145 -F 0x400 $outfolder/regions_tmp.bam > $outfolder/tmp_bam5.bam
echo "tmp5 created"
samtools view -h -q 30 -f 147 -F 0x400 $outfolder/regions_tmp.bam > $outfolder/tmp_bam6.bam
echo "tmp6 created"
samtools view -h -q 30 -f 161 -F 0x400 $outfolder/regions_tmp.bam > $outfolder/tmp_bam7.bam
echo "tmp7 created"
samtools view -h -q 30 -f 163 -F 0x400 $outfolder/regions_tmp.bam > $outfolder/tmp_bam8.bam
echo "tmp8 created"

samtools merge $outfolder/$id.filtered_fixed.bam $outfolder/tmp_bam1.bam \
$outfolder/tmp_bam2.bam $outfolder/tmp_bam3.bam $outfolder/tmp_bam4.bam \
$outfolder/tmp_bam5.bam $outfolder/tmp_bam6.bam $outfolder/tmp_bam7.bam \
$outfolder/tmp_bam8.bam
rm $outfolder/*tmp*

samtools index $outfolder/$id.filtered_fixed.bam

touch $outfolder/$id.fragments.tsv
echo "sample_id"$'\t'"gene"$'\t'"length" > $outfolder/$id.fragments.tsv

        samtools view -f 1 -F 0x400 $outfolder/$id.filtered_fixed.bam "chr3:187459936-187463582" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo $id$'\t'"BCL6"$'\t'$line >> $outfolder/$id.fragments.tsv
        done

        samtools view -f 1 -F 0x400 $outfolder/$id.filtered_fixed.bam "chr2:136873818-136875814" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo $id$'\t'"CXCR4"$'\t'$line >> $outfolder/$id.fragments.tsv
        done

        samtools view -f 1 -F 0x400 $outfolder/$id.filtered_fixed.bam "chr18:60981723-60986853" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo $id$'\t'"BCL2"$'\t'$line >> $outfolder/$id.fragments.tsv
        done

        samtools view -f 1 -F 0x400 $outfolder/$id.filtered_fixed.bam "chr12:122459690-122461321" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo $id$'\t'"BCL7A"$'\t'$line >> $outfolder/$id.fragments.tsv
        done

        samtools view -f 1 -F 0x400 $outfolder/$id.filtered_fixed.bam "chrX:12993125-12994482" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo $id$'\t'"TMSB4X"$'\t'$line >> $outfolder/$id.fragments.tsv
        done

        samtools view -f 1 -F 0x400 $outfolder/$id.filtered_fixed.bam "chr6:37137828-37140123" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo $id$'\t'"PIM1"$'\t'$line >> $outfolder/$id.fragments.tsv
        done

        samtools view -f 1 -F 0x400 $outfolder/$id.filtered_fixed.bam "chr16:11348988-11350159" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo $id$'\t'"SOCS1"$'\t'$line >> $outfolder/$id.fragments.tsv
        done

        samtools view -f 1 -F 0x400 $outfolder/$id.filtered_fixed.bam "chr14:106322324-106330947" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo $id$'\t'"IGHM"$'\t'$line >> $outfolder/$id.fragments.tsv
        done

        samtools view -f 1 -F 0x400 $outfolder/$id.filtered_fixed.bam "chr8:128748093-128752594" | cut -f 9 | awk '{print sqrt($0^2)}' | while read -r line
        do
        echo $id$'\t'"MYC"$'\t'$line >> $outfolder/$id.fragments.tsv
        done

echo "$id finished"

done

done
module unload samtools
module unload sambamba
