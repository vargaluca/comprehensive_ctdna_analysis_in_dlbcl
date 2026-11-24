suppressPackageStartupMessages(suppressWarnings(library(readr)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))

args <- commandArgs(trailingOnly = T)
table <- args[1]
outfile <- args[2]
id <- args[3]

# table <- read.csv("/Volumes/Extreme_SSD/comprehensive_ctdna_analysis_in_dlbcl/data/TIS-DL001-501_variants.csv", header = T, sep = ",")
table <- read.csv(table, header = T, sep = ",")

## read maf file, filter for SNPs
## 1: writing a loop: sort all variants by genomic coordinates and compare variants with next one. 

suppressWarnings(table <- table[order(as.numeric(gsub("chr", "", table[,"chromosome"])), table[,"position"]),]) #order variants by genomic coordinate

## looping through variants

table[,"distance"] <- NA
table[,"mutation_pattern"] <- NA
table[,"md_tag"] <- NA
table[,"flanking_tumor"] <- table[,"flanking_bps"]

table <- table %>%
  mutate(vaf = t_alt_count/t_depth * 100)

## loop contents: iterate through lines, grab position and substract it from next
## in line variant's position. If the length of this number is smaller or equal to 150 (it's 170
## in the article), put both lines into a new dataframe. Substract three, it'll be the number of flanking
## bases needed. Remove duplicate lines.
## in the end. 1-1 flanking bases for each variant (trinucleotides). 
## Combine using grep wildcards (ideally: trinucleotide1 - '.' wildcards in place of 
## bases inbetween - trinucleotide2), write these in a new column. Transform file
## with variant1, position, variant2, position, etc info. Save file.

position <- 0
next_row <- 0

df <- table
df <- df[F,]
df2 <- df

## loop through variants to get potential phased pairs into two dataframes if it's not emtpy

for (row in 1:nrow(table)) { #iterate for every line in dataframe
  position <- table[row, "position"]
  next_row <- row + 1
  position2 <- table[next_row, "position"]
  distance <- position2 - position #distance between two variants
  table[row, "distance"] <- distance 
  if (next_row <= nrow(table)) { #this was necessary so it doesn't throw an error message at the last row
    if ((table[row, "chromosome"] == table[next_row, "chromosome"]) && (distance <= 170) && distance >= 3) { #if on the same chromosome and over 3 and under 170 distance
      df <- rbind(df, table[row,]) #first variants
      df2 <- rbind(df2, table[next_row,]) #second variants 
    }
  }
}

if(nrow(df) != 0){
  print("variants found, processing")

## replace variant in flanking bases to tumor variant
for (row in 1:nrow(df)) {
  substr(df[row, "flanking_tumor"], 2, 2) <- df[row, "tumor_allele"]
}

for (row in 1:nrow(df2)) {
  substr(df2[row, "flanking_tumor"], 2, 2) <- df2[row, "tumor_allele"]
}

df2 <- rename_with(df2, ~paste0("2_", .), 4:length(df2))

df2.2 <- df2 %>%  select("2_position", "2_reference_allele", "2_tumor_allele", 
         "2_t_depth", "2_t_ref_count", "2_t_alt_count",
         "2_flanking_bps", "2_flanking_tumor", "2_vaf")

phased_data <- cbind(df, df2.2) #combine two tables
phased_data <- phased_data %>%
  select("chromosome", "gene", "entrez_id", "position", "2_position",    
         "reference_allele", "2_reference_allele", "tumor_allele", "2_tumor_allele",
        "t_depth", "2_t_depth", "t_ref_count", "2_t_ref_count",
         "t_alt_count", "2_t_alt_count", "vaf", "2_vaf", "flanking_bps", 
         "2_flanking_bps", "flanking_tumor", "2_flanking_tumor", "distance", 
         "mutation_pattern", "md_tag")

# creating mutation pattern
dots_number <- 0
for (row in 1:nrow(phased_data)) { #create patterns to search in bam files
  dots_number <- phased_data[row, "distance"] - 3
    phased_data[row, "mutation_pattern"] <- paste0(rep(".", times = dots_number), collapse = "")
  phased_data[row, "mutation_pattern"] <- paste0(phased_data[row, "flanking_tumor"], phased_data[row, "mutation_pattern"], phased_data[row, "2_flanking_tumor"])
}

# creating md-tag pattern
number <- 0
for (row in 1:nrow(phased_data)) { #create patterns to search in bam files
  number <- phased_data[row, "distance"] - 1
  ref1 <- phased_data[row, "reference_allele"]
  ref2 <- phased_data[row, "2_reference_allele"]
  phased_data[row, "md_tag"] <- paste0(ref1, number, ref2, collapse = "")
}

id_col <- rep(id, times = nrow(phased_data))
phased_data$sample_id <- id_col
phased_data <- phased_data %>%
  relocate(sample_id, .before = chromosome)
              
phased_data[phased_data == ""] <- NA
write_csv(phased_data, file = outfile, na = "NA")

}else{
  print("no variants found")
}


## Filtering step: samtools view and grep, moving lines to a new file if grep is
## able to find the pattern.
## samtools view recal_aligned_reads/redo_calls_maxevents10/TIS-DL142-015_1-0X_S11.recalibrated.bam 
## "chr1:144901100-144901200" | grep 'ACG....AAC'

