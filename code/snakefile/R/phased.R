suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(maftools))
suppressPackageStartupMessages(library(janitor))

args <- commandArgs(trailingOnly = T)
# maffile <- "/disk/work/diagnostics/DLBC/lv1/ncbi_data/maf_files/all_files/subsets_local/LS1065_maftools.maf"
maffile <- args[1]
outfile <- args[2]
# outfile <- "/disk/work/diagnostics/DLBC/lv1/DLBCL_liquid_K21/tables/phased_variants/lane_03/TIS-DL001-501.unfiltered_phased.tsv"
id <- args[3]
# id <- "TIS-DL001-501"

## read maf file, filter for SNPs

suppressWarnings(maf <- read.maf(maffile, verbose = F))

# maf@data$gnomadAF_AF[is.na(maf@data$gnomadAF_AF)] <- 0
# maf@data$AF[is.na(maf@data$AF)] <- 0
# maf@maf.silent$gnomadAF_AF[is.na(maf@data$gnomadAF_AF)] <- 0
# maf@maf.silent$AF[is.na(maf@data$AF)] <- 0


#maf@data <- maf@data %>%
# select(!c("Start_Position", "End_Position")) %>%
# dplyr::rename("Start_Position" = "i.Start_Position", "End_Position" = "i.End_Position")

#maf@maf.silent <- maf@maf.silent %>%
# select(!c("Start_Position", "End_Position")) %>%
# dplyr::rename("Start_Position" = "i.Start_Position", "End_Position" = "i.End_Position")

maf_data <- maf@data
maf_silent <- maf@maf.silent
maf_all <- rbind(maf_data, maf_silent, fill = T) %>%
  filter(Variant_Type == "SNP")
maf_all$gnomadAF_AF[is.na(maf_all$gnomadAF_AF)] <- 0
maf_all$AF[is.na(maf_all$AF)] <- 0

maf_all <- data.frame(lapply(maf_all, function(x) unlist(x)))
maf_all$AF <- as.double(maf_all$AF)
maf_all$gnomadAF_AF <- as.double(maf_all$gnomadAF_AF)

## 1: writing a loop: sort all variants by genomic coordinates and compare variants with next one. 

suppressWarnings(maf_all <- maf_all[order(as.numeric(gsub("chr", "", maf_all[,"Chromosome"])), maf_all[,"Start_Position"]),]) #order variants by genomic coordinate

## soft filtering variants

# maf_all_data <- remove_empty(maf_all, "cols") %>%
#   mutate(flanking_tumor = flanking_bps) %>%
#   mutate(VAF = t_alt_count/t_depth * 100) %>%
#   filter(t_alt_count >= 5 & t_depth >= 50) %>%
#   filter(!str_detect(FILTER, "contamination|normal_artifact|strand_bias")) %>%
#   select(Hugo_Symbol, Entrez_Gene_Id, Chromosome, Start_Position, Variant_Type,
#          Reference_Allele, Tumor_Seq_Allele2, dbSNP_RS, t_depth, t_ref_count, t_alt_count,
#          Consequence, flanking_bps, flanking_tumor, gnomadAF_AC, gnomadAF_AF, AF, VAF, FILTER) #filter for at least 5 reads and 50 depth and selected flags

## hard filtering variants

maf_all_data <- remove_empty(maf_all, "cols") %>%
 filter(FILTER == "PASS") %>%
 select(Hugo_Symbol, Entrez_Gene_Id, Chromosome, Start_Position, Variant_Type,
        Reference_Allele, Tumor_Seq_Allele2, dbSNP_RS, t_depth, t_ref_count, t_alt_count,
        Consequence, flanking_bps, gnomadAF_AC, gnomadAF_AF, AF, FILTER) %>%
 mutate(flanking_tumor = flanking_bps) %>%
 mutate(VAF = t_alt_count/t_depth * 100)
 # filter(t_alt_count >= 5 & t_depth >= 50)

# normally t_alt_count >= 5 & t_depth >= 50, changed it for ncbi data

# maf_all_data$gnomadAF_AF[is.na(maf_all_data$gnomadAF_AF)] <- 0
# maf_all_data$AF[is.na(maf_all_data$AF)] <- 0

maf_all_data <- maf_all_data %>%
  filter(gnomadAF_AF < 0.0001) %>%
  filter(AF < 0.01) #filter out gnomad<0.01 and ensembl MAF<0.0001
 
maf_all_data[,"distance"] <- NA
maf_all_data[,"mutation_pattern"] <- NA
maf_all_data[,"md_tag"] <- NA


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

df <- maf_all_data
df <- df[F,]
df2 <- df

## loop through variants to get potential phased pairs into two dataframes if it's not emtpy

for (row in 1:nrow(maf_all_data)) { #iterate for every line in dataframe
  position <- maf_all_data[row, "Start_Position"]
  next_row <- row + 1
  position2 <- maf_all_data[next_row, "Start_Position"]
  distance <- position2 - position #distance between two variants
  maf_all_data[row, "distance"] <- distance 
  if (next_row <= nrow(maf_all_data)) { #this was necessary so it doesn't throw an error message at the last row
    if ((maf_all_data[row, "Chromosome"] == maf_all_data[next_row, "Chromosome"]) && (distance <= 150) && distance >= 3) { #if on the same chromosome and over 3 and under 150 distance
      df <- rbind(df, maf_all_data[row,]) #first variants
      df2 <- rbind(df2, maf_all_data[next_row,]) #second variants 
    }
  }
}

if(nrow(df) != 0){
  print("variants found, processing")
  
## replace variant in flanking bases to tumor variant
for (row in 1:nrow(df)) {
  substr(df[row, "flanking_tumor"], 2, 2) <- df[row, "Tumor_Seq_Allele2"]
}

for (row in 1:nrow(df2)) {
  substr(df2[row, "flanking_tumor"], 2, 2) <- df2[row, "Tumor_Seq_Allele2"]
}

df2 <- rename_with(df2, ~paste0("2_", .), 3:length(df2))

df2.2 <- df2 %>%  select("2_Start_Position", "2_Reference_Allele", "2_Tumor_Seq_Allele2", 
         "2_dbSNP_RS", "2_t_depth", "2_t_ref_count", "2_t_alt_count", "2_Consequence",
         "2_flanking_bps", "2_flanking_tumor", "2_VAF", "2_FILTER")

phased_data <- cbind(df, df2.2) #combine two tables
phased_data <- phased_data %>%
  select("Hugo_Symbol", "Entrez_Gene_Id", "Chromosome", "Start_Position", "2_Start_Position",    
         "Reference_Allele", "2_Reference_Allele", "Tumor_Seq_Allele2", "2_Tumor_Seq_Allele2",
         "dbSNP_RS", "2_dbSNP_RS", "t_depth", "2_t_depth", "t_ref_count", "2_t_ref_count",
         "t_alt_count", "2_t_alt_count", "VAF", "2_VAF", "Consequence", "2_Consequence", "flanking_bps", 
         "2_flanking_bps", "flanking_tumor", "2_flanking_tumor", "FILTER", "2_FILTER", "distance", 
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
  ref1 <- phased_data[row, "Reference_Allele"]
  ref2 <- phased_data[row, "2_Reference_Allele"]
  phased_data[row, "md_tag"] <- paste0(ref1, number, ref2, collapse = "")
}

id_col <- rep(id, times = nrow(phased_data))
phased_data$sample_id <- id_col
phased_data <- phased_data %>%
  relocate(sample_id, .before = Hugo_Symbol)
              
phased_data[phased_data == ""] <- NA
write_delim(phased_data, file = outfile, na = "NA")

}else{
  print("no variants found")
}


## Filtering step: samtools view and grep, moving lines to a new file if grep is
## able to find the pattern.
## samtools view recal_aligned_reads/redo_calls_maxevents10/TIS-DL142-015_1-0X_S11.recalibrated.bam 
## "chr1:144901100-144901200" | grep 'ACG....AAC'


#filtered_phased <- read_delim("/disk/work/diagnostics/DLBC/lv1/DLBCL_genclass_K21/tables/phased_variants/TIS-DL142-015_1-0X_S11.filtered_phased.tsv")

