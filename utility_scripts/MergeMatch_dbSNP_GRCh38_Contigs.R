#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

# Load required libraries
library(tidyverse)

# Load the two files
DBSNP <- read_tsv(paste0("dbSNP_b",args[1],"_MetaData.txt"))
ASSEMBLY <- read_tsv("reference_fasta_metadata.txt")

# Merge the two files by the Assembly_GI and contig length as a cross check
MERGE <- full_join(ASSEMBLY, DBSNP, by = c("ASSEMBLY_GI", "SLEN"))

# Write out the full file for reference
write_tsv(MERGE, paste0("dbSNP",args[1],"_GRCh38_Full_CrossMap.txt"))

## Write out a list of contigs that only exist in dbSNP that can be removed from the annotation VCF
# Extract only thos lines in the dbSNP file with no matching contig in the assembly used for alignment
DBSNP_ONLY <- MERGE %>% filter(is.na(Contig)) %>% select(CONTIG)
write_tsv(DBSNP_ONLY, paste0("contigs_2_remove_from_dbSNP",args[1],".txt"), col_names = FALSE)
# Create a bed file version as the chromosome list does not seem to work
DBSNP_ONLY <- MERGE %>% filter(is.na(Contig)) %>% select(CONTIG, SLEN) %>% add_column(Start = 0, .after = "CONTIG")
write_tsv(DBSNP_ONLY, paste0("contigs_2_remove_from_dbSNP",args[1],".bed"), col_names = FALSE)

## Generate a bcftools compatible rename key
# Extract the lines with matching contigs in dbSNP file and the GRCh38 assembly used for alignment
RENAME_KEY <- MERGE %>% filter(!is.na(Contig) & !is.na(CONTIG)) %>% select(CONTIG, Contig)
write_delim(RENAME_KEY, paste0("GRCh38_dbSNP",args[1],"_2_UCSC_Contigs.txt"), delim = " ", col_names = FALSE)
