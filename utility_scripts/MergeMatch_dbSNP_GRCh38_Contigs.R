#!/usr/bin/Rscript

# Load required libraries
library(tidyverse)

# Load the two files
DBSNP <- read_tsv("dbSNP_b152_MetaData.txt")
ASSEMBLY <- read_tsv("reference_fasta_metadata.txt")

# Merge the two files by the Assembly_GI and contig length as a cross check
MERGE <- full_join(ASSEMBLY, DBSNP, by = c("ASSEMBLY_GI", "SLEN"))

# Write out the full file for reference
write_tsv(MERGE, "dbSNP152_GRCh38_Full_CrossMap.txt")

## Write out a list of contigs that only exist in dbSNP that can be removed from the annotation VCF
# Extract only thos lines in the dbSNP file with no matching contig in the assembly used for alignment
DBSNP_ONLY <- MERGE %>% filter(is.na(Contig)) %>% select(CONTIG)
write_tsv(DBSNP_ONLY, "contigs_2_remove_from_dbSNP152.txt", col_names = FALSE)
# Create a bed file version as the chromosome list does not seem to work
DBSNP_ONLY <- MERGE %>% filter(is.na(Contig)) %>% select(CONTIG, SLEN) %>% add_column(Start = 0, .after = "CONTIG")
write_tsv(DBSNP_ONLY, "contigs_2_remove_from_dbSNP152.bed", col_names = FALSE)

## Generate a bcftools compatible rename key
# Extract the lines with matching contigs in dbSNP file and the GRCh38 assembly used for alignment
RENAME_KEY <- MERGE %>% filter(!is.na(Contig) & !is.na(CONTIG)) %>% select(CONTIG, Contig)
write_delim(RENAME_KEY, "GRCh38_dbSNP152_2_UCSC_Contigs.txt", delim = " ", col_names = FALSE)
