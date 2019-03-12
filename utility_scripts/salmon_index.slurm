#!/usr/bin/env bash

##############
## Extract Fasta
## - use gffread submodule in cufflinks
## - processes plus and negative strand genes correctly, producing fasta file with nucleotides in transcript order

module load cufflinks/2.2.1

gffread -h

gffread Homo_sapiens.GRCh38.95.ucsc.gtf -g ${GENOME} -w Homo_sapiens.GRCh38.95.ucsc.fasta

# Check the number of unique transcripts produced to ensure entire GTF was converted
grep "^>" Homo_sapiens.GRCh38.95.ucsc.fasta | cut -d" " -f1 | sort | uniq | wc -l
# 206601

cut -f9 Homo_sapiens.GRCh38.95.ucsc.gtf | grep "transcript_id" | sed 's/; transcript_id "/\'$'\t''/g' | sed 's/"; transcript_version/\'$'\t''/g' | cut -f2 | sort | uniq | wc -l
# 206601

### The number of transcripts match between the input GTF and output fasta file

##############
## Make Salmon Index files
##
##

module load salmon/0.12.0
module load gcc/8.2.0

# Create index
salmon -help

salmon index -t Homo_sapiens.GRCh38.95.ucsc.fasta -i salmonIndex_Homo_sapiens.GRCh38.95.ucsc --type quasi -k 31