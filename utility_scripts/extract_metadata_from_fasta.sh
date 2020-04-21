#!/usr/bin/env bash

# Goal: Extract the meta-data from the GRCh38 reference genome fasta files downloaded from NCBI

# Usage: extract_metadata_from_fasta.sh <INPUT_FASTA>

### Extract Metadata fasta header lines
## This is designed to work with NCBI downloads with the following metadata flags in the fasta header lines
#AC: sequence accession.version
#gi: sequence gi
#LN: sequence length
#rg: region
#    - chromosome to which unlocalized scaffolds are assigned,
#       e.g. chr1
#    - region on chromosome within which alt-scaffolds or patch
#      scaffolds are placed, e.g. chr6:28696604-33335493
#    - not present for chromosomes, other replicons, or unplaced
#      scaffolds
#    - coordinates are 1-based
#rl: role of the sequence in the assembly
#    - possible values are: Chromosome, Mitochondrion, unlocalized,
#      unplaced, alt-scaffold fix-patch, novel-patch, decoy
#M5: md5 checksum of the sequence as a single string of uppercase
#    letters without line breaks (as produced by Samtools or Picard)
#AS: assembly-name
#hm: hard-masked regions, either a single span, two spans separated by
#    a comma, or "multiple" if more than two spans were hard-masked
#    - coordinates are 1-based
#tp: topology
#    - circular for chrM and chrEBV
#    - not present for linear chromosomes and scaffolds
## Note: The EBV contig has an SP: tag that is not expect and does not exist in others, we convert this to the AS: tag
## Process:
### Find header lines as they start with ">"
### Remove the ">" character
### Replace the single SP: tag with AS: tag
### For lines missing the "rg:" tag replace the "rl:" that other wise follows it with "rg:.  rl:"
### For lines missing the "hm:" add the "hm:" to the end of the row with two space separations
### For lines missing the "tp:" tag replace the "hm:" that other wise follows it with "tp:.  hm:"
### sed artistry thanks to Christophe Legendre
zcat ${1} \
    | \
    grep "^>" \
    | \
    sed 's/>//' \
    | \
    sed 's/ SP:/ AS:/' \
    | \
    sed '/rg:/! s/rl:/rg:.  rl:/' \
    | \
    sed '/hm:/! s/$/  hm:./' \
    | \
    sed '/tp:/! s/ hm:/ tp:.  hm:/' \
    > temp_fasta_metadata.txt

## Create a header for the final file
echo -e Contig"\t"ASSEMBLY_ACC"\t"ASSEMBLY_GI"\t"SLEN"\t"REGION"\t"ROLE"\t"MD5_CHECKSUM"\t"ASSEMBLY_NAME"\t"TOPOLOGY"\t"MASKED_REGIONS > reference_fasta_metadata.txt

# Remove Tag keys from each data line and convert to tab separated
sed 's/AC://g' temp_fasta_metadata.txt \
    | \
    sed 's/gi://g' \
    | \
    sed 's/LN://g' \
    | \
    sed 's/rg://g' \
    | \
    sed 's/rl://g' \
    | \
    sed 's/M5://g' \
    | \
    sed 's/AS://g' \
    | \
    sed 's/tp://g' \
    | \
    sed 's/hm://g' \
    | \
    awk '{IFS="  " ; OFS="\t" ; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' >> reference_fasta_metadata.txt
