#!/usr/bin/env bash

# Usage: ./build_centromere_and_heterochromatin_bed_files.sh  <resources.ini>

module load BEDTools/2.29.0-GCC-8.2.0-2.31.1

## This file is created for usage in jetstream_build_package

# Check if resources.ini was provided on the command line
if [ -n "$1" ]; then
  echo "Required ini file detected"
else
  echo "Input INI file not provided, exiting due to missing requirement"
  exit 1
fi

# Read required variables from configuration file
. ${1}

# navigate to parent directory
cd ${PARENT_DIR}

# make a folder for bcftools
mkdir -p public_databases/ncbi

# enter created directory
cd public_databases/ncbi

# Download the Modeled centromeres and heterochromatin regions
wget https://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/38/Modeled_regions_for_GRCh38.tsv

# Make a bed file of all of the Modeled_regions_for_GRCh38
gawk -F'\t' 'NR > 1 { OFS = "\t" ; print "chr"$2,$3,$4,$5,$1}' Modeled_regions_for_GRCh38.tsv | sed -e "s///" > Modeled_regions_for_GRCh38.bed

# Make a centromere only bed file
grep CEN Modeled_regions_for_GRCh38.bed > Modeled_regions_for_GRCh38_centromere.bed

# Make heterochromitin only bed file
grep HET Modeled_regions_for_GRCh38.bed > Modeled_regions_for_GRCh38_heterochromitin.bed

# Get the regions of the black list that intersect with the centromeres
bedtools intersect -wb -a Modeled_regions_for_GRCh38_centromere.bed \
	-b /home/tgenref/homo_sapiens/grch38_hg38/public_databases/encode/Blacklist-2.0/lists/hg38-blacklist.v2.bed.gz | \
	cut -f6- | \
	sort -k1,1V -k2,2n -k3,3n > Modeled_regions_for_GRCh38_centromere_blacklist_intersect.bed

# Merge the centromere regions and the overlaping encode blacklist regions together
cat Modeled_regions_for_GRCh38_centromere.bed \
	Modeled_regions_for_GRCh38_centromere_blacklist_intersect.bed | \
	cut -f1-4 | \
	sort -k1,1V -k2,2n -k3,3n | \
	bedtools merge > Modeled_regions_for_GRCh38_centromere_blacklist_merged.bed

# Write this full document as a README
cat $0 > README

