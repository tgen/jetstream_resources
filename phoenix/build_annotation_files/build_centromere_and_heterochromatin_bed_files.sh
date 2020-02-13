#!/usr/bin/env bash

# Usage: ./build_centromere_and_heterochromatin_bed_files.sh  <resources.ini>

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

# Write this full document as a README
cat $0 > README

