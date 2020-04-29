#!/usr/bin/env bash

# Usage: ./build_centromere_and_heterochromatin_bed_files.sh  <resources.ini>

## This file is created for usage in jetstream_build_package

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -ue

# Check resources.ini was provided on the command line
if [ -n "$1" ]
then
  echo "Required ini file detected"
else
  echo "Input INI file not provided, exiting due to missing requirement"
  exit 1
fi

# Read required variables from configuration file
. ${1}

####################################
## Load Required Tools
###################################
if [ ${ENVIRONMENT} == "TGen" ]
then
  module load BEDTools/2.29.0-GCC-8.2.0-2.31.1
elif [ ${ENVIRONMENT} == "LOCAL" ]
then
  echo
  echo "Assuming required tools are available in $PATH"
  echo
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Enviroment Variable"
  echo "Only TGen or LOCAL are supported"
  exit 1
fi

####################################
## Create Expected Folder Structure
###################################

# Make top level directory if not available
if [ -e ${PARENT_DIR} ]
then
    echo "Parent directory: ${PARENT_DIR} exists, moving into it"
    cd ${PARENT_DIR}
else
    echo "Parent directory NOT fount, creating and moving into it now"
    mkdir -p ${PARENT_DIR}
    cd ${PARENT_DIR}
fi

# Make public_databases folder if not available
if [ -e public_databases ]
then
    echo "Public Databases folder exists, moving into it"
    cd public_databases
else
    echo "Public Databases folder NOT fount, creating and moving into it now"
    mkdir -p public_databases
    cd public_databases
fi

# Make ncbi folder if not available
if [ -e ncbi ]
then
    echo "ncbi folder exists, moving into it"
    cd ncbi
else
    echo "ncbi folder NOT fount, creating and moving into it now"
    mkdir -p ncbi
    cd ncbi
fi

####################################
## Download and Manipulate the centromere file
###################################

# Initialize a bcftools index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Download the Modeled centromeres and heterochromatin regions
wget https://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/38/Modeled_regions_for_GRCh38.tsv
fc -ln -1 >> README
echo >> README

# Make a bed file of all of the Modeled_regions_for_GRCh38
gawk -F'\t' 'NR > 1 { OFS = "\t" ; print "chr"$2,$3,$4,$5,$1}' Modeled_regions_for_GRCh38.tsv | sed -e "s/^M//" > Modeled_regions_for_GRCh38.bed
fc -ln -1 >> README
echo >> README

# Make a centromere only bed file
grep CEN Modeled_regions_for_GRCh38.bed > Modeled_regions_for_GRCh38_centromere.bed
fc -ln -1 >> README
echo >> README

# Make heterochromitin only bed file
grep HET Modeled_regions_for_GRCh38.bed > Modeled_regions_for_GRCh38_heterochromitin.bed
fc -ln -1 >> README
echo >> README

# Get the regions of the black list that intersect with the centromeres
bedtools intersect -wb -a Modeled_regions_for_GRCh38_centromere.bed \
	-b ${PARENT_DIR}/public_databases/encode/Blacklist-2.0/lists/hg38-blacklist.v2.bed.gz | \
	cut -f6- | \
	sort -k1,1V -k2,2n -k3,3n > Modeled_regions_for_GRCh38_centromere_blacklist_intersect.bed
fc -ln -1 >> README
echo >> README

# Merge the centromere regions and the overlaping encode blacklist regions together
cat Modeled_regions_for_GRCh38_centromere.bed Modeled_regions_for_GRCh38_centromere_blacklist_intersect.bed | \
	cut -f1-4 | \
	sort -k1,1V -k2,2n -k3,3n | \
	bedtools merge > Modeled_regions_for_GRCh38_centromere_blacklist_merged.bed
fc -ln -1 >> README
echo >> README

echo
echo "Process Complete"
echo
