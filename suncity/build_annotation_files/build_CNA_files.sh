#!/usr/bin/env bash

# Usage: ./build_CNA_files.sh  <resources.ini>

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
  module load GATK/4.1.7.0-GCCcore-8.3.0-Java-1.8
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
if [ -e ucsc ]
then
    echo "ucsc folder exists, moving into it"
    cd ucsc
else
    echo "ncbi folder NOT fount, creating and moving into it now"
    mkdir -p ucsc
    cd ucsc
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
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
fc -ln -1 >> README
echo >> README

# Make a bed file of all of the Modeled_regions_for_GRCh37
zcat gap.txt.gz | gawk -F'\t' 'NR > 1 { OFS = "\t" ; print $2,$3,$4,$8}' | sed -e "s/^M//" > Modeled_regions_for_GRCh37.bed
fc -ln -1 >> README
echo >> README

# Make a centromere only bed file
grep cent Modeled_regions_for_GRCh37.bed > Modeled_regions_for_GRCh37_centromere.bed
fc -ln -1 >> README
echo >> README

# Make heterochromitin only bed file
grep hete Modeled_regions_for_GRCh37.bed > Modeled_regions_for_GRCh37_heterochromitin.bed
fc -ln -1 >> README
echo >> README

# Get the regions of the black list that intersect with the centromeres
bedtools intersect -wb -a Modeled_regions_for_GRCh37_centromere.bed \
  -b /home/tgenref/homo_sapiens/grch37_hg19/hs37d5_suncity/public_databases/encode/Blacklist-2.0/lists/hg19-blacklist.v2.bed |\
  cut -f5- |\
  sort -k1,1V -k2,2n -k3,3n > Modeled_regions_for_GRCh37_centromere_blacklist_intersect.bed
fc -ln -1 >> README
echo >> README

# Merge the centromere regions and the overlaping encode blacklist regions together
cat Modeled_regions_for_GRCh37_centromere.bed Modeled_regions_for_GRCh37_centromere_blacklist_intersect.bed | \
  cut -f1-4 | \
  sort -k1,1V -k2,2n -k3,3n | \
  bedtools merge > Modeled_regions_for_GRCh37_centromere_blacklist_merged.bed
fc -ln -1 >> README
echo >> README

cat /home/tgenref/homo_sapiens/grch37_hg19/hs37d5_suncity/public_databases/encode/Blacklist-2.0/lists/hg19-blacklist.v2.bed Modeled_regions_for_GRCh37_centromere_blacklist_merged.bed | cut -f1-3 | sort -k1,1V -k2,2n > Encode_deny_list_with_ucsc_centromere.bed
fc -ln -1 >> README
echo >> README

bedtools merge -i Encode_deny_list_with_ucsc_centromere.bed | cut -c 4- > Encode_deny_list_with_ucsc_centromere.merged.bed
fc -ln -1 >> README
echo >> README

#Prep centromere for suncity usage
cat Modeled_regions_for_GRCh37_centromere_blacklist_merged.bed | cut -c 4- > Modeled_regions_for_GRCh37_centromere_blacklist_merged_final.bed

####################################
## Download and Manipulate the Bismap mappability file
###################################

cd ${PARENT_DIR}/public_databases

# Make bismap folder if not available
if [ -e bismap ]
then
    echo "bismap folder exists, moving into it"
    cd bismap
else
    echo "bismap folder NOT fount, creating and moving into it now"
    mkdir -p bismap
    cd bismap
fi

# Initialize a bismap README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Download the Umap k100 hg19 Single-read mappability file
wget https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bed.gz
fc -ln -1 >> README
echo >> README

# Unzip the bed file and remove the header
zcat k100.umap.bed.gz | awk 'NR >1' > k100.umap.no_header.bed
fc -ln -1 >> README
echo >> README

# Index the mappability file
gatk IndexFeatureFile \
     --input k100.umap.no_header.bed
fc -ln -1 >> README
echo >> README

echo
echo "Process Complete"
echo
