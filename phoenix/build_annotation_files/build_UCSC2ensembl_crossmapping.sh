#!/usr/bin/env bash

# Usage: ./build_UCSC2ensembl_crossmapping.sh  <resources.ini>

## This file is created for usage with bcftools annotate to update contig names in the snpSniffer
## output vcf from "chr1" to "1" as the tools database does not support chr in the contig string

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
## Navigate Directory Structure
###################################

# Check top level directory if not available
if [ -e ${TOPLEVEL_DIR} ]
then
    echo "Top level directory: ${TOPLEVEL_DIR} exists, moving into it"
    cd ${TOPLEVEL_DIR}
else
    echo "Top level directory NOT found, IT IS REQUIRED, EXITING"
    exit 1
fi

# Create directory for tool resources
if [ -e tool_resources ]
then
    echo "tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "tool_resources directory NOT found, creating and moving into it now"
    mkdir tool_resources
    cd tool_resources
fi

# Create directory for tool
if [ -e "bcftools" ]
then
    echo "The bcftools directory exists, exiting to prevent overwriting existing index"
    exit 2
else
    echo "The bcftools directory was NOT found, creating and moving into it now"
    mkdir bcftools
    cd bcftools
fi

# Initialize a bcftools index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

echo -e chr1"\t"1 > GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr2"\t"2 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr3"\t"3 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr4"\t"4 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr5"\t"5 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr6"\t"6 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr7"\t"7 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr8"\t"8 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr9"\t"9 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr10"\t"10 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr11"\t"11 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr12"\t"12 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr13"\t"13 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr14"\t"14 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr15"\t"15 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr16"\t"16 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr17"\t"17 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr18"\t"18 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr19"\t"19 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr20"\t"20 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr21"\t"21 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chr22"\t"22 >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chrX"\t"X >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
echo -e chrY"\t"Y >> GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
