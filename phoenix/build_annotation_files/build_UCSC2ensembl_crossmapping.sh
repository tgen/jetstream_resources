#!/usr/bin/env bash

# Usage: ./build_UCSC2ensembl_crossmapping.sh  <resources.ini>

## This file is created for usage with bcftools annotate to update contig names in the snpSniffer
## output vcf from "chr1" to "1" as the tools database does not support chr in the contig string

# Check if resources.ini was provided on the command line
if [ -n "$1" ]; then
  echo "Required ini file detected"
else
  echo "Input INI file not provided, exiting due to missing requirement"
  exit 1
fi

# Read required variables from configuration file
. ${1}

# navigate to gnomad directory
cd ${TOPLEVEL_DIR}/tool_resources

# make a folder for bcftools
mkdir -p bcftools

# enter folder
cd bcftools

# Write this full document as a README
cat $0 > README

# Added user information and timestamp to README
USER=`whoami`
DATE=`date`
echo "Downloaded and Processed by:  ${USER}" >> README
echo ${DATE} >> README

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
