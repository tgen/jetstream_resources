#!/usr/bin/env bash

# Usage: create_genome_reference.sh <Config.ini>

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
if [ $ENVIRONMENT == "TGen"]
then
  module load SAMtools/1.10-GCC-8.2.0-2.31.1
else
  echo
  echo "Assuming required tools are available in $PATH"
  echo
fi

####################################
## Configure and make Directory Structure
###################################

# Make top level directory if not available
if [ -e ${TOPLEVEL_DIR} ]
then
    echo "Top level directory: ${TOPLEVEL_DIR} exists, moving into it"
    cd ${TOPLEVEL_DIR}
else
    echo "Top level directory NOT found, creating and moving into it now"
    mkdir -p ${TOPLEVEL_DIR}
    cd ${TOPLEVEL_DIR}
fi

# Initialize a top level README
touch README
echo >> README
echo "Reference Genome and related files required for JetStream ${WORKFLOW_NAME} Workflow" >> README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
echo >> README
echo "
## The GRCh38 reference genome is represented in different public locations and not all have defined source information
## This file tracks the source and generation of the reference and annotation files used in the Jetstream phoenix workflow

## NOTE:
### - Because BWA is alt aware and STAR does not support alternative contigs we are creating two versions
### - Because BWAKIT and the Broad bundle lack sufficient source information we are creating references independently
### ---- However, the BWAKIT version is identical to the broad bundle version, albeit with different base file names
### - Also, review of what is used by GDC showed they use an identical base genome (1-22, X, Y, M, supercontigs, decoys, EBV) but with out alternative contigs
### ---- To this they added a series of HPV and other viral genomes
### ---- They use the same reference for BWA and STAR alignments
"  >> README

####################################
## Create reference genomes from known sources
####################################

# Make reference_genome directory if not available
if [ -e genome_reference ]
then
    echo "Genome Reference directory exists, moving into it"
    cd genome_reference
else
    echo "Genome Reference directory NOT found, creating and moving into it now"
    mkdir genome_reference
    cd genome_reference
fi

# Initialize a reference_genome README
touch README
echo >> README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Download GRC/NCBI README
echo "Download README from NCBI" >> README
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt
fc -ln -1 >> README
echo >> README

####################################
## Download BWA REFERENCE GENOME
####################################

echo "## Download reference genome fasta for BWA " >> README
echo >> README

# Create a downloads directory to store downloads after processing
if [ -e downloads ]
then
    echo "Downloads directory exists, doing nothing"
else
    echo "Downloads directory DOES NOT exist, creating it now"
    mkdir downloads
fi

# Download the primary assembly with decoy sequences:
echo "Download primary assembly with decoy sequences:" >> README
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
fc -ln -1 >> README
# Archive a copy of the download file
echo "Archive a copy of the download file" >> README
cp GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz downloads/
fc -ln -1 >> README
# Decompress the reference
gunzip GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
fc -ln -1 >> README
echo >> README

# Download bwakit, to get the HLA contigs used by bwakit and the production GRCh38 pipeline used at Broad
echo "Download bwakit, to get the HLA contigs used by bwakit and the production GRCh38 pipeline used at Broad:" >> README
wget https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2
fc -ln -1 >> README
# Decompress the bwakit download
tar xvjf bwakit-0.7.15_x64-linux.tar.bz2
fc -ln -1 >> README
echo >> README

### Create a fasta file with just the HLA from BWAKIT
# Find the line with the first HLA region
echo "The bwakit hs38DH-extra.fa file contains decoy and HLA contigs, need to extract just the HLA ones" >> README
echo "Hand curated the name of the first HLA contig, found line number with awk as follows:" >> README
awk '/HLA-A\*01:01:01:01/{print FNR}' bwa.kit/resource-GRCh38/hs38DH-extra.fa >> README
fc -ln -1 >> README
# Identified line
# 75967
# Print all lines including this line and afterwards
echo "Print all lines after and including first HLA line" >> README
awk 'NR>=75967' bwa.kit/resource-GRCh38/hs38DH-extra.fa > bwakit_HLA.fasta
fc -ln -1 >> README
echo >> README

### Concatenate the fasta files to make final reference genome fasta
echo "Concatenate the fasta files to make final reference genome fasta to be used by BWA" >> README
cat GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna bwakit_HLA.fasta > GRCh38tgen_decoy_alts_hla.fa
fc -ln -1 >> README
echo >> README

# Add symbolic link to indicate which FASTA is used by BWA
ln -s GRCh38tgen_decoy_alts_hla.fa BWA_FASTA

# Create faidx and dict files
echo "Create BWA faidx index using samtools" >> README

samtools faidx GRCh38tgen_decoy_alts_hla.fa
fc -ln -1 >> README
echo >> README

echo "Create BWA dictionary file using samtools" >> README
samtools dict --assembly GRCh38 \
    --species "Homo sapiens" \
    --uri "downloads/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz" \
    --output GRCh38tgen_decoy_alts_hla.dict \
    GRCh38tgen_decoy_alts_hla.fa
fc -ln -1 >> README
echo >> README

# Clean up the directory to store the downloads
mv bwakit-0.7.15_x64-linux.tar.bz2 downloads
mv bwa.kit downloads
mv README_analysis_sets.txt downloads
mv bwakit_HLA.fasta downloads
rm GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna

# Add flag to top level to indicate process is complete
touch ${TOPLEVEL_DIR}/GENOME_FASTA_GENERATION_COMPLETE

####################################
## Download STAR REFERENCE GENOME
####################################

echo >> README
echo "## Download reference genome fasta for STAR " >> README
echo >> README

# Since STAR does not support alternative contigs, download version wtih out
echo "Download fasta file for STAR, version without alternate contigs and no HLA alleles added" >> README
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
fc -ln -1 >> README
echo >> README

# Archive a copy in the downloads section then decompress
echo "Archive compressed STAR fasta file download then decompress" >> README
cp GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz downloads/
fc -ln -1 >> README
gunzip GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
fc -ln -1 >> README
echo >> README

# Rename fastq file
echo "Rename the downloaded file to decrease filename length and make consistent with bwa fasta" >> README
mv GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna GRCh38tgen_decoy.fa
fc -ln -1 >> README
echo >> README

# Add symbolic link to indicate which FASTA is used by STAR
ln -s GRCh38tgen_decoy.fa STAR_FASTA

# Create faidx and dict files
echo "Create STAR faidx index using samtools" >> README
samtools faidx GRCh38tgen_decoy.fa
fc -ln -1 >> README
echo >> README

echo "Create STAR dictionary file using samtools" >> README
samtools dict --assembly GRCh38 \
    --species "Homo sapiens" \
    --uri "downloads/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz" \
    --output GRCh38tgen_decoy.dict \
    GRCh38tgen_decoy.fa
fc -ln -1 >> README
echo >> README

# Add flag to top level to indicate process is complete
touch ${TOPLEVEL_DIR}/RNA_FASTA_GENERATION_COMPLETE