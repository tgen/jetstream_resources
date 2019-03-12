#!/usr/bin/env bash

## The GRCh38 reference genome is represented in different public locations and not all have defined source information
## This file tracks the source and generation of the reference and annotation files used in the Jetstream phoenix workflow

## NOTE:
### - Because BWA is alt aware and STAR does not support alternative contigs we are creating two versions
### - Because BWAKIT and the Broad bundle lack sufficient source information we are creating references independently
### ---- However, the BWAKIT version is identical to the broad bundle version, albeit with different base file names
### - Also, review of what is used by GDC showed they use an identical base genome (1-22, X, Y, M, supercontigs, decoys, EBV) but with out alternative contigs
### ---- To this they added a series of HPV and other viral genomes
### ---- They use the same reference for BWA and STAR alignments

# BROAD_BUNDLE: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/?pli=1
# BWAKIT: https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2
# GDC: https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files

# Confirm BWAKIT and Broad Bundle are identical in sequence
#diff /home/tgenref/homo_sapiens/grch38_hg38/hs38dh/genome_references/hs38DH.fa /home/tgenref/homo_sapiens/grch38_hg38/broad_resource_bundle/Homo_sapiens_assembly38.fasta
# A massive number of events are seen BUT this seems to relate to the differences in the FASTA line lengths as the broad bundle lengths are longer 101 letters versus 71 (which we get with NCBI download as well)

#tr -d '\n' < /home/tgenref/homo_sapiens/grch38_hg38/hs38dh/genome_references/hs38DH.fa > hs38DH.fa
#tr -d '\n' < /home/tgenref/homo_sapiens/grch38_hg38/broad_resource_bundle/Homo_sapiens_assembly38.fasta > Homo_sapiens_assembly38.fasta
#diff Homo_sapiens_assembly38.fasta hs38DH.fa | wc-l
# 0
# No differences detected SO THE BWAKIT AND BROAD BUNDLE FASTA FILES ARE IDENTICAL!!

# Check contigs
#grep "^>" /home/tgenref/homo_sapiens/grch38_hg38/hs38dh/genome_references/hs38DH.fa > hs38DH_contigs.txt
#grep "^>" /home/tgenref/homo_sapiens/grch38_hg38/broad_resource_bundle/Homo_sapiens_assembly38.fasta > Homo_sapiens_assembly38_contigs.txt

#diff hs38DH_contigs.txt Homo_sapiens_assembly38_contigs.txt
# No contig differences

####################################
## Configure and make Directory Structure
####################################

PARENT_DIR="/scratch/jkeats"
TOPLEVEL_DIR="hs38d1"
CREATOR="Jonathan Keats"

# Change to parent directory
cd ${PARENT_DIR}

# Make top level directory if not available
if [ -e ${TOPLEVEL_DIR} ]
then
    echo "Top level directory: ${TOPLEVEL_DIR} exists, moving into it"
    cd ${TOPLEVEL_DIR}
else
    echo "Top level directory NOT fount, creating and moving into it now"
    mkdir ${TOPLEVEL_DIR}
    cd ${TOPLEVEL_DIR}
fi

# Initialize a top level README
touch README
echo "Reference Genome and related files required for JetStream Phoenix Workflow" >> README


####################################
## Create reference genomes from known sources
####################################

# Make reference_genome directory if not available
if [ -e genome_reference ]
then
    echo "Genome Reference directory exists, moving into it"
    cd genome_reference
else
    echo "Genome Reference directory NOT fount, creating and moving into it now"
    mkdir genome_reference
    cd genome_reference
fi

# Initialize a reference_genome README
touch README_TGen
echo >> README_TGen
echo "Created and downloaded by ${CREATOR}" >> README_TGen
date >> README_TGen
echo >> README_TGen

# Download GRC/NCBI README
echo "Download README from NCBI" >> README_TGen
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt
fc -ln -1 >> README_TGen

exit 1
####################################
## BUILD BWA REFERENCE GENOME
####################################

# Download the primary assembly with decoy sequences:
echo "Download primary assembly with decoy sequences:" >> README_TGen
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
fc -ln -1 >> README_TGen
# Decompress the reference
gunzip GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
fc -ln -1 >> README_TGen
echo >> README_TGen

# Download bwakit, to get the HLA contigs used by bwakit and the production GRCh38 pipeline used at Broad
echo "Download bwakit, to get the HLA contigs used by bwakit and the production GRCh38 pipeline used at Broad:" >> README_TGen
wget https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2
fc -ln -1 >> README_TGen
# Decompress the bwakit download
tar xvjf bwakit-0.7.15_x64-linux.tar.bz2
fc -ln -1 >> README_TGen
echo >> README_TGen

### Create a fasta file with just the HLA from BWAKIT

# Find the line with the first HLA region
echo "The bwakit hs38DH-extra.fa file contains decoy and HLA contigs, need to extract just the HLA ones" >> README_TGen
awk '/HLA-A\*01:01:01:01/{print FNR}' bwa.kit/resource-GRCh38/hs38DH-extra.fa >> README_TGen
# Identified line
# 75967
# Print all lines including this line and afterwards
awk 'NR>=75967' bwa.kit/resource-GRCh38/hs38DH-extra.fa > bwakit_HLA.fasta
fc -ln -1 >> README_TGen
echo >> README_TGen

### Concatenate the fasta files to make final reference genome fasta
echo "Concatenate the fasta files to make final reference genome fasta to be used by BWA" >> README_TGen
cat GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna bwakit_HLA.fasta > GRCh38_hs38d1_Alts_HLA.fa
fc -ln -1 >> README_TGen

# Add symbolic link to indicate which FASTA is used by BWA
ln -s GRCh38_hs38d1_Alts_HLA.fa BWA_FASTA

# Clean up the directory to store the downloads
if [ -e downloads ]
then
    echo "Downloads directory exists, doing nothing"
else
    echo "Downloads directory DOES NOT exist, creating it now"
    mkdir downloads
fi

mv bwa.kit downloads
mv GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna downloads
mv bwakit_HLA.fasta downloads

exit 1
# Get a contig list for comparisons
grep "^>" GRCh38_hs38d1_Alts_HLA.fa > GRCh38_hs38d1_Alts_HLA_contigs.txt
diff GRCh38_hs38d1_Alts_HLA_contigs.txt hs38DH_contigs.txt
diff GRCh38_hs38d1_Alts_HLA_contigs.txt Homo_sapiens_assembly38_contigs.txt
# Lots of differences in contigs, seem to be related to the rl: annotation element
diff GRCh38_hs38d1_Alts_HLA_contigs.txt hs38DH_contigs.txt | grep "725020268"
#< >chrUn_JTFH01001998v1_decoy  AC:JTFH01001998.1  gi:725020268  LN:2001  rl:decoy  M5:35916d4135a2a9db7bc0749955d7161a  AS:hs38d1
#> >chrUn_JTFH01001998v1_decoy  AC:JTFH01001998.1  gi:725020268  LN:2001  rl:unplaced  M5:35916d4135a2a9db7bc0749955d7161a  AS:hs38d1
diff GRCh38_hs38d1_Alts_HLA_contigs.txt Homo_sapiens_assembly38_contigs.txt | grep "725020268"
#< >chrUn_JTFH01001998v1_decoy  AC:JTFH01001998.1  gi:725020268  LN:2001  rl:decoy  M5:35916d4135a2a9db7bc0749955d7161a  AS:hs38d1
#> >chrUn_JTFH01001998v1_decoy  AC:JTFH01001998.1  gi:725020268  LN:2001  rl:unplaced  M5:35916d4135a2a9db7bc0749955d7161a  AS:hs38d1

### Check versus HS38DH from BWAKIT
diff GRCh38_hs38d1_Alts_HLA.fa /home/tgenref/homo_sapiens/grch38_hg38/hs38dh/genome_references/hs38DH.fa
# Differences with repeat masked lowecase and uppercase AGCT on the decoy contigs
diff GRCh38_hs38d1_Alts_HLA.fa /home/tgenref/homo_sapiens/grch38_hg38/broad_resource_bundle/Homo_sapiens_assembly38.fasta

# Add symbolic link to indicate which FASTA is used by BWA

####################################
## BUILD STAR REFERENCE GENOME
####################################

# Since start does not support alternative contigs, download version wtih out
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz

# Rename fastq file
mv GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna GRCh38_hs38d1.fa

# Add symbolic link to indicate which FASTA is used by STAR