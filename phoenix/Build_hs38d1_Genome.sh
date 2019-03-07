#!/usr/bin/env bash

## The GRCh38 reference genome is represented in different public locations and not all have defined source information
## This file tracks the source and generation of the reference and annotation files used in the Jetstream phoenix workflow

## NOTE:
### - Because BWA is alt aware and STAR does not support alternative contigs we are creating two versions
### - Because BWAKIT and the Broad bundle lack sufficient source information we are creating references independently
### ---- However, the BWA version is identical to the one bwakit version which is identical to the broad bundle version, albeit with different names
### - Also, review of what is used by GDC showed they use an identical base genome (1-22, X, Y, M, supercontigs, decoys, EBV) but with out alternative contigs
### ---- To this they added a series of HPV and other viral genomes
### ---- They use the same reference for BWA and STAR alignments

# BROAD_BUNDLE: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/?pli=1
# BWAKIT: https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2
# GDC: https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files

# Confirm BWAKIT and Broad Bundle are identical in sequence
diff /home/tgenref/homo_sapiens/grch38_hg38/hs38dh/genome_references/hs38DH.fa /home/tgenref/homo_sapiens/grch38_hg38/broad_resource_bundle/Homo_sapiens_assembly38.fasta
# A massive number of events are seen BUT this seems to relate to the differences in the FASTA line lengths as the broad bundle lengths are longer 101 letters versus 71 (which we get with NCBI download as well)

tr -d '\n' < /home/tgenref/homo_sapiens/grch38_hg38/hs38dh/genome_references/hs38DH.fa > hs38DH.fa
tr -d '\n' < /home/tgenref/homo_sapiens/grch38_hg38/broad_resource_bundle/Homo_sapiens_assembly38.fasta > Homo_sapiens_assembly38.fasta
diff Homo_sapiens_assembly38.fasta hs38DH.fa | wc-l
# 0
# No differences detected SO THE BWAKIT AND BROAD BUNDLE FASTA FILES ARE IDENTICAL!!

# Check contigs
grep "^>" /home/tgenref/homo_sapiens/grch38_hg38/hs38dh/genome_references/hs38DH.fa > hs38DH_contigs.txt
grep "^>" /home/tgenref/homo_sapiens/grch38_hg38/broad_resource_bundle/Homo_sapiens_assembly38.fasta > Homo_sapiens_assembly38_contigs.txt

diff hs38DH_contigs.txt Homo_sapiens_assembly38_contigs.txt
# No contig differences

# Download GRC/NCBI README
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt

####################################
## BUILD BWA REFERENCE GENOME
####################################

# Download GRC/NCBI README
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt

# Download the primary assembly with decoy sequences:
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
# Decompress the reference
gunzip GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz

# Download bwakit, to get the HLA contigs used by bwakit and the production GRCh38 pipeline used at Broad
wget https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2
# Decompress the bwakit download
tar xvjf bwakit-0.7.15_x64-linux.tar.bz2


### Create a fasta file with just the HLA from BWAKIT

# Find the line with the first HLA region
awk '/HLA-A\*01:01:01:01/{print FNR}' bwa.kit/resource-GRCh38/hs38DH-extra.fa
# Identified line
# 75967
# Print all lines including this line and afterwards
awk 'NR>=75967' bwa.kit/resource-GRCh38/hs38DH-extra.fa > bwakit_HLA.fasta

### Concatenate the fasta files to make final reference genome fasta
cat GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna bwakit_HLA.fasta > GRCh38_hs38d1_Alts_HLA.fa

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