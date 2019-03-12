#!/bin/bash -i

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history

## The GRCh38 reference genome is represented in different public locations and not all have defined source information
## This file tracks the source and generation of the reference and annotation files used in the Jetstream phoenix workflow

## NOTE:
### - Because BWA is alt aware and STAR does not support alternative contigs we are creating two versions
### - Because BWAKIT and the Broad bundle lack sufficient source information we are creating references independently
### ---- However, the BWAKIT version is identical to the broad bundle version, albeit with different base file names
### - Also, review of what is used by GDC showed they use an identical base genome (1-22, X, Y, M, supercontigs, decoys, EBV) but with out alternative contigs
### ---- To this they added a series of HPV and other viral genomes
### ---- They use the same reference for BWA and STAR alignments

####################################
## Configure and make Directory Structure
####################################

PATH_TO_REPO="/home/jkeats/git_repositories/jetstream_resources/"
PARENT_DIR="/scratch/jkeats"
TOPLEVEL_DIR="hg38tgen"
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
echo >> README
echo "Reference Genome and related files required for JetStream Phoenix Workflow" >> README
echo >> README
echo "For details on file creation see the associated github repository:"
echo "https://github.com/tgen/jetstream_resources"

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
echo >> README
echo "For details on file creation see the associated github repository:" >> README_TGen
echo "https://github.com/tgen/jetstream_resources" >> README_TGen
echo "Created and downloaded by ${CREATOR}" >> README_TGen
date >> README_TGen
echo >> README_TGen

# Download GRC/NCBI README
echo "Download README from NCBI" >> README_TGen
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt
fc -ln -1 >> README_TGen

####################################
## Download BWA REFERENCE GENOME
####################################

echo "####################################" >> README_TGen
echo "## Download reference genome fasta for BWA " >> README_TGen
echo "####################################" >> README_TGen
echo >> README_TGen

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
echo "Hand curated the name of the first HLA contig, found line number with awk as follows:" >> README_TGen
awk '/HLA-A\*01:01:01:01/{print FNR}' bwa.kit/resource-GRCh38/hs38DH-extra.fa >> README_TGen
fc -ln -1 >> README_TGen
# Identified line
# 75967
# Print all lines including this line and afterwards
echo "Print all lines after and including first HLA line" >> README_TGen
awk 'NR>=75967' bwa.kit/resource-GRCh38/hs38DH-extra.fa > bwakit_HLA.fasta
fc -ln -1 >> README_TGen
echo >> README_TGen

### Concatenate the fasta files to make final reference genome fasta
echo "Concatenate the fasta files to make final reference genome fasta to be used by BWA" >> README_TGen
cat GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna bwakit_HLA.fasta > GRCh38tgen_decoy_alts_hla.fa
fc -ln -1 >> README_TGen
echo >> README_TGen

# Add symbolic link to indicate which FASTA is used by BWA
ln -s GRCh38tgen_decoy_alts_hla.fa BWA_FASTA

# Clean up the directory to store the downloads
if [ -e downloads ]
then
    echo "Downloads directory exists, doing nothing"
else
    echo "Downloads directory DOES NOT exist, creating it now"
    mkdir downloads
fi

mv bwakit-0.7.15_x64-linux.tar.bz2 downloads
mv bwa.kit downloads
mv GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna downloads
mv bwakit_HLA.fasta downloads

# Create faidx and dict files
echo "Create faidx index using samtools" >> README_TGen
module load samtools/1.9
fc -ln -1 >> README_TGen
samtools faidx GRCh38tgen_decoy_alts_hla.fa
fc -ln -1 >> README_TGen
echo >> README_TGen

echo "Create dictionary file using samtools" >> README_TGen
samtools dict --assembly GRCh38 \
    --species "Homo sapiens" \
    --uri "downloads/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna" \
    --output GRCh38_hs38d1_Alts_HLA.dict \
    GRCh38tgen_decoy_alts_hla.fa
fc -ln -1 >> README_TGen
echo >> README_TGen

####################################
## INDEX BWA REFERENCE GENOME
####################################
cd ..

if [ -e tool_specific_resources ]
then
    echo "tool_specific_resources directory exists, moving into it"
    cd tool_specific_resources
else
    echo "tool_specific_resources directory NOT fount, creating and moving into it now"
    mkdir tool_specific_resources
    cd tool_specific_resources
fi

if [ -e bwa ]
then
    echo "The BWA directory exists, moving into it"
    cd bwa
else
    echo "The BWA directory NOT fount, creating and moving into it now"
    mkdir bwa
    cd bwa
fi

# Initialize a bwa index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "BWA Index creation details:" >> README
echo >> README
echo "There are no details on how the alt-aware required .alt file is created" >> README
echo "Because of this we are copying and renaming the file from BWAKIT" >> README
echo >> README

# Copy in the .alt file from bwa.kit
echo "Copy in the .alt file from bwa.kit" >> README
cp ../../genome_reference/downloads/bwa.kit/resource-GRCh38/hs38DH.fa.alt GRCh38tgen_decoy_alts_hla.fa.alt
fc -ln -1 >> README
echo >> README

# Create a symbolic link to the reference genome
ln -s ../../genome_reference/GRCh38tgen_decoy_alts_hla.fa GRCh38tgen_decoy_alts_hla.fa

# Create bwa index files using bwa utility script
echo "Create bwa index as follows:" >> README
sbatch --export ALL,FASTA='GRCh38tgen_decoy_alts_hla.fa' ${PATH_TO_REPO}/utility_scripts/bwa_index.slurm
fc -ln -1 >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/bwa_index.slurm >> README
echo >> README
echo >> README

####################################
## Download STAR REFERENCE GENOME
####################################

# move back to the genome_reference directory
cd ../../genome_reference

echo "####################################" >> README_TGen
echo "## Download reference genome fasta for STAR " >> README_TGen
echo "####################################" >> README_TGen
echo >> README_TGen

# Since STAR does not support alternative contigs, download version wtih out
echo "Download fasta file for STAR, version without alternate contigs and no HLA alleles added" >> README_TGen
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
fc -ln -1 >> README_TGen
gunzip GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
fc -ln -1 >> README_TGen
echo >> README_TGen

# Rename fastq file
echo "Rename the downloaded file to decrease filename length and make consistent with bwa fasta" >> README_TGen
mv GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna GRCh38tgen_decoy.fa
fc -ln -1 >> README_TGen
echo >> README_TGen

# Add symbolic link to indicate which FASTA is used by STAR
ln -s GRCh38tgen_decoy.fa STAR_FASTA

####################################
## Download Gene Model File
####################################

# move back to top level directory
cd ${PARENT_DIR}/${TOPLEVEL_DIR}

# Make gene_model directory if not available
if [ -e gene_model ]
then
    echo "Gene Model directory exists, moving into it"
    cd gene_model
else
    echo "Gene Model directory NOT fount, creating and moving into it now"
    mkdir gene_model
    cd gene_model
fi

# Make ensembl_v95 directory if not available
if [ -e ensembl_v95 ]
then
    echo "Gene Model directory exists, moving into it"
    cd ensembl_v95
else
    echo "Gene Model directory NOT fount, creating and moving into it now"
    mkdir ensembl_v95
    cd ensembl_v95
fi

# Initialize a gene_model specific README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "Ensembl version 95 download and manipulations to align with reference genomes" >> README
echo >> README

# Download the gene model gtf from ensembl
echo "Download the gene model gtf from ensembl" >> README
wget ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz
fc -ln -1 >> README
gunzip Homo_sapiens.GRCh38.95.gtf.gz
fc -ln -1 >> README

####################################
## Generate Gene Model Specific - tool_specific_resources
####################################

# Make gene_model specific tool_specific_resources directory if not available
if [ -e tool_specific_resources ]
then
    echo "Gene Model tool_specific_resources directory exists, moving into it"
    cd tool_specific_resources
else
    echo "Gene Model tool_specific_resources directory NOT fount, creating and moving into it now"
    mkdir tool_specific_resources
    cd tool_specific_resources
fi


####################################
## Generate Salmon Index
####################################


exit 0

####################################
## Genome Build Validations
####################################

### Compare BWAKIT and Broad Bundle versions

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


### Test TGen created fasta file versus other source fasta files

# Get a contig list for comparisons
grep "^>" GRCh38tgen_decoy_alts_hla.fa > GRCh38_hs38d1_Alts_HLA_contigs.txt
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
diff GRCh38tgen_decoy_alts_hla.fa /home/tgenref/homo_sapiens/grch38_hg38/hs38dh/genome_references/hs38DH.fa
# Differences with repeat masked lowecase and uppercase AGCT on the decoy contigs
diff GRCh38tgen_decoy_alts_hla.fa /home/tgenref/homo_sapiens/grch38_hg38/broad_resource_bundle/Homo_sapiens_assembly38.fasta