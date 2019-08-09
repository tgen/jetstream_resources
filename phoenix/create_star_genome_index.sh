#!/usr/bin/env bash

# Usage: create_star_genome_index.sh <Config.ini>

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -ue

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

# Check that the reference genome for RNA was created successfully
if [ -e RNA_FASTA_GENERATION_COMPLETE ]
then
    echo "RNA fasta exists, moving forward"
else
    echo "RNA fasta generation complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
fi

# Check gene_model directory if not available
if [ -e gene_model ]
then
    echo "Gene Model directory exists, moving into it"
    cd gene_model
else
    echo "Gene Model directory NOT found, IT IS REQUIRED, exiting as it should be created by a previous step"
    exit 1
fi

# Check specific gene model directory
if [ -e ${GENE_MODEL_NAME} ]
then
    echo "Specific Gene Model directory exists, moving into it"
    cd ${GENE_MODEL_NAME}
else
    echo "Specific Gene Model directory NOT found, IT IS REQUIRED, exiting as it should be created by a previous step"
    exit 1
fi

# Check that the required GTF was created successfully
if [ -e GENE_MODEL_GTF_GENERATION_COMPLETE ]
then
    echo "Required gene model GTF exists, moving forward"
else
    echo "Required gene model GTF DOES NOT exist, exiting"
    exit 2

# Make gene_model specific tool_resources directory if not available
if [ -e tool_resources ]
then
    echo "Gene Model tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "Gene Model tool_resources directory NOT found, creating and entering it now"
    mkdir tool_resources
    cd tool_resources
fi

####################################
## Generate STAR Index
####################################

# Make star index directory if not available
if [ -e "star_${STAR_VERSION}" ]
then
    echo "STAR directory exists, moving into it"
    cd star_${STAR_VERSION}
else
    echo "STAR directory NOT fount, creating and moving into it now"
    mkdir star_${STAR_VERSION}
    cd star_${STAR_VERSION}
fi

# Initialize a star specific README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/phoenix" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Create the STAR index files
echo "Create STAR index files for 75bp read length" >> README
sbatch --export ALL,GTF="${STAR_GTF}",FASTA="${STAR_GENOME}",SJDB_OVERHANG="74",INDEX_DIR="75bpReads" ${PATH_TO_REPO}/utility_scripts/star_index.sh
fc -ln -1 >> README
echo >> README

echo "Create STAR index files for 100bp read length" >> README
sbatch --export ALL,GTF="${STAR_GTF}",FASTA="${STAR_GENOME}",SJDB_OVERHANG="99",INDEX_DIR="100bpReads" ${PATH_TO_REPO}/utility_scripts/star_index.sh
fc -ln -1 >> README
echo >> README

echo "Create STAR index files for 150bp read length" >> README
sbatch --export ALL,GTF="${STAR_GTF}",FASTA="${STAR_GENOME}",SJDB_OVERHANG="149",INDEX_DIR="150bpReads" ${PATH_TO_REPO}/utility_scripts/star_index.sh
fc -ln -1 >> README
echo >> README

echo >> README
echo "Specific script code as follows:" >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/star_index.sh >> README
echo >> README