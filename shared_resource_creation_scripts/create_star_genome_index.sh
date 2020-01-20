#!/usr/bin/env bash

# Usage: create_star_genome_index.sh <Config.ini> <star_index_lengths.csv>

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

# Check STAR Index Lengths CSV was provided on the command line
if [ -n "$2" ]
then
  echo "Required csv file detected"
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

# Check that the reference genome for RNA was created successfully
if [ -e GENOME_FASTA_GENERATION_COMPLETE ]
then
    echo "Genome fasta exists, moving forward"
else
    echo "Genome fasta generation complete flag NOT found"
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
fi

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
echo >> README
echo "------------------------------------------------------" >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

####################################
## Determine required variables
####################################

# Determine the reference genome fasta full path
echo "Determine the full path filename of the reference genome fasta" >> README
echo "REFERENCE_GENOME_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_LINK} ".gz"`" >> README
REFERENCE_GENOME_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_LINK} ".gz"`
echo "REFERENCE_GENOME_FASTA=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_GENOME_FILENAME}" >> README
REFERENCE_GENOME_FASTA=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_GENOME_FILENAME}
echo >> README

# Determine the GTF file full path
echo "Determine the full path filename of the GTF file" >> README
echo "GTF_FILE=`basename ${GENE_MODEL_DOWNLOAD_LINK} ".gz"`" >> README
GTF_FILE=`basename ${GENE_MODEL_DOWNLOAD_LINK} ".gz"`
echo "REFERENCE_GENOME_FASTA=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_GENOME_FILENAME}" >> README
GENE_MODEL_GTF=${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/${GTF_FILE}
echo >> README

####################################
## Create the STAR index files
####################################

for line in `cat $2`
do

    OVERHANG=`echo ${line} | cut -d"," -f1`
    DIR=`echo ${line} | cut -d"," -f2`
    echo "Create STAR index files for ${DIR}" >> README
    sbatch --export ALL,STAR_VERSION="${STAR_VERSION}",GTF="${GENE_MODEL_GTF}",FASTA="${REFERENCE_GENOME_FASTA}",SJDB_OVERHANG="${OVERHANG}",INDEX_DIR="${DIR}" ${PATH_TO_REPO}/utility_scripts/star_index.sh

done
fc -ln -1 >> README
echo >> README
echo "Index Generation Script: ${PATH_TO_REPO}/utility_scripts/star_index.sh" >> README

echo "------------------------------------------------------" >> README
echo  >> README