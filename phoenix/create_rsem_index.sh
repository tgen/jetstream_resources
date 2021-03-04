#!/usr/bin/env bash

# Usage: create_rsem_index.sh <Config.ini>

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
    echo "Gene Model directory NOT found, IT IS REQUIRED, EXITING"
    exit 1
fi

# Check specific gene model directory
if [ -e ${GENE_MODEL_NAME} ]
then
    echo "Specific Gene Model directory exists, moving into it"
    cd ${GENE_MODEL_NAME}
else
    echo "Specific Gene Model directory NOT found, IT IS REQUIRED, EXITING"
    exit 1
fi

# Check to ensure that transcriptome fasta was created
if [ -e TRANSCRIPTOME_FASTA_GENERATION_COMPLETE ]
then
    echo "Transcriptome fasta exists, moving forward"
else
    echo "Transcriptome fasta does not exist or was not created succesfully"
    echo "Please try again or check gene model output"
    exit 1
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

# Make rsem index directory if not available
if [ -e "rsem_${RSEM_VERSION}" ]
then
    echo "RSEM directory exists, moving into it"
    cd rsem_${RSEM_VERSION}
else
    echo "RSEM directory NOT found, creating and moving into it now"
    mkdir rsem_${RSEM_VERSION}
    cd rsem_${RSEM_VERSION}
fi

####################################
## Generate RSEM Index
####################################

# Initialize a rsem specific README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Create the Salmon index
if [ ${ENVIRONMENT} == "TGen" ]
then
  # Load the expected rsem module
  module load RSEM/${RSEM_VERSION}
  fc -ln -1 >> README
elif [ ${ENVIRONMENT} == "LOCAL" ]
then
  echo
  echo "RSEM Index will be created on the local compute"
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Enviroment Variable"
  touch FAILED_RSEM_INDEX
  echo "FAILED_RSEM_INDEX" >> README
  exit 1
fi

# Generate RSEM Index Files
rsem-prepare-reference \
	--gtf ${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/${GENE_MODEL_FILENAME} \
	${TOPLEVEL_DIR}/genome_reference/${REFERENCE_RNA_GENOME_NAME} \
	hg38tgen

# Error Capture
if [ "$?" = "0" ]
then
    echo "PASSED_RSEM_INDEX" >> README
else
    touch FAILED_RSEM_INDEX
    echo "FAILED_RSEM_INDEX" >> README
    exit 1
fi

echo >> README
