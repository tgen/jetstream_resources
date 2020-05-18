#!/usr/bin/env bash

# Usage: create_salmon_index.sh <Config.ini>

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

# Make salmon index directory if not available
if [ -e "salmon_${SALMON_VERSION}" ]
then
    echo "Salmon directory exists, moving into it"
    cd salmon_${SALMON_VERSION}
else
    echo "Salmon directory NOT found, creating and moving into it now"
    mkdir salmon_${SALMON_VERSION}
    cd salmon_${SALMON_VERSION}
fi

####################################
## Generate Salmon Index
####################################

# Initialize a salmon specific README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Determine the fullpath to the transcriptome fasta file
GENE_MODEL_BASENAME=`basename ${GENE_MODEL_FILENAME} ".gtf"`
GENE_MODEL_TRANSCRIPTOME_FASTA=${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/${GENE_MODEL_BASENAME}.transcriptome.fasta

# Create the Salmon index
if [ $ENVIRONMENT == "TGen"]
then
  # Submit index generation job to the slurm scheduler
  sbatch --export ALL,SALMON_VERSION="${SALMON_VERSION}",TRANSCRIPTOME_FASTA="${GENE_MODEL_TRANSCRIPTOME_FASTA}" ${PATH_TO_REPO}/utility_scripts/salmon_index.sh
  fc -ln -1 >> README
elif [ $ENVIRONMENT == "LOCAL"]
then
  echo
  echo "SALMON Index will be created on the local compute"

  # Generate BWA Index Files
  salmon index -t ${GENE_MODEL_TRANSCRIPTOME_FASTA} -i salmon_quasi_75merPlus --type quasi -k 31

  # Error Capture
  if [ "$?" = "0" ]
  then
      echo "PASSED_SALMON_INDEX" >> README
  else
      touch FAILED_SALMON_INDEX
      echo "FAILED_SALMON_INDEX" >> README
      exit 1
  fi
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Enviroment Variable"
  touch FAILED_SALMON_INDEX
  echo "FAILED_SALMON_INDEX" >> README
  exit 1
fi

echo "Create salmon index to support typical paired-end seqeuncing with read lengths >=75bp" >> README
echo >> README
echo "Specific script code as follows:" >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/salmon_index.sh >> README
echo >> README