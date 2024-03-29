#!/usr/bin/env bash

# Usage: create_bwa_genome_index.sh <Config.ini>

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

# Check that the reference genome was created successfully
if [ -e GENOME_FASTA_GENERATION_COMPLETE ]
then
    echo "Genome fasta exists, moving forward"
    cd genome_reference
else
    echo "Genome fasta generation complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
fi

####################################
## Generate BWA index
####################################

# Initialize a bwa index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "BWA Index creation details:" >> README
echo >> README

FASTA=$REFERENCE_DNA_GENOME_NAME

# Create bwa index files
# shellcheck disable=SC1020
if [ ${ENVIRONMENT} == "TGen" ]
then
  # Submit sequenza resource generation job to the slurm scheduler
  sbatch --export ALL,FASTA="${FASTA}" ${PATH_TO_REPO}/utility_scripts/sequenza_gc_wiggle.sh
  fc -ln -1 >> README
elif [ ${ENVIRONMENT} == "LOCAL" ]
then
  echo
  echo "Sequenza resources will be created on the local compute"
  sequenza−utils gc_wiggle −w 50 --fasta ${FASTA} -o ${FASTA::-3}.gc50Base.wig.gz
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Enviroment Variable"
  exit 1
fi
