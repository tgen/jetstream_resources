#!/usr/bin/env bash

# Usage: create_bowtie2_genome_index.sh <Config.ini>

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
else
    echo "Genome fasta generation complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
fi

####################################
## Generate BWA index
####################################

if [ -e tool_resources ]
then
    echo "tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "tool_resources directory NOT found, creating and moving into it now"
    mkdir tool_resources
    cd tool_resources
fi

if [ -e "bowtie2_${BOWTIE2_VERSION}" ]
then
    echo "The Bowtie2 directory exists, exiting to prevent overwriting existing index"
    exit 2
else
    echo "The Bowtie2 directory was NOT found, creating and moving into it now"
    mkdir bowtie2_${BOWTIE2_VERSION}
    cd bowtie2_${BOWTIE2_VERSION}
fi

# Initialize a bowtie2 index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "Bowtie2 Index creation details:" >> README
echo >> README

# Determine the expected FASTA sequence filename based on the download link
FASTA_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_LINK} ".gz"`

# Determine the basename of the fasta to name the bowtie index outputs
BOWTIE_BASE=`basename ${FASTA_FILENAME} ".fa"`

# Create a symbolic link to the reference genome
ln -s ../../genome_reference/${FASTA_FILENAME} ${FASTA_FILENAME}

# Create bowtie2 index files using bowtie utility script
echo "Create bowtie index as follows:" >> README
sbatch --export ALL,FASTA="${FASTA_FILENAME}",BOWTIE2_MODULE="${BOWTIE2_MODULE}",BOWTIE_BASE="${BOWTIE_BASE}" ${PATH_TO_REPO}/utility_scripts/bowtie2_index.sh
fc -ln -1 >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/bowtie2_index.sh >> README
echo >> README
echo >> README