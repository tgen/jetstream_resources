#!/usr/bin/env bash

# Usage: create_snpSniffer_reference.sh <Config.ini>

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

# Create tool resources directory if needed
if [ -e tool_resources ]
then
    echo "tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "tool_resources directory NOT found, creating and moving into it now"
    mkdir tool_resources
    cd tool_resources
fi

# Create snpSniffer directory or exit if existing
if [ -e "snpSniffer" ]
then
    echo "The snpSniffer directory exists, exiting to prevent overwriting existing index"
    exit 2
else
    echo "The snpSniffer directory was NOT found, creating and moving into it now"
    mkdir snpSniffer
    cd snpSniffer
fi

####################################
## Download Required snpSniffer files
####################################

# Initialize a README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "BWA Index creation details:" >> README
echo >> README

# Download snpSniffer files
echo "Downloading snpSniffer files as follows:" >> README
wget ${SNPSNIFFER_POSITION_TABLE}
fc -ln -1 >> README
echo >> README

wget ${SNPSNIFFER_DATABASE}
fc -ln -1 >> README
echo >> README
echo >> README