#!/usr/bin/env bash

# Usage: create_disease_specific_resources.sh <Config.ini>

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
    exit 2
fi

# Check gene_model directory if not available
if [ -e gene_model ]
then
    echo "Gene Model directory exists, moving into it"
    cd gene_model
else
    echo "Gene Model directory NOT found, IT IS REQUIRED, EXITING"
    exit 2
fi

# Make disease_specific_resources directory
if [ -e disease_specific_resources ]
then
    echo "disease_specific_resources directory exists, exiting to prevent overwriting"
    exit 2
else
    echo "disease_specific_resources directory NOT found, creating it now and moving into it"
    mkdir disease_specific_resources
    cd disease_specific_resources
fi

# Make specific disease directory
if [ -e ${DISEASE_NAME} ]
then
    echo "Specific disease directory exists, exiting to prevent overwriting"
    exit 2
else
    echo "Specific disease directory NOT found, creating it now and moving into it"
    mkdir ${DISEASE_NAME}
    cd ${DISEASE_NAME}
fi

echo "Downloading disease resource files"
echo "Downloading disease resource files" >> README
echo >> README
echo "wget https://raw.githubusercontent.com/tgen/MMRF_CoMMpass/master/immunoglobulin_transcripts/Homo_sapiens_GRCh38_103_ucsc_ig_ENST_to_filter_out.tsv" >> README

wget https://raw.githubusercontent.com/tgen/MMRF_CoMMpass/master/immunoglobulin_transcripts/Homo_sapiens_GRCh38_103_ucsc_ig_ENST_to_filter_out.tsv

echo >> README
echo "Done..."