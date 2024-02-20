#!/usr/bin/env bash

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -ue

# Automated Script to download and configure VEP database for JetStream workflow

# Usage: ./create_vep_database.sh <resources.ini>
# Sbatch Usage: sbatch --nodes=1 --cpus-per-task=10 --time=0-72:00:00 create_vep_database.sh <resources.ini>

# Check that resources.ini was provided on the command line
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
## Load Required Tools
###################################

if [ ${ENVIRONMENT} == "TGen" ]
then
  # This release will use containers for processing beyond base unix tasks
  module load singularity
elif [ ${ENVIRONMENT} == "LOCAL" ]
then
  echo
  echo "Assuming singularity is available in $PATH"
  echo
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Enviroment Variable"
  echo "Only TGen or LOCAL are supported"
  exit 1
fi

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

# Check that the reference genome for was created successfully
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

#######

# Check if version specific directory exists
if [ -e vep ]
then
    echo "VEP directory already exists, exiting to prevent overwriting"
    echo "ERROR - CHECK YOUR CONFIGURATION"
    exit 1
else
    echo "VEP directory not found, creating it and moving into is"
    mkdir vep
    cd vep
fi

# Initialize a process specific README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

echo "Download version specific VEP database" >> README
echo "wget ${VEP_CACHE_DOWNLOAD_LINK}" >> README
wget ${VEP_CACHE_DOWNLOAD_LINK}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: VEP Cache Download"
else
    touch FAILED_VEP_CACHE_DOWNLOAD
    echo "FAILED: VEP Cache Download" >> README
    exit 1
fi
echo >> README


## WARNING - No provided MD5sum available


# Determine the downloaded VEP cache filename
echo "## Determine the downloaded VEP Cache filename" >> README
echo "    VEP_CACHE_DOWNLOAD_FILENAME=`basename ${VEP_CACHE_DOWNLOAD_LINK}`" >> README
VEP_CACHE_DOWNLOAD_FILENAME=`basename ${VEP_CACHE_DOWNLOAD_LINK}`
echo >> README

# Decompress the provided cache package
echo "Decompress Package for usage" >> README
tar xzf ${VEP_CACHE_DOWNLOAD_FILENAME}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: Decompression of VEP Cache"
else
    touch FAILED_VEP_CACHE_DECOMPRESSION
    echo "FAILED: Decompression of VEP Cache" >> README
    exit 1
fi
echo >> README
