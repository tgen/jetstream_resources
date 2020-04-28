#!/usr/bin/env bash

# Usage: ./build_encode_blacklist.sh  <resources.ini>

## This file is created for usage with bcftools annotate to update contig names in the snpSniffer
## output vcf from "chr1" to "1" as the tools database does not support chr in the contig string

# Check if resources.ini was provided on the command line
if [ -n "$1" ]; then
  echo "Required ini file detected"
else
  echo "Input INI file not provided, exiting due to missing requirement"
  exit 1
fi

# Read required variables from configuration file
. ${1}

####################################
## Create Expected Folder Structure
###################################

# Make top level directory if not available
if [ -e ${PARENT_DIR} ]
then
    echo "Parent directory: ${PARENT_DIR} exists, moving into it"
    cd ${PARENT_DIR}
else
    echo "Parent directory NOT fount, creating and moving into it now"
    mkdir -p ${PARENT_DIR}
    cd ${PARENT_DIR}
fi

# Make public_databases folder if not available
if [ -e public_databases ]
then
    echo "Public Databases folder exists, moving into it"
    cd public_databases
else
    echo "Public Databases folder NOT fount, creating and moving into it now"
    mkdir -p public_databases
    cd public_databases
fi

# Make encode folder if not available
if [ -e encode ]
then
    echo "encode folder exists, moving into it"
    cd encode
else
    echo "encode folder NOT fount, creating and moving into it now"
    mkdir -p encode
    cd encode
fi

####################################
## Create resource file
####################################

# Download the V2 blacklist bundle
wget https://github.com/Boyle-Lab/Blacklist/archive/v2.0.tar.gz

# Unpack the download
tar xvzf v2.0.tar.gz

# enter unpacked directory
cd Blacklist-2.0/

# Initialize a reference_genome README
touch README
echo >> README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Write this full document as a README
cat $0 >> README
