#!/usr/bin/env bash

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -ue

# Automated Script to download and configure VEP database for Phoenix workflow

# Usage: ./create_vep_database.sh <resources.ini>
# Sbatch Usage: sbatch --nodes=1 --cpus-per-task=10 --time=0-72:00:00 create_vep_database.sh <resources.ini>

# Check to resources.ini was provided on the command line
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
cd ${PARENT_DIR}

# Check if tool specific resources directory exists
if [ -e tool_specific_resources ]
then
    echo "Tool Specific Resources directory exists, moving into it"
    cd tool_specific_resources
else
    echo "Tool Specific Resources directory not found, creating it and moving into is"
    mkdir tool_specific_resources
    cd tool_specific_resources
fi

# Check if VEP directory exists
if [ -e vep ]
then
    echo "VEP directory exists, moving into it"
    cd vep
else
    echo "VEP directory not found, creating it and moving into is"
    mkdir vep
    cd vep
fi

# Check if version specific directory exists
if [ -e v${ENSEMBL_VERSION} ]
then
    echo "VEP version ${ENSEMBL_VERSION} directory already exists, exiting to prevent overwriting"
    echo "ERROR - CHECK YOUR CONFIGURATION"
    exit 1
else
    echo "VEP version ${ENSEMBL_VERSION} directory not found, creating it and moving into is"
    mkdir v${ENSEMBL_VERSION}
    cd v${ENSEMBL_VERSION}
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
wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/variation/indexed_vep_cache/homo_sapiens_vep_${ENSEMBL_VERSION}_GRCh38.tar.gz
fc -ln -1 >> README
echo >> README

echo "Download version specific VEP database file checksums" >> README
wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/variation/indexed_vep_cache/CHECKSUMS
fc -ln -1 >> README
echo >> README

echo "Compare checksums" >> README
# extract checksum for the vep version as the file contains checksums for multiple species
# Ensembl now uses "sum" for check sum validation
# Extract the provided checksum and number of 512bit blocks
PROVIDED_CHECKSUM=`grep "homo_sapiens_vep_${ENSEMBL_VERSION}_GRCh38.tar.gz" CHECKSUMS | cut -d" " -f1`
PROVIDED_512bitBLOCKS=`grep "homo_sapiens_vep_${ENSEMBL_VERSION}_GRCh38.tar.gz" CHECKSUMS | cut -d" " -f2`
# Calculate the checksum of the downlaoded file
VALIDATION_SUM=`sum homo_sapiens_vep_${ENSEMBL_VERSION}_GRCh38.tar.gz`
VALIDATION_CHECKSUM=`echo ${VALIDATION_SUM} | cut -d" " -f1`
VALIDATION_512bitBLOCKS=`echo ${VALIDATION_SUM} | cut -d" " -f2`
# Validate Checksum
if [ ${PROVIDED_CHECKSUM} -eq ${VALIDATION_CHECKSUM} ]
then
  echo "Complete: checksum validation"
else
  echo "FAILED: checksum validation"
  touch FAILED_CHECKSUM_VALIDATION
  exit 1
fi
# Validate 512 bit blocks
if [ ${PROVIDED_512bitBLOCKS} -eq ${VALIDATION_512bitBLOCKS} ]
then
  echo "Complete: checksum 512bit blocks validation"
else
  echo "FAILED: checksum 512bit blocks validation"
  touch FAILED_CHECKSUM_512bitBLOCK_VALIDATION
  exit 1
fi

echo >> README

echo "Decompress Package for usage" >> README
tar xzf homo_sapiens_vep_${ENSEMBL_VERSION}_GRCh38.tar.gz
fc -ln -1 >> README
echo >> README

############################
###
### Download and Processing Notes
###
############################