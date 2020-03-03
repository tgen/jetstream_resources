#!/usr/bin/env bash

# Automated Script to download and build dbSNP VCF files for usage in Suncity workflow

# Usage: ./build_dbSNP_b138_broadBundle.sh <config.ini>

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
## Load Required Tools
###################################

# Load required modules
module load BCFtools/1.10.1-foss-2019a
module load R/3.6.1-phoenix

####################################
## Create Expected Folder Structure
###################################

# Make top level directory if not available
if [ -e ${TOPLEVEL_DIR} ]
then
    echo "Top level directory: ${TOPLEVEL_DIR} exists, moving into it"
    cd ${TOPLEVEL_DIR}
else
    echo "Top level directory NOT fount, creating and moving into it now"
    mkdir -p ${TOPLEVEL_DIR}
    cd ${TOPLEVEL_DIR}
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

# Make dbsnp folder if not available
if [ -e dbsnp ]
then
    echo "dbSNP folder exists, moving into it"
    cd dbsnp
else
    echo "dbSNP folder NOT fount, creating and moving into it now"
    mkdir -p dbsnp
    cd dbsnp
fi

# Make dbSNP release version folder if not available
if [ -e ${DBSNP_RELEASE_VERSION} ]
then
    echo "dbSNP ${DBSNP_RELEASE_VERSION} folder exists, exiting to prevent overwrite"
    exit 1
else
    echo "dbSNP ${DBSNP_RELEASE_VERSION} folder NOT fount, creating and moving into it now"
    mkdir -p ${DBSNP_RELEASE_VERSION}
    cd ${DBSNP_RELEASE_VERSION}
fi

####################################
## Download and Process Files as Needed
###################################

# Initialize a gnomAD specific README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Download the files
wget ${BROAD_BUNDLE_DBSNP_DOWNLOAD_LINK}
wget ${BROAD_BUNDLE_DBSNP_DOWNLOAD_MD5SUM_LINK}

# Update check sum file to remove Broad Insitute path
sed -i 's/\/humgen\/gsa-scr1\/pub\/bundle\/2.8\/b37\///g' dbsnp_138.b37.vcf.gz.md5

# Check MD5 checksums - All downloads pass
md5sum --check dbsnp_138.b37.vcf.gz.md5

# Convert downloaded file to BCF
bcftools view --threads 4 --output-type b --output-file dbsnp_138.b37.bcf dbsnp_138.b37.vcf.gz

# Index file
bcftools index --thread 4 dbsnp_138.b37.bcf
