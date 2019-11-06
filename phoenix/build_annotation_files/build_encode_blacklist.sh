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

# navigate to parent directory
cd ${PARENT_DIR}

# make a folder for bcftools
mkdir -p public_databases/encode

# enter created directory
cd public_databases/encode

# Download the V2 blacklist bundle
wget https://github.com/Boyle-Lab/Blacklist/archive/v2.0.tar.gz

# Unpack the download
tar xvzf v2.0.tar.gz

# enter unpacked directory
cd Blacklist-2.0/

# Write this full document as a README
cat $0 > README