#!/usr/bin/env bash

## Dynamic download and Compile of any version of STAR
# Usage compile_bcftools.sh <version (ie. "1.10.1")>

## Load required modules
module load gcc/8.2.0

############################
## Parameterized Code
############################
echo
echo

# Capture version and ensure one was provided
if [ -z "$1" ]
  then
    echo "ERROR - No argument supplied"
    echo "Usage compile_bcftools.sh <version (ie. "1.10.1")>"
    echo
    exit 1
  else
    VERSION=$1
fi

# Make the code Dynamic to any user
USER=`whoami`

# Navigate to users $HOME directory
cd $HOME

# Check if downloads directory exists, if not create it
if [ -e downloads ]
then
    echo "Found an existing downloads directory"
else
    echo "No downloads directory found"
    echo "Creating a downloads directory"
    mkdir downloads
fi

# Navigate into the downloads directory
cd downloads

# Download bcftools version (should be bundled with htslib)
wget https://github.com/samtools/bcftools/releases/download/1.10.1/bcftools-${VERSION}.tar.bz2

# Decompress package
tar xvjf bcftools-${VERSION}.tar.bz2

# Enter decompressed folder to complile
cd bcftools-${VERSION}

make

echo
echo ##########################
echo
echo "All Done"
echo

exit 0

## IMPORTANT
# In order to use the BCFtools plugins, this environment variable must be set and point to the correct location:

export BCFTOOLS_PLUGINS=/path/to/bcftools/plugins
export BCFTOOLS_PLUGINS=/home/jkeats/downloads/bcftools-1.10.1/plugins/

