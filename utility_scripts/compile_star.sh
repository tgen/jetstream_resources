#!/usr/bin/env bash

## Dynamic download and Compile of any version of STAR
# Usage compile_star.sh <version (ie. "2.7.0f")>

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
    echo "Usage compile_star.sh <version (ie. "2.7.0f")>"
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

# Download Star version
wget https://github.com/alexdobin/STAR/archive/${VERSION}.tar.gz

# Decompress package
tar xvzf ${VERSION}.tar.gz

# Enter decompressed folder to complile
cd STAR-${VERSION}/source/

make STAR

echo
echo ##########################
echo
echo "All Done"
echo

exit 0
