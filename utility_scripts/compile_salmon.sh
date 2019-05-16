#!/usr/bin/env bash

## Dynamic download and Compile of any version of Salmon
# Usage compile_salmon.sh <version (ie. "2.7.0f")>

## Load required modules
module load gcc/8.2.0
module load boost/1.69.0
module load cmake/3.9.3

############################
## Parameterized Code
############################
echo
echo

# Capture version and ensure one was provided
if [ -z "$1" ]
  then
    echo "ERROR - No argument supplied"
    echo "Usage compile_salmon.sh <version (ie. "0.12.0")>"
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

# Download Salmon version (the downloads and tar archives have "v" infront of the VERSION)
wget https://github.com/COMBINE-lab/salmon/archive/v${VERSION}.tar.gz

# Decompress package
tar xvzf v${VERSION}.tar.gz

# Enter decompressed folder, create build directory and enter it
cd salmon-${VERSION}
mkdir build
cd build

# Needed to point to path for Intel TBB library otherwise got an error about g++ version not being > 5.2
cmake -DTBB_INSTALL_DIR=/usr/lib64/ ../

# Compile
make

# Finalize the install
make install

echo
echo ##########################
echo
echo "All Done"
echo

exit 0
