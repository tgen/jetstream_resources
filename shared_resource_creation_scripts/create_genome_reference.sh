#!/usr/bin/env bash

# Usage: create_genome_reference.sh <Config.ini>

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

## Load required modules to ensure needed tools are in your path
module load ${SAMTOOLS_VERSION}

####################################
## Configure and make Directory Structure
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

# Initialize a top level README
touch README
echo >> README
echo "Reference Genome and related files required for JetStream ${WORKFLOW_NAME} Workflow" >> README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "
Genome downloaded from ${GENOME_SOURCE}
Gene models downloaded from ${GENEMODEL_SOURCE}
"  >> README

####################################
## Create reference genomes from known sources
####################################

# Make reference_genome directory if not available
if [ -e genome_reference ]
then
    echo "Genome Reference directory exists, moving into it"
    cd genome_reference
else
    echo "Genome Reference directory NOT fount, creating and moving into it now"
    mkdir genome_reference
    cd genome_reference
fi

# Initialize a reference_genome README
touch README
echo >> README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README


####################################
## Download BWA REFERENCE GENOME
####################################

echo "## Download reference fasta from ${GENOME_SOURCE}" >> README
echo >> README
echo "wget ${GENOME_FASTA_DOWNLOAD_LINK}" >> README
wget ${GENOME_FASTA_DOWNLOAD_LINK}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: download fasta checksum"
else
    touch FAILED_DOWNLOAD_FASTA_CHECKSUM
    echo "FAILED: download fasta checksum" >> README
    exit 1
fi

echo "## Download reference fasta checksum from ${GENOME_SOURCE}" >> README
echo >> README
echo "wget ${GENOME_FASTA_MD5_DOWNLOAD_LINK}" >> README
wget ${GENOME_FASTA_MD5_DOWNLOAD_LINK}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: download fasta"
else
    touch FAILED_DOWNLOAD_FASTA
    echo "FAILED: download fasta" >> README
    exit 1
fi

# Check MD5SUM

# Decompressed the downloaded reference fasta
echo "## Determine the downloaded fasta filename" >> README
echo "    GENOME_FASTA_DOWNLOAD_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_LINK}`" >> README
GENOME_FASTA_DOWNLOAD_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_LINK}`
echo >> README

echo "## Decompress the Downloaded FASTA file" >> README
echo "    gunzip ${GENOME_FASTA_DOWNLOAD_FILENAME}" >> README
gunzip ${GENOME_FASTA_DOWNLOAD_FILENAME}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: gunzip fasta"
else
    touch FAILED_GUNZIP_FASTA
    echo "FAILED: gunzip fasta" >> README
    exit 1
fi
echo >> README

# Create faidx and dict files
echo "## Determine the decompressed FASTA filename" >> README
echo "    GENOME_FASTA_DECOMPRESSED_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_FILENAME} ".gz"`" >> README
GENOME_FASTA_DECOMPRESSED_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_FILENAME} ".gz"`
echo >> README

echo "## Create faidx index using samtools" >> README
echo "    samtools faidx ${GENOME_FASTA_DECOMPRESSED_FILENAME}" >> README
samtools faidx ${GENOME_FASTA_DECOMPRESSED_FILENAME}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: samtools faidx"
else
    touch FAILED_SAMTOOLS_FAIDX
    echo "FAILED: samtools faidx" >> README
    exit 1
fi
echo >> README

echo "## Determine the decompressed FASTA basename" >> README
echo "    GENOME_FASTA_DECOMPRESSED_BASENAME=`basename ${GENOME_FASTA_DECOMPRESSED_FILENAME} ".fa"`" >> README
GENOME_FASTA_DECOMPRESSED_BASENAME=`basename ${GENOME_FASTA_DECOMPRESSED_FILENAME} ".fa"`
echo >> README

echo "## Create BWA dictionary file using samtools" >> README
echo "    samtools dict --assembly ${GENOME_ASSEMBLY_NAME} --species "${SPECIES}" --output ${GENOME_FASTA_DECOMPRESSED_BASENAME}.dict ${GENOME_FASTA_DECOMPRESSED_FILENAME}" >> README
samtools dict --assembly ${GENOME_ASSEMBLY_NAME} \
    --species "${SPECIES}" \
    --output ${GENOME_FASTA_DECOMPRESSED_BASENAME}.dict \
    ${GENOME_FASTA_DECOMPRESSED_FILENAME}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: samtools dict"
else
    touch FAILED_SAMTOOLS_DICT
    echo "FAILED: samtools dict" >> README
    exit 1
fi
echo >> README

# Add flag to top level to indicate process is complete
echo "## Add process completion flag to top level direcory" >> README
echo "    touch ${TOPLEVEL_DIR}/GENOME_FASTA_GENERATION_COMPLETE" >> README
touch ${TOPLEVEL_DIR}/GENOME_FASTA_GENERATION_COMPLETE

