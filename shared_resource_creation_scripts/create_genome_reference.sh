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
module load samtools/${SAMTOOLS_VERSION}

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
echo "wget ${GENOME_FASTA_DOWNLOAD_LINK}"
wget ${GENOME_FASTA_DOWNLOAD_LINK}

echo "## Download reference fasta checksum from ${GENOME_SOURCE}" >> README
echo >> README
echo "wget ${GENOME_FASTA_MD5_DOWNLOAD_LINK}"
wget ${GENOME_FASTA_MD5_DOWNLOAD_LINK}

# Check MD5SUM

# Decompressed the downloaded reference fasta
GENOME_FASTA_DOWNLOAD_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_LINK}`

gunzip ${GENOME_FASTA_DOWNLOAD_FILENAME}

# Add symbolic link to indicate which FASTA is used by BWA
GENOME_FASTA_DECOMPRESSED_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_FILENAME} ".gz"`
ln -s ${GENOME_FASTA_DECOMPRESSED_FILENAME} BWA_FASTA

# Create faidx and dict files
echo "Create faidx index using samtools" >> README
echo "samtools faidx ${GENOME_FASTA_NAME}"
samtools faidx ${GENOME_FASTA_NAME}
fc -ln -1 >> README
echo >> README

echo "Create BWA dictionary file using samtools" >> README
GENOME_FASTA_DECOMPRESSED_BASENAME=`basename ${GENOME_FASTA_NAME} ".fa"`
samtools dict --assembly ${GENOME_ASSEMBLY_NAME} \
    --species ${SPECIES} \
    --uri ${GENOME_FASTA_DECOMPRESSED_FILENAME} \
    --output ${GENOME_FASTA_DECOMPRESSED_BASENAME}.dict \
    ${GENOME_FASTA_DECOMPRESSED_FILENAME}
fc -ln -1 >> README
echo >> README


# Add flag to top level to indicate process is complete
touch ${TOPLEVEL_DIR}/GENOME_FASTA_GENERATION_COMPLETE

