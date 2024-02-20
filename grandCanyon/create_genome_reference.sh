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
## Configure and make Directory Structure
###################################

# Make top level directory if not available
if [ -e ${TOPLEVEL_DIR} ]
then
    echo "Top level directory: ${TOPLEVEL_DIR} exists, moving into it"
    cd ${TOPLEVEL_DIR}
else
    echo "Top level directory NOT found, creating and moving into it now"
    mkdir -p ${TOPLEVEL_DIR}
    cd ${TOPLEVEL_DIR}
fi

# Check if reference_genome directory exists or not, exit if so to prevent creating errors
# Make reference_genome directory if not available
if [ -e genome_reference ]
then
    echo "Genome Reference directory exists, exiting to prevent overwrite"
    echo "----WARNING----"
    exit 1
else
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
    Gene models downloaded from ${GENE_MODEL_SOURCE}
    "  >> README

    echo "Genome Reference directory NOT found, creating and moving into it now"
    mkdir genome_reference
    cd genome_reference
fi


####################################
## Create reference genomes from known sources
####################################

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
## Download REFERENCE GENOME Files
####################################

# Download bgzip compressed reference genome
echo "## Download reference fasta from ${GENOME_SOURCE}" >> README
echo "    aws s3 --no-sign-request cp ${GENOME_FASTA_DOWNLOAD_LINK} ." >> README
singularity exec -B ${REF_DIR} ${AWS_CLI_SIF} aws s3 --no-sign-request cp ${GENOME_FASTA_DOWNLOAD_LINK} .
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: download fasta"
else
    touch FAILED_DOWNLOAD_FASTA
    echo "FAILED: download fasta" >> README
    exit 1
fi
echo >> README

# Determine the downloaded fasta filename
echo "## Determine the downloaded fasta filename" >> README
echo "    GENOME_FASTA_DOWNLOAD_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_LINK}`" >> README
GENOME_FASTA_DOWNLOAD_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_LINK}`
echo >> README

#T2T consortium provides the MD5 checksums as a string on the website
# Calculate the checksum of the downloaded file
VALIDATION_SUM=`md5sum $GENOME_FASTA_DOWNLOAD_FILENAME`
# Extract the resulting hash
VALIDATION_SUM_RESULT=`echo $VALIDATION_SUM | cut -d" " -f1`
# Validate Checksum
if [ ${VALIDATION_SUM_RESULT} == ${GENOME_FASTA_MD5} ]
then
  echo "Complete: checksum validation"
else
  echo "FAILED: checksum validation"
  touch FAILED_CHECKSUM_VALIDATION
  exit 1
fi

# Download bgzip compressed reference genome index
echo "## Download reference fasta bgzip index from ${GENOME_SOURCE}" >> README
echo "    aws s3 --no-sign-request cp ${GENOME_FASTA_INDEX_DOWNLOAD_LINK} ." >> README
singularity exec -B ${REF_DIR} ${AWS_CLI_SIF} aws s3 --no-sign-request cp ${GENOME_FASTA_INDEX_DOWNLOAD_LINK} .
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: download bgzip fasta index"
else
    touch FAILED_DOWNLOAD_FASTA
    echo "FAILED: download bgzip fasta index" >> README
    exit 1
fi
echo >> README

####################################
## Process reference genome to produce needed files
###################################

# Determine the expected decompressed fasta filename
echo "## Determine the decompressed FASTA filename" >> README
echo "    GENOME_FASTA_DECOMPRESSED_FILENAME=${REFERENCE_DNA_GENOME_NAME}" >> README
GENOME_FASTA_DECOMPRESSED_FILENAME=${REFERENCE_DNA_GENOME_NAME}
echo >> README

# Decompressed the downloaded reference fasta
echo "## Decompress the Downloaded FASTA file" >> README
echo "    gunzip -c ${GENOME_FASTA_DOWNLOAD_FILENAME} > ${GENOME_FASTA_DECOMPRESSED_FILENAME}" >> README
gunzip -c ${GENOME_FASTA_DOWNLOAD_FILENAME} > ${GENOME_FASTA_DECOMPRESSED_FILENAME}
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

# Create faidx file
echo "## Create faidx index using samtools" >> README
echo "    samtools faidx ${GENOME_FASTA_DECOMPRESSED_FILENAME}" >> README
singularity exec -B ${REF_DIR} ${SAMTOOLS_SIF} samtools faidx ${GENOME_FASTA_DECOMPRESSED_FILENAME}
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


# Create dict file
echo "## Determine the decompressed FASTA basename" >> README
echo "    GENOME_FASTA_DECOMPRESSED_BASENAME=`basename ${GENOME_FASTA_DECOMPRESSED_FILENAME} ".fa"`" >> README
GENOME_FASTA_DECOMPRESSED_BASENAME=`basename ${GENOME_FASTA_DECOMPRESSED_FILENAME} ".fa"`
echo >> README

echo "## Create dictionary file using samtools" >> README
echo "    samtools dict --assembly ${GENOME_ASSEMBLY_NAME} --species "${SPECIES}" --output ${GENOME_FASTA_DECOMPRESSED_BASENAME}.dict ${GENOME_FASTA_DECOMPRESSED_FILENAME}" >> README
singularity exec -B ${REF_DIR} ${SAMTOOLS_SIF} samtools dict \
  --assembly ${GENOME_ASSEMBLY_NAME} \
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

####################################
## Create chunk/scatter intervals
###################################

# NOT NEEDED WITH CHM13, will be by chromosome as there are no gaps to leverage

####################################
## Cleanup
###################################

# Add flag to top level to indicate process is complete
echo "## Add process completion flag to top level direcory" >> README
echo "    touch ${TOPLEVEL_DIR}/GENOME_FASTA_GENERATION_COMPLETE" >> README
touch ${TOPLEVEL_DIR}/GENOME_FASTA_GENERATION_COMPLETE

# If the fastas are the same for DNA and RNA, RNA FASTA GENERATION is also complete
if [[ "${REFERENCE_DNA_GENOME_NAME}" == "${REFERENCE_RNA_GENOME_NAME}" ]]
then
  touch ${TOPLEVEL_DIR}/RNA_FASTA_GENERATION_COMPLETE
fi