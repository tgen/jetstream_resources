#!/usr/bin/env bash

# Usage: create_bowtie2_genome_index.sh <Config.ini>

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
## Navigate Directory Structure
###################################

# Check top level directory if not available
if [ -e ${TOPLEVEL_DIR} ]
then
    echo "Top level directory: ${TOPLEVEL_DIR} exists, moving into it"
    cd ${TOPLEVEL_DIR}
else
    echo "Top level directory NOT found, IT IS REQUIRED, EXITING"
    exit 1
fi

# Check that the reference genome was created successfully
# Since bowtie is not alt-aware we will use the same reference as used by STAR for RNA alignments
if [ -e RNA_FASTA_GENERATION_COMPLETE ]
then
    echo "Genome fasta, used for RNA alignments exists, moving forward"
else
    echo "Genome fasta, used for RNA alignments, generation complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
fi

if [ -e tool_resources ]
then
    echo "tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "tool_resources directory NOT found, creating and moving into it now"
    mkdir tool_resources
    cd tool_resources
fi

if [ -e "bowtie2_${BOWTIE2_VERSION}" ]
then
    echo "The Bowtie2 directory exists, exiting to prevent overwriting existing index"
    exit 2
else
    echo "The Bowtie2 directory was NOT found, creating and moving into it now"
    mkdir bowtie2_${BOWTIE2_VERSION}
    cd bowtie2_${BOWTIE2_VERSION}
fi

####################################
## Generate BOWTIE index
####################################

# Initialize a bowtie2 index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "Bowtie2 Index creation details:" >> README
echo >> README

# Determine the basename of the fasta to name the bowtie index outputs
BOWTIE_BASE=`basename ${REFERENCE_RNA_GENOME_NAME} ".fa"`

# Create a symbolic link to the reference genome
ln -s ../../genome_reference/${REFERENCE_RNA_GENOME_NAME} ${REFERENCE_RNA_GENOME_NAME}

# Create index files
echo "Create bowtie index as follows:" >> README
if [ $ENVIRONMENT == "TGen"]
then
  # Create bowtie2 index files using bowtie utility script
  sbatch --export ALL,FASTA="${REFERENCE_RNA_GENOME_NAME}",BOWTIE2_MODULE="${BOWTIE2_MODULE}",BOWTIE_BASE="${BOWTIE_BASE}" ${PATH_TO_REPO}/utility_scripts/bowtie2_index.sh
  fc -ln -1 >> README
elif [ $ENVIRONMENT == "LOCAL"]
then
  echo
  echo "Bowtie2 Index will be created on the local compute"

  # Generate Bowtie Index Files
  bowtie2-build --threads ${LOCAL_COMPUTE_CORES} ${REFERENCE_RNA_GENOME_NAME} ${BOWTIE_BASE}

  # Error Capture
  if [ "$?" = "0" ]
  then
      echo "PASSED_BOWTIE_INDEX" >> README
  else
      touch FAILED_BOWTIE_INDEX
      echo "FAILED_BOWTIE_INDEX" >> README
      exit 1
  fi
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Enviroment Variable"
  touch FAILED_BOWTIE_INDEX
  echo "FAILED_BOWTIE_INDEX" >> README
  exit 1
fi

echo >> README
cat ${PATH_TO_REPO}/utility_scripts/bowtie2_index.sh >> README
echo >> README
echo >> README