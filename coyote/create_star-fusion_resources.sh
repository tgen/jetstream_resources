#!/usr/bin/env bash

# This is a custom script for the coyote pipeline as it does not have a plug-n-play
# ctat library available from broad. As such we are building from scratch.
# Usage: create_star-fusion_resources.sh <Config.ini>

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

# Check gene_model directory if not available
if [ -e gene_model ]
then
    echo "Gene Model directory exists, moving into it"
    cd gene_model
else
    echo "Gene Model directory NOT found, IT IS REQUIRED, exiting as it should be created by a previous step"
    exit 1
fi

# Check specific gene model directory
if [ -e ${GENE_MODEL_NAME} ]
then
    echo "Specific Gene Model directory exists, moving into it"
    cd ${GENE_MODEL_NAME}
else
    echo "Specific Gene Model directory NOT found, IT IS REQUIRED, exiting as it should be created by a previous step"
    exit 1
fi

# Check that the required GTF was created successfully
if [ -e GENE_MODEL_GTF_GENERATION_COMPLETE ]
then
    echo "Required gene model GTF exists, moving forward"
else
    echo "Required gene model GTF DOES NOT exist, exiting"
    exit 2
fi

# Make gene_model specific tool_resources directory if not available
if [ -e tool_resources ]
then
    echo "Gene Model tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "Gene Model tool_resources directory NOT found, creating and entering it now"
    mkdir tool_resources
    cd tool_resources
fi

# Make star-fusion directory if not available
if [ -e "starFusion" ]
then
    echo "starFusion directory exists, moving into it"
    cd starFusion
else
    echo "starFusion directory NOT found, creating and moving into it now"
    mkdir starFusion
    cd starFusion
fi

# Make a build folder
mkdir -p starFusion_Resources
cd starFusion_Resources

# Create symbolic link
REFERENCE_RNA_GENOME_FASTA=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_RNA_GENOME_NAME}
GENE_MODEL_PATH=${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/${GENE_MODEL_FILENAME}

####################################
## Build Annotations from Source
###################################

# Since this is a custom build we need the entire Dfam.hmm database file
if [ ! -e Dfam.hmm ]; then
  wget --no-check-certificate https://dfam.org/releases/current/families/Dfam.hmm.gz
  gunzip Dfam.hmm.gz
fi

wget --no-check-certificate https://www.dfam.org/releases/Dfam_3.1/infrastructure/dfamscan.pl.gz
gunzip dfamscan.pl.gz
chmod +x dfamscan.pl

# Add current DIR to path
CURRENT_DIR=`pwd`
export PATH=$PATH:$CURRENT_DIR

if [ ${ENVIRONMENT} == "TGen" ]
then
  # Load required modules
  module load STAR-Fusion/1.8.1-GCC-8.2.0-2.31.1-Perl-5.28.1-Python-3.7.2
  module load blast/2.7.1
  module load hmmer/3.2.1

  # We need to hmmpress the Dfam.hmm for quicker operation
  if [ ! -e Dfam.hmm.h3m ]; then
    hmmpress Dfam.hmm
  fi

  # Use provided starFusion build script
  /packages/easybuild/software/STAR-Fusion/1.8.1-GCC-8.2.0-2.31.1-Perl-5.28.1-Python-3.7.2/ctat-genome-lib-builder/prep_genome_lib.pl \
  --CPU 20 \
  --max_readlength 150 \
  --genome_fa ${REFERENCE_RNA_GENOME_FASTA} \
  --gtf ${GENE_MODEL_PATH} \
  --pfam_db current \
  --dfam_db ${CURRENT_DIR}/Dfam.hmm
elif [ ${ENVIRONMENT} == "LOCAL" ]
then
  echo
  echo "Assuming required tools are available in $PATH"
  # Ensure the starFusion repository is available on your system and the subfolder "ctat-genome-lib-builder" is available in the $PATH

  # We need to hmmpress the Dfam.hmm for quicker operation
  if [ ! -e Dfam.hmm.h3m ]; then
    hmmpress Dfam.hmm
  fi

  prep_genome_lib.pl \
  --CPU 20 \
  --max_readlength 150 \
  --genome_fa ${REFERENCE_RNA_GENOME_FASTA} \
  --gtf ${GENE_MODEL_PATH} \
  --pfam_db current \
  --dfam_db ${CURRENT_DIR}/Dfam.hmm
  echo
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Environment Variable"
  echo "Only TGen or LOCAL are supported"
  exit 1
fi
