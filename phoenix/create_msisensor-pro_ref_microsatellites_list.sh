#!/usr/bin/env bash

# Usage: create_msisensor-pro_ref_microsatellites_list.sh <Config.ini>

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
if [ -e GENOME_FASTA_GENERATION_COMPLETE ]
then
    echo "Genome fasta exists, moving forward"
else
    echo "Genome fasta generation complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
fi

# Create directory for tool resources
if [ -e tool_resources ]
then
    echo "tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "tool_resources directory NOT found, creating and moving into it now"
    mkdir tool_resources
    cd tool_resources
fi

# Create directory for MSIsensor-pro version
if [ -e "msisensor_pro_v${MSISENSOR_PRO_VERSION}" ]
then
    echo "The MSIsensor-pro directory exists, exiting to prevent overwriting existing microsatellites list"
    exit 2
else
    echo "The MSIsensor-pro directory was NOT found, creating and moving into it now"
    mkdir msisensor_pro_v${MSISENSOR_PRO_VERSION}
    cd msisensor_pro_v${MSISENSOR_PRO_VERSION}
fi

####################################
## Generate MSIsensor-pro reference genome microsatellites list
####################################

# Initialize a MSIsensor scan README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "MSIsensor-pro reference genome microsatellites list creation details:" >> README

echo >> README

# Determine required input variable fullpath
REFERENCE_DNA_GENOME_FASTA="${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}"
REFERENCE_DNA_GENOME_MICROSATELLITES_LIST="${REFERENCE_DNA_GENOME_NAME/.fa/.microsatellites.list}"

# Create MSIsensor-pro ref microsatellites file
if [ ${ENVIRONMENT} == "TGen" ]
then
  # Submit index generation job to the slurm scheduler
  sbatch --export ALL,FASTA="${REFERENCE_DNA_GENOME_FASTA}",MICROSATELLITES_LIST="${REFERENCE_DNA_GENOME_MICROSATELLITES_LIST}",MSISENSOR_PRO_MODULE="${MSISENSOR_PRO_MODULE}" ${PATH_TO_REPO}/utility_scripts/create_msisensor-pro_ref_microsatellites_list.sh
  fc -ln -1 >> README
elif [ ${ENVIRONMENT} == "LOCAL" ]
then
  echo
  echo "MSIsensor-pro will be run on the local compute"

  # Generate list of homopolymers and miscrosatelites from reference genome
  msisensor-pro scan -d "$FASTA" -o "$REFERENCE_DNA_GENOME_MICROSATELLITES_LIST"

  # Error Capture
  if [ "$?" = "0" ]
  then
    echo "PASSED_MSISENSOR_PRO_SCAN" >> README
  else
    touch FAILED_MSISENSOR_PRO_SCAN
    echo "FAILED_MSISENSOR_PRO_SCAN" >> README
    exit 1
  fi
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Enviroment Variable"
  touch FAILED_MSISENSOR_PRO_SCAN
  echo "FAILED_MSISENSOR_PRO_SCAN" >> README
  exit 1
fi

echo "Create msisensor-pro reference microsatellites list as follows:" >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/create_msisensor-pro_ref_microsatellites_list.sh >> README
echo >> README
echo >> README
