#!/usr/bin/env bash

# Usage: create_exome_capture_resources.sh <Config.ini> <capture_kits.csv>

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -ue

# Check to resources.ini was provided on the command line
if [ -n "$1" ]
then
  echo "Required ini file detected"
else
  echo "Input INI file not provided, exiting due to missing requirement"
  exit 1
fi

# Read required variables from configuration file
. ${1}

# Check to resources.ini was provided on the command line
if [ -n "$2" ]
then
  echo "Required csv List of capture kits detected"
  echo "Provided kits are:"
  cat $2
  echo
else
  echo "Required capture kit list NOT provided, exiting due to missing requirement"
  exit 1
fi

####################################
## Navigate Directory Structure
###################################

# Check top level directory if not available exit
if [ -e ${TOPLEVEL_DIR} ]
then
    echo "Top level directory: ${TOPLEVEL_DIR} exists, moving into it"
    cd ${TOPLEVEL_DIR}
else
    echo "Top level directory NOT found, IT IS REQUIRED, EXITING"
    exit 1
fi

# Check that reference genome is created as the dict file is required
# Check that the reference genome was created successfully
if [ -e GENOME_FASTA_GENERATION_COMPLETE ]
then
    echo "Genome fasta exists, moving forward"
else
    echo "Genome fasta generation complete flag NOT found"
    echo "Try again later as this is required"
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

# Make capturetool_resources directory if not available
if [ -e capture_kits ]
then
    echo "Capture Kits directory exists, moving into it"
    cd capture_kits
else
    echo "Capture Kits directory NOT found, creating and entering it now"
    mkdir capture_kits
    cd capture_kits
fi

####################################
## Generate Required Capture Kit Resoureces
####################################

# load required tools
module load GATK/4.1.4.0-GCCcore-8.2.0-Java-1.8
module load BEDTools/2.29.0-GCC-8.2.0-2.31.1
module load Python/3.7.2-foss-2019a

# Initialize a capture kit specific README
echo >> README
echo "------------------------------------------------------" >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/phoenix" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Create the capture kit files
for line in `cat ${2}`
do

    KIT_CODE=`echo ${line} | cut -d"," -f1`
    KIT_NAME=`echo ${line} | cut -d"," -f2`
    TARGETS_BED=`echo ${line} | cut -d"," -f3`
    BAITS_BED=`echo ${line} | cut -d"," -f4`

    # Check if capture kit directory already exists in this location
    if [ -e $KIT_CODE ]
    then
        echo "The $KIT_CODE capture kit directory already exists it will not be created to prevent overwriting existing data"
        echo "The $KIT_CODE capture kit directory already exists it will not be created to prevent overwriting existing data" >> README
    else
        echo "The $KIT_CODE capture kit directory was NOT found, creating and entering it now"
        echo "The $KIT_CODE capture kit directory was NOT found, creating and entering it now" >> README
        # Touch a file to create a kit code to product link
        touch ${KIT_CODE}_${KIT_NAME}
        mkdir $KIT_CODE
        cd $KIT_CODE

        # Run script to create needed files from expected inputs
        python3 ${PATH_TO_REPO}/utility_scripts/make_exome_refpack.py \
            -t ${PARENT_DIR}/capture_kits/${KIT_CODE}/source_files_ucsc/${TARGETS_BED} \
            -b ${PARENT_DIR}/capture_kits/${KIT_CODE}/source_files_ucsc/${BAITS_BED} \
            -r ${REFERENCE_DNA_GENOME_BASENAME} \
            -g ${GENE_MODEL_GTF} \
            -o ${KIT_CODE}_${GENOME_SUBVERSION_NAME}_${GENE_MODEL_NAME}


        # Back up to main directory
        cd ..
    fi

    ## Exit after first as a test
    exit 0

done

echo "------------------------------------------------------" >> README