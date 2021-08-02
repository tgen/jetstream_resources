#!/usr/bin/env bash

# Usage: create_snpEff_db.sh <Config.ini>

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

# Check that the reference genome for RNA was created successfully
if [ -e RNA_FASTA_GENERATION_COMPLETE ]
then
    echo "RNA fasta exists, moving forward"
else
    echo "RNA fasta generation complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
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

# Make snpEff directory if not available
if [ -e "snpEff_${SNPEFF_VERSION}" ]
then
    echo "snpEff directory exists, moving into it"
    cd snpEff_${SNPEFF_VERSION}
else
    echo "snpEff directory NOT found, creating and moving into it now"
    mkdir snpEff_${SNPEFF_VERSION}
    cd snpEff_${SNPEFF_VERSION}
fi

####################################
## Generate Custom snpEff Database
####################################

# Initialize a snpEff specific README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Copy the snpEff config file to current directory
echo "Copy snpEff config file with any required modifications to directory" >> README
cp ${PATH_TO_REPO}/utility_files/snpEff.config .
fc -ln -1 >> README
echo >> README

# Capture the full path to the config file to pass to snpEff
echo "Capture the full path of the snpEff config file" >> README
SNPEFF_CONFIG_PATH=`realpath snpEff.config`
fc -ln -1 >> README
echo >> README

# Make snpEff data directory if not available
if [ -e data ]
then
    echo "snpEff data directory exists, moving into it"
    cd data
else
    echo "snpEff data directory NOT found, creating and moving into it now"
    mkdir data
    cd data
fi

# Create directoy for this new annotation set
if [ -e ${SNPEFF_DB_NAME} ]
then
    echo "snpEff ${SNPEFF_DB_NAME} directory exists, exiting to prevent errors"
    exit 2
else
    echo "snpEff ${SNPEFF_DB_NAME} directory NOT found, creating and moving into it now"
    mkdir ${SNPEFF_DB_NAME}
    cd ${SNPEFF_DB_NAME}
fi

# Determine the reference genome fasta full path
echo "Determine the full path filename of the reference genome fasta" >> ../../README
echo "REFERENCE_DNA_GENOME_FASTA=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}" >> ../../README
REFERENCE_DNA_GENOME_FASTA=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}
echo >> ../../README

# Determine the GTF file full path
echo "Determine the full path filename of the GTF file" >> ../..README
echo "GENE_MODEL_GTF=${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/${GENE_MODEL_FILENAME}" >> ../..README
GENE_MODEL_GTF=${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/${GENE_MODEL_FILENAME}
echo >> ../../README

# Copy GTF to new location, update name, and compress
echo "Copy GTF to snpEff folder" >> ../../README
cp ${GENE_MODEL_GTF} genes.gtf
fc -ln -1 >> ../../README
echo >> ../../README
gzip genes.gtf
fc -ln -1 >> ../../README
echo >> ../../README

## Copy genome fasta, update name, and compress
cd ..
if [ -e genomes ]
then
    echo "snpEff genomes directory already exists, entering it now"
    cd genomes
else
    echo "snpEff genomes directory DOES NOT exits, creating it now"
    mkdir genomes
    cd genomes
fi

# Check that the planned reference fasta does not exist
if [ -e ${SNPEFF_DB_NAME}.fa.gz ]
then
    echo "Expected reference already exists in genomes directory"
else
    echo "Genome reference fasta DOES NOT exits in genomes directory"
    cp ${REFERENCE_DNA_GENOME_FASTA} ${SNPEFF_DB_NAME}.fa
    fc -ln -1 >> ../../README
    echo >> ../../README
    gzip ${SNPEFF_DB_NAME}.fa
    fc -ln -1 >> ../../README
    echo >> ../../README
fi

# Navigate to main snpEff folder and launch creation script
cd ../..

# Create snpEff database
if [ $ENVIRONMENT == "TGen" ]
then
  sbatch --export ALL,ENVIRONMENT="${ENVIRONMENT}",SNPEFF_VERSION="${SNPEFF_VERSION}",SNPEFF_DB_NAME="${SNPEFF_DB_NAME}",SNPEFF_CONFIG_PATH="${SNPEFF_CONFIG_PATH}" ${PATH_TO_REPO}/utility_scripts/build_snpEff_db.sh
  fc -ln -1 >> README
elif [ $ENVIRONMENT == "LOCAL" ]
then
  echo "Assuming required tools are available in $PATH or defiend in the .ini file"
  java -jar ${SNPEFF} build -gtf22 -v ${SNPEFF_DB_NAME} -c ${SNPEFF_CONFIG_PATH}

  touch CREATED_SNPEFF_${SNPEFF_DB_NAME}_DATABASE
  echo "CREATED_SNPEFF_${SNPEFF_DB_NAME}_DATABASE" >> README
  echo
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Enviroment Variable"
  echo "Only TGen or LOCAL are supported"
  exit 1
fi
