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

####################################
## Generate Custom snpEff Database
####################################

# Make snpEff directory if not available
if [ -e "snpEff_${SNPEFF_VERSION}" ]
then
    echo "snpEff directory exists, exiting to prevent overwrite"
    exit 1
else
    echo "snpEff directory NOT fount, creating and moving into it now"
    mkdir snpEff_${SNPEFF_VERSION}
    cd snpEff_${SNPEFF_VERSION}
fi

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
mkdir data
cd data

# Create directoy for this new annotation set
mkdir ${SNPEFF_DB_NAME}
cd ${SNPEFF_DB_NAME}

# Determine the reference genome fasta full path
echo "Determine the full path filename of the reference genome fasta" >> ../../README
echo "REFERENCE_GENOME_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_LINK} ".gz"`" >> ../../README
REFERENCE_GENOME_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_LINK} ".gz"`
echo "REFERENCE_GENOME_FASTA=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_GENOME_FILENAME}" >> ../../README
REFERENCE_GENOME_FASTA=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_GENOME_FILENAME}
echo >> ../../README

# Determine the GTF file full path
echo "Determine the full path filename of the GTF file" >> ../..README
echo "GTF_FILE=`basename ${GENE_MODEL_DOWNLOAD_LINK} ".gz"`" >> ../..README
GTF_FILE=`basename ${GENE_MODEL_DOWNLOAD_LINK} ".gz"`
echo "GENE_MODEL_GTF=${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/${GTF_FILE}" >> ../..README
GENE_MODEL_GTF=${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/${GTF_FILE}
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
    cp ${REFERENCE_GENOME_FASTA} ${SNPEFF_DB_NAME}.fa
    fc -ln -1 >> ../../README
    echo >> ../../README
    gzip ${SNPEFF_DB_NAME}.fa
    fc -ln -1 >> ../../README
    echo >> ../../README
fi

# Navigate to main snpEff folder and launch creation script
cd ../..
sbatch --export ALL,SNPEFF_VERSION="${SNPEFF_VERSION}",SNPEFF_DB_NAME="${SNPEFF_DB_NAME}",SNPEFF_CONFIG_PATH="${SNPEFF_CONFIG_PATH}" ${PATH_TO_REPO}/utility_scripts/build_snpEff_db.sh
fc -ln -1 >> README
echo >> README
