#!/usr/bin/env bash


# Usage: create_snpEff_db_index.sh <Config.ini>

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -ue

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

####################################
## Generate Custom snpEff Database
####################################

# Make snpEff directory if not available
if [ -e "snpEff_${SNPEFF_VERSION}" ]
then
    echo "snpEff directory exists, moving into it"
    cd snpEff_${SNPEFF_VERSION}
else
    echo "snpEff directory NOT fount, creating and moving into it now"
    mkdir snpEff_${SNPEFF_VERSION}
    cd snpEff_${SNPEFF_VERSION}
fi

# Initialize a snpEff specific README
echo >> README_${SNPEFF_DB_NAME}
echo "For details on file creation see the associated github repository:" >> README_${SNPEFF_DB_NAME}
echo "https://github.com/tgen/jetstream_resources/phoenix" >> README_${SNPEFF_DB_NAME}
echo "Created and downloaded by ${CREATOR}" >> README_${SNPEFF_DB_NAME}
date >> README_${SNPEFF_DB_NAME}
echo >> README_${SNPEFF_DB_NAME}

# Copy the snpEff config file to current directory
echo "Copy snpEff config file with any required modifications to directory" >> README_${SNPEFF_DB_NAME}
cp ${PATH_TO_REPO}/utility_files/snpEff.config .
fc -ln -1 >> README_${SNPEFF_DB_NAME}
echo >> README_${SNPEFF_DB_NAME}

# Capture the full path to the config file to pass to snpEff
echo "Capture the full path of the snpEff config file" >> README_${SNPEFF_DB_NAME}
SNPEFF_CONFIG_PATH=`realpath snpEff.config`
fc -ln -1 >> README_${SNPEFF_DB_NAME}
echo >> README_${SNPEFF_DB_NAME}

# Make snpEff data directory if not available
if [ -e data ]
then
    echo "snpEff data directory exists, moving into it"
    cd data
else
    echo "snpEff data directory NOT fount, creating and moving into it now"
    mkdir data
    cd data
fi

# Create directoy for this new annotation set
if [ -e ${SNPEFF_DB_NAME} ]
then
    echo "snpEff ${SNPEFF_DB_NAME} directory exists, exiting to prevent errors"
    exit 2
else
    echo "snpEff ${SNPEFF_DB_NAME} directory NOT fount, creating and moving into it now"
    mkdir ${SNPEFF_DB_NAME}
    cd ${SNPEFF_DB_NAME}
fi

# Copy GTF to new location, update name, and compress
echo "Copy GTF to snpEff folder" >> ../../README_${SNPEFF_DB_NAME}
cp ${GENE_MODEL_GTF} genes.gtf
fc -ln -1 >> ../../README_${SNPEFF_DB_NAME}
echo >> ../../README_${SNPEFF_DB_NAME}
gzip genes.gtf
fc -ln -1 >> ../../README_${SNPEFF_DB_NAME}
echo >> ../../README_${SNPEFF_DB_NAME}

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
    fc -ln -1 >> ../../README_${SNPEFF_DB_NAME}
    echo >> ../../README_${SNPEFF_DB_NAME}
    gzip ${SNPEFF_DB_NAME}.fa
    fc -ln -1 >> ../../README_${SNPEFF_DB_NAME}
    echo >> ../../README_${SNPEFF_DB_NAME}
fi

# Navigate to main snpEff folder and launch creation script
cd ../..
sbatch --export ALL,SNPEFF_VERSION="${SNPEFF_VERSION}",SNPEFF_DB_NAME="${SNPEFF_DB_NAME}",SNPEFF_CONFIG_PATH="${SNPEFF_CONFIG_PATH}" ${PATH_TO_REPO}/utility_scripts/build_snpEff_db.sh
fc -ln -1 >> README_${SNPEFF_DB_NAME}
echo >> README_${SNPEFF_DB_NAME}
