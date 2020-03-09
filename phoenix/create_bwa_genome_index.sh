#!/usr/bin/env bash

# Usage: create_bwa_genome_index.sh <Config.ini>

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

# Check that the reference genome was created successfully
if [ -e GENOME_FASTA_GENERATION_COMPLETE ]
then
    echo "Genome fasta exists, moving forward"
else
    echo "Genome fasta generation complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
fi

####################################
## Generate BWA index
####################################

if [ -e tool_resources ]
then
    echo "tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "tool_resources directory NOT found, creating and moving into it now"
    mkdir tool_resources
    cd tool_resources
fi

if [ -e "bwa_${BWA_VERSION}" ]
then
    echo "The BWA directory exists, exiting to prevent overwriting existing index"
    exit 2
else
    echo "The BWA directory was NOT found, creating and moving into it now"
    mkdir bwa_${BWA_VERSION}
    cd bwa_${BWA_VERSION}
fi

# Initialize a bwa index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/phoenix" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "BWA Index creation details:" >> README
echo >> README
echo "There are no details on how the alt-aware required .alt file is created" >> README
echo "Because of this we are copying and renaming the file from BWAKIT" >> README
echo >> README

# Copy in the .alt file from bwa.kit
echo "Copy in the .alt file from bwa.kit" >> README
cp ../../genome_reference/downloads/bwa.kit/resource-GRCh38/hs38DH.fa.alt GRCh38tgen_decoy_alts_hla.fa.alt
fc -ln -1 >> README
echo >> README

# Create a symbolic link to the reference genome
ln -s ../../genome_reference/GRCh38tgen_decoy_alts_hla.fa GRCh38tgen_decoy_alts_hla.fa

# Create bwa index files using bwa utility script
echo "Create bwa index as follows:" >> README
sbatch --export ALL,FASTA="GRCh38tgen_decoy_alts_hla.fa",BWA_VERSION="${BWA_VERSION}" ${PATH_TO_REPO}/utility_scripts/bwa_index.sh
fc -ln -1 >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/bwa_index.sh >> README
echo >> README
echo >> README