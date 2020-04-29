#!/usr/bin/env bash

# Usage: ./build_lymphocyte_count_windows.sh <Config.ini>

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

# Create directory for tool
if [ -e "tgen_lymphocyteReceptor_counts" ]
then
    echo "The tgen_lymphocyteReceptor_counts directory exists, exiting to prevent overwriting existing index"
    exit 2
else
    echo "The tgen_lymphocyteReceptor_counts directory was NOT found, creating and moving into it now"
    mkdir tgen_lymphocyteReceptor_counts
    cd tgen_lymphocyteReceptor_counts
fi

# Initialize a bwa index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

## This will write out a bed file with the human b38 locations for the B-cell and T-cell receptor loci
## These are used to get a crude count to determine if a sample is B-cell or T-cell enriched

## Loci extracted from the ensembl v97 GTF by visual identification by Jonathan Keats

## IGK (IGKC - IGKV2-40
# chr2:88857161-89330429

## IGH (IGHA2 - IGHVIII-82)
# chr14:105583731-106879812

## IGL (IGLV4-69 - IGLC7)
# chr22:22030934-22923034

## TCR-alpha (TRAV1-1 - TRAC)
# chr14:21621838-22552156

## TCR-beta (TRBV1 - TRBC2) there is one V past the C2 but didn't include it
# chr7:142299177-142802748

### Print a bed file with the regions of interest (0-base corrected positions)
echo -e chr2"\t"88857160"\t"90367699"\t"IGK > lymphocyteReceptor_loci.bed
echo -e chr14"\t"105583730"\t"106879812"\t"IGH >> lymphocyteReceptor_loci.bed
echo -e chr22"\t"21991098"\t"22923034"\t"IGL >> lymphocyteReceptor_loci.bed
echo -e chr14"\t"21621837"\t"22552156"\t"TRA >> lymphocyteReceptor_loci.bed
echo -e chr7"\t"142299176"\t"142802748"\t"TRB >> lymphocyteReceptor_loci.bed
