#!/usr/bin/env bash

# Usage: ./build_Myeloma_FISH_Probe_Locations.sh  <resources.ini>

## This file is created for usage in the build package process

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
if [ -e "seq_fish" ]
then
    echo "The seq_fish directory exists, exiting to prevent overwriting existing index"
    exit 2
else
    echo "The seq_fish directory was NOT found, creating and moving into it now"
    mkdir seq_fish
    cd seq_fish
fi

# Initialize a seq_fish index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

wget https://raw.githubusercontent.com/tgen/MMRF_CoMMpass/master/myeloma_FISH_probe_locations/GRCh38_Myeloma_FISH_Probe_Locations.bed

echo
echo "Process Complete"
echo