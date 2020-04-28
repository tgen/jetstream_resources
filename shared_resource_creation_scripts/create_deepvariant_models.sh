#!/usr/bin/env bash

# Usage: create_deepvariant_models.sh <Config.ini>

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

if [ -e tool_resources ]
then
    echo "tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "tool_resources directory NOT found, creating and moving into it now"
    mkdir tool_resources
    cd tool_resources
fi

if [ -e "deepvariant_${DEEPVARIANT_VERSION}" ]
then
    echo "The deepvariant directory exists, exiting to prevent overwriting existing index"
    exit 2
else
    echo "The deepvariant directory was NOT found, creating and moving into it now"
    mkdir deepvariant_${DEEPVARIANT_VERSION}
    cd deepvariant_${DEEPVARIANT_VERSION}
fi

####################################
## Download deepvariant models
####################################

# Initialize a deepvariant model README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "Deepvariant creation details:" >> README
echo >> README

# Download and cleanup deepvariant release zip
echo "Create deepvariant model as follows:" >> README
wget https://github.com/google/deepvariant/releases/download/${DEEPVARIANT_VERSION}/deepvariant.zip
fc -ln -1 >> README
echo >> README
unzip -qq deepvariant.zip
fc -ln -1 >> README
# Need to output to /dev/null because find likes to report the directories as missing after they have been moved already
find . -name "${DEEPVARIANT_VERSION}+data" -exec mv {} . \; 2> /dev/null
fc -ln -1 >> README
echo >> README
rm -rf deepvariant.zip deepvariant
fc -ln -1 >> README
echo >> README
