#!/usr/bin/env bash

# Usage: create_canfam3toRosCfam_liftover_chain.sh <Config.ini>

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

####################################
## Navigate Directory Structure
###################################

# Check liftover directory if not available
LIFTOVER_DIR=$(dirname ${PARENT_DIR})/liftover_files
if [ -e ${LIFTOVER_DIR} ]
then
    echo "Liftover directory: ${LIFTOVER_DIR} exists, moving into it"
    cd ${LIFTOVER_DIR}
else
    echo "Liftover directory NOT found, creating and moving into it"
    mkdir -p ${LIFTOVER_DIR}
    cd ${LIFTOVER_DIR}
fi

####################################
## Generate liftover chain file(s)
####################################

# Initialize README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "${GENOME_SUBVERSION_NAME} liftover chain file creation details:" >> README
echo >> README

# Download chain file
echo "Downloading ${CANFAM3_TO_ROS_CHAIN_DOWNLOAD_LINK}" >> README
wget --no-check-certificate ${CANFAM3_TO_ROS_CHAIN_DOWNLOAD_LINK}
fc -ln -1 >> README
echo >> README

# Download aliases so that we can fix the chain file
echo "Downloading ${ROS_CHAIN_ALIASES_DOWNLOAD_LINK}" >> README
wget --no-check-certificate ${ROS_CHAIN_ALIASES_DOWNLOAD_LINK}
fc -ln -1 >> README
echo >> README

CANFAM3_TO_ROS_CHAIN=`basename ${CANFAM3_TO_ROS_CHAIN_DOWNLOAD_LINK}`
ROS_CHAIN_ALIASES=`basename ${ROS_CHAIN_ALIASES_DOWNLOAD_LINK}`

# The chain file is labelled for the assembly, not genome reference. Changing that and uncompressing to make parsing easier
gunzip -c ${CANFAM3_TO_ROS_CHAIN} > canFam3To${GENOME_ASSEMBLY_NAME}.chain

# We know that they have some text issues for lines 40 and 41, so we swap them here, also printing the first two columns only
awk -v OFS='\t' '{ if(NR==40 || NR==41) print $1, $3 ; else print $1, $2 }' ${ROS_CHAIN_ALIASES} | tail -n+2 > ${ROS_CHAIN_ALIASES}_fixed.txt
mv ${ROS_CHAIN_ALIASES}_fixed.txt ${ROS_CHAIN_ALIASES}

# Update the chain file based on the aliases
while read line; do 
  name=$(echo ${line} | cut -d' ' -f1)
  alias=$(echo ${line} | cut -d' ' -f2)
  sed -i "s/${name}/${alias}/g" canFam3To${GENOME_ASSEMBLY_NAME}.chain
done < ${ROS_CHAIN_ALIASES}

# Fixing the canFam3 source name to match our expections - no chr and chrUn_* contigs have .1 at the end
sed 's/chr//g' canFam3To${GENOME_ASSEMBLY_NAME}.chain | sed 's/Un_//g' | awk '{ if ($3 ~ /(^JH|^A)/) $3 = $3 ".1" }1' > ${FINAL_CHAIN_NAME::-3}

gzip ${FINAL_CHAIN_NAME::-3}

# Remove temp files
rm canFam3To${GENOME_ASSEMBLY_NAME}.chain 

echo
echo "Process Complete"
echo
