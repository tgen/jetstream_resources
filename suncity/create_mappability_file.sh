#!/usr/bin/env bash

# Usage generic : $0 <resources.ini>
# Usage example suncity: $0 <suncity_resources.ini>

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

## Get executables in PATH
module load ${GATK_MODULE}

### make target directories:
if [[ -d ${TOPLEVEL_DIR} && -w ${TOPLEVEL_DIR} ]] ;
then
  mkdir -p ${TOPLEVEL_DIR}/public_databases/bismap/ && cd ${TOPLEVEL_DIR}/public_databases/bismap/ || ( echo -e "ERROR entering directory: ${TOPLEVEL_DIR}/public_databases/bismap/" && exit 1 )
else
  echo -e "ERROR: DIR NOT FOUND << ${TOPLEVEL_DIR}>>; This directory MUST be created before running this script << $0 >>";
  exit 1
fi

echo -e ${PWD}

## check if file already exist in TARGET directory
if [[ -e k100.umap.bed.gz ]]
then
  rm k100.umap.bed.gz &>/dev/null
fi

## Downloading compressed BED file
wget ${BISMAP_UMPA_DOWNLOAD_LINK}

if [[ $? -eq 0 ]]
then
  zcat k100.umap.bed.gz | awk 'NR >1' > k100.umap.no_header.bed
  gatk IndexFeatureFile --feature-file k100.umap.no_header.bed
  if [[ $? -eq 0 ]]
  then
      echo -e "Mappability Files Creation Completed SUCCESSFULLY"
  else
    echo -e "ERROR: GATK issue;\nMappability Files Creation: FAILED"
    exit 1
  fi
else
  echo -e "ERROR: Download Issue;\nMappability Files Creation: FAILED"
  exit 1
fi
