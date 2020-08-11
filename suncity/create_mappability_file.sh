#!/usr/bin/env bash

# Usage : $0 <resources.ini>


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

## Check if executables in PATH
if [[ $(type -P "gatk") ]]
then
  echo "SUCCESS: gatk found in PATH" ;
else
  echo "ERROR: gatk NOT in PATH; Need gatk installed; At least version 4.1.4.0" ;
  exit 2
fi

### make target directories:
if [[ -d ${TOPLEVEL_DIR} && -w ${TOPLEVEL_DIR} ]] ;
then
  mkdir -p ${TOPLEVEL_DIR}/public_databases/bismap/
  cd ${TOPLEVEL_DIR}/public_databases/bismap/
  if [[ $? -ne 0 ]]
  then
    echo -e "ERROR entering directory: ${TOPLEVEL_DIR}/public_databases/bismap/"
    exit 1
  fi
else
  echo -e "ERROR: DIR NOT FOUND << ${TOPLEVEL_DIR}>>; This directory MUST be created before running this script << $0 >>";
  exit 1
fi

echo -e "curdir: ${PWD}"

# Initialize a samtools_stats index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/suncity" >> README
echo "and ">> README
echo "https://bismap.hoffmanlab.org/" >> README
echo >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "$0 <resources.ini>  --> Usage details:" >> README
echo "The input is an ini file with specific defined variables; Example of ini file can be found in the pipeline folders; such as << suncity_resources.ini >>" >> README
echo >> README
echo "Why these files?" >> README
echo "Mappability Files are used with GATK CNV tool to get better accuracy in copy number calls." > README
echo >> README
echo >> README



## check if file already exist in TARGET directory
if [[ -e k100.umap.bed.gz ]]
then
  echo -e "File << k100.umap.bed.gz >> already exists. Exiting to prevent overwriting it."
  exit 2
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

