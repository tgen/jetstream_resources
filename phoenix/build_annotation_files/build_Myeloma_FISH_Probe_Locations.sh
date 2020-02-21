#!/usr/bin/env bash

# Usage: ./build_Myeloma_FISH_Probe_Locations.sh  <resources.ini>

## This file is created for usage in the build package process

# Check if resources.ini was provided on the command line
if [ -n "$1" ]; then
  echo "Required ini file detected"
else
  echo "Input INI file not provided, exiting due to missing requirement"
  exit 1
fi

# Read required variables from configuration file
. ${1}

# navigate to parent directory
cd ${PARENT_DIR}

# make a folder for bcftools
mkdir -p hg38tgen/tool_resources/seq_fish 

# enter created directory
cd hg38tgen/tool_resources/seq_fish

module load liftOver

wget -O GRCh37_Myeloma_FISH_Probe_Locations.bed https://raw.githubusercontent.com/tgen/MMRF_CoMMpass/master/myeloma_FISH_probe_locations/Myeloma_FISH_Probe_Locations.txt

awk -F'\t' '{ $1 = "chr"$1 ; OFS = "\t" ; print $0}' GRCh37_Myeloma_FISH_Probe_Locations.bed | sort -k1,1V -k2,2n> GRCh37_Myeloma_FISH_Probe_Locations_with_chr.bed

liftOver -minMatch=0.5 GRCh37_Myeloma_FISH_Probe_Locations_with_chr.bed /home/tgenref/homo_sapiens/liftover_files/hg19ToHg38.over.chain.gz GRCh38_Myeloma_FISH_Probe_Locations.bed unMapped

# Write this full document as a README
cat $0 > README
