#!/usr/bin/env bash

# Usage: ./build_delly_annotations_e98.sh  <resources.ini>

## This file is created for usage in jetstream_build_package

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -eu

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
## Load Required Tools
###################################
if [ ${ENVIRONMENT} == "TGen" ]
then
  # This version of python3 with pandas is required
  module load python/3.6.0
elif [ ${ENVIRONMENT} == "LOCAL" ]
then
  echo
  echo "Assuming required tools are available in $PATH"
  echo
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Enviroment Variable"
  echo "Only TGen or LOCAL are supported"
  exit 1
fi

####################################
## Define Functions
###################################

function make_file_list_biotypes_to_keep(){
local FILE_OUT=$1
## List of Biotypes to KEEP in the Final GTF
echo -e "IG_C_gene
IG_C_pseudogene
IG_D_gene
IG_J_gene
IG_J_pseudogene
IG_V_gene
IG_V_pseudogene
lincRNA
miRNA
Mt_rRNA
Mt_tRNA
processed_transcript
protein_coding
rRNA
snoRNA
snRNA
TR_C_gene
TR_D_gene
TR_J_gene
TR_J_pseudogene
TR_V_gene
TR_V_pseudogene" > ${FILE_OUT}
}

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

# Make delly specific tool_resources directory if not available
if [ -e delly ]
then
    echo "Gene Model delly directory exists, moving into it"
    cd delly
else
    echo "Gene Model delly directory NOT found, creating and entering it now"
    mkdir delly
    cd delly
fi


####################################
## Script Code
###################################

# Dynamically assign variables used in script code
URL_FILE_CONTIG_EXCL="https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/human.hg38.excl.tsv"
EXPECTED_GTF_FILE=${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/${GENE_MODEL_FILENAME}
DIR_OUT=`pwd`

FILE_LIST_BIOTYPES=${DIR_OUT}/list_biotypes_to_keep_for_grep.txt
make_file_list_biotypes_to_keep ${FILE_LIST_BIOTYPES}
INTERMEDIATE_GTF=${DIR_OUT}/$(basename ${EXPECTED_GTF_FILE%.*})_gene.gtf
BED_ANNOFILE_ALL_BIOTYPES=${DIR_OUT}/delly_anno_$(basename ${EXPECTED_GTF_FILE%.*}).all_biotypes.bed
BED_ANNOFILE_KEPT_BIOTYPES=${DIR_OUT}/delly_anno_$(basename ${EXPECTED_GTF_FILE%.*}).bed

echo -e "Extracting lines from GTF with gene information only ..."
cat ${EXPECTED_GTF_FILE} | awk -F"\t" ' $3=="gene" { OFS="\t" ; print} ' > ${INTERMEDIATE_GTF}
if [[ $? -ne 0 ]] ;	then echo -e "ERROR: Subset GTF by GENE only FAILED; Aborting." ; exit 1 ; fi


echo -e "parsing GTF ... and making all_biotype bed file ..."
python $(dirname $(readlink -f $0))/../../utility_scripts/parse_gtf.py --gtf ${INTERMEDIATE_GTF} --out ${BED_ANNOFILE_ALL_BIOTYPES}
echo -e "subsetting by specific biotypes ..."
cat ${BED_ANNOFILE_ALL_BIOTYPES} | grep -wf ${FILE_LIST_BIOTYPES}  > ${BED_ANNOFILE_KEPT_BIOTYPES}
if [[ $? -ne 0 ]] ;	then echo -e "ERROR: grep FAILED ;Aborting ;" ; exit 1 ; fi

echo -e "Downloading the GRCh38 contig exclusion file provided by Delly's Author ... (if file exists; otherwise Error and Manually intervention is required) "
echo -e "running command: << wget ${URL_FILE_CONTIG_EXCL} >>"
wget ${URL_FILE_CONTIG_EXCL}
if [[ $? -ne 0 ]] ;	then echo -e "ERROR: wget FAILED for file ${URL_FILE_CONTIG_EXCL} ;Aborting ;" ; exit 1 ; fi
mv human.hg38.excl.tsv hg38.excl


echo -e "Cleaning temporary files ..."
rm ${BED_ANNOFILE_ALL_BIOTYPES} ${INTERMEDIATE_GTF} ${FILE_LIST_BIOTYPES} 

echo -e "Build Anno File for Delly: DONE"
exit

