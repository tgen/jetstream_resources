#!/usr/bin/env bash

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

# Check gene_model directory if not available
if [ -e gene_model ]
then
    echo "Gene Model directory exists, moving into it"
    cd gene_model
else
    echo "Gene Model directory NOT found, IT IS REQUIRED, EXITING"
    exit 1
fi

# Check specific gene model directory
if [ -e ${GENE_MODEL_NAME} ]
then
    echo "Specific Gene Model directory exists, moving into it"
    cd ${GENE_MODEL_NAME}
else
    echo "Specific Gene Model directory NOT found, IT IS REQUIRED, EXITING"
    exit 1
fi

## module to load
module load HTSlib/1.10-GCC-8.2.0-2.31.1
module load BEDTools/2.29.0-GCC-8.2.0-2.31.1
module load python/3.6.2
module load parallel/20171222

# --------------------------------------------------
# REQUIRED SCRIPTS
# MODIFIED APPROPRIATELY
# --------------------------------------------------
BASH_SCRIPT_CDS_PADDING="${PATH_TO_REPO}/utility_scripts/padCDS/CDS_padding_Xbp.sh"
DIR_REQUIRED_SCRIPTS=${PATH_TO_REPO}/utility_scripts/padCDS/scripts
export PATH=${DIR_REQUIRED_SCRIPTS}:${PATH}
SNPEFF_JAR_FP=${SNPEFF}
SNPEFF_CONFIG_FILE_FP=${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/tool_resources/snpEff_${SNPEFF_VERSION}/snpEff.config

# --------------------------------------------------
## INPUTS
# --------------------------------------------------
TOOLNAME="snpeff"
DIR_WORK=${PWD}
CPUS=10
GENOME_VERSION=${SNPEFF_DB_NAME}
GTF=${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/${GENE_MODEL_FILENAME}
NBP=2 ## Number of bases to pad to the CDS intervals [default is 2]
##OUTPUTS
FILE_SNPEFF_DEBUG_OUTPUT=${TOOLNAME}_output_with_listed_canonical_Transcripts.bed
FILE_LIST_SNPEFF_CANONICAL_TRANSCRIPTS=${TOOLNAME}_canonical_Transcripts_1col.txt


# --------------------------------------------------
## Some Checkings ...
# --------------------------------------------------
if [[ ! -w $(dirname ${GTF}) ]]
then
  echo -e "ERROR: This Directory Where the GTF is currently residing MUST be writable: $(dirname ${GTF}) \nChange DIR permission or move GTF to another writable directory; Aborting; "
  exit 1
fi

if [[ ! -w ${PWD} ]]
then
  echo -e "ERROR: CURRENT DIR MUST be writable: ${PWD} \nChange DIR permission or move to a writable directory; Aborting; "
  exit 1
fi



##@@@@@@@@@@@@@@@@@
## MAIN
##@@@@@@@@@@@@@@@@@
##@@@@@@@@@@@@@@@@@
## SNPEFF
##@@@@@@@@@@@@@@@@@
TOOLNAME="snpeff"
echo -e "capturing the list of ${TOOLNAME} canonical transcripts"
java -Xmx8g -jar ${SNPEFF_JAR_FP} -d -v -canon -config ${SNPEFF_CONFIG_FILE_FP} ${GENOME_VERSION} ${GTF} &> ${FILE_SNPEFF_DEBUG_OUTPUT}

# --------------------------------------------------
## parsing the ${TOOLNAME} output to only get the list of CANONICAL TRANSCRIPTS
# --------------------------------------------------
echo -e "parsing ${FILE_SNPEFF_DEBUG_OUTPUT} "
cat  ${FILE_SNPEFF_DEBUG_OUTPUT} | awk '/Canonical transcripts:/,/done/' | grep -vE "Canonical transcripts:|done." | cut -f3- | awk 'NR>1' | cut -f3 > ${FILE_LIST_SNPEFF_CANONICAL_TRANSCRIPTS}

echo -e "count of canonical transcripts captured: $(cat ${FILE_LIST_SNPEFF_CANONICAL_TRANSCRIPTS} | wc -l)"

# --------------------------------------------------
## FILTERING
# --------------------------------------------------
echo -e "Filtering by canonical transcripts the GTF file '${GTF}'"
FLT_GTF_FILTERED=${GTF/.gtf/.flt_canonical_${TOOLNAME}.gtf}
grep -wFf ${FILE_LIST_SNPEFF_CANONICAL_TRANSCRIPTS} ${GTF} > ${FLT_GTF_FILTERED}

# --------------------------------------------------
## PADDING CDS
# --------------------------------------------------
OUTFILE_SNPEFF=${FLT_GTF_FILTERED/.gtf/.cds_padded.bed}
echo -e "padding CDS ... "
bash ${BASH_SCRIPT_CDS_PADDING} \
--dir_work ${DIR_WORK} \
--gtf ${FLT_GTF_FILTERED} \
--outfilename ${OUTFILE_SNPEFF} \
--cpus ${CPUS} \
--nbp ${NBP} \
--cleanup


##@@@@@@@@@@@@@@@@
## VEP
##@@@@@@@@@@@@@@@@
TOOLNAME="vep"
echo -e "running mySQL ... "
module load DBD-mysql/4.050-foss-2019a-Perl-5.28.1 ;
mysql -u anonymous -h ensembldb.ensembl.org -P ${ENSEMBL_PORT} -e "SELECT transcript.stable_id FROM ${ENSEMBL_DATABASE}.transcript JOIN ${ENSEMBL_DATABASE}.gene ON ${ENSEMBL_DATABASE}.transcript.transcript_id=${ENSEMBL_DATABASE}.gene.canonical_transcript_id WHERE transcript.stable_id LIKE 'ENS%'" > ${GENOME_VERSION}_canonical_transcript_ids_from_mysql_command.txt;

echo -e "removing header ..."
awk 'NR>1' ${GENOME_VERSION}_canonical_transcript_ids_from_mysql_command.txt > ${GENOME_VERSION}_canonical_transcript_ids_from_mysql_command.noHeader.txt

echo -e "count of canonical transcripts captured: $(cat ${GENOME_VERSION}_canonical_transcript_ids_from_mysql_command.noHeader.txt | wc -l)"

# --------------------------------------------------
# FILTERING
# --------------------------------------------------
echo -e "filtering by ${TOOLNAME} canonical_transcript_id the GTF file ... using fgrep ... "
FLT_GTF_FILTERED=${GTF/.gtf/.flt_canonical_${TOOLNAME}.gtf}
fgrep -wf ${GENOME_VERSION}_canonical_transcript_ids_from_mysql_command.noHeader.txt ${GTF} > ${FLT_GTF_FILTERED}


# --------------------------------------------------
## PADDING CDS
# --------------------------------------------------
OUTFILE_VEP=${FLT_GTF_FILTERED/.gtf/.cds_padded.bed}
echo -e "padding CDS ... "
bash ${BASH_SCRIPT_CDS_PADDING} \
--dir_work ${DIR_WORK} \
--gtf ${FLT_GTF_FILTERED} \
--outfilename ${OUTFILE_VEP} \
--cpus ${CPUS} \
--nbp ${NBP} \
--cleanup


## INTERSECTING with CAPTURE  : PERFORMED in the PIPELINE not Here
#KIT_BED_EXTENDED=/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/capture_kits/S5U/S5U_hg38tgen_ensembl_v98.extended.bed
#bedtools intersect -a ${OUTFILE_VEP} -b ${KIT_BED_EXTENDED} > ${OUTFILE_VEP/.bed/.isec_S5U.bed}
