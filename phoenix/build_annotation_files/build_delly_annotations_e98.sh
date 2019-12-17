#!/usr/bin/env bash

# This version of python3 with pandas is required
module load python/3.6.0

set -euo pipefail

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

ENSEMBL_VERSION=98
URL_FILE_CONTIG_EXCL="https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/human.hg38.excl.tsv"

DIR_IN=/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v${ENSEMBL_VERSION}
EXPECTED_GTF_FILE=${DIR_IN}/Homo_sapiens.GRCh38.${ENSEMBL_VERSION}.ucsc.gtf
DIR_OUT=/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v${ENSEMBL_VERSION}/tool_resources/delly
#DIR_OUT="."


COUNT_GTF_FILES_FOUND=$(find $(dirname ${EXPECTED_GTF_FILE}) -type f -name "$(basename ${EXPECTED_GTF_FILE})" | wc -l)
if [[ ${COUNT_GTF_FILES_FOUND} -eq 1 ]] ; 
then 
	GTF_FILE=$(find $(dirname ${EXPECTED_GTF_FILE}) -type f -name "$(basename ${EXPECTED_GTF_FILE})")
else
	echo -e "ERROR: FNF:  Expected file << ${EXPECTED_GTF_FILE} >>  NOT FOUND" 
	exit 1
fi

if [[ ! -e $(dirname $(readlink -f $0))/../../utility_scripts/parse_gtf.py && ! -e $(dirname $(readlink -f $0))/../../utility_scripts/GTF.py ]]
then
	echo -e "ERROR: Missing Python scripts; Scripts parse_gtf and GTF.py scripts MUST be in << utility_scripts >> directory "
	exit 1
fi

## definition of the output files here:
mkdir -p ${DIR_OUT}
FILE_LIST_BIOTYPES=${DIR_OUT}/list_biotypes_to_keep_for_grep.txt
make_file_list_biotypes_to_keep ${FILE_LIST_BIOTYPES}
INTERMEDIATE_GTF=${DIR_OUT}/$(basename ${GTF_FILE%.*})_gene.gtf
BED_ANNOFILE_ALL_BIOTYPES=${DIR_OUT}/delly_anno_$(basename ${GTF_FILE%.*}).all_biotypes.bed
BED_ANNOFILE_KEPT_BIOTYPES=${DIR_OUT}/delly_anno_$(basename ${GTF_FILE%.*}).bed



echo -e "Extracting lines from GTF with gene information only ..."
cat ${GTF_FILE} | awk -F"\t" ' $3=="gene" { OFS="\t" ; print} ' > ${INTERMEDIATE_GTF}
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

