#!/usr/bin/env bash

# Automated Script to download and build CLINVAR VCF file for usage in Phoenix workflow

# Usage: ./build_clinvar_20190715.sh

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
## Load Required Tools
###################################
if [ ${ENVIRONMENT} == "TGen" ]
then
  module load BCFtools/1.10.1-foss-2019a
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
## Create Expected Folder Structure
###################################

# Make top level directory if not available
if [ -e ${PARENT_DIR} ]
then
    echo "Parent directory: ${PARENT_DIR} exists, moving into it"
    cd ${PARENT_DIR}
else
    echo "Parent directory NOT fount, creating and moving into it now"
    mkdir -p ${PARENT_DIR}
    cd ${PARENT_DIR}
fi

# Make public_databases folder if not available
if [ -e public_databases ]
then
    echo "Public Databases folder exists, moving into it"
    cd public_databases
else
    echo "Public Databases folder NOT fount, creating and moving into it now"
    mkdir -p public_databases
    cd public_databases
fi

# Make clinvar folder if not available
if [ -e clinvar ]
then
    echo "clinvar folder exists, moving into it"
    cd clinvar
else
    echo "clinvar folder NOT fount, creating and moving into it now"
    mkdir -p clinvar
    cd clinvar
fi

# Make clinvar release version folder if not available
if [ -e ${CLINVAR_VERSION} ]
then
    echo "clinvar ${CLINVAR_VERSION} folder exists, exiting to prevent overwrite"
    exit 1
else
    echo "clinvar ${CLINVAR_VERSION} folder NOT fount, creating and moving into it now"
    mkdir -p ${CLINVAR_VERSION}
    cd ${CLINVAR_VERSION}
fi

####################################
## Parameterized Code
####################################

# Added user information and timestamp to README
USER=`whoami`
DATE=`date`
echo "Downloaded and Processed by:  ${USER}" > README
echo ${DATE} >> README

# Download the current clinvar files
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2019/clinvar_20190715.vcf.g*

# Remove local NIH path from MD5 file
sed -i 's/\/panfs\/pan1\/clintest\/ftp_test_prod\/vcf\/vcf_GRCh38\///g' clinvar_20190715.vcf.gz.md5

# Check MD5 checksums
CHECKSUM_STATUS=`md5sum --check clinvar_20190715.vcf.gz.md5 | cut -d" " -f2`

if [ $CHECKSUM_STATUS == "OK" ]
then
  echo "Checksum Validation Passed"
elif [ $CHECKSUM_STATUS == "FAILED" ]
then
  echo
  md5sum clinvar_20190715.vcf.gz
  cat clinvar_20190715.vcf.gz.md5
  echo "Checksum Validation FAILED"
  echo
  exit 1
else
  echo "Unexpected Validation Event"
  exit 1
fi

# Check contig sequences to determine if they match our hg38tgen reference genome with UCSC contif names
zcat clinvar_20190715.vcf.gz | cut -f1 | grep -v "#" | sort | uniq
# Found 1-22,X,Y,MT, and NW_009646201.1, which is a GRCh38.p13 patch for the ABO*A1.02 allele

## Remove the NW_009646201.1 patch and convert the contig names in clinvar_20190715.vcf.gz
bcftools filter \
	--threads 4 \
	--targets ^NW_009646201.1  \
	--output-type b \
	--output temp_droppedContigs.bcf \
	clinvar_20190715.vcf.gz
bcftools index --threads 4 temp_droppedContigs.bcf

## Update the remaining contigs to the UCSC Naming convention
bcftools annotate \
	--threads 4 \
	--rename-chrs ${PATH_TO_REPO}/utility_files/GRCh38_Cosmic90_2_UCSC_Contigs.txt \
	--output-type b \
	--output temp_renamed.bcf \
	temp_droppedContigs.bcf
bcftools index --threads 4 temp_renamed.bcf

# Determine the full name and path of the DNA genome.fa.fai index file
echo "Determine the full path filename of the DNA reference genome.fa.fai file" >> README
echo "REFERENCE_DNA_GENOME_FAI=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}.fai" >> README
REFERENCE_DNA_GENOME_FAI=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}.fai
echo >> README

## Update the contig header so lengths are included
bcftools reheader \
	--fai ${REFERENCE_DNA_GENOME_FAI} \
	temp_renamed.bcf \
	| \
	bcftools view \
	--output-type b \
	--output-file clinvar_20190715_hg38tgen.bcf
bcftools index --threads 4 clinvar_20190715_hg38tgen.bcf

# Remove all the temp files
rm temp_*
