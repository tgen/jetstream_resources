#!/usr/bin/env bash

# Automated Script to download and build COSMIC coding and non-coding mutation VCF files for usage in Phoenix workflow

# Usage: ./build_cosmic.sh

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
    echo "Parent directory NOT found, creating and moving into it now"
    mkdir -p ${PARENT_DIR}
    cd ${PARENT_DIR}
fi

# Make public_databases folder if not available
if [ -e public_databases ]
then
    echo "Public Databases folder exists, moving into it"
    cd public_databases
else
    echo "Public Databases folder NOT found, creating and moving into it now"
    mkdir -p public_databases
    cd public_databases
fi

# Make cosmic folder if not available
if [ -e cosmic ]
then
    echo "COSMIC folder exists, moving into it"
    cd cosmic
else
    echo "COSMIC folder NOT found, creating and moving into it now"
    mkdir -p cosmic
    cd cosmic
fi

# Make COSMIC release version folder if not available
if [ -e v${COSMIC_VERSION} ]
then
    echo "COSMIC v${COSMIC_VERSION} folder exists, exiting to prevent overwrite"
    exit 1
else
    echo "COSMIC v${COSMIC_VERSION} folder NOT found, creating and moving into it now"
    mkdir -p v${COSMIC_VERSION}
    cd v${COSMIC_VERSION}
fi

###################################
## Download and process the coding mutation VCF file
###################################

# Added user information and timestamp to README
USER=`whoami`
DATE=`date`
echo "Downloaded and Processed by:  ${USER}" >> README
echo ${DATE} >> README

# Generate authentication token
#EMAIL=jkeats@tgen.org
#PASSWORD=NotMyRealPassword123
echo "Please enter your COSMIC account email"
read COSMIC_ACCOUNT_EMAIL

echo "Please enter your COSMIC website authentication password"
read PASSWORD

TOKEN=`echo "${COSMIC_ACCOUNT_EMAIL}:${PASSWORD}" | base64`

## DOWNLOAD Coding Mutations file
# Get download link, this generates a custom url and access key for the next step, returned as a JSON, which is parsed into the download key
RESPONSE=$(curl -k -H "Authorization: Basic ${TOKEN}" "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v${COSMIC_VERSION}/VCF/CosmicCodingMuts.vcf.gz")
DOWNLOAD_KEY=$(echo $RESPONSE | sed 's/"}//g' | sed 's/{.*.vcf.gz//g')
# Download Cosmic Coding Mutations file (Authentication token expires after 1 hour)
curl -k "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v${COSMIC_VERSION}/VCF/CosmicCodingMuts.vcf.gz${DOWNLOAD_KEY}" > CosmicCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.gz

# The downloaded file is not bgzip so it needs to be decompressed with gunzip before processing
gunzip CosmicCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.gz

# Check the contigs in the file, for matching the contigs need to match and coming from ensemble they will be 1 not chr1, etx..
echo
echo "Unique contig list in Cosmic Coding Mutations File:"
bcftools view -H CosmicCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf | cut -f1 | sort | uniq
echo
### As expected only contigs in the VCF are 1-22, X, Y, MT in both files

# Convert to block compressed and create .tbi index, seems needed, not sure why
bcftools view \
    --threads 4 \
    --output-type z \
    --output-file CosmicCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.bgz \
    CosmicCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf
bcftools index --tbi CosmicCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.bgz

# Update the contig names using bcftools
bcftools annotate \
	--threads 4 \
	--rename-chrs ${PATH_TO_REPO}/utility_files/GRCh38_Cosmic90_2_UCSC_Contigs.txt \
	--output-type b \
	--output temp_renamed.bcf \
	CosmicCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.bgz
bcftools index --threads 4 temp_renamed.bcf

# Determine the full name and path of the DNA genome.fa.fai index file
echo "Determine the full path filename of the DNA reference genome.fa.fai file" >> README
echo "REFERENCE_DNA_GENOME_FAI=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}.fai" >> README
REFERENCE_DNA_GENOME_FAI=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}.fai
echo >> README

# Add updated contig header lines and add missing ones, required current bcftools version
bcftools reheader \
	--fai ${REFERENCE_DNA_GENOME_FAI} \
	temp_renamed.bcf \
	| \
	bcftools view \
	--output-type b \
	--output-file CosmicCodingMuts_v${COSMIC_VERSION}_hg38_tempe.bcf
bcftools index --threads 4 CosmicCodingMuts_v${COSMIC_VERSION}_hg38_tempe.bcf

rm temp_*
rm CosmicCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.bgz.tbi
rm CosmicCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.bgz
rm CosmicCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf

############################
###
### Download and process the non-coding VCF file
###
############################

# Download Cosmic Non-Coding Mutations file (Authentication token expires after 1 hour)
RESPONSE=$(curl -k -H "Authorization: Basic ${TOKEN}" "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v${COSMIC_VERSION}/VCF/CosmicNonCodingVariants.vcf.gz")
DOWNLOAD_KEY=$(echo $RESPONSE | sed 's/"}//g' | sed 's/{.*.vcf.gz//g')
curl -k "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v${COSMIC_VERSION}/VCF/CosmicNonCodingVariants.vcf.gz${DOWNLOAD_KEY}" > CosmicNonCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.gz

# The downloaded file is not bgzip so it needs to be decompressed with gunzip before processing
gunzip CosmicNonCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.gz

# Check the contigs in the file, for matching the contigs need to match and coming from ensemble they will be 1 not chr1, etx..
echo
echo "Unique contig list in Cosmic Coding Mutations File:"
bcftools view -H CosmicNonCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf | cut -f1 | sort | uniq
echo
### As expected only contigs in the VCF are 1-22, X, Y, MT in both files

# Convert to block compressed and create .tbi index, seems needed, not sure why
bcftools view \
    --threads 4 \
    --output-type z \
    --output-file CosmicNonCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.bgz \
    CosmicNonCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf
bcftools index --tbi CosmicNonCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.bgz

# Update the contig names using bcftools
bcftools annotate \
	--threads 4 \
	--rename-chrs ${PATH_TO_REPO}/utility_files/GRCh38_Cosmic90_2_UCSC_Contigs.txt \
	--output-type b \
	--output temp_renamed.bcf \
	CosmicNonCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.bgz
bcftools index --threads 4 temp_renamed.bcf

# Add updated contig header lines and add missing ones, required current bcftools version
bcftools reheader \
	--fai ${REFERENCE_DNA_GENOME_FAI} \
	temp_renamed.bcf \
	| \
	bcftools view \
	--output-type b \
	--output-file CosmicNonCodingMuts_v${COSMIC_VERSION}_hg38_tempe.bcf
bcftools index --threads 4 CosmicNonCodingMuts_v${COSMIC_VERSION}_hg38_tempe.bcf

rm temp_*
rm CosmicNonCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.bgz.tbi
rm CosmicNonCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf.bgz
rm CosmicNonCodingMuts_v${COSMIC_VERSION}_GRCh38.vcf



