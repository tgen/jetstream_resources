#!/bin/sh

# Automated Script to download and build dbSNP build 154 VCF file for usage in Phoenix workflow

# Usage: ./build_dbSNP_b154.sh

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
  module load R/3.6.1-phoenix
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

# Make dbSNP folder if not available
if [ -e dbsnp ]
then
    echo "dbSNP folder exists, moving into it"
    cd dbsnp
else
    echo "dbSNP folder NOT fount, creating and moving into it now"
    mkdir -p dbsnp
    cd dbsnp
fi

# Make dbSNP release version folder if not available
if [ -e ${DBSNP_VERSION} ]
then
    echo "dbSNP ${DBSNP_VERSION} folder exists, exiting to prevent overwrite"
    exit 1
else
    echo "dbSNP ${DBSNP_VERSION} folder NOT fount, creating and moving into it now"
    mkdir -p ${DBSNP_VERSION}
    cd ${DBSNP_VERSION}
fi

####################################
## Download and Manipulate the dbSNP File
###################################

# Initialize a bcftools index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Determine the full name and path of the original DNA reference genome download
DOWNLOADED_FASTA_GZ_NAME=`basename ${REFERENCE_DNA_BASE_DOWNLOAD}`
DOWNLOADED_FASTA_GZ_FULLPATH=${TOPLEVEL_DIR}/genome_reference/downloads/${DOWNLOADED_FASTA_GZ_NAME}

# Determine the full name and path of the DNA genome.fa.fai index file
echo "Determine the full path filename of the DNA reference genome.fa.fai file" >> README
echo "REFERENCE_DNA_GENOME_FAI=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}.fai" >> README
REFERENCE_DNA_GENOME_FAI=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}.fai
echo >> README

# Download the files
wget ftp.ncbi.nih.gov/snp/archive/b154/release_notes.txt
wget ftp.ncbi.nih.gov/snp/archive/b154/VCF/GCF_000001405.38.bgz
wget ftp.ncbi.nih.gov/snp/archive/b154/VCF/GCF_000001405.38.bgz.md5
wget ftp.ncbi.nih.gov/snp/archive/b154/VCF/GCF_000001405.38.bgz.tbi
wget ftp.ncbi.nih.gov/snp/archive/b154/VCF/GCF_000001405.38.bgz.tbi.md5

# Check MD5 checksums
CHECKSUM_STATUS=`md5sum --check GCF_000001405.38.bgz.md5 | cut -d" " -f2`

if [ $CHECKSUM_STATUS == "OK" ]
then
  echo "GCF_000001405.38.bgz.md5 Checksum Validation Passed"
elif [ $CHECKSUM_STATUS == "FAILED" ]
then
  echo
  cat GCF_000001405.38.bgz.md5
  echo "GCF_000001405.38.bgz.md5 Checksum Validation FAILED"
  echo
  exit 1
else
  echo "Unexpected Validation Event"
  exit 1
fi

CHECKSUM_STATUS=`md5sum --check GCF_000001405.38.bgz.tbi.md5 | cut -d" " -f2`

if [ $CHECKSUM_STATUS == "OK" ]
then
  echo "GCF_000001405.38.bgz.tbi.md5 Checksum Validation Passed"
elif [ $CHECKSUM_STATUS == "FAILED" ]
then
  echo
  cat GCF_000001405.38.bgz.tbi.md5
  echo "GCF_000001405.38.bgz.tbi.md5 Checksum Validation FAILED"
  echo
  exit 1
else
  echo "Unexpected Validation Event"
  exit 1
fi

# Check contig seqeunces to determine if they match our hg38tgen reference genome with UCSC contif names
bcftools view -h GCF_000001405.38.bgz | cut -f1 | grep -v "#" | sort | uniq
# This showed the file uses standard NCBI nucleotide contigs like NC_000001.11 versus CM000663.2 versus chr1

## Need to build a rename key
# Step 1 - get the meta-data for the contigs in the dbSNP VCF
${PATH_TO_REPO}/utility_scripts/get_dbSNPvcf_contig_mappings.sh ${1} GCF_000001405.38.bgz b154
# Step 2 - Get the meta-data from the Reference Genome downloaded from NBCI
${PATH_TO_REPO}/utility_scripts/extract_metadata_from_fasta.sh ${DOWNLOADED_FASTA_GZ_FULLPATH}
# Step 3 - Merge output files and generate list of contigs to remove the dbSNP vcf as they are not in the assembly and the rename key
Rscript ${PATH_TO_REPO}/utility_scripts/MergeMatch_dbSNP_GRCh38_Contigs.R

# Now remove contigs that are not wanted in the dbSNP vcf as they don't exist in our refence genome (p1 versus p13 issues)
bcftools filter \
	--threads 8 \
	--targets-file ^contigs_2_remove_from_dbSNP154.bed \
	--output-type b \
	--output temp_droppedContigs.bcf \
	GCF_000001405.38.bgz
bcftools index --threads 8 temp_droppedContigs.bcf

# Now rename the contigs in the processed BCF file
bcftools annotate \
	--threads 8 \
	--rename-chrs GRCh38_dbSNP154_2_UCSC_Contigs.txt \
	--output-type b \
	--output temp_renamed.bcf \
	temp_droppedContigs.bcf
bcftools index --threads 8 temp_renamed.bcf

# Fix header to match our reference genome (in case it matters)
bcftools reheader \
	--threads 4 \
	--fai ${REFERENCE_DNA_GENOME_FAI} \
	temp_renamed.bcf \
	| \
	bcftools view \
	--threads 4 \
	--output-type b \
	--output-file dbSNP_b154_hg38tgen.bcf
bcftools index --threads 4 dbSNP_b154_hg38tgen.bcf

# Get stats to confirm all processes worked
bcftools index --threads 4 --stats dbSNP_b154_hg38tgen.bcf

# Remove temp files
rm temp_*

echo
echo "Process Complete"
echo