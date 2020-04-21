#!/usr/bin/env bash

# Usage: get_dbSNPvcf_contig_mappings.sh <RESOURCE.ini> <INPUT.VCF> <VERSION>

# Takes a dbSNP vcf with "NC_xxxxxx.x" contig headers to matching table

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
if [ $ENVIRONMENT == "TGen"]
then
  ## NCBI eUtils required Perl module
  module load perl/5.24.1
  EUTILS_PATH=/home/jkeats/downloads/edirect/
  module load BCFtools/1.10.1-foss-2019a
  JSON_AWK=/home/jkeats/downloads/JSON.awk-1.3/JSON.awk
else
  echo
  echo "Assuming required tools are available in $PATH"
  echo
fi

#################################################
## Paramaterized Code

# Get the list of contigs in the provided VCF file
bcftools view -h ${2} | grep "##contig" | sed 's/##contig=<ID=//g' | sed 's/>//g' > contig_list.txt

# Determine the number of contigs to process
CONTIG_COUNT=`wc -l contig_list.txt | cut -f1`
echo "Found ${CONTIG_COUNT} contigs in the ${2} file"

# Create table with header
echo -e CONTIG"\t"GI"\t"SLEN"\t"GENOME_TYPE"\t"GENOME_SUBTYPE"\t"GENOME_SUBNAME"\t"ASSEMBLY_GI"\t"ASSEMBLY_ACC"\t"EXTRA_INFO > dbSNP_${3}_MetaData.txt

echo
# Parse out matching information from NCBI nucleotide database
for CONTIG in `cat contig_list.txt`
do
    # Define current contig
    echo "--------------------------------------------"
    echo "Working on contig: ${CONTIG}"

    # Get a JSON record for the contig from NCBI nucleotide database
    ${EUTILS_PATH}esearch -db nucleotide -query "${CONTIG}" | ${EUTILS_PATH}efetch -format docsum -mode json > ${CONTIG}.json

    ### Flatten the JSON record using JSON.awk and extract need variables
    # Extract the matching UID, which is also the GI for the query
    GI=`awk -f ${JSON_AWK} ${CONTIG}.json | grep "result" | grep -w "uid" | cut -f2 | sed 's/"//g'`
    # Extract the length of the contig
    SLEN=`awk -f ${JSON_AWK} ${CONTIG}.json | grep "result" | grep -w "slen" | cut -f2`
    # Extract the genome type
    GENOME_TYPE=`awk -f ${JSON_AWK} ${CONTIG}.json | grep "result" | grep -w "genome" | cut -f2 | sed 's/"//g'`
    # Extract the genome subtype
    GENOME_SUBTYPE=`awk -f ${JSON_AWK} ${CONTIG}.json | grep "result" | grep -w "subtype" | grep -v "statistics" | cut -f2 | sed 's/"//g'`
    # Extract the genome subname
    GENOME_SUBNAME=`awk -f ${JSON_AWK} ${CONTIG}.json | grep "result" | grep -w "subname" | cut -f2 | sed 's/"//g'`
    # Extract the matching genome assembly GI, for matching with GRCh38 contigs in BAM files
    ASSEMBLY_GI=`awk -f ${JSON_AWK} ${CONTIG}.json | grep "result" | grep -w "assemblygi" | cut -f2 | sed 's/"//g'`
    # Extract the matching genome assembly accession number, for matching with GRCh38 contigs in BAM files
    ASSEMBLY_ACC=`awk -f ${JSON_AWK} ${CONTIG}.json | grep "result" | grep -w "assemblyacc" | cut -f2 | sed 's/"//g'`
    # Extract the matching extra annotations (pipe|separated|string|of|information)
    EXTRA_INFO=`awk -f ${JSON_AWK} ${CONTIG}.json | grep "result" | grep -w "extra" | cut -f2 | sed 's/"//g'`

    ### Send extracted results to new table
    echo -e ${CONTIG}"\t"${GI}"\t"${SLEN}"\t"${GENOME_TYPE}"\t"${GENOME_SUBTYPE}"\t"${GENOME_SUBNAME}"\t"${ASSEMBLY_GI}"\t"${ASSEMBLY_ACC}"\t"${EXTRA_INFO} >> dbSNP_${3}_MetaData.txt

    # Remove JSON record
    rm ${CONTIG}.json

done

