#!/usr/bin/env bash

# Usage: create_genome_reference.sh <Config.ini>

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
  module load ${SAMTOOLS_VERSION}
  module load ${GATK_MODULE}
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
## Configure and make Directory Structure
###################################

# Make top level directory if not available
if [ -e ${TOPLEVEL_DIR} ]
then
    echo "Top level directory: ${TOPLEVEL_DIR} exists, moving into it"
    cd ${TOPLEVEL_DIR}
else
    echo "Top level directory NOT found, creating and moving into it now"
    mkdir -p ${TOPLEVEL_DIR}
    cd ${TOPLEVEL_DIR}
fi

# Check if reference_genome directory exists or not, exit if so to prevent creating errors

# Make reference_genome directory if not available
if [ -e genome_reference ]
then
    echo "Genome Reference directory exists, exiting to prevent overwrite"
    echo "----WARNING----"
    exit 1
else
    # Initialize a top level README
    touch README
    echo >> README
    echo "Reference Genome and related files required for JetStream ${WORKFLOW_NAME} Workflow" >> README
    echo >> README
    echo "For details on file creation see the associated github repository:" >> README
    echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
    echo "Created and downloaded by ${CREATOR}" >> README
    date >> README
    echo >> README
    echo "
    Genome downloaded from ${GENOME_SOURCE}
    Gene models downloaded from ${GENE_MODEL_SOURCE}
    "  >> README

    echo "Genome Reference directory NOT found, creating and moving into it now"
    mkdir genome_reference
    cd genome_reference
fi


####################################
## Create reference genomes from known sources
####################################

# Initialize a reference_genome README
touch README
echo >> README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README


####################################
## Download BWA REFERENCE GENOME
####################################

echo "## Download reference fasta from ${GENOME_SOURCE}" >> README
echo "    wget ${GENOME_FASTA_DOWNLOAD_LINK}" >> README
wget ${GENOME_FASTA_DOWNLOAD_LINK}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: download fasta"
else
    touch FAILED_DOWNLOAD_FASTA
    echo "FAILED: download fasta" >> README
    exit 1
fi
echo >> README

if [ ${GENOME_FASTA_MD5_DOWNLOAD_LINK} != "NA" ]
then
  echo "## Download reference fasta checksum from ${GENOME_SOURCE}" >> README
  echo "    wget ${GENOME_FASTA_MD5_DOWNLOAD_LINK}" >> README
  wget ${GENOME_FASTA_MD5_DOWNLOAD_LINK}
  # Error Capture
  if [ "$?" = "0" ]
  then
      echo "Completed: download fasta checksum"
  else
      touch FAILED_DOWNLOAD_FASTA_CHECKSUM
      echo "FAILED: download fasta checksum" >> README
      exit 1
  fi
  echo >> README
fi

# Determine the downloaded fasta filename
echo "## Determine the downloaded fasta filename" >> README
echo "    GENOME_FASTA_DOWNLOAD_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_LINK}`" >> README
GENOME_FASTA_DOWNLOAD_FILENAME=`basename ${GENOME_FASTA_DOWNLOAD_LINK}`
echo >> README

# Check MD5SUM
if [ ${GENOME_SOURCE} == "ensembl" ]
then
  echo "ENSEMBL is supported"
  # Ensembl now uses "sum" for check sum validation
  # Extract the provided checksum and number of 512bit blocks
  PROVIDED_CHECKSUM=`grep ${GENOME_FASTA_DOWNLOAD_FILENAME} CHECKSUMS | cut -d" " -f1`
  PROVIDED_512bitBLOCKS=`grep ${GENOME_FASTA_DOWNLOAD_FILENAME} CHECKSUMS | cut -d" " -f2`
  # Calculate the checksum of the downlaoded file
  VALIDATION_SUM=`sum ${GENOME_FASTA_DOWNLOAD_FILENAME}`
  VALIDATION_CHECKSUM=`echo ${VALIDATION_SUM} | cut -d" " -f1`
  VALIDATION_512bitBLOCKS=`echo ${VALIDATION_SUM} | cut -d" " -f2`
  # Validate Checksum
  if [ ${PROVIDED_CHECKSUM} -eq ${VALIDATION_CHECKSUM} ]
  then
    echo "Complete: checksum validation"
  else
    echo "FAILED: checksum validation"
    touch FAILED_CHECKSUM_VALIDATION
    exit 1
  fi
  # Validate 512 bit blocks
  if [ ${PROVIDED_512bitBLOCKS} -eq ${VALIDATION_512bitBLOCKS} ]
  then
    echo "Complete: checksum 512bit blocks validation"
  else
    echo "FAILED: checksum 512bit blocks validation"
    touch FAILED_CHECKSUM_512bitBLOCK_VALIDATION
    exit 1
  fi
elif [ ${GENOME_SOURCE} == "NCBI" ]
then
  echo "NCBI is supported"
  #NCBI provides a single file with the MD5 checksums for each file in the directory
  grep $GENOME_FASTA_DOWNLOAD_FILENAME md5checksums.txt > ${GENOME_FASTA_DOWNLOAD_FILENAME}.md5
  # Calculate the checksum of the downloaded file
  VALIDATION_SUM=`md5sum --check ${GENOME_FASTA_DOWNLOAD_FILENAME}.md5`
  # Extract the result, expect "OK"
  VALIDATION_SUM_RESULT=`echo $VALIDATION_SUM | cut -d" " -f2`
  # Validate Checksum
  if [ ${VALIDATION_SUM_RESULT} == "OK" ]
  then
    echo "Complete: checksum validation"
  else
    echo "FAILED: checksum validation"
    touch FAILED_CHECKSUM_VALIDATION
    exit 1
  fi
elif [ ${GENOME_SOURCE} == "1000G" ]
then
  echo "1000G is supported"
  echo "WARNING - No checksum provided the file CANNOT be validated"
else
  echo "Current Genome Source is NOT SUPPORTED"
  exit 1
fi

# Determine the expected decompressed fasta filename
echo "## Determine the decompressed FASTA filename" >> README
echo "    GENOME_FASTA_DECOMPRESSED_FILENAME=${REFERENCE_DNA_GENOME_NAME}" >> README
GENOME_FASTA_DECOMPRESSED_FILENAME=${REFERENCE_DNA_GENOME_NAME}
echo >> README

# Decompressed the downloaded reference fasta
echo "## Decompress the Downloaded FASTA file" >> README
echo "    gunzip -c ${GENOME_FASTA_DOWNLOAD_FILENAME} > ${GENOME_FASTA_DECOMPRESSED_FILENAME}" >> README
gunzip -c ${GENOME_FASTA_DOWNLOAD_FILENAME} > ${GENOME_FASTA_DECOMPRESSED_FILENAME}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: gunzip fasta"
else
    touch FAILED_GUNZIP_FASTA
    echo "FAILED: gunzip fasta" >> README
    exit 1
fi
echo >> README

# Create faidx and dict files
echo "## Create faidx index using samtools" >> README
echo "    samtools faidx ${GENOME_FASTA_DECOMPRESSED_FILENAME}" >> README
samtools faidx ${GENOME_FASTA_DECOMPRESSED_FILENAME}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: samtools faidx"
else
    touch FAILED_SAMTOOLS_FAIDX
    echo "FAILED: samtools faidx" >> README
    exit 1
fi
echo >> README

echo "## Determine the decompressed FASTA basename" >> README
echo "    GENOME_FASTA_DECOMPRESSED_BASENAME=`basename ${GENOME_FASTA_DECOMPRESSED_FILENAME} ".fa"`" >> README
GENOME_FASTA_DECOMPRESSED_BASENAME=`basename ${GENOME_FASTA_DECOMPRESSED_FILENAME} ".fa"`
echo >> README

echo "## Create BWA dictionary file using samtools" >> README
echo "    samtools dict --assembly ${GENOME_ASSEMBLY_NAME} --species "${SPECIES}" --output ${GENOME_FASTA_DECOMPRESSED_BASENAME}.dict ${GENOME_FASTA_DECOMPRESSED_FILENAME}" >> README
samtools dict --assembly ${GENOME_ASSEMBLY_NAME} \
    --species "${SPECIES}" \
    --output ${GENOME_FASTA_DECOMPRESSED_BASENAME}.dict \
    ${GENOME_FASTA_DECOMPRESSED_FILENAME}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: samtools dict"
else
    touch FAILED_SAMTOOLS_DICT
    echo "FAILED: samtools dict" >> README
    exit 1
fi
echo >> README

echo "Create 2bit genome reference for CHIP analysis steps" >> README
echo "    ${FATOTWOBIT} ${REFERENCE_DNA_GENOME_NAME} ${REFERENCE_DNA_GENOME_NAME%.fa}.2bit" >> README
${FATOTWOBIT} ${REFERENCE_DNA_GENOME_NAME} ${REFERENCE_DNA_GENOME_NAME%.fa}.2bit
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: faToTwoBit file creation"
else
    touch FAILED_faToTwoBit_FILE
    echo "FAILED: faToTwoBit file" >> README
    exit 1
fi
echo >> README

# Create chunk/scatter intervals

# Generate calling interval windows
echo "## Create sequence analysis grouping chunks" >> README
echo "    gatk ScatterIntervalsByNs --OUTPUT ${GENOME_FASTA_DECOMPRESSED_BASENAME}_NorACGT.interval_list --REFERENCE ${GENOME_FASTA_DECOMPRESSED_FILENAME}" >> README
gatk ScatterIntervalsByNs \
  --OUTPUT ${GENOME_FASTA_DECOMPRESSED_BASENAME}_NorACGT.interval_list \
  --REFERENCE ${GENOME_FASTA_DECOMPRESSED_FILENAME}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: gatk ScatterInvervalsByNs"
else
    touch FAILED_GATK_SCATTERINVERVALSbyNS
    echo "FAILED: gatk ScatterIntervalsByNs" >> README
    exit 1
fi
echo >> README

echo "    grep -w -v "Nmer" ${GENOME_FASTA_DECOMPRESSED_BASENAME}_NorACGT.interval_list > ${GENOME_FASTA_DECOMPRESSED_BASENAME}_ACGT.interval_list" >> README
grep -w -v "Nmer" ${GENOME_FASTA_DECOMPRESSED_BASENAME}_NorACGT.interval_list > ${GENOME_FASTA_DECOMPRESSED_BASENAME}_ACGT.interval_list
echo >> README

echo "    mkdir chunk_intervals" >> README
mkdir chunk_intervals
echo >> README

echo "    gatk IntervalListTools --SCATTER_COUNT 42 --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW --UNIQUE true --SORT true --BREAK_BANDS_AT_MULTIPLES_OF 0 --INPUT ${GENOME_FASTA_DECOMPRESSED_BASENAME}_ACGT.interval_list --OUTPUT chunk_intervals"
gatk IntervalListTools \
  --SCATTER_COUNT 50 \
  --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
  --UNIQUE true \
  --SORT true \
  --BREAK_BANDS_AT_MULTIPLES_OF 0 \
  --INPUT ${GENOME_FASTA_DECOMPRESSED_BASENAME}_ACGT.interval_list \
  --OUTPUT chunk_intervals
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: gatk IntervalListTools"
else
    touch FAILED_GATK_IntervalListTools
    echo "FAILED: gatk IntervalListTools" >> README
    exit 1
fi
echo >> README

# Create template strings of contigs and calling intervals for Jetstream YAML document

# Calling contigs: - {"contig": chr1, "length": 248956422}
cut -f2,3 ${GENOME_FASTA_DECOMPRESSED_BASENAME}.dict \
  | \
  grep -v "VN" \
  | \
  sed 's/SN://g' \
  | \
  sed 's/LN://g' \
  | \
  awk '{OFS="" ; print "- {\"contig\": ", $1, ", \"length\": ", $2, "}"}' > ${GENOME_FASTA_DECOMPRESSED_BASENAME}_calling_contigs.txt

# Calling intervals: - {"contig": "chr1", "length": 197665, "start": 10001, "stop": 207666}
grep -w "ACGTmer" ${GENOME_FASTA_DECOMPRESSED_BASENAME}_ACGT.interval_list \
  | \
  cut -f1,2,3 \
  | \
  awk '{OFS="" ; print "- {\"contig\": \"", $1, "\", \"length\": ", $3-$2, ", \"start\": ", $2,  ", \"stop\": ", $3, "}"}' > ${GENOME_FASTA_DECOMPRESSED_BASENAME}_calling_intervals.txt
# Sequence Groupings: - ['chr1'] OR - ['chr18', 'chr19', 'chr20']
cut -f2,3 ${GENOME_FASTA_DECOMPRESSED_BASENAME}.dict \
  | \
  grep -v "VN" \
  | \
  sed 's/SN://g' \
  | \
  sed 's/LN://g' \
  | \
  awk '{OFS="" ; print "- [\x27", $1, "\x27]"}' > ${GENOME_FASTA_DECOMPRESSED_BASENAME}_sequence_groupings.txt

# Add flag to top level to indicate process is completeq
echo "## Add process completion flag to top level direcory" >> README
echo "    touch ${TOPLEVEL_DIR}/GENOME_FASTA_GENERATION_COMPLETE" >> README
touch ${TOPLEVEL_DIR}/GENOME_FASTA_GENERATION_COMPLETE

# If the fastas are the same for DNA and RNA, RNA FASTA GENERATION is also complete
if [[ "${REFERENCE_DNA_GENOME_NAME}" == "${REFERENCE_RNA_GENOME_NAME}" ]]
then
  touch ${TOPLEVEL_DIR}/RNA_FASTA_GENERATION_COMPLETE
fi
