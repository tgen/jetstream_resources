#!/usr/bin/env bash

# Usage: create_gene_model.sh <Config.ini>

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
    exit 2
fi

# Check that the reference genome was created successfully
if [ -e GENOME_FASTA_GENERATION_COMPLETE ]
then
    echo "Genome fasta exists, moving forward"
else
    echo "Genome fasta generation complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
fi

# Check gene_model directory if not available
if [ -e gene_model ]
then
    echo "Gene Model directory exists, moving into it"
    cd gene_model
else
    echo "Gene Model directory NOT found, creating and moving into it now"
    mkdir gene_model
    cd gene_model
fi

# Make specific gene model directory
if [ -e ${GENE_MODEL_NAME} ]
then
    echo "Specific Gene Model directory exists, exiting to prevent overwriting"
    exit 2
else
    echo "Specific Gene Model directory NOT found, creating it now and moving into it"
    mkdir ${GENE_MODEL_NAME}
    cd ${GENE_MODEL_NAME}
fi

####################################
## Download Gene Model File
####################################

# Initialize a gene_model specific README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "${GENE_MODEL_NAME} download and manipulated to align with reference genomes" >> README
echo >> README

# Download the gene model gtf
echo "Download the gene model gtf"
echo "Download the gene model gtf" >> README
echo "wget ${GENE_MODEL_DOWNLOAD_LINK}" >> README
wget --no-check-certificate ${GENE_MODEL_DOWNLOAD_LINK}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: download GTF"
else
    touch FAILED_DOWNLOAD_GTF
    echo "FAILED: download GTF" >> README
    exit 1
fi
echo >> README

# Determine GTF filename
echo "Capture the downloaded GTF filename"
echo "Capture the downloaded GTF filename" >> README
echo "GTF_FILE_GZ=`basename ${GENE_MODEL_DOWNLOAD_LINK}`" >> README
GTF_FILE_GZ=`basename ${GENE_MODEL_DOWNLOAD_LINK}`

# Determine if GTF checksum was provided and if so download it
if [ ${GENE_MODEL_MD5_DOWNLOAD_LINK} != "NA" ]
then
  echo "## Download GTF checksum from ${GENE_MODEL_SOURCE}" >> README
  echo "    wget ${GENE_MODEL_MD5_DOWNLOAD_LINK}" >> README
  wget --no-check-certificate ${GENE_MODEL_MD5_DOWNLOAD_LINK}
  # Error Capture
  if [ "$?" = "0" ]
  then
      echo "Completed: download GTF checksum"
  else
      touch FAILED_DOWNLOAD_GTF_CHECKSUM
      echo "FAILED: download GTF checksum" >> README
      exit 1
  fi
  echo >> README
fi


# Check MD5SUM
if [ ${GENE_MODEL_SOURCE} == "ensembl" ]
then
  echo "ENSEMBL is supported"
  # Ensembl now uses "sum" for check sum validation
  # Extract the provided checksum and number of 512bit blocks
  PROVIDED_CHECKSUM=`grep ${GTF_FILE_GZ} CHECKSUMS | cut -d" " -f1`
  PROVIDED_512bitBLOCKS=`grep ${GTF_FILE_GZ} CHECKSUMS | cut -d" " -f2`
  # Calculate the checksum of the downlaoded file
  VALIDATION_SUM=`sum ${GTF_FILE_GZ}`
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
elif [ ${GENOME_SOURCE} == "OTHER" ]
then
  echo "OTHER is supported"
  echo "WARNING - No checksum provided the file CANNOT be validated"
  exit 2
else
  echo "Current Genome Source is NOT SUPPORTED"
  exit 1
fi

# Determine expected decompressed GTF filename
echo "Capture the expected decompressed GTF filename"
echo "Capture the expected decompressed GTF filename" >> README
echo "GTF_FILE=`basename ${GENE_MODEL_DOWNLOAD_LINK} ".gz"`" >> README
GTF_FILE=`basename ${GENE_MODEL_DOWNLOAD_LINK} ".gz"`

# Decompress GTF
echo "Decompressing GTF"
echo "Decompressing GTF" >> README
echo "gunzip -c ${GTF_FILE_GZ} > ${GTF_FILE}" >> README
gunzip -c ${GTF_FILE_GZ} > ${GTF_FILE}
# Error Capture
if [ "$?" = "0" ]
then
    echo "Completed: gunzip GTF"
else
    touch FAILED_GUNZIP_GTF
    echo "FAILED: gunzip GTF" >> README
    exit 1
fi
echo >> README

# Capture the decompressed filename basename
echo "Capture the decompressed filename basename"
echo "Capture the decompressed filename basename" >> README
echo "GTF_FILE_BASE=`basename ${GTF_FILE} ".gtf"`" >> README
GTF_FILE_BASE=`basename ${GTF_FILE} ".gtf"`
echo >> README

####################################
## Determine required variables
####################################

# Determine the reference genome fasta full path
echo "Determine the full path filename of the reference genome fasta" >> README
echo "REFERENCE_GENOME_FASTA=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}" >> README
REFERENCE_GENOME_FASTA=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}
echo >> README

# Determine the full path to the reference fasta.fa.fai file
echo "Determine full path filename of the reference genome fasta.fa.fai index file" >> README
echo "REFERENCE_GENOME_FAI=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}.fai" >> README
REFERENCE_GENOME_FAI=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}.fai
echo >> README

# Determine the full path to the reference dict file
echo "Determine the basename of the reference genome fasta to infer the .dict filename" >> README
echo "REFERENCE_GENOME_DICT_BASENAME=${REFERENCE_DNA_GENOME_NAME%.fa}" >> README
REFERENCE_GENOME_DICT_BASENAME=${REFERENCE_DNA_GENOME_NAME%.fa}
echo >> README

####################################
## Check to ensure the GTF contigs and Genome Contigs Match
####################################

# Get a list of the unique contigs in the GTF
cut -f1 ${GTF_FILE} | grep -v "#" | sort | uniq > temp_gtf_contigs.txt

# Determine if any GTF contig DOES NOT exist in the reference genome (genome.fa.fai is a unique contig list)
echo "Checking that all GTF contigs exist in reference genome"
for contig in `cat temp_gtf_contigs.txt`
do
  CONTIG_CHECK=`awk -v var1="$contig" '{if($1 == var1) print $0}' ${REFERENCE_GENOME_FAI} | wc -l`
  if [ ${CONTIG_CHECK} == 1 ]
  then
    echo "GTF Contig: ${contig} exists in reference genome fasta"
  elif [ ${CONTIG_CHECK} == 0 ]
  then
    echo
    echo "ERROR GTF Contig: ${contig} DOES NOT exist in reference genome fasta"
    echo
    exit 1
  else
    echo "WARNING - UNEXPECTED EVENT"
    echo "Contig count in reference genome is not 0 or 1 as expected"
    echo "Contig ${contig} was seen ${CONTIG_CHECK} times in the reference.fa.fai"
    exit 2
  fi
done

####################################
## Create required files produced from the GTF
####################################
# Create a bed file for the start and stop for each gene
awk -F '[\t"]' '$1 !~ /^#/ { if (a[$10] == "" ) { a[$10] = $1 ; b[$10] = $4 ; c[$10] = $5 ; next } ;
        if ($4 < b[$10]) { b[$10] = $4 } ;
        if ($5 > c[$10]) { c[$10] = $5 }
} END {
for (i in a) {
        OFS = "\t" ; print a[i], b[i], c[i], i
}
}' ${GTF_FILE} | sort -k1,1V -k2,2n -k3,3n > ${GTF_FILE_BASE}.gene.bed

# Create a bed file for the start and stop of each exon for each gene
awk -F '[\t"]' '$1 !~ /^#/ { if ($3 == "exon") { OFS = "\t" ; print $1, $4, $5, $10 }}' ${GTF_FILE} | sort -k1,1V -k2,2n -k3,3n > ${GTF_FILE_BASE}.exon.bed

# Create reflat file from GTF for Picard RNAseqMetrics
# Uses gtfToGenePred from UCSC
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
echo "Creating RefFlat for picrard rnaSeqMetrics"
echo "Create refflat file from GTF for picard rnaSeqMetrics" >> README
#echo ${GTFTOGENEPRED_BINARY} -genePredExt -ignoreGroupsWithoutExons ${GTF_FILE} /dev/stdout | awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${GTF_FILE_BASE}.refFlat.txt >> README
${GTFTOGENEPRED_BINARY} -genePredExt \
    -ignoreGroupsWithoutExons \
    ${GTF_FILE} \
    /dev/stdout \
    | \
    awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${GTF_FILE_BASE}.refFlat.txt
fc -ln -1 >> README
echo >> README

# Extract the ribosomal RNA locations and create file for usage with Picard RNAseqMetrics
echo "Creating ribosomal RNA locations file for Picar rnaSeqMetrics"
echo "Create ribosomal interval file from GTF for picard rnaSeqMetrics" >> README
#echo "grep -w "rRNA" ${GTF_FILE} | cut -f1,4,5,7,9 | sed 's/gene_id "//g' | sed 's/"; transcript_id "/\'$'\t''/g' | cut -f1-5 > temp_RibosomalLocations.txt" >> README
grep -w "rRNA" ${GTF_FILE} \
    | \
    cut -f1,4,5,7,9 \
    | \
    sed 's/gene_id "//g' \
    | \
    sed 's/"; transcript_id "/\'$'\t''/g' \
    | \
    cut -f1-5 > temp_RibosomalLocations.txt
fc -ln -1 >> README
echo >> README

echo "Created final ribosome interval list file" >> README
echo "cat ${TOPLEVEL_DIR}/genome_reference/${REFERENCE_GENOME_DICT_BASENAME}.dict temp_RibosomalLocations.txt > ${GTF_FILE_BASE}.ribo.interval_list" >> README
cat ${TOPLEVEL_DIR}/genome_reference/${REFERENCE_GENOME_DICT_BASENAME}.dict temp_RibosomalLocations.txt > ${GTF_FILE_BASE}.ribo.interval_list
echo >> README

# Create refflat file from GTF fro IGV gene model track with HUGO IDs
echo "Creating IGV RefFlat File"
echo "Create refflat file from GTF for IGV genemodel tracks with HUGO IDs" >> README
#echo "${GTFTOGENEPRED_BINARY} -genePredExt -ignoreGroupsWithoutExons -geneNameAsName2 ${GTF_FILE} /dev/stdout | awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${GTF_FILE_BASE}.refFlat.hugoID.txt" >> README
${GTFTOGENEPRED_BINARY} -genePredExt \
    -ignoreGroupsWithoutExons \
    -geneNameAsName2 \
    ${GTF_FILE} \
    /dev/stdout \
    | \
    awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${GTF_FILE_BASE}.refFlat.hugoID.txt
fc -ln -1 >> README
echo >> README

# Create transcriptome fasta file derived from processed GTF
echo "Submitting transcriptome fasta generation process"
echo "Create transcriptome fasta file from the processed GTF for tools like Salmon" >> README
#echo "sbatch --parsable --export ALL,GENOME="${REFERENCE_GENOME_FASTA}",GTF="${GTF_FILE}",OUTPUT="${GTF_FILE_BASE}.transcriptome.fasta" ${PATH_TO_REPO}/utility_scripts/create_transcript_fasta.sh" >> README
sbatch --parsable --export ALL,GENOME="${REFERENCE_GENOME_FASTA}",GTF="${GTF_FILE}",OUTPUT="${GTF_FILE_BASE}.transcriptome.fasta" ${PATH_TO_REPO}/utility_scripts/create_transcript_fasta.sh
fc -ln -1 >> README
echo >> README
echo >> README
echo "----------------------------------------------------------------------------">> README
echo "Specific script code as follows:" >> README
echo >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/create_transcript_fasta.sh >> README
echo >> README

# Indicate GTF was created successfully
touch GENE_MODEL_GTF_GENERATION_COMPLETE

# Clean-up directory
rm temp_*

# Indicate completed
echo "All Processes Completed"
