#!/usr/bin/env bash

# Usage: create_gene_model.sh <Config.ini>

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -ue

# Read required variables from configuration file
. ${1}

####################################
## Download Gene Model File
####################################

# move to top level directory
mkdir -p ${TOPLEVEL_DIR}
cd ${TOPLEVEL_DIR}

# Make gene_model directory if not available
if [ -e gene_model ]
then
    echo "Gene Model directory exists, moving into it"
    cd gene_model
else
    echo "Gene Model directory NOT fount, creating and moving into it now"
    mkdir gene_model
    cd gene_model
fi

# Make ensembl_v95 directory if not available
if [ -e ${GENE_MODEL_NAME} ]
then
    echo "Specific Gene Model directory exists, moving into it"
    cd ${GENE_MODEL_NAME}
else
    echo "Specific Gene Model directory NOT fount, creating and moving into it now"
    mkdir ${GENE_MODEL_NAME}
    cd ${GENE_MODEL_NAME}
fi

# Initialize a gene_model specific README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "${GENE_MODEL_NAME} download and manipulated to align with reference genomes" >> README
echo >> README

# Download the gene model gtf
echo "Download the gene model gtf" >> README
wget ${GENE_MODEL_DOWNLOAD_LINK}
fc -ln -1 >> README
# Capture GTF file name
echo "Capture the downloaded GTF filename" >> README
GTF_FILE_GZ=`basename ${GENE_MODEL_DOWNLOAD_LINK}`
fc -ln -1 >> README
gunzip ${GTF_FILE_GZ}
fc -ln -1 >> README
echo >> README
# Capture the decompressed filename
echo "Capture the decompressed filename" >> README
GTF_FILE_FLAT=`basename ${GTF_FILE_GZ} ".gz"`
fc -ln -1 >> README
echo >> README
# Capture the decompressed filename basename
echo "Capture the decompressed filename basename" >> README
GTF_FILE_BASE=`basename ${GTF_FILE_FLAT} ".gtf"`
fc -ln -1 >> README
echo >> README

####################################
## Create updated gene model with ucsc style contig names - JK Method, generates results identical to RR and AC
####################################

## FILE ACCOUNTING
echo "Record the number of lines in original file" >> README
wc -l ${GTF_FILE_FLAT} >> README
fc -ln -1 >> README
echo >> README
# 2737564 Homo_sapiens.GRCh38.95.gtf
# Set variable for testing
INPUT_GTF_LINES=`cat ${GTF_FILE_FLAT} | wc -l`

echo "For processing the file needs to be split into 3 parts" >> README
echo "--- 1) header section" >> README
echo "--- 2) after header column 1, this will be used to update contig names" >> README
echo "--- 2) after header columns 2-end" >> README
echo >> README

# Create a header file
echo "Extract header line to temp file" >> README
grep "^#" ${GTF_FILE_FLAT} > temp_header
fc -ln -1 >> README
echo >> README

# Create first column file
echo "Extract column 1 without the header lines to temp file" >> README
grep -v "^#" ${GTF_FILE_FLAT} | cut -f1 > temp_c1.txt
fc -ln -1 >> README
echo >> README


# Create 2-final column file
echo "Extract column 2 and subsequent without the header lines to temp file" >> README
grep -v "^#" ${GTF_FILE_FLAT} | cut -f2- > temp_c2plus.txt
fc -ln -1 >> README
echo >> README


# Update column 1 file with new contig names
echo "Update column 1 file with new contig names" >> README
for line in `cat ${PATH_TO_REPO}/utility_files/ensembl95_ucsc_mappings_keats.csv`
do
    OLD=`echo ${line} | cut -d, -f1`
    NEW=`echo ${line} | cut -d, -f2`
    echo "Changing:  $OLD  to  $NEW" >> README

    # update in place
    sed -i "s/\b${OLD}/${NEW}/g" temp_c1.txt
done
fc -ln -1 >> README
echo >> README


#Create final file
echo "Create final updated GTF file with UCSC contig names" >> README
paste temp_c1.txt temp_c2plus.txt > temp_new.gtf
fc -ln -1 >> README
cat temp_header temp_new.gtf > ${GTF_FILE_BASE}.ucsc.gtf
fc -ln -1 >> README
echo >> README

# Check final GTF line count
echo "Check Final GTF file length" >> README
wc -l ${GTF_FILE_BASE}.ucsc.gtf >> README
fc -ln -1 >> README
echo >> README
# 2737564 Homo_sapiens.GRCh38.95.ucsc.gtf  ## Matches starting line count!!

# Confirm the starting and final GTF have the same line count
echo "Confirm the starting and final GTF have the same line count" >> README
FINAL_GTF_LINES=`cat ${GTF_FILE_BASE}.ucsc.gtf | wc -l`
fc -ln -1 >> README
echo >> README

if [ ${INPUT_GTF_LINES} -eq ${FINAL_GTF_LINES} ]
then
    echo "Input and Final GTF files match"
    echo "Input and Final GTF files match"  >> README
else
    echo
    echo
    echo "####### ERROR #########"
    echo "Input and Final file lengths DO NOT match"
    echo "Input = ${INPUT_GTF_LINES}"
    echo "Final = ${FINAL_GTF_LINES}"
    echo
    echo
    exit 1
fi
echo >> README

# Clean-up directory
rm temp_*
mkdir downloads
mv ${GTF_FILE_FLAT} downloads

# Create transcriptome fasta file derived from processed GTF
# This submission records the jobID so the next step does not start until this is complete
echo "Create transcriptome fasta file from the processed GTF for tools like Salmon" >> README
GTF_FASTA_JOBID=$(sbatch --parsable --export ALL,GENOME='${REFERENCE_RNA_GENOME_FASTA}',GTF='${GTF_FILE_BASE}.ucsc.gtf',OUTPUT='${GTF_FILE_BASE}.ucsc.transcriptome.fasta' ${PATH_TO_REPO}/utility_scripts/create_transcript_fasta.sh)
fc -ln -1 >> README
echo >> README
echo "Specific script code as follows:" >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/create_transcript_fasta.sh >> README
echo >> README
