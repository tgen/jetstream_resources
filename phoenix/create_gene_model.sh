#!/usr/bin/env bash

# Usage: create_gene_model.sh <Config.ini>

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -ue

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

# Check that the reference genome for RNA was created successfully
if [ -e RNA_FASTA_GENERATION_COMPLETE ]
then
    echo "RNA fasta exists, moving forward"
else
    echo "RNA fasta generation complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
fi

# Check gene_model directory if not available
if [ -e gene_model ]
then
    echo "Gene Model directory exists, moving into it"
    cd gene_model
else
    echo "Gene Model directory NOT fount, creating and moving into it now"
    mkdir gene_model
    cd gene_model
fi

# Make specific gene model directory
if [ -e ${GENE_MODEL_NAME} ]
then
    echo "Specific Gene Model directory exists, exiting to prevent overwriting"
    exit 2
else
    echo "Specific Gene Model directory NOT fount, creating it now and moving into it"
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
echo "https://github.com/tgen/jetstream_resources/phoenix" >> README
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

# Check that all contigs exist in renaming key
cut -f1 ${GTF_FILE_FLAT} | grep -v "#" | sort | uniq > temp_gtf_unique_contig_list
GTF_CONTIG_NUMBER=`wc -l temp_gtf_unique_contig_list | awk '{print $1}'`
RENAME_CONTIG_NUMBER=`wc -l ${PATH_TO_REPO}/utility_files/ensembl97_ucsc_mappings_keats.csv | awk '{print $1}'`
# Test if the rename key and GTF have the same number of contigs
if [ ${GTF_CONTIG_NUMBER} -eq ${RENAME_CONTIG_NUMBER} ]
then
    echo "Number of unique contigs lines match in GTF and rename key"
else
    echo
    echo "WARNING - number of contigs DOES NOT match between GTF and rename key"
    echo "Exiting, this MUST BE FIXED"
    echo
    exit 2
fi
#Just because the number of contigs matches doesn't mean they are the same ones
#Make a list of unique rename original ensembl contig names
cut -d"," -f1 ${PATH_TO_REPO}/utility_files/ensembl97_ucsc_mappings_keats.csv > temp_rename_key_contig_list
#Merge the two lists and test if they all match using unique count as they should all be 2
MATCHING_CONTIG_NUMBER=`cat temp_gtf_unique_contig_list temp_rename_key_contig_list | sort | uniq -c | awk '{print $1}' | grep "2" | wc -l | awk '{print $1}'`
#Test if the number of matching contigs is correct
if [ ${GTF_CONTIG_NUMBER} -eq ${MATCHING_CONTIG_NUMBER} ]
then
    echo "All contigs names match in the GTF and rename key"
else
    echo
    echo "WARNING - number of unique contig names DOES NOT match between GTF and rename key"
    echo "Exiting, this MUST BE FIXED"
    echo
    exit 2
fi

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
for line in `cat ${PATH_TO_REPO}/utility_files/ensembl97_ucsc_mappings_keats.csv`
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

# Create a bed file for the start and stop for each gene
awk -F '[\t"]' '$1 !~ /^#/ { if (a[$10] == "" ) { a[$10] = $1 ; b[$10] = $4 ; c[$10] = $5 ; next } ;
        if ($4 < b[$10]) { b[$10] = $4 } ;
        if ($5 > c[$10]) { c[$10] = $5 }
} END {
for (i in a) {
        OFS = "\t" ; print a[i], b[i], c[i], i 
}
}' ${GTF_FILE_BASE}.ucsc.gtf | sort -k1,1V -k2,2n -k3,3n > ${GTF_FILE_BASE}.ucsc.gene.bed

# Create a bed file for the start and stop of each exon for each gene
awk -F '[\t"]' '$1 !~ /^#/ { if ($3 == "exon") { OFS = "\t" ; print $1, $4, $5, $10 }}' ${GTF_FILE_BASE}.ucsc.gtf | sort -k1,1V -k2,2n -k3,3n > ${GTF_FILE_BASE}.ucsc.exon.bed 

# Create reflat file from GTF for Picard RNAseqMetrics
# Uses gtfToGenePred from UCSC
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
echo "Create refflat file from GTF for picard rnaSeqMetrics" >> README
${GTFTOGENEPRED_BINARY} -genePredExt \
    -ignoreGroupsWithoutExons \
    ${GTF_FILE_BASE}.ucsc.gtf \
    /dev/stdout \
    | \
    awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${GTF_FILE_BASE}.ucsc.refFlat.txt
fc -ln -1 >> README
echo >> README

# Extract the ribosomal RNA locations and create file for usage with Picard RNAseqMetrics
echo "Create ribosomal interval file from GTF for picard rnaSeqMetrics" >> README
grep -w "rRNA" ${GTF_FILE_BASE}.ucsc.gtf \
    | \
    cut -f1,4,5,7,9 \
    | \
    sed 's/gene_id "//g' \
    | \
    sed 's/"; transcript_id "/\'$'\t''/g' \
    | \
    cut -f1-5 > temp_RibosomalLocations.txt
fc -ln -1 >> README
cat ${REFERENCE_RNA_GENOME_DICT} temp_RibosomalLocations.txt > ${GTF_FILE_BASE}.ucsc.ribo.interval_list
fc -ln -1 >> README
echo >> README

echo "Create refflat file from GTF for IGV genemodel tracks with HUGO IDs" >> README
${GTFTOGENEPRED_BINARY} -genePredExt \
    -ignoreGroupsWithoutExons \
    -geneNameAsName2 \
    ${GTF_FILE_BASE}.ucsc.gtf \
    /dev/stdout \
    | \
    awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${GTF_FILE_BASE}.ucsc.refFlat.hugoID.txt
fc -ln -1 >> README
echo >> README


# Clean-up directory
rm temp_*
mkdir downloads
mv ${GTF_FILE_FLAT} downloads

# Create transcriptome fasta file derived from processed GTF
# This submission records the jobID so the next step does not start until this is complete
echo "Create transcriptome fasta file from the processed GTF for tools like Salmon" >> README
sbatch --parsable --export ALL,GENOME="${REFERENCE_RNA_GENOME_FASTA}",GTF="${GTF_FILE_BASE}.ucsc.gtf",OUTPUT="${GTF_FILE_BASE}.ucsc.transcriptome.fasta" ${PATH_TO_REPO}/utility_scripts/create_transcript_fasta.sh
fc -ln -1 >> README
echo >> README
echo "Specific script code as follows:" >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/create_transcript_fasta.sh >> README
echo >> README

# Indicate GTF was created successfully
touch GENE_MODEL_GTF_GENERATION_COMPLETE









