#!/usr/bin/env bash

# Usage: create_genderCheck_SNP_list.sh <Config.ini>

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -ue

# Read required variables from configuration file
. ${1}

# Load required modules
module load samtools/1.9
module load BEDTools/2.26.0

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

#ToDo: Fix needed checks like for dbSNP VCF and exon filter file
## Needed inputs
# 1) dbSNP list
# 2) exon bed file

# Check that the gene model GTF exists
if [ -e ${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/GENE_MODEL_GTF_GENERATION_COMPLETE ]
then
    echo "Gene model exists, moving forward"
else
    echo "Gene model complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
fi

####################################
## Generate Gender Check SNP list
####################################

if [ -e tool_resources ]
then
    echo "tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "tool_resources directory NOT found, creating and moving into it now"
    mkdir tool_resources
    cd tool_resources
fi

if [ -e "tgen_gender_check" ]
then
    echo "The tgen_gender_check directory exists, exiting to prevent overwriting existing index"
    exit 2
else
    echo "The tgen_gender_check directory was NOT found, creating and moving into it now"
    mkdir tgen_gender_check
    cd tgen_gender_check
fi

# Initialize a tgen gender check README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/phoenix" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "TGen Gender Check Resources creation details:" >> README
echo >> README

# Filter the dbSNP vcf to get a set of variants on chrX for gender check process
echo "Extract just the common SNV on chrX from the dbSNP vcf" >> README
echo "Remove any variants in the pseudoautosomal regions" >> README
bcftools filter \
    --threads 8 \
    --regions chrX \
    --targets ^chrX:10001-2781479,chrX:155701383-156030895 \
    --include 'INFO/VC == "SNV" & INFO/COMMON == 1' \
    --output-type b \
    --output temp_chrx_common_snv.bcf \
    /home/tgenref/homo_sapiens/grch38_hg38/public_databases/dbsnp/b152/dbSNP_b152_hg38tgen.bcf
fc -ln -1 >> README
echo >> README

# need to filter to snps --types snps \

# Create a list of exon coordinates from the GTF file, and correct start to bed format
echo "Create an BED file with the exons in GTF" >> README
grep "exon" ${TOPLEVEL_DIR}/gene_model/${GENE_MODEL_NAME}/Homo_sapiens.GRCh38.97.ucsc.gtf | cut -f1,4,5 | awk '{OFS="\t" ; print $1, $2-1, $3}' > temp_exons_all.bed
fc -ln -1 >> README
echo >> README

# Many exons overlap so use bedtools to collapse the overlaps
echo "Sort the BED file to prepare for collapsing" >> README
bedtools sort -i temp_exons_all.bed > temp_exons_all_sorted.bed
fc -ln -1 >> README
echo >> README

echo "Collapse the BED file exons to removed overlaps and prevent selection of the same SNP multiple times" >> README
bedtools merge -i temp_exons_all_sorted.bed > temp_exons_collapsed.bed
fc -ln -1 >> README
echo >> README

# Filter the list to known exon coordiantes so the list is useful for genomes, exomes, and likely even RNAseq
echo "Filter the chrX snps to just those in the exon regions but gender determination in genomes and exomes" >> README
bcftools filter \
    --targets-file temp_exons_collapsed.bed \
    --output-type z \
    --output chrx_common_dbSNPv152_snv_exons.vcf.gz \
    temp_chrx_common_snv.bcf
fc -ln -1 >> README
echo >> README

# The vcf needs to be tbi indexed to work with freebayes
echo "Create an index file for usage with freebayes" >> README
bcftools index \
    --tbi \
    chrx_common_dbSNPv152_snv_exons.vcf.gz
fc -ln -1 >> README
echo >> README

# Create bed file for usage with freebayes to run genotyping
echo "Create a BED file of the positions in question for usage with freebayes" >> README
bcftools view -H chrx_common_dbSNPv152_snv_exons.vcf.gz | cut -f1,2 | awk '{OFS="\t" ; print $1, $2-1, $2}' > chrx_common_dbSNPv152_snv_exons.bed
fc -ln -1 >> README
echo >> README

rm temp_*

echo >> README
echo >> README