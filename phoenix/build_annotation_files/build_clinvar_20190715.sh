#!/usr/bin/env bash

# Automated Script to download and build CLINVAR VCF file for usage in Phoenix workflow

# Usage: ./build_clinvar_20190715.sh

# Write this full document as a README
cat $0 > README

# Added user information and timestamp to README
USER=`whoami`
DATE=`date`
echo "Downloaded and Processed by:  ${USER}" >> README
echo ${DATE} >> README

# Load required modules
module load samtools/1.9

# Define any needed variables
FAI=/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa.fai

####################################
##
## Parameterized Code
##
####################################

# Download the current clinvar files
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2019/clinvar_20190715.vcf.g*
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2019/clinvar_20190715_papu.vcf.g*

# Check MD5 checksums - All downloads pass
md5sum clinvar_20190715_papu.vcf.gz
cat clinvar_20190715_papu.vcf.gz.md5
md5sum clinvar_20190715.vcf.gz
cat clinvar_20190715.vcf.gz.md5

# Check contig seqeunces to determine if they match our hg38tgen reference genome with UCSC contif names
zcat clinvar_20190715.vcf.gz | cut -f1 | grep -v "#" | sort | uniq
# Found 1-22,X,Y,MT, and NW_009646201.1, which is a GRCh38.p13 patch for the ABO*A1.02 allele
zcat clinvar_20190715_papu.vcf.gz | cut -f1 | grep -v "#" | sort | uniq
# Found only Y

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
	--rename-chrs /home/jkeats/git_repositories/GRCh38_CrossMapping/GRCh38_Cosmic89_2_UCSC_Contigs.txt \
	--output-type b \
	--output temp_renamed.bcf \
	temp_droppedContigs.bcf
bcftools index --threads 4 temp_renamed.bcf

## Update the contig header (need to use an unreleased version with this new feature) so lengths are included
/home/jkeats/downloads/bcftools/bcftools reheader \
	--fai ${FAI} \
	temp_renamed.bcf \
	| \
	bcftools view \
	--output-type b \
	--output-file clinvar_20190715_hg38tgen.bcf
bcftools index --threads 4 clinvar_20190715_hg38tgen.bcf

# Remove all the temp files
rm temp_*

############################
###
### Download and Processing Notes
###
############################