#!/usr/bin/env bash

# Automated Script to download and build COSMIC coding and non-coding mutation VCF files for usage in Phoenix workflow

# Usage: ./build_cosmic_v89.sh

# Write this full document as a README
cat $0 > README

# Added user information and timestamp to README
USER=`whoami`
DATE=`date`
echo "Downloaded and Processed by:  ${USER}" >> README
echo ${DATE} >> README

# Modules needed
module load samtools/1.9

# Define FAI file
FAI=/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa.fai

# Generate authentication token
EMAIL=jkeats@tgen.org
#PASSWORD=NotMyRealPassword123
echo "Please enter your COSMIC website authentication password"
read PASSWORD

TOKEN=`echo "${EMAIL}:${PASSWORD}" | base64`

############################
###
### Download and process the coding mutation VCF file
###
############################

## DOWNLOAD Coding Mutations file
# Get download link, this generates a custom url and access key for the next step, returned as a JSON, which is parsed into the download key
RESPONSE=`curl -H "Authorization: Basic ${TOKEN}" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v89/VCF/CosmicCodingMuts.vcf.gz`
DOWNLOAD_KEY=`echo $RESPONSE | sed 's/"}//g' | sed 's/{"url":"https:\/\/cog.sanger.ac.uk\/cosmic\/GRCh38\/cosmic\/v89\/VCF\/CosmicCodingMuts.vcf.gz//g'`
# Download Cosmic Coding Mutations file (Authentication token expires after 1 hour)
curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v89/VCF/CosmicCodingMuts.vcf.gz${DOWNLOAD_KEY}" > CosmicCodingMuts_v89_GRCh38.vcf.gz

# The downloaded file is not bgzip so it needs to be decompressed with gunzip before processing
gunzip CosmicCodingMuts_v89_GRCh38.vcf.gz

# Check the contigs in the file, for matching the contigs need to match and coming from ensemble they will be 1 not chr1, etx..
echo
echo "Unique contig list in Cosmic Coding Mutations File:"
bcftools view -H CosmicCodingMuts_v89_GRCh38.vcf | cut -f1 | sort | uniq
echo
### As expected only contigs in the VCF are 1-22, X, Y, MT in both files

# Convert to block compressed and create .tbi index, seems needed, not sure why
bcftools view \
    --threads 4 \
    --output-type z \
    --output-file CosmicCodingMuts_v89_GRCh38.vcf.bgz \
    CosmicCodingMuts_v89_GRCh38.vcf
bcftools index --tbi CosmicCodingMuts_v89_GRCh38.vcf.bgz

# Update the contig names using bcftools
bcftools annotate \
	--threads 4 \
	--rename-chrs /home/jkeats/git_repositories/GRCh38_CrossMapping/GRCh38_Cosmic89_2_UCSC_Contigs.txt \
	--output-type b \
	--output temp_renamed.bcf \
	CosmicCodingMuts_v89_GRCh38.vcf.bgz
bcftools index --threads 4 temp_renamed.bcf

# Add updated contig header lines and add missing ones, required current bcftools version
/home/jkeats/downloads/bcftools/bcftools reheader \
	--fai ${FAI} \
	temp_renamed.bcf \
	| \
	bcftools view \
	--output-type b \
	--output-file CosmicCodingMuts_v89_hg38tgen.bcf
bcftools index --threads 4 CosmicCodingMuts_v89_hg38tgen.bcf

rm temp_*

############################
###
### Download and process the non-coding VCF file
###
############################

# Download Cosmic Non-Coding Mutations file (Authentication token expires after 1 hour)
RESPONSE=`curl -H "Authorization: Basic ${TOKEN}" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v89/VCF/CosmicNonCodingVariants.vcf.gz`
DOWNLOAD_KEY=`echo $RESPONSE | sed 's/"}//g' | sed 's/{"url":"https:\/\/cog.sanger.ac.uk\/cosmic\/GRCh38\/cosmic\/v89\/VCF\/CosmicNonCodingVariants.vcf.gz//g'`
curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v89/VCF/CosmicNonCodingVariants.vcf.gz${DOWNLOAD_KEY}" > CosmicNonCodingMuts_v89_GRCh38.vcf.gz

# The downloaded file is not bgzip so it needs to be decompressed with gunzip before processing
gunzip CosmicNonCodingMuts_v89_GRCh38.vcf.gz

# Check the contigs in the file, for matching the contigs need to match and coming from ensemble they will be 1 not chr1, etx..
echo
echo "Unique contig list in Cosmic Coding Mutations File:"
bcftools view -H CosmicNonCodingMuts_v89_GRCh38.vcf | cut -f1 | sort | uniq
echo
### As expected only contigs in the VCF are 1-22, X, Y, MT in both files

# Convert to block compressed and create .tbi index, seems needed, not sure why
bcftools view \
    --threads 4 \
    --output-type z \
    --output-file CosmicNonCodingMuts_v89_GRCh38.vcf.bgz \
    CosmicNonCodingMuts_v89_GRCh38.vcf
bcftools index --tbi CosmicNonCodingMuts_v89_GRCh38.vcf.bgz

# Update the contig names using bcftools
bcftools annotate \
	--threads 4 \
	--rename-chrs /home/jkeats/git_repositories/GRCh38_CrossMapping/GRCh38_Cosmic89_2_UCSC_Contigs.txt \
	--output-type b \
	--output temp_renamed.bcf \
	CosmicNonCodingMuts_v89_GRCh38.vcf.bgz
bcftools index --threads 4 temp_renamed.bcf

# Add updated contig header lines and add missing ones, required current bcftools version
/home/jkeats/downloads/bcftools/bcftools reheader \
	--fai ${FAI} \
	temp_renamed.bcf \
	| \
	bcftools view \
	--output-type b \
	--output-file CosmicNonCodingMuts_v89_hg38tgen.bcf
bcftools index --threads 4 CosmicNonCodingMuts_v89_hg38tgen.bcf

rm temp_*

############################
###
### Download and Processing Notes
###
############################

