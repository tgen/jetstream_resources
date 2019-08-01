#!/usr/bin/env bash

# Automated Script to download and build GnomAD release 2.1.1 VCF files for usage in Phoenix workflow

# Usage: ./build_gnomAD_r2.1.1.sh

# Write this full document as a README
cat $0 > README

# Added user information and timestamp to README
USER=`whoami`
DATE=`date`
echo "Downloaded and Processed by:  ${USER}" >> README
echo ${DATE} >> README

# Load required modules
module load samtools/1.9

############################
###
### Download and Processing Exome File
###
############################
echo
echo

# Download the current exome VCF
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
# 85.31 GiB, MD5: cff8d0cfed50adc9211d1feaed2d4ca7

# Check the md5sum
echo "cff8d0cfed50adc9211d1feaed2d4ca7  gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz" > gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.md5
md5sum --check gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.md5 | tee --append README
echo

## Using these massive VCF files we discovered there a significant benefits to using bcf files
bcftools view \
    --threads 8 \
    --output-type b \
    --output-file gnomad.exomes.r2.1.1.sites.liftover_grch38.bcf \
    gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
bcftools index --threads 4 gnomad.exomes.r2.1.1.sites.liftover_grch38.bcf

## Also using these files for annotation removing the massive amount of data in the INFO column significantly improves runtime
bcftools annotate \
    --threads 8 \
    --remove INFO \
    --output-type b \
    --output gnomad.exomes.r2.1.1.sites.liftover_grch38_NoINFO.bcf
    gnomad.exomes.r2.1.1.sites.liftover_grch38.bcf
bcftools index --threads 4 gnomad.exomes.r2.1.1.sites.liftover_grch38_NoINFO.bcf

rm gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
rm gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.md5

############################
###
### Download and Processing Genome File
###
############################
echo
echo

# Download the current genome VCF
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz
# 743.06 GiB, MD5: 83de3d5b52669f714e810d4fcf047c18

# Check the md5sum
echo "83de3d5b52669f714e810d4fcf047c18  gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz" > gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz.md5
md5sum --check gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz.md5 | tee --append README
echo

## Using these massive VCF files we discovered there a significant benefits to using bcf files
bcftools view \
    --threads 8 \
    --output-type b \
    --output-file gnomad.genomes.r2.1.1.sites.liftover_grch38.bcf \
    gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz
bcftools index --threads 4 gnomad.genomes.r2.1.1.sites.liftover_grch38.bcf

## Also using these files for annotation removing the massive amount of data in the INFO column significantly improves runtime
bcftools annotate \
    --threads 8 \
    --remove INFO \
    --output-type b \
    --output gnomad.genomes.r2.1.1.sites.liftover_grch38_NoINFO.bcf
    gnomad.genomes.r2.1.1.sites.liftover_grch38.bcf
bcftools index --threads 4 gnomad.genomes.r2.1.1.sites.liftover_grch38_NoINFO.bcf

rm gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz
rm gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz.md5


############################
###
### Download and Processing Notes
###
############################

