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
module load gatk/4.1.3.0

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
    --output gnomad.exomes.r2.1.1.sites.liftover_grch38_NoINFO.bcf \
    gnomad.exomes.r2.1.1.sites.liftover_grch38.bcf
bcftools index --threads 4 gnomad.exomes.r2.1.1.sites.liftover_grch38_NoINFO.bcf

bcftools stats --threads 8 gnomad.exomes.r2.1.1.sites.liftover_grch38.bcf > gnomad.exomes.r2.1.1.sites.liftover_grch38.bcf.stats
plot-vcfstats --no-PDF --title "GnomAD Exomes" -p plots_vcfstats_exomes gnomad.exomes.r2.1.1.sites.liftover_grch38.bcf.stats

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
    --output gnomad.genomes.r2.1.1.sites.liftover_grch38_NoINFO.bcf \
    gnomad.genomes.r2.1.1.sites.liftover_grch38.bcf
bcftools index --threads 4 gnomad.genomes.r2.1.1.sites.liftover_grch38_NoINFO.bcf

bcftools stats --threads 8 gnomad.genomes.r2.1.1.sites.liftover_grch38.bcf > gnomad.genomes.r2.1.1.sites.liftover_grch38.bcf.stats
plot-vcfstats --no-PDF --title "GnomAD Genomes" -p plots_vcfstats_genomes gnomad.genomes.r2.1.1.sites.liftover_grch38.bcf.stats

rm gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz
rm gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz.md5

############################
###
### Create Genome and Exome versions for MuTect Input
###
############################

# Make Allele Frequency Only VCF
# Clear ID and QUAL and delete all INFO fields other than AF
bcftools annotate \
    --threads 8 \
    --remove ID,QUAL,^INFO/AF \
    --output-type z \
    --output gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz \
    gnomad.exomes.r2.1.1.sites.liftover_grch38.bcf
# Create TBI index
bcftools index --threads 4 --tbi gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz

bcftools view \
    --threads 8 \
    --apply-filters PASS \
    --output-type z \
    --output-file gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutect2.vcf.gz \
    gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz

mv gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutect2.vcf.gz gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz
rm gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz.tbi
bcftools index --threads 4 --tbi gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz


# Select Common Bi-Allelic SNPs got Mutects Contamination Testing
bcftools view \
    --threads 8 \
    --apply-filters PASS \
    --include 'AF>=0.05' \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --output-type z \
    --output-file gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutectContamination.vcf.gz \
    gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz
# Create TBI index
bcftools index --threads 4 --tbi gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutectContamination.vcf.gz


# Make Allele Frequency Only VCF
# Clear ID and QUAL and delete all INFO fields other than AF
bcftools annotate \
    --threads 8 \
    --remove ID,QUAL,^INFO/AF \
    --output-type z \
    --output gnomad.genomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz \
    gnomad.genomes.r2.1.1.sites.liftover_grch38.bcf
# Create TBI index
bcftools index --threads 4 --tbi gnomad.genomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz

bcftools view \
    --threads 8 \
    --apply-filters PASS \
    --output-type z \
    --output-file gnomad.genomes.r2.1.1.sites.liftover_grch38_ForMutect2.vcf.gz \
    gnomad.genomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz

mv gnomad.genomes.r2.1.1.sites.liftover_grch38_ForMutect2.vcf.gz gnomad.genomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz

rm gnomad.genomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz.tbi
bcftools index --threads 4 --tbi gnomad.genomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz

# Select Common Bi-Allelic SNPs got Mutects Contamination Testing
bcftools view \
    --threads 8 \
    --apply-filters PASS \
    --include 'AF>=0.05' \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --output-type z \
    --output-file gnomad.genomes.r2.1.1.sites.liftover_grch38_ForMutectContamination.vcf.gz \
    gnomad.genomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz
# Create TBI index
bcftools index --threads 4 --tbi gnomad.genomes.r2.1.1.sites.liftover_grch38_ForMutectContamination.vcf.gz

## BROAD VERSION OF THE SAME PROCESSES (https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect_resources.wdl)
# Clear ID and QUAL and delete all INFO fields other than AF
# Save off the header for later:
# grep '^#' ${input_vcf} > header &
# Get all lines in the file except the header:
# Preserve all fields before INFO, Grab only the AF annotation from the INFO Field
# replace ID (3rd) and QUAL (6th) columns with '.' (empty):
# grep -v "^#" ${input_vcf} | sed -e 's#\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t.*;AF=\([0-9]*\.[e0-9+-]*\).*#\1\t\2\t.\t\4\t\5\t.\t\7\tAF=\8#g' > simplified_body &
# Wait for background processes to finish:
# wait
# Consolidate files:
# cat header simplified_body > simplified.vcf
# Zip the VCF:
# bgzip simplified.vcf -O ${output_name}.vcf.gz
# Index output file:
# gatk --java-options "-Xmx${command_mem}g" IndexFeatureFile -F ${output_name}.vcf.gz
# Cleanup:
# rm -f header body simplified_info simplified_body simplified.vcf simplified.vcf.idx
# Select Common Bi-Allelic SNPs
# gatk --java-options "-Xmx${command_mem}g" SelectVariants \
#            -V ${input_vcf} \
#            -select-type SNP -restrict-alleles-to BIALLELIC \
#            -select "AF > ${minimum_allele_frequency}" \
#            -O ${output_name}.vcf.gz \
#            --lenient



############################
###
### Download and Processing Notes
###
############################

