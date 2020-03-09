#!/usr/bin/env bash

# Automated Script to download and build GnomAD release 2.1.1 VCF files for usage in Phoenix workflow

# Usage: ./build_gnomAD_r2.1.1.sh <config.ini>

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

# Load required modules
module load samtools/1.9
module load gatk/4.1.3.0

####################################
## Create Expected Folder Structure
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

# Make public_databases folder if not available
if [ -e public_databases ]
then
    echo "Public Databases folder exists, moving into it"
    cd public_databases
else
    echo "Public Databases folder NOT found, creating and moving into it now"
    mkdir -p public_databases
    cd public_databases
fi

# Make gnomad folder if not available
if [ -e gnomad ]
then
    echo "gnomAD folder exists, moving into it"
    cd gnomad
else
    echo "gnomAD folder NOT found, creating and moving into it now"
    mkdir -p gnomad
    cd gnomad
fi

# Make gnomad release version folder if not available
if [ -e ${GNOMAD_RELEASE_VERSION} ]
then
    echo "gnomAD ${GNOMAD_RELEASE_VERSION} folder exists, exiting to prevent overwrite"
    exit 1
else
    echo "gnomAD ${GNOMAD_RELEASE_VERSION} folder NOT found, creating and moving into it now"
    mkdir -p ${GNOMAD_RELEASE_VERSION}
    cd ${GNOMAD_RELEASE_VERSION}
fi

####################################
## Download and Process Files as Needed
###################################

# Initialize a gnomAD specific README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

############################
###
### Download and Processing Exome File
###
############################
echo
echo

# Download the current exome VCF
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz


# Check the md5sum
echo "f034173bf6e57fbb5e8ce680e95134f2  gnomad.exomes.r2.1.1.sites.vcf.bgz" > gnomad.exomes.r2.1.1.sites.vcf.bgz.md5
md5sum --check gnomad.exomes.r2.1.1.sites.vcf.bgz.md5 | tee --append README
echo

bcftools index --tbi --threads 4 gnomad.exomes.r2.1.1.sites.vcf.bgz

## Using these massive VCF files we discovered there a significant benefits to using bcf files
bcftools view \
    --threads 8 \
    --output-type b \
    --output-file gnomad.exomes.r2.1.1.sites.bcf \
    gnomad.exomes.r2.1.1.sites.vcf.bgz
bcftools index --threads 4 gnomad.exomes.r2.1.1.sites.bcf

bcftools stats --threads 8 gnomad.exomes.r2.1.1.sites.bcf > gnomad.exomes.r2.1.1.sites.bcf.stats
#plot-vcfstats --no-PDF --title "GnomAD Exomes" -p plots_vcfstats_exomes gnomad.exomes.r2.1.1.sites.bcf.stats

## Prep for Mutect2 steps

# Make Allele Frequency Only VCF
# Clear ID and QUAL and delete all INFO fields other than AF
bcftools annotate \
    --threads 8 \
    --remove ID,QUAL,^INFO/AF \
    --output-type b \
    --output temp_ForMutect.bcf \
    gnomad.exomes.r2.1.1.sites.bcf
# Create TBI index
bcftools index --threads 4 temp_ForMutect.bcf

bcftools view \
    --threads 8 \
    --apply-filters PASS \
    --output-type z \
    --output-file gnomad.exomes.r2.1.1.sites_ForMutect.vcf.gz \
    temp_ForMutect.bcf
# Create TBI index
bcftools index --threads 4 --tbi gnomad.exomes.r2.1.1.sites_ForMutect.vcf.gz


# Select Common Bi-Allelic SNPs got Mutects Contamination Testing
bcftools view \
    --threads 8 \
    --apply-filters PASS \
    --include 'AF>=0.05' \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --output-type z \
    --output-file gnomad.exomes.r2.1.1.sites_ForMutectContamination.vcf.gz \
    gnomad.exomes.r2.1.1.sites_ForMutect.vcf.gz
# Create TBI index
bcftools index --threads 4 --tbi gnomad.exomes.r2.1.1.sites_ForMutectContamination.vcf.gz

# clean up
rm temp*

############################
###
### Download and Processing Genome File
###
############################
echo
echo

# Download the current genome VCF
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz


# Check the md5sum
echo "e6eadf5ac7b2821b40f350da6e1279a2  gnomad.genomes.r2.1.1.sites.vcf.bgz" > gnomad.genomes.r2.1.1.sites.vcf.bgz.md5
md5sum --check gnomad.genomes.r2.1.1.sites.vcf.bgz.md5 | tee --append README
echo

bcftools index --tbi --threads 4 gnomad.genomes.r2.1.1.sites.vcf.bgz

## Using these massive VCF files we discovered there a significant benefits to using bcf files
bcftools view \
    --threads 8 \
    --output-type b \
    --output-file gnomad.genomes.r2.1.1.sites.bcf \
    gnomad.genomes.r2.1.1.sites.vcf.bgz
bcftools index --threads 4 gnomad.genomes.r2.1.1.sites.bcf


bcftools stats --threads 8 gnomad.genomes.r2.1.1.sites.bcf > gnomad.genomes.r2.1.1.sites.bcf.stats
#plot-vcfstats --no-PDF --title "GnomAD Genomes" -p plots_vcfstats_genomes gnomad.genomes.r2.1.1.sites.bcf.stats


## Prep for Mutect2 steps

# Make Allele Frequency Only VCF
# Clear ID and QUAL and delete all INFO fields other than AF
bcftools annotate \
    --threads 8 \
    --remove ID,QUAL,^INFO/AF \
    --output-type b \
    --output temp_ForMutect.bcf \
    gnomad.genomes.r2.1.1.sites.bcf
# Create TBI index
bcftools index --threads 4 temp_ForMutect.bcf

bcftools view \
    --threads 8 \
    --apply-filters PASS \
    --output-type z \
    --output-file gnomad.genomes.r2.1.1.sites_ForMutect.vcf.gz \
    temp_ForMutect.bcf
# Create TBI index
bcftools index --threads 4 --tbi gnomad.genomes.r2.1.1.sites_ForMutect.vcf.gz


# Select Common Bi-Allelic SNPs got Mutects Contamination Testing
bcftools view \
    --threads 8 \
    --apply-filters PASS \
    --include 'AF>=0.05' \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --output-type z \
    --output-file gnomad.genomes.r2.1.1.sites_ForMutectContamination.vcf.gz \
    gnomad.genomes.r2.1.1.sites_ForMutect.vcf.gz
# Create TBI index
bcftools index --threads 4 --tbi gnomad.genomes.r2.1.1.sites_ForMutectContamination.vcf.gz

# clean up
rm temp*

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

# Write this full document as a README
cat $0 > README