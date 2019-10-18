#!/usr/bin/env bash

# Automated Script to download and build GnomAD release 3.0 genome VCF file for usage in Phoenix workflow

# Usage: ./build_gnomAD_r3.0.sh <resources.ini>
# Sbatch Usage: sbatch --nodes=1 --cpus-per-task=10 --time=0-72:00:00 build_gnomAD_r3.0.sh <resources.ini>

# Check to resources.ini was provided on the command line
if [ -n "$1" ]; then
  echo "Required ini file detected"
else
  echo "Input INI file not provided, exiting due to missing requirement"
  exit 1
fi

# Read required variables from configuration file
. ${1}

# navigate to gnomad directory
cd ${PARENT_DIR}/public_databases/gnomad

# make version specific directory
mkdir -p r3.0

# Move into the version specific directory
cd r3.0

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
### Download and Processing Genome File
###
############################
echo
echo

# Download the current exome VCF
wget https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz
#235.68 GiB, MD5: f3501102192975da34b5d2c32f7c0791

# Check the md5sum
echo "f3501102192975da34b5d2c32f7c0791  gnomad.genomes.r3.0.sites.vcf.bgz" > gnomad.genomes.r3.0.sites.vcf.bgz.md5
md5sum --check gnomad.genomes.r3.0.sites.vcf.bgz.md5 | tee --append README
echo

## Using these massive VCF files we discovered there a significant benefits to using bcf files
## Also for internal usage we limit the VCF to just the PASS variants
bcftools view \
    --threads 8 \
    --apply-filters PASS \
    --output-type b \
    --output-file gnomad.genomes.r3.0.sites.pass.bcf \
    gnomad.genomes.r3.0.sites.vcf.bgz
bcftools index --threads 4 gnomad.genomes.r3.0.sites.pass.bcf

bcftools stats --threads 8 gnomad.genomes.r3.0.sites.pass.bcf > gnomad.genomes.r3.0.sites.pass.bcf.stats
plot-vcfstats --no-PDF --title "GnomADr3.0 Genomes" -p plots_vcfstats_genomes gnomad.genomes.r3.0.sites.pass.bcf.stats

## Also using these files for annotation removing the massive amount of data in the INFO column significantly improves runtime
## For utility with mutect2 make simple version with minimal data for annotation yes/no flags and a germline reference
bcftools annotate \
    --threads 8 \
    --remove ID,QUAL,^INFO/AF \
    --output-type z \
    --output gnomad.genomes.r3.0.sites.pass.AnnotationReference.bcf \
    gnomad.genomes.r3.0.sites.pass.bcf
# Create CSI index (This will be file used for bcftools database annotation purposes
bcftools index --threads 4 gnomad.genomes.r3.0.sites.pass.AnnotationReference.bcf

# Create bgzip vcf.gz version for MuTect2 germline reference resource
bcftools view \
    --threads 8 \
    --output-type z \
    --output-file gnomad.genomes.r3.0.sites.pass.AnnotationReference.vcf.gz \
    gnomad.genomes.r3.0.sites.pass.AnnotationReference.bcf
# Create the needed TBI index (hope two indexes will not cause and issue)
bcftools index --threads 4 --tbi gnomad.genomes.r3.0.sites.pass.AnnotationReference.vcf.gz

# Select Common Bi-Allelic SNPs got Mutects Contamination Testing
bcftools view \
    --threads 8 \
    --include 'AF>=0.05' \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --output-type z \
    --output-file gnomad.genomes.r3.0.sites.pass.ForMutectContamination.vcf.gz \
    gnomad.genomes.r3.0.sites.pass.AnnotationReference.vcf.gz
# Create TBI index
bcftools index --threads 4 --tbi gnomad.genomes.r3.0.sites.pass.ForMutectContamination.vcf.gz

## BROAD VERSION OF THE SAME PROCESSES (https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect_resources.wdl)

############################
###
### Download and Processing Notes
###
############################

