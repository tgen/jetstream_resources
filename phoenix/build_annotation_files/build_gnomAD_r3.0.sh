#!/usr/bin/env bash

# Automated Script to download and build GnomAD release 3.0 genome VCF file for usage in Phoenix workflow

# Usage: ./build_gnomAD_r3.0.sh <resources.ini>
# Sbatch Usage: sbatch --nodes=1 --cpus-per-task=10 --time=0-72:00:00 build_gnomAD_r3.0.sh <resources.ini>

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
if [ ${ENVIRONMENT} == "TGen" ]
then
  module load BCFtools/1.10.1-foss-2019a
  module load Python/3.7.2-foss-2019a
elif [ ${ENVIRONMENT} == "LOCAL" ]
then
  echo
  echo "Assuming required tools are available in $PATH"
  echo
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Enviroment Variable"
  echo "Only TGen or LOCAL are supported"
  exit 1
fi

####################################
## Create Expected Folder Structure
###################################

# Make top level directory if not available
if [ -e ${PARENT_DIR} ]
then
    echo "Parent directory: ${PARENT_DIR} exists, moving into it"
    cd ${PARENT_DIR}
else
    echo "Parent directory NOT fount, creating and moving into it now"
    mkdir -p ${PARENT_DIR}
    cd ${PARENT_DIR}
fi

# Make public_databases folder if not available
if [ -e public_databases ]
then
    echo "Public Databases folder exists, moving into it"
    cd public_databases
else
    echo "Public Databases folder NOT fount, creating and moving into it now"
    mkdir -p public_databases
    cd public_databases
fi

# Make dbSNP folder if not available
if [ -e gnomad ]
then
    echo "gnomad folder exists, moving into it"
    cd gnomad
else
    echo "gnomad folder NOT fount, creating and moving into it now"
    mkdir -p gnomad
    cd gnomad
fi

# Make r3.0 release version folder if not available
if [ -e r3.0 ]
then
    echo "r3.0 folder exists, exiting to prevent overwrite"
    exit 1
else
    echo "r3.0 folder NOT fount, creating and moving into it now"
    mkdir -p r3.0
    cd r3.0
fi

####################################
## Download and Manipulate the gnomAD file
###################################

# Initialize a gnomeAD r3.0 index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

############################
###
### Download and Processing Genome File
###
############################
echo
echo

# Download the current genome VCF
wget https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz
#235.68 GiB, MD5: f3501102192975da34b5d2c32f7c0791
fc -ln -1 >> README
echo >> README

# Check the md5sum
echo "f3501102192975da34b5d2c32f7c0791  gnomad.genomes.r3.0.sites.vcf.bgz" > gnomad.genomes.r3.0.sites.vcf.bgz.md5
fc -ln -1 >> README
echo >> README
md5sum --check gnomad.genomes.r3.0.sites.vcf.bgz.md5 | tee --append README
fc -ln -1 >> README
echo >> README
echo

## Using these massive VCF files we discovered there a significant benefits to using bcf files
## Also for internal usage we limit the VCF to just the PASS variants
bcftools view \
    --threads 8 \
    --apply-filters PASS \
    --output-type b \
    --output-file gnomad.genomes.r3.0.sites.pass.bcf \
    gnomad.genomes.r3.0.sites.vcf.bgz
fc -ln -1 >> README
echo >> README

bcftools index --threads 4 gnomad.genomes.r3.0.sites.pass.bcf
fc -ln -1 >> README
echo >> README

bcftools stats --threads 8 gnomad.genomes.r3.0.sites.pass.bcf > gnomad.genomes.r3.0.sites.pass.bcf.stats
fc -ln -1 >> README
echo >> README

plot-vcfstats --no-PDF --title "GnomADr3.0 Genomes" -p plots_vcfstats_genomes gnomad.genomes.r3.0.sites.pass.bcf.stats
fc -ln -1 >> README
echo >> README

## Also using these files for annotation removing the massive amount of data in the INFO column significantly improves runtime
## For utility with mutect2 make simple version with minimal data for annotation yes/no flags and a germline reference
bcftools annotate \
    --threads 8 \
    --remove ID,QUAL,^INFO/AF \
    --output-type b \
    --output gnomad.genomes.r3.0.sites.pass.AnnotationReference.bcf \
    gnomad.genomes.r3.0.sites.pass.bcf
fc -ln -1 >> README
echo >> README

# Create CSI index (This will be file used for bcftools database annotation purposes
bcftools index --threads 4 gnomad.genomes.r3.0.sites.pass.AnnotationReference.bcf
fc -ln -1 >> README
echo >> README

# Create bgzip vcf.gz version for MuTect2 germline reference resource
bcftools view \
    --threads 8 \
    --output-type z \
    --output-file gnomad.genomes.r3.0.sites.pass.AnnotationReference.vcf.gz \
    gnomad.genomes.r3.0.sites.pass.AnnotationReference.bcf
fc -ln -1 >> README
echo >> README

# Create the needed TBI index (hope two indexes will not cause and issue)
bcftools index --threads 4 --tbi gnomad.genomes.r3.0.sites.pass.AnnotationReference.vcf.gz
fc -ln -1 >> README
echo >> README

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
fc -ln -1 >> README
echo >> README

# Create TBI index
bcftools index --threads 4 --tbi gnomad.genomes.r3.0.sites.pass.ForMutectContamination.vcf.gz
fc -ln -1 >> README
echo >> README

## BROAD VERSION OF THE SAME PROCESSES (https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect_resources.wdl)

############################
###
### Download and Processing Notes
###
############################

