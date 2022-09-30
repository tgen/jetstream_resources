#!/usr/bin/env bash

# Automated Script to download and build GnomAD genome VCF file for usage in Phoenix workflow

# Usage: ./build_gnomAD.sh <resources.ini>
# Sbatch Usage: sbatch --nodes=1 --cpus-per-task=10 --time=0-72:00:00 build_gnomAD.sh <resources.ini>

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

# Make r${GNOMAD_GENOME_VERSION} release version folder if not available
if [ -e r${GNOMAD_GENOME_VERSION} ]
then
    echo "r${GNOMAD_GENOME_VERSION} folder exists, exiting to prevent overwrite"
    exit 1
else
    echo "r${GNOMAD_GENOME_VERSION} folder NOT fount, creating and moving into it now"
    mkdir -p r${GNOMAD_GENOME_VERSION}
    cd r${GNOMAD_GENOME_VERSION}
fi

####################################
## Download and Manipulate the gnomAD file
###################################

# Initialize a gnomeAD r${GNOMAD_GENOME_VERSION} index README
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

for chr in chr{1..22} chrX chrY; do
  sbatch -o /dev/null -c 4 --wait -p overflow --wrap="module load BCFtools/1.10.1-foss-2019a singularity/3.7.1-phoenix ;
singularity exec docker://google/cloud-sdk gsutil cp gs://gcp-public-data--gnomad/release/${GNOMAD_GENOME_VERSION}/vcf/genomes/gnomad.genomes.v${GNOMAD_GENOME_VERSION}.sites.${chr}.vcf.bgz - |
bcftools view --apply-filters PASS --output-type u --output-file gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.${chr}.bcf -" > /dev/null 2>&1 &
done
fc -ln -1 >> README
echo >> README

# Wait for the download and bcftools conversions to finish
wait

# Check the md5sum
# echo "f3501102192975da34b5d2c32f7c0791  gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.vcf.bgz" > gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.vcf.bgz.md5
# fc -ln -1 >> README
# echo >> README
# md5sum --check gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.vcf.bgz.md5 | tee --append README
# fc -ln -1 >> README
# echo >> README
# echo
#### Expectation is that 1. gsutil is more reliable (rsync-ish) than wget, 2. bcftools would probably complain if something weird happens

## Also using these files for annotation removing the massive amount of data in the INFO column significantly improves runtime
## For utility with mutect2 make simple version with minimal data for annotation yes/no flags and a germline reference
for chr in chr{1..22} chrX chrY; do
  sbatch -o /dev/null -c 8 --wait -p defq,overflow --wrap="module load BCFtools/1.10.1-foss-2019a ;
bcftools annotate --threads 8 --remove ID,QUAL,^INFO/AC,^INFO/AN,^INFO/AF,^INFO/n_alt_alleles --output-type u --output gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.reduced.${chr}.bcf gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.${chr}.bcf" > /dev/null 2>&1 &
done
fc -ln -1 >> README
echo >> README

# Wait for reduced bcfs
wait

## Concatentate the bcfs for a monolithic file
concat_cmd="bcftools concat --output-type b --output gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.bcf"
for chr in chr{1..22} chrX chrY; do
  concat_cmd+=" gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.${chr}.bcf"
done
# We are creating this for users other than us, we use the reduced and annotationReference bcfs in the pipeline
sbatch -o concat_pass_sites.slurm.out -c 4 --wrap="module load BCFtools/1.10.1-foss-2019a; $concat_cmd"
fc -ln -1 >> README
echo >> README

bcftools index --threads 8 gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.bcf
fc -ln -1 >> README
echo >> README

## Concatentate the bcfs for a monolithic file
concat_cmd="bcftools concat --output-type b --output gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.reduced.bcf"
for chr in chr{1..22} chrX chrY; do
  concat_cmd+=" gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.reduced.${chr}.bcf"
done
$concat_cmd
fc -ln -1 >> README
echo >> README

# Create CSI index (This will be file used for bcftools database annotation purposes
bcftools index --threads 8 gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.reduced.bcf
fc -ln -1 >> README
echo >> README

# AnnotationReference
bcftools annotate \
    --threads 8 \
    --remove ID,QUAL,^INFO/AF \
    --output-type b \
    --output gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.AnnotationReference.bcf \
    gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.reduced.bcf
fc -ln -1 >> README
echo >> README

# Create CSI index (This will be file used for bcftools database annotation purposes
bcftools index --threads 8 gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.AnnotationReference.bcf
fc -ln -1 >> README
echo >> README

# Create bgzip vcf.gz version for MuTect2 germline reference resource
bcftools view \
    --threads 8 \
    --output-type z \
    --output-file gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.AnnotationReference.vcf.gz \
    gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.AnnotationReference.bcf
fc -ln -1 >> README
echo >> README

# Create the needed TBI index (hope two indexes will not cause and issue)
bcftools index --threads 8 --tbi gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.AnnotationReference.vcf.gz
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
    --output-file gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.ForMutectContamination.vcf.gz \
    gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.AnnotationReference.vcf.gz
fc -ln -1 >> README
echo >> README

# Create TBI index
bcftools index --threads 4 --tbi gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.ForMutectContamination.vcf.gz
fc -ln -1 >> README
echo >> README

# Wait for pending slurm jobs (Should only be the one for creating the monolithic pass bcf)
wait

# Clean up temp chr bcfs
for chr in chr{1..22} chrX chrY; do
  rm gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.${chr}.bcf
  rm gnomad.genomes.r${GNOMAD_GENOME_VERSION}.sites.pass.reduced.${chr}.bcf
done

## BROAD VERSION OF THE SAME PROCESSES (https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect_resources.wdl)
