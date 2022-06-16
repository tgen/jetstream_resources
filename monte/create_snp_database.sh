#!/usr/bin/env bash

# Usage: create_snp_database.sh <Config.ini>

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

# Check if the tool resources directory exists
if [ -e public_databases ]
then
    echo "public_databases directory exists, moving into it"
    cd public_databases
else
    echo "public_databases directory NOT found, creating and moving into it now"
    mkdir public_databases
    cd public_databases
fi

####################################
## Prepping common variables
####################################

REF_PATH=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_DNA_GENOME_NAME}
CHAIN_PATH=$(dirname ${PARENT_DIR})/liftover_files/${FINAL_CHAIN_NAME}


####################################
## Loading necessary modules
####################################

module load BCFtools/1.10.1-foss-2019a
module load GATK/4.1.8.0-GCCcore-8.3.0-Java-1.8

####################################
## Download dogSD files
####################################

if [ -e "dog_sd" ]
then
    echo "The dog_sd directory exists, moving into it"
    cd dog_sd
else
    echo "The dog_sd directory was NOT found, creating and moving into it now"
    mkdir dog_sd ; chmod 777 dog_sd
    cd dog_sd
fi

touch README
echo >> README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

echo "wget ${SNP_FILE_DOWNLOAD_LINK}" >> README
echo "wget ${SNP_FILE_MD5SUM_DOWNLOAD_LINK}" >> README
echo >> README

dogSNP=`basename ${SNP_FILE_DOWNLOAD_LINK}`
dogMD5=`basename ${SNP_FILE_MD5SUM_DOWNLOAD_LINK}`

if [ -e ${dogSNP} ] || [ -e ${dogSNP::-4} ]
then
    echo "${dogSNP} or ${dogSNP::-4} exists, moving on"
else
    wget ${SNP_FILE_DOWNLOAD_LINK}
    wget ${SNP_FILE_MD5SUM_DOWNLOAD_LINK}

if [ `md5sum ${dogSNP} | cut -d " " -f 1` != `grep -m 1 ${dogSNP} ${dogMD5} | cut -d " " -f 1` ]; then
    echo "The md5s do not match for ${dogSNP}, please run the script to again."
    exit 1
else
    echo "md5sum ${dogSNP} | cut -d " " -f 1" >> README
    md5sum ${dogSNP} | cut -d " " -f 1 >> README
    md5sum ${dogSNP} | cut -d " " -f 1
    echo "grep -m 1 ${dogSNP} ${dogMD5} | cut -d " " -f 1" >> README
    grep -m 1 ${dogSNP} ${dogMD5} | cut -d " " -f 1 >> README
    grep -m 1 ${dogSNP} ${dogMD5} | cut -d " " -f 1
fi
fi
echo >> README

############################
###
### Create Genome and Exome versions for MuTect Input
###
############################


if [ -e ${dogSNP::-4} ]
then
    echo "${dogSNP::-4} exists, moving on"
else
    bzip2 -d ${dogSNP}
fi

# We need to rename the contigs to match chain file
for contig in $(bcftools view -h ${dogSNP::-4} | grep contig | cut -d'=' -f3 | cut -d',' -f1); do 
  new_contig=$(echo ${contig} | sed 's/chr//g' | sed 's/Un_//g' | awk '{ if ($1 ~ /(^JH|^A)/) $1 = $1 ".1" }1')
  echo "${contig} ${new_contig}" >> chr_rename.txt
done

# Make Allele Frequency Only VCF
# Clear ID and QUAL and delete all INFO fields other than AF
bcftools annotate \
    --threads 8 \
    --remove ID,QUAL,^INFO/AF \
    --rename-chrs chr_rename.txt \
    --output-type z \
    --output ${dogSNP::-8}_ForMutect.vcf.gz \
    ${dogSNP::-4}

# Create TBI index
bcftools index --threads 4 --tbi ${dogSNP::-8}_ForMutect.vcf.gz

bcftools view \
    --threads 8 \
    --apply-filters PASS \
    --output-type z \
    --output-file ${dogSNP::-8}_ForMutect2.vcf.gz \
    ${dogSNP::-8}_ForMutect.vcf.gz

mv ${dogSNP::-8}_ForMutect2.vcf.gz ${dogSNP::-8}_ForMutect.vcf.gz
rm ${dogSNP::-8}_ForMutect.vcf.gz.tbi
bcftools index --threads 4 --tbi ${dogSNP::-8}_ForMutect.vcf.gz

# Select Common Bi-Allelic SNPs got Mutects Contamination Testing
bcftools view \
    --threads 8 \
    --apply-filters PASS \
    --include 'AF>=0.05' \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --output-type z \
    --output-file ${dogSNP::-8}_ForMutectContamination.vcf.gz \
    ${dogSNP::-8}_ForMutect.vcf.gz
# Create TBI index
bcftools index --threads 4 --tbi ${dogSNP::-8}_ForMutectContamination.vcf.gz

echo "Processing liftover from canfam3.1 to ros_cfam1.0 on ${dogSNP}"

gatk LiftoverVcf \
  --INPUT ${TOPLEVEL_DIR}/public_databases/dog_sd/${dogSNP::-8}_ForMutectContamination.vcf.gz \
  --OUTPUT ${TOPLEVEL_DIR}/public_databases/dog_sd/${dogSNP::-8}_ForMutectContamination.${GENOME_SUBVERSION_NAME}.vcf.gz \
  --CHAIN ${CHAIN_PATH} \
  --REFERENCE_SEQUENCE ${REF_PATH} \
  --REJECT ${TOPLEVEL_DIR}/public_databases/dog_sd/${dogSNP::-8}.${GENOME_SUBVERSION_NAME}.failed.vcf.gz\
  --LIFTOVER_MIN_MATCH 0.9 \
  --RECOVER_SWAPPED_REF_ALT \
  --WARN_ON_MISSING_CONTIG

cd ${TOPLEVEL_DIR}/public_databases


####################################
## Download EVA snp files
####################################

if [ -e "eva" ]
then
    echo "The eva directory exists, exiting to prevent overwriting existing index"
    exit 2
else
    echo "The eva directory was NOT found, creating and moving into it now"
    mkdir eva
    cd eva
fi

touch README
echo >> README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

echo "wget ${SNP_VCF_FILE_DOWNLOAD_LINK}*" >> README
echo >> README
wget ${SNP_VCF_FILE_DOWNLOAD_LINK}*

echo "wget ${SNP_VCF_MD5SUM_DOWNLOAD_LINK}" >> README
echo >> README
wget ${SNP_VCF_MD5SUM_DOWNLOAD_LINK}

evaSNP=`basename ${SNP_VCF_FILE_DOWNLOAD_LINK}`
evaMD5=`basename ${SNP_VCF_MD5SUM_DOWNLOAD_LINK}`
if [ `md5sum ${evaSNP} | cut -d " " -f 1` != `grep -m 1 ${evaSNP} ${evaMD5} | cut -d " " -f 1` ]; then
    echo "The md5s do not match for ${evaSNP}, please run the script to again."
    exit 1
else
    echo "md5sum ${evaSNP} | cut -d " " -f 1" >> README
    md5sum ${evaSNP} | cut -d " " -f 1 >> README
    md5sum ${evaSNP} | cut -d " " -f 1
    echo "grep -m 1 ${evaSNP} ${evaMD5} | cut -d " " -f 1" >> README
    grep -m 1 ${evaSNP} ${evaMD5} | cut -d " " -f 1 >> README
    grep -m 1 ${evaSNP} ${evaMD5} | cut -d " " -f 1
fi
echo >> README

# We need to rename the contigs to match chain file
for contig in $(bcftools view -h ${evaSNP} | grep contig | cut -d'=' -f3 | cut -d',' -f1); do 
  new_contig=$(echo ${contig} | sed 's/chr0//g' | sed 's/chr//g')
  echo "${contig} ${new_contig}" >> chr_rename.txt
done

bcftools annotate \
    --threads 8 \
    --rename-chrs chr_rename.txt \
    --output-type z \
    --output ${evaSNP::-8}.chr_renamed.vcf.gz \
    ${evaSNP}

mv ${evaSNP::-8}.chr_renamed.vcf.gz ${evaSNP}

echo "Processing liftover from canfam3.1 to ros_cfam1.0 on ${evaSNP}"

gatk LiftoverVcf \
  --INPUT ${TOPLEVEL_DIR}/public_databases/eva/${evaSNP} \
  --OUTPUT ${TOPLEVEL_DIR}/public_databases/eva/${evaSNP::-8}.${GENOME_SUBVERSION_NAME}.vcf.gz \
  --CHAIN ${CHAIN_PATH} \
  --REFERENCE_SEQUENCE ${REF_PATH} \
  --REJECT ${TOPLEVEL_DIR}/public_databases/eva/${evaSNP::-8}.${GENOME_SUBVERSION_NAME}.failed.vcf.gz\
  --LIFTOVER_MIN_MATCH 0.9 \
  --RECOVER_SWAPPED_REF_ALT \
  --WARN_ON_MISSING_CONTIG
