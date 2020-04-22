#!/usr/bin/env bash

#SBATCH --job-name="create_snpEff_DB"
#SBATCH --time=0-24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8
#SBATCH --partition defq,overflow

# Usage: sbatch --export ALL,SNPEFF_VERSION="v4_3t",SNPEFF_DB_NAME="grch38.97",SNPEFF_CONFIG_PATH="/home/snpEff.config" build_snpEff_db.sh

if [ $ENVIRONMENT == "TGen" ]
then
  # Load required module
  module load snpEff/${SNPEFF_VERSION}
else
  # Set variables from inputs
  ENVIRONMENT=${1}
  SNPEFF_VERSION=${2}
  SNPEFF_DB_NAME=${3}
  SNPEFF_CONFIG_PATH=${4}
  echo "Assuming required tools are available in $PATH"
fi

# Call snpEff build, normally this would be java -jar snpEff but the module load points to a wrapper script with module load snpEff/v4_3t
# snpEff == java -jar snpEff.jar

snpEff build -gtf22 -v ${SNPEFF_DB_NAME} -c ${SNPEFF_CONFIG_PATH}

touch CREATED_SNPEFF_${SNPEFF_DB_NAME}_DATABASE

echo "CREATED_SNPEFF_${SNPEFF_DB_NAME}_DATABASE" >> README