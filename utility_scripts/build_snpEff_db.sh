#!/usr/bin/env bash

#SBATCH --job-name="create_snpEff_DB"
#SBATCH --time=0-24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8

# Usage: sbatch --export ALL,SNPEFF_VERSION="v4_3t",SNPEFF_DB_NAME="grch38.97" build_snpEff_db.sh


module load snpEff/${SNPEFF_VERSION}

# Call snpEff build, normally this would be java -jar snpEff but the module load points to a wrapper script
# With module load snpEff/v4_3t
$ snpEff == java -jar snpEff.jar
snpEff build -gtf22 -v ${SNPEFF_DB_NAME}

touch CREATED_SNPEFF_${SNPEFF_DB_NAME}_DATABASE

echo "CREATED_SNPEFF_${SNPEFF_DB_NAME}_DATABASE" >> README_${SNPEFF_DB_NAME}