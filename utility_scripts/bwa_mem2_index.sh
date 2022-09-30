#!/bin/bash
#SBATCH --job-name="bwa_mem2_index"
#SBATCH --time=0-24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 16G
#SBATCH --partition defq

## Script to create bwa index files

## Usage:
## sbatch --export ALL,FASTA="genome.fa",BWA_MEM2_VERSION="2.2.1" bwa_index.sh

module load bwa-mem2/${BWA_MEM2_VERSION}

bwa-mem2 index ${FASTA}

# Error Capture
if [ "$?" = "0" ]
then
    echo "PASSED_BWA_INDEX" >> README
else
    touch FAILED_BWA_INDEX
    echo "FAILED_BWA_INDEX" >> README
    exit 1
fi
