#!/bin/bash
#SBATCH --job-name="bwa_index"
#SBATCH --time=0-24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8

## Script to create bwa index files

## Usage:
## sbatch --export ALL,FASTA='genome.fa' bwa_index.slurm

module load bwa/0.7.17

bwa index ${FASTA}

# Error Capture
if [ "$?" = "0" ]
then
    echo "PASSED_BWA_INDEX" >> README
else
    touch FAILED_BWA_INDEX
    echo "FAILED_BWA_INDEX" >> README
    exit 1
fi