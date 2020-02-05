#!/bin/bash
#SBATCH --job-name="bowtie2_index"
#SBATCH --time=0-24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 10
#SBATCH --partition defq,overflow

## Script to create bowtie2 index files

## Usage:
## sbatch --export ALL,FASTA="genome.fa",BOWTIE2_MODULE="Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1",BOWTIE_BASE="genome" bowtie2_index.sh

module load ${BOWTIE2_MODULE}

bowtie2-build --threads 10 ${FASTA} ${BOWTIE_BASE}

# Error Capture
if [ "$?" = "0" ]
then
    echo "PASSED_BWA_INDEX" >> README
else
    touch FAILED_BWA_INDEX
    echo "FAILED_BWA_INDEX" >> README
    exit 1
fi