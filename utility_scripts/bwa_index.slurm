#!/bin/bash
#SBATCH --job-name="bwa_index"
#SBATCH --time=0-24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8

module load bwa/0.7.17

echo ${FASTA}

bwa index ${FASTA}