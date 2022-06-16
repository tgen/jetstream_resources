#!/bin/bash
#SBATCH --job-name="bwa_sequenza_gc_wiggle"
#SBATCH --time=0-24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 2
#SBATCH --partition defq

## Script to create bwa index files

## Usage:
## sbatch --export ALL,FASTA="genome.fa" sequenza_gc_wiggle.sh

module load sequenza/3.0.0

sequenza−utils gc_wiggle −w 50 --fasta ${FASTA} -o ${FASTA::-3}.gc50Base.wig.gz
