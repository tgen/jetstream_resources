#!/bin/bash
#SBATCH --job-name="star_index"
#SBATCH --time=0-24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 20

## Script to create star index files

## Usage:
## sbatch --export ALL,GTF='transcriptome.gtf',FASTA='genome.fa',SJDB_OVERHANG='74',INDEX_DIR='75bpReads' star_index.slurm

set -ux
module load STAR/2.6.1d

#GTF=../../Homo_sapiens.GRCh38.95.ucsc.gtf
#FASTA=../../../../Homo_sapiens.GRCh38.noalt.fasta
#SJDB_OVERHANG="${1}"
#INDEX_DIR=./GRCh38.95_${1}bp/

if [ -d "${INDEX_DIR}" ]; then
  echo "Index already exists: ${INDEX_DIR}"
  exit 1
else
  mkdir -p "${INDEX_DIR}"
fi

STAR \
  --runMode genomeGenerate \
  --genomeDir "${INDEX_DIR}" \
  --runThreadN 19 \
  --sjdbOverhang "${SJDB_OVERHANG}" \
  --genomeFastaFiles "${FASTA}" \
  --sjdbGTFfile "${GTF}"

# Error Capture
if [ "$?" = "0" ]
then
    echo "PASSED_STAR_INDEX" >> README
else
    touch FAILED_STAR_INDEX
    echo "FAILED_STAR_INDEX" >> README
    exit 1
fi