#!/bin/bash
#SBATCH --job-name="create_transcriptome_fasta"
#SBATCH --time=0-24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8

## Script to create a transcriptome fasta file from an input GTF and Genome.fa files

## Usage:
## sbatch --export ALL,GENOME='genome.fa',GTF='transcriptome.gtf',OUTPUT='output.fasta' create_transcript_fasta.slurm

## Extract Transcriptome Fasta
## - use gffread submodule in cufflinks
## - processes plus and negative strand genes correctly, producing fasta file with nucleotides in transcript order

module load cufflinks/2.2.1

gffread ${GTF} -g ${GENOME} -w ${OUTPUT}

# Error Capture
if [ "$?" = "0" ]
then
    echo "PASSED_CREATE_TRANSCRIPTOME_FASTA" >> README
else
    touch FAILED_CREATE_TRANSCRIPTOME_FASTA
    echo "FAILED_CREATE_TRANSCRIPTOME_FASTA" >> README
    exit 1
fi

# Check the number of unique transcripts produced to ensure entire GTF was converted
FASTA_TRANSCRIPTS=`grep "^>" ${OUTPUT} | cut -d" " -f1 | sort | uniq | wc -l`
GTF_TRANSCRIPTS=`cut -f9 ${GTF} | grep "transcript_id" | sed 's/; transcript_id "/\'$'\t''/g' | sed 's/"; transcript_version/\'$'\t''/g' | cut -f2 | sort | uniq | wc -l`

if [ ${GTF_TRANSCRIPTS} -eq ${FASTA_TRANSCRIPTS} ]
then
    echo "Output transcripts match as expected"
    echo "Output transcripts match as expected" >> README
    exit 0
else
    echo "Output transcripts DO NOT match as expected"
    echo "Output transcripts DO NOT match as expected" >> README
    echo "INPUT GTF TRANSCRIPT COUNT = ${GTF_TRANSCRIPTS}" >> README
    echo "OUTPUT FASTA TRANSCRIPT COUNT = ${FASTA_TRANSCRIPTS}" >> README
    exit 1
fi