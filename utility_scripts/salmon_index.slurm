#!/bin/bash
#SBATCH --job-name="salmon_index"
#SBATCH --time=0-24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8

##############
## Make Salmon Index files
##
##

module load salmon/0.12.0

# NOTES: This will build the quasi-mapping-based index, using an auxiliary k-mer hash over k-mers of length 31. While quasi-mapping will make used of arbitrarily long matches between the query and reference, the k size selected here will act as the minimum acceptable length for a valid match. Thus, a smaller value of k may slightly improve sensitivty. We find that a k of 31 seems to work well for reads of 75bp or longer, but you might consider a smaller k if you plan to deal with shorter reads.

salmon index -t ${TRANSCRIPTOME_FASTA} -i salmon_quasi_75merPlus --type quasi -k 31

# Error Capture
if [ "$?" = "0" ]
then
    echo "PASSED_SALMON_INDEX" >> README
else
    touch FAILED_SALMON_INDEX
    echo "FAILED_SALMON_INDEX" >> README
    exit 1
fi