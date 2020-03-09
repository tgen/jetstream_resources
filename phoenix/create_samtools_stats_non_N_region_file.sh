#!/usr/bin/env bash

# Usage: create_bwa_genome_index.sh <Config.ini>

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -ue

# Read required variables from configuration file
. ${1}

####################################
## Navigate Directory Structure
###################################

# Check top level directory if not available
if [ -e ${TOPLEVEL_DIR} ]
then
    echo "Top level directory: ${TOPLEVEL_DIR} exists, moving into it"
    cd ${TOPLEVEL_DIR}
else
    echo "Top level directory NOT found, IT IS REQUIRED, EXITING"
    exit 1
fi

# Check that the reference genome was created successfully
if [ -e GENOME_FASTA_GENERATION_COMPLETE ]
then
    echo "Genome fasta exists, moving forward"
else
    echo "Genome fasta generation complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
fi

####################################
## Generate BWA index
####################################

if [ -e tool_resources ]
then
    echo "tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "tool_resources directory NOT found, creating and moving into it now"
    mkdir tool_resources
    cd tool_resources
fi

if [ -e "samtools_stats" ]
then
    echo "The samtools_stats directory exists, exiting to prevent overwriting existing non_N_region_file."
    exit 2
else
    echo "The samtools_stats directory was NOT found, creating and moving into it now"
    mkdir samtools_stats
    cd samtools_stats
fi

# Initialize a samtools_stats index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/phoenix" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "samtools stats --target-regions non_N_region_file creation details:" >> README
echo >> README
echo "By default samtools stats consideres all bases, including N's, in the reference genome when calculating stats." > README
echo "The 1 based regions file generated by make_samtools_stats_non_N_region_file_from_fasta.awk excludes all N bases" > README 
echo "and is passed to the --target-regions option of samtools stats to ensure coverage and mapping statistics are " >> README
echo "calculated based on non N genome space." >> README
echo >> README

# Create bwa index files using bwa utility script
echo "Create samtools stats non N region file of primary contigs with chrX and chrY removed as follows:" >> README
sbatch --export ALL ${PATH_TO_REPO}/utility_scripts/make_samtools_stats_non_N_region_file_from_fasta.awk ${REFERENCE_DNA_GENOME_FASTA} GRCh38tgen_decoy_alts_hla_samstats_no_N_1based_primary_contigs_no_chrX_chrY.txt chr1CONTIG_SEPchr2CONTIG_SEPchr3CONTIG_SEPchr4CONTIG_SEPchr5CONTIG_SEPchr6CONTIG_SEPchr7CONTIG_SEPchr8CONTIG_SEPchr9CONTIG_SEPchr10CONTIG_SEPchr11CONTIG_SEPchr12CONTIG_SEPchr13CONTIG_SEPchr14CONTIG_SEPchr15CONTIG_SEPchr16CONTIG_SEPchr17CONTIG_SEPchr18CONTIG_SEPchr19CONTIG_SEPchr20CONTIG_SEPchr21CONTIG_SEPchr22
fc -ln -1 >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/make_samtools_stats_non_N_region_file_from_fasta.awk >> README
echo >> README
echo >> README
