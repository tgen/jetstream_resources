#!/bin/bash

# This script is used to generate a 1 based regions file of non N genomic space of provided contigs to be used with the samtools stats --target-regions option.

# Usage: make_samtools_stats_non_N_region_file_from_fasta.awk /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa GRCh38tgen_decoy_alts_hla_samstats_no_N_1based.txt chr1CONTIG_SEPchr2 
#		 make_samtools_stats_non_N_region_file_from_fasta.awk /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa GRCh38tgen_decoy_alts_hla_samstats_no_N_1based.txt chr1CONTIG_SEPchr2CONTIG_SEPchr3CONTIG_SEPchr4CONTIG_SEPchr5CONTIG_SEPchr6CONTIG_SEPchr7CONTIG_SEPchr8CONTIG_SEPchr9CONTIG_SEPchr10CONTIG_SEPchr11CONTIG_SEPchr12CONTIG_SEPchr13CONTIG_SEPchr14CONTIG_SEPchr15CONTIG_SEPchr16CONTIG_SEPchr17CONTIG_SEPchr18CONTIG_SEPchr19CONTIG_SEPchr20CONTIG_SEPchr21CONTIG_SEPchr22

# Full or relative path to reference fasta
FASTA=$1

# Relative path to output file 
OUTPUT=$2

# "CONTIG_SEP" delimited list of contigs to include in the OUTPUT
CONTIG_LIST=$3

gawk -v CONTIG_LIST=${CONTIG_LIST} 'BEGIN {
	START=0
	CONTIG="beginning"
	split(CONTIG_LIST, CONTIGS, "CONTIG_SEP")
	for (i in CONTIGS) CONTIGS[CONTIGS[i]] = ""
	SKIP_CONTIG="No"
};
$1 ~ /^>/ {
	if (CONTIG != "beginning" && b[length(b)] != "N" && SKIP_CONTIG == "No") {
		OFS = "\t"
		print CONTIG,TOTAL-START,TOTAL-1
	}
	split($1, a, ">")
	TOTAL=1
	START=0
	RESET="No"
	split(a[2], c, " ")
	CONTIG=c[1]
	if (CONTIG in CONTIGS == 0) {
		SKIP_CONTIG="Yes"
	}
	next
}
{
	if ( SKIP_CONTIG == "No" ) {
		split($0, b, "")
		for (i=1; i <= length($0); i++) {
			
			if (b[i] == "N") {
				if ( RESET == "Yes" ) {
					OFS = "\t"
					print CONTIG,TOTAL-START,TOTAL-1
					RESET="No"
				}
				TOTAL+=1
				START=0
			} else {
				RESET="Yes"
				START+=1
				TOTAL+=1
			}
		}
	}
} END {
	if (b[length(b)] != "N" && SKIP_CONTIG == "No") {
		OFS = "\t"
		print CONTIG,TOTAL-START,TOTAL-1
	}
}' ${FASTA} > ${OUTPUT}



