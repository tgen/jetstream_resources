#!/bin/sh

# dbSNP build 152 Downloaded July 27, 2019 by Jonathan Keats

# Write this full document as a README
cat $0 > README

# Load required modules
module load samtools/1.9
module load R/3.4.4

# Define any needed variables
DOWNLOADED_FASTA_GZ=/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/downloads/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
FAI=/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa.fai 

# Download the files
wget ftp.ncbi.nih.gov/snp/redesign/latest_release/release_notes.txt
wget ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/CHECKSUMS
wget ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.38.bgz
wget ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.38.bgz.md5
wget ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.38.bgz.tbi
wget ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.38.bgz.tbi.md5

# Check MD5 checksums - All downloads pass
md5sum --check CHECKSUMS

# Check contig seqeunces to determine if they match our hg38tgen reference genome with UCSC contif names
bcftools view -h GCF_000001405.38.bgz | cut -f1 | grep -v "#" | sort | uniq
# This showed the file uses standard NCBI nucleotide contigs like NC_000001.11 versus CM000663.2 versus chr1

## Need to build a rename key
# Step 1 - get the meta-data for the contigs in the dbSNP VCF
/home/jkeats/git_repositories/GRCh38_CrossMapping/utility_scripts/get_dbSNPvcf_contig_mappings.sh GCF_000001405.38.bgz b152
# Step 2 - Get the meta-data from the Reference Genome downloaded from NBCI
/home/jkeats/git_repositories/GRCh38_CrossMapping/utility_scripts/extract_metadata_from_fasta.sh ${DOWNLOADED_FASTA_GZ}
# Step 3 - Merge output files and generate list of contigs to remove the dbSNP vcf as they are not in the assembly and the rename key
Rscript /home/jkeats/git_repositories/GRCh38_CrossMapping/utility_scripts/MergeMatch_dbSNP_GRCh38_Contigs.R

# Now remove contigs that are not wanted in the dbSNP vcf as they don't exist in our refence genome (p1 versus p13 issues)
bcftools filter \
	--threads 8 \
	--targets-file ^/home/jkeats/git_repositories/GRCh38_CrossMapping/contigs_2_remove_from_dbSNP152.bed \
	--output-type b \
	--output temp_droppedContigs.bcf \
	GCF_000001405.38.bgz
bcftools index --threads 8 temp_droppedContigs.bcf

# Now rename the contigs in the processed BCF file
bcftools annotate \
	--threads 8 \
	--rename-chrs /home/jkeats/git_repositories/GRCh38_CrossMapping/GRCh38_dbSNP152_2_UCSC_Contigs.txt \
	--output-type b \
	--output temp_renamed.bcf \
	temp_droppedContigs.bcf
bcftools index --threads 8 temp_renamed.bcf

# Fix header to match our reference genome (in case it matters)
/home/jkeats/downloads/bcftools/bcftools reheader \
	--threads 4 \
	--fai ${FAI} \
	temp_renamed.bcf \
	| \
	bcftools view \
	--threads 4 \
	--output-type b \
	--output-file dbSNP_b152_hg38tgen.bcf
bcftools index --threads 4 dbSNP_b152_hg38tgen.bcf

# Get stats to confirm all processes worked
bcftools index --threads 4 --stats dbSNP_b152_hg38tgen.bcf

# Remove temp files
rm temp_*


