#!/bin/bash -i
HELP="Build_hg38tgen.sh
usage: Build_hg38tgen.sh <PARENT_DIR> <REFDIRNAME(optional)> <CREATOR(optional)> 
"
### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -ue

PARENT_DIR=${1:-$(pwd)}
TOPLEVEL_DIR=${2:-hg38tgen}
CREATOR=${3:-${USER}}

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PATH_TO_REPO="$(dirname "${SCRIPT_DIR}")"


####################################
## Configure and make Directory Structure
####################################

# Change to parent directory
cd ${PARENT_DIR}

# Make top level directory if not available
if [ -e ${TOPLEVEL_DIR} ]
then
    echo "Top level directory: ${TOPLEVEL_DIR} exists, moving into it"
    cd ${TOPLEVEL_DIR}
else
    echo "Top level directory NOT fount, creating and moving into it now"
    mkdir ${TOPLEVEL_DIR}
    cd ${TOPLEVEL_DIR}
fi

# Initialize a top level README
touch README
echo >> README
echo "Reference Genome and related files required for JetStream Phoenix Workflow" >> README
echo >> README
echo "For details on file creation see the associated github repository:"
echo "https://github.com/tgen/jetstream_resources"
echo >> README
echo "
## The GRCh38 reference genome is represented in different public locations and not all have defined source information
## This file tracks the source and generation of the reference and annotation files used in the Jetstream phoenix workflow

## NOTE:
### - Because BWA is alt aware and STAR does not support alternative contigs we are creating two versions
### - Because BWAKIT and the Broad bundle lack sufficient source information we are creating references independently
### ---- However, the BWAKIT version is identical to the broad bundle version, albeit with different base file names
### - Also, review of what is used by GDC showed they use an identical base genome (1-22, X, Y, M, supercontigs, decoys, EBV) but with out alternative contigs
### ---- To this they added a series of HPV and other viral genomes
### ---- They use the same reference for BWA and STAR alignments
"  >> README

####################################
## Create reference genomes from known sources
####################################

# Make reference_genome directory if not available
if [ -e genome_reference ]
then
    echo "Genome Reference directory exists, moving into it"
    cd genome_reference
else
    echo "Genome Reference directory NOT fount, creating and moving into it now"
    mkdir genome_reference
    cd genome_reference
fi

# Initialize a reference_genome README
touch README
echo >> README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Download GRC/NCBI README
echo "Download README from NCBI" >> README
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt
fc -ln -1 >> README
echo >> README

####################################
## Download BWA REFERENCE GENOME
####################################

echo "####################################" >> README
echo "## Download reference genome fasta for BWA " >> README
echo "####################################" >> README
echo >> README

# Create a downloads directory to store downloads after processing
if [ -e downloads ]
then
    echo "Downloads directory exists, doing nothing"
else
    echo "Downloads directory DOES NOT exist, creating it now"
    mkdir downloads
fi

# Download the primary assembly with decoy sequences:
echo "Download primary assembly with decoy sequences:" >> README
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
fc -ln -1 >> README
# Archive a copy of the download file
echo "Archive a copy of the download file" >> README
cp GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz downloads/
fc -ln -1 >> README
# Decompress the reference
gunzip GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
fc -ln -1 >> README
echo >> README

# Download bwakit, to get the HLA contigs used by bwakit and the production GRCh38 pipeline used at Broad
echo "Download bwakit, to get the HLA contigs used by bwakit and the production GRCh38 pipeline used at Broad:" >> README
wget https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2
fc -ln -1 >> README
# Decompress the bwakit download
tar xvjf bwakit-0.7.15_x64-linux.tar.bz2
fc -ln -1 >> README
echo >> README

### Create a fasta file with just the HLA from BWAKIT
# Find the line with the first HLA region
echo "The bwakit hs38DH-extra.fa file contains decoy and HLA contigs, need to extract just the HLA ones" >> README
echo "Hand curated the name of the first HLA contig, found line number with awk as follows:" >> README
awk '/HLA-A\*01:01:01:01/{print FNR}' bwa.kit/resource-GRCh38/hs38DH-extra.fa >> README
fc -ln -1 >> README
# Identified line
# 75967
# Print all lines including this line and afterwards
echo "Print all lines after and including first HLA line" >> README
awk 'NR>=75967' bwa.kit/resource-GRCh38/hs38DH-extra.fa > bwakit_HLA.fasta
fc -ln -1 >> README
echo >> README

### Concatenate the fasta files to make final reference genome fasta
echo "Concatenate the fasta files to make final reference genome fasta to be used by BWA" >> README
cat GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna bwakit_HLA.fasta > GRCh38tgen_decoy_alts_hla.fa
fc -ln -1 >> README
echo >> README

# Add symbolic link to indicate which FASTA is used by BWA
ln -s GRCh38tgen_decoy_alts_hla.fa BWA_FASTA

# Create faidx and dict files
echo "Create BWA faidx index using samtools" >> README
module load samtools/1.9
fc -ln -1 >> README
samtools faidx GRCh38tgen_decoy_alts_hla.fa
fc -ln -1 >> README
echo >> README

echo "Create BWA dictionary file using samtools" >> README
samtools dict --assembly GRCh38 \
--species "Homo sapiens" \
--uri "downloads/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz" \
--output GRCh38tgen_decoy_alts_hla.dict \
GRCh38tgen_decoy_alts_hla.fa
fc -ln -1 >> README
echo >> README

# Clean up the directory to store the downloads
mv bwakit-0.7.15_x64-linux.tar.bz2 downloads
mv bwa.kit downloads
mv README_analysis_sets.txt downloads
mv bwakit_HLA.fasta downloads
rm GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna

####################################
## INDEX BWA REFERENCE GENOME
####################################
cd ..

if [ -e tool_resources ]
then
    echo "tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "tool_resources directory NOT fount, creating and moving into it now"
    mkdir tool_resources
    cd tool_resources
fi

if [ -e "bwa_0.7.17" ]
then
    echo "The BWA directory exists, moving into it"
    cd "bwa_0.7.17"
else
    echo "The BWA directory NOT fount, creating and moving into it now"
    mkdir "bwa_0.7.17"
    cd "bwa_0.7.17"
fi

# Initialize a bwa index README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "BWA Index creation details:" >> README
echo >> README
echo "There are no details on how the alt-aware required .alt file is created" >> README
echo "Because of this we are copying and renaming the file from BWAKIT" >> README
echo >> README

# Copy in the .alt file from bwa.kit
echo "Copy in the .alt file from bwa.kit" >> README
cp ../../genome_reference/downloads/bwa.kit/resource-GRCh38/hs38DH.fa.alt GRCh38tgen_decoy_alts_hla.fa.alt
fc -ln -1 >> README
echo >> README

# Create a symbolic link to the reference genome
ln -s ../../genome_reference/GRCh38tgen_decoy_alts_hla.fa GRCh38tgen_decoy_alts_hla.fa

# Create bwa index files using bwa utility script
echo "Create bwa index as follows:" >> README
sbatch --export ALL,FASTA='GRCh38tgen_decoy_alts_hla.fa' ${PATH_TO_REPO}/utility_scripts/bwa_index.slurm
fc -ln -1 >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/bwa_index.slurm >> README
echo >> README
echo >> README

####################################
## Download STAR REFERENCE GENOME
####################################

# move back to the genome_reference directory
cd ../../genome_reference

echo "####################################" >> README
echo "## Download reference genome fasta for STAR " >> README
echo "####################################" >> README
echo >> README

# Since STAR does not support alternative contigs, download version wtih out
echo "Download fasta file for STAR, version without alternate contigs and no HLA alleles added" >> README
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
fc -ln -1 >> README
echo >> README

# Archive a copy in the downloads section then decompress
echo "Archive compressed STAR fasta file download then decompress" >> README
cp GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz downloads/
fc -ln -1 >> README
gunzip GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
fc -ln -1 >> README
echo >> README

# Rename fastq file
echo "Rename the downloaded file to decrease filename length and make consistent with bwa fasta" >> README
mv GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna GRCh38tgen_decoy.fa
fc -ln -1 >> README
echo >> README

# Add symbolic link to indicate which FASTA is used by STAR
ln -s GRCh38tgen_decoy.fa STAR_FASTA

# Create faidx and dict files
echo "Create STAR faidx index using samtools" >> README
samtools faidx GRCh38tgen_decoy.fa
fc -ln -1 >> README
echo >> README

echo "Create STAR dictionary file using samtools" >> README
samtools dict --assembly GRCh38 \
--species "Homo sapiens" \
--uri "downloads/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz" \
--output GRCh38tgen_decoy.dict \
GRCh38tgen_decoy.fa
fc -ln -1 >> README
echo >> README

####################################
## Download Gene Model File
####################################

# move back to top level directory
cd ${PARENT_DIR}/${TOPLEVEL_DIR}

# Make gene_model directory if not available
if [ -e gene_model ]
then
    echo "Gene Model directory exists, moving into it"
    cd gene_model
else
    echo "Gene Model directory NOT fount, creating and moving into it now"
    mkdir gene_model
    cd gene_model
fi

# Make ensembl_v95 directory if not available
if [ -e ensembl_v95 ]
then
    echo "Gene Model directory exists, moving into it"
    cd ensembl_v95
else
    echo "Gene Model directory NOT fount, creating and moving into it now"
    mkdir ensembl_v95
    cd ensembl_v95
fi

# Initialize a gene_model specific README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "Ensembl version 95 download and manipulations to align with reference genomes" >> README
echo >> README

# Download the gene model gtf from ensembl
echo "Download the gene model gtf from ensembl" >> README
wget ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz
fc -ln -1 >> README
gunzip Homo_sapiens.GRCh38.95.gtf.gz
fc -ln -1 >> README
echo >> README

####################################
## Create updated gene model with ucsc style contig names - JK Method, generates results identical to RR and AC
####################################

## FILE ACCOUNTING
echo "Record the number of lines in original file" >> README
wc -l Homo_sapiens.GRCh38.95.gtf >> README
fc -ln -1 >> README
echo >> README
# 2737564 Homo_sapiens.GRCh38.95.gtf
# Set variable for testing
INPUT_GTF_LINES=`cat Homo_sapiens.GRCh38.95.gtf | wc -l`

echo "For processing the file needs to be split into 3 parts" >> README
echo "--- 1) header section" >> README
echo "--- 2) after header column 1, this will be used to update contig names" >> README
echo "--- 2) after header columns 2-end" >> README
echo >> README

# Create a header file
echo "Extract header line to temp file" >> README
grep "^#" Homo_sapiens.GRCh38.95.gtf > temp_header
fc -ln -1 >> README
echo >> README

# Create first column file
echo "Extract column 1 without the header lines to temp file" >> README
grep -v "^#" Homo_sapiens.GRCh38.95.gtf | cut -f1 > temp_c1.txt
fc -ln -1 >> README
echo >> README


# Create 2-final column file
echo "Extract column 2 and subsequent without the header lines to temp file" >> README
grep -v "^#" Homo_sapiens.GRCh38.95.gtf | cut -f2- > temp_c2plus.txt
fc -ln -1 >> README
echo >> README


# Update column 1 file with new contig names
echo "Update column 1 file with new contig names" >> README
for line in `cat ${PATH_TO_REPO}/utility_files/ensembl95_ucsc_mappings_keats.csv`
do
    OLD=`echo ${line} | cut -d, -f1`
    NEW=`echo ${line} | cut -d, -f2`
    echo "Changing:  $OLD  to  $NEW" >> README

    # update in place
    sed -i "s/\b${OLD}/${NEW}/g" temp_c1.txt
done
fc -ln -1 >> README
echo >> README


#Create final file
echo "Create final updated GTF file with UCSC contig names" >> README
paste temp_c1.txt temp_c2plus.txt > temp_new.gtf
fc -ln -1 >> README
cat temp_header temp_new.gtf > Homo_sapiens.GRCh38.95.ucsc.gtf
fc -ln -1 >> README
echo >> README

# Check final GTF line count
echo "Check Final GTF file length" >> README
wc -l Homo_sapiens.GRCh38.95.ucsc.gtf >> README
fc -ln -1 >> README
echo >> README
# 2737564 Homo_sapiens.GRCh38.95.ucsc.gtf  ## Matches starting line count!!

# Confirm the starting and final GTF have the same line count
echo "Confirm the starting and final GTF have the same line count" >> README
FINAL_GTF_LINES=`cat Homo_sapiens.GRCh38.95.ucsc.gtf | wc -l`
fc -ln -1 >> README
echo >> README

if [ ${INPUT_GTF_LINES} -eq ${FINAL_GTF_LINES} ]
then
    echo "Input and Final GTF files match"
    echo "Input and Final GTF files match"  >> README
else
    echo
    echo
    echo "####### ERROR #########"
    echo "Input and Final file lengths DO NOT match"
    echo "Input = ${INPUT_GTF_LINES}"
    echo "Final = ${FINAL_GTF_LINES}"
    echo
    echo
    exit 1
fi
echo >> README

# Clean-up directory
rm temp_*
mkdir downloads
mv Homo_sapiens.GRCh38.95.gtf downloads

# Create transcriptome fasta file derived from processed GTF
# This submission records the jobID so the next step does not start until this is complete
echo "Create transcriptome fasta file from the processed GTF for tools like Salmon" >> README
GTF_FASTA_JOBID=$(sbatch --parsable --export ALL,GENOME='../../genome_reference/GRCh38tgen_decoy.fa',GTF='Homo_sapiens.GRCh38.95.ucsc.gtf',OUTPUT='Homo_sapiens.GRCh38.95.ucsc.transcriptome.fasta' ${PATH_TO_REPO}/utility_scripts/create_transcript_fasta.slurm)
fc -ln -1 >> README
echo >> README
echo "Specific script code as follows:" >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/create_transcript_fasta.slurm >> README
echo >> README

####################################
## Generate Gene Model Specific - tool_resources
####################################

# Make gene_model specific tool_resources directory if not available
if [ -e tool_resources ]
then
    echo "Gene Model tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "Gene Model tool_resources directory NOT fount, creating and moving into it now"
    mkdir tool_resources
    cd tool_resources
fi


####################################
## Generate Salmon Index
####################################

# Make salmon index directory if not available
if [ -e "salmon_0.12.0" ]
then
    echo "Salmon directory exists, moving into it"
    cd "salmon_0.12.0"
else
    echo "Salmon directory NOT fount, creating and moving into it now"
    mkdir "salmon_0.12.0"
    cd "salmon_0.12.0"
fi

# Initialize a salmon specific README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Create the Salmon index
echo "Create salmon index to support typical paired-end seqeuncing with read lengths >=75bp" >> README
sbatch --dependency=afterok:${GTF_FASTA_JOBID} --export ALL,TRANSCRIPTOME_FASTA='../../Homo_sapiens.GRCh38.95.ucsc.transcriptome.fasta' ${PATH_TO_REPO}/utility_scripts/salmon_index.slurm
fc -ln -1 >> README
echo >> README
echo "Specific script code as follows:" >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/salmon_index.slurm >> README
echo >> README

####################################
## Generate Salmon Index
####################################

cd ..
# Make star index directory if not available
if [ -e "star_2.6.1d" ]
then
    echo "STAR directory exists, moving into it"
    cd "star_2.6.1d"
else
    echo "STAR directory NOT fount, creating and moving into it now"
    mkdir "star_2.6.1d"
    cd "star_2.6.1d"
fi

# Initialize a star specific README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Create the STAR index files
echo "Create STAR index files for 75bp read length" >> README
sbatch --export ALL,GTF='../../Homo_sapiens.GRCh38.95.ucsc.gtf',FASTA='../../../../genome_reference/GRCh38tgen_decoy.fa',SJDB_OVERHANG='74',INDEX_DIR='75bpReads' ${PATH_TO_REPO}/utility_scripts/star_index.slurm
fc -ln -1 >> README
echo >> README

echo "Create STAR index files for 100bp read length" >> README
sbatch --export ALL,GTF='../../Homo_sapiens.GRCh38.95.ucsc.gtf',FASTA='../../../../genome_reference/GRCh38tgen_decoy.fa',SJDB_OVERHANG='99',INDEX_DIR='100bpReads' ${PATH_TO_REPO}/utility_scripts/star_index.slurm
fc -ln -1 >> README
echo >> README

echo "Create STAR index files for 150bp read length" >> README
sbatch --export ALL,GTF='../../Homo_sapiens.GRCh38.95.ucsc.gtf',FASTA='../../../../genome_reference/GRCh38tgen_decoy.fa',SJDB_OVERHANG='149',INDEX_DIR='150bpReads' ${PATH_TO_REPO}/utility_scripts/star_index.slurm
fc -ln -1 >> README
echo >> README

echo >> README
echo "Specific script code as follows:" >> README
echo >> README
cat ${PATH_TO_REPO}/utility_scripts/star_index.slurm >> README
echo >> README
exit 0

####################################
## Genome Build Validations
####################################

### Compare BWAKIT and Broad Bundle versions

# BROAD_BUNDLE: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/?pli=1
# BWAKIT: https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2
# GDC: https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files

# Confirm BWAKIT and Broad Bundle are identical in sequence
diff /home/tgenref/homo_sapiens/grch38_hg38/hs38dh/genome_references/hs38DH.fa /home/tgenref/homo_sapiens/grch38_hg38/broad_resource_bundle/Homo_sapiens_assembly38.fasta
# A massive number of events are seen BUT this seems to relate to the differences in the FASTA line lengths as the broad bundle lengths are longer 101 letters versus 71 (which we get with NCBI download as well)

tr -d '\n' < /home/tgenref/homo_sapiens/grch38_hg38/hs38dh/genome_references/hs38DH.fa > hs38DH.fa
tr -d '\n' < /home/tgenref/homo_sapiens/grch38_hg38/broad_resource_bundle/Homo_sapiens_assembly38.fasta > Homo_sapiens_assembly38.fasta
diff Homo_sapiens_assembly38.fasta hs38DH.fa | wc-l
# 0
# No differences detected SO THE BWAKIT AND BROAD BUNDLE FASTA FILES ARE IDENTICAL!!

# Check contigs
grep "^>" /home/tgenref/homo_sapiens/grch38_hg38/hs38dh/genome_references/hs38DH.fa > hs38DH_contigs.txt
grep "^>" /home/tgenref/homo_sapiens/grch38_hg38/broad_resource_bundle/Homo_sapiens_assembly38.fasta > Homo_sapiens_assembly38_contigs.txt

diff hs38DH_contigs.txt Homo_sapiens_assembly38_contigs.txt
# No contig differences


### Test TGen created fasta file versus other source fasta files

# Get a contig list for comparisons
grep "^>" GRCh38tgen_decoy_alts_hla.fa > GRCh38_hs38d1_Alts_HLA_contigs.txt
diff GRCh38_hs38d1_Alts_HLA_contigs.txt hs38DH_contigs.txt
diff GRCh38_hs38d1_Alts_HLA_contigs.txt Homo_sapiens_assembly38_contigs.txt
# Lots of differences in contigs, seem to be related to the rl: annotation element
diff GRCh38_hs38d1_Alts_HLA_contigs.txt hs38DH_contigs.txt | grep "725020268"
#< >chrUn_JTFH01001998v1_decoy  AC:JTFH01001998.1  gi:725020268  LN:2001  rl:decoy  M5:35916d4135a2a9db7bc0749955d7161a  AS:hs38d1
#> >chrUn_JTFH01001998v1_decoy  AC:JTFH01001998.1  gi:725020268  LN:2001  rl:unplaced  M5:35916d4135a2a9db7bc0749955d7161a  AS:hs38d1
diff GRCh38_hs38d1_Alts_HLA_contigs.txt Homo_sapiens_assembly38_contigs.txt | grep "725020268"
#< >chrUn_JTFH01001998v1_decoy  AC:JTFH01001998.1  gi:725020268  LN:2001  rl:decoy  M5:35916d4135a2a9db7bc0749955d7161a  AS:hs38d1
#> >chrUn_JTFH01001998v1_decoy  AC:JTFH01001998.1  gi:725020268  LN:2001  rl:unplaced  M5:35916d4135a2a9db7bc0749955d7161a  AS:hs38d1

### Check versus HS38DH from BWAKIT
diff GRCh38tgen_decoy_alts_hla.fa /home/tgenref/homo_sapiens/grch38_hg38/hs38dh/genome_references/hs38DH.fa
# Differences with repeat masked lowecase and uppercase AGCT on the decoy contigs
diff GRCh38tgen_decoy_alts_hla.fa /home/tgenref/homo_sapiens/grch38_hg38/broad_resource_bundle/Homo_sapiens_assembly38.fasta
