#!/usr/bin/env bash


# Usage: create_snpEff_db_index.sh <Config.ini>

### Setting as an interactive BASH session and forcing history to capture commands to a log/README file
HISTFILE=~/.bash_history
set -o history
set -ue

# Check resources.ini was provided on the command line
if [ -n "$1" ]
then
  echo "Required ini file detected"
else
  echo "Input INI file not provided, exiting due to missing requirement"
  exit 1
fi

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

# Check that the reference genome for RNA was created successfully
if [ -e RNA_FASTA_GENERATION_COMPLETE ]
then
    echo "RNA fasta exists, moving forward"
else
    echo "RNA fasta generation complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
fi

# Check gene_model directory if not available
if [ -e gene_model ]
then
    echo "Gene Model directory exists, moving into it"
    cd gene_model
else
    echo "Gene Model directory NOT found, IT IS REQUIRED, exiting as it should be created by a previous step"
    exit 1
fi

# Check specific gene model directory
if [ -e ${GENE_MODEL_NAME} ]
then
    echo "Specific Gene Model directory exists, moving into it"
    cd ${GENE_MODEL_NAME}
else
    echo "Specific Gene Model directory NOT found, IT IS REQUIRED, exiting as it should be created by a previous step"
    exit 1
fi

# Check that the required GTF was created successfully
if [ -e GENE_MODEL_GTF_GENERATION_COMPLETE ]
then
    echo "Required gene model GTF exists, moving forward"
else
    echo "Required gene model GTF DOES NOT exist, exiting"
    exit 2
fi

# Make gene_model specific tool_resources directory if not available
if [ -e tool_resources ]
then
    echo "Gene Model tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "Gene Model tool_resources directory NOT found, creating and entering it now"
    mkdir tool_resources
    cd tool_resources
fi

# Make star-fusion directory if not available
if [ -e "starFusion_${STAR_FUSION_SOURCE_VERSION}" ]
then
    echo "starFusion directory exists, moving into it"
    cd starFusion_${STAR_FUSION_SOURCE_VERSION}
else
    echo "starFusion directory NOT found, creating and moving into it now"
    mkdir starFusion_${STAR_FUSION_SOURCE_VERSION}
    cd starFusion_${STAR_FUSION_SOURCE_VERSION}
fi

####################################
## Download Star-Fusion Plug and play Resource File
####################################

# Initialize a snpEff specific README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Download the full plug-n-play resourece bundle
echo "Download STAR-Fusion Source MD5" >> README
wget ${STAR_FUSION_SOURCE_MD5}
fc -ln -1 >> README
echo >> README

echo "Download STAR-Fusion source bundle" >> README
wget ${STAR_FUSION_SOURCE}
fc -ln -1 >> README
echo >> README

# Capture the md5sum filename
echo "Capture the md5sum filename" >> README
MD5_FILENAME=`basename ${STAR_FUSION_SOURCE_MD5}`
fc -ln -1 >> README
echo >> README

# Check the MD5 checksum matches and capture error if it doesn't
echo "Check the MD5 checksum matches" >> README
echo "     md5sum --check ${MD5_FILENAME}" >> README
md5sum --check ${MD5_FILENAME}

if [ $? -eq 0 ]
then
   echo "Checksum validated" >> README
   touch CHECKSUM_VALIDATED
else
    echo "FAILED Checksum Validation" >> README
    touch FAILED_CHECKSUM_VALIDATION
    exit 1
fi

# Capture the buncle tar filename
echo "Capture the md5sum filename" >> README
BUNDLE_FILENAME=`basename ${STAR_FUSION_SOURCE}`
fc -ln -1 >> README
echo >> README

# Decompress the downloaded plug-n-play file
echo "Decompress the downloaded bundle" >> README
echo "    tar xvzf ${BUNDLE_FILENAME}" >> README
tar xvzf ${BUNDLE_FILENAME}

if [ $? -eq 0 ]
then
   echo "Tar Archive Extraction Complete" >> README
   touch TAR_ARCHIVE_EXTRACTED
else
    echo "FAILED Tar Archive Extraction" >> README
    touch FAILED_CHECKSUM_VALIDATION
    exit 1
fi

# enter the decompressed folder
BUNDLE_FOLDER=`basename ${BUNDLE_FILENAME} ".tar.gz"`
cd ${BUNDLE_FOLDER}

# Create symbolic link
REFERENCE_RNA_GENOME_FASTA=${TOPLEVEL_DIR}/genome_reference/${REFERENCE_RNA_GENOME_NAME}
ln -s ${REFERENCE_RNA_GENOME_FASTA} ${REFERENCE_RNA_GENOME_NAME}

# We will use the GTF provided in the bundle, see below
'''
Comparing the gtf downloaded in the "GRCh38_gencode_v32_CTAT_lib_Dec062019.source/gencode.v32.annotation.gtf" the GTF
name (gencode.v32.annotation.gtf), number of lines, and diff matches the file downloaded via (ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz)
This file content is one of 3 in the Content = Comprehensive Gene Annotation
Description= It contains the comprehensive gene annotations on the reference chromsomes only, This is the "Main Annotation File" for most users
The star-fusion download includes a fasta "GRCh38.primary_assembly.genome.fa" with chr1-22, X,Y, M, and GL/KI supercontigs, total 194 contigs
  - GL format --> GL000008.2
  - KI format --> KI270302.1

Our RNA reference genome is chr1-22, X, Y, M, GL/KI supercontigs, KN/JTFH decoys, and chrEBV, total xxx contigs. 2580 contigs
  - GL format --> chrUn_GL000195v1, chr4_GL000008v2_random
  - KI format --> chrUn_KI270302v1, chr1_KI270706v1_random
  - KN format --> chrUn_KN707606v1_decoy
  - JTFH format --> chrUn_JTFH01000001v1_decoy

gencode.v32.annotation.gtf
  - chr1-22, X, Y, M
gencode.v32.primary_assembly.annotation.gtf
  - chr1-22, X, Y, M, GL000009.2, KI270442.1

Therefore, we can use the same file as Star-fusion recommends and provides in the source download with our geneome as
all contigs match
'''

####################################
## Build Annotations from Source
###################################

# A specific dfamscan perl script is needed for the build and needs to be added to the $PATH
# Get dfamscan.pl
wget https://www.dfam.org/releases/current/infrastructure/dfamscan.pl.gz
gunzip dfamscan.pl.gz
chmod +x dfamscan.pl

# Add current DIR to path
CURRENT_DIR=`pwd`
export PATH=$PATH:$CURRENT_DIR

if [ ${ENVIRONMENT} == "TGen" ]
then
  # Load required modules
  module load STAR-Fusion/1.8.1-GCC-8.2.0-2.31.1-Perl-5.28.1-Python-3.7.2
  module load blast/2.7.1
  module load hmmer/3.2.1

  # Use provided starFusion build script
  /packages/easybuild/software/STAR-Fusion/1.8.1-GCC-8.2.0-2.31.1-Perl-5.28.1-Python-3.7.2/ctat-genome-lib-builder/prep_genome_lib.pl \
  --CPU 20 \
  --max_readlength 150 \
  --genome_fa GRCh38tgen_decoy.fa \
  --gtf gencode.v32.annotation.gtf \
  --fusion_annot_lib fusion_lib.*.dat.gz \
  --annot_filter_rule AnnotFilterRule.pm \
  --pfam_db current \
  --dfam_db human \
  --human_gencode_filter
elif [ ${ENVIRONMENT} == "LOCAL" ]
then
  echo
  echo "Assuming required tools are available in $PATH"
  # Ensure the starFusion repository is available on your system and the subfolder "ctat-genome-lib-builder" is available in the $PATH
  prep_genome_lib.pl \
  --CPU 20 \
  --max_readlength 150 \
  --genome_fa GRCh38tgen_decoy.fa \
  --gtf gencode.v32.annotation.gtf \
  --fusion_annot_lib fusion_lib.*.dat.gz \
  --annot_filter_rule AnnotFilterRule.pm \
  --pfam_db current \
  --dfam_db human \
  --human_gencode_filter
  echo
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Enviroment Variable"
  echo "Only TGen or LOCAL are supported"
  exit 1
fi
