#!/usr/bin/env bash


# Usage: create_snpEff_db_index.sh <Config.ini>

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

####################################
## Download Star-Fusion Plug and play Resource File
####################################

# Make star-fusion directory if not available
if [ -e "starFusion_${STARFUSION_PnP_ANNOTATION_VERSION}" ]
then
    echo "starFusion directory exists, moving into it"
    cd starFusion_${STARFUSION_PnP_ANNOTATION_VERSION}
else
    echo "starFusion directory NOT fount, creating and moving into it now"
    mkdir starFusion_${STARFUSION_PnP_ANNOTATION_VERSION}
    cd starFusion_${STARFUSION_PnP_ANNOTATION_VERSION}
fi

# Initialize a snpEff specific README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/phoenix" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Download the full plug-n-play resourece bundle
echo "Download STAR-Fusion plug-n-play MD5" >> README
wget ${STARFUSION_PnP_MD5_DOWNLOAD_LINK}
fc -ln -1 >> README
echo >> README

echo "Download STAR-Fusion plug-n-play bundle" >> README
wget ${STARFUSION_PnP_DOWNLOAD_LINK}
fc -ln -1 >> README
echo >> README

# Capture the md5sum filename
echo "Capture the md5sum filename" >> README
MD5_FILENAME=`basename ${STARFUSION_PnP_MD5_DOWNLOAD_LINK}`
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
BUNDLE_FILENAME=`basename ${STARFUSION_PnP_DOWNLOAD_LINK}`
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