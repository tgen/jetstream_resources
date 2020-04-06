#!/usr/bin/env bash

# Usage: create_broad_resource_bundle.sh <Config.ini>

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

# Check if the tool resources directory exists
if [ -e public_databases ]
then
    echo "public_databases directory exists, moving into it"
    cd public_databases
else
    echo "public_databases directory NOT found, creating and moving into it now"
    mkdir public_databases
    cd public_databases
fi

if [ -e "broad_resource_bundle" ]
then
    echo "The broad_resource_bundle directory exists, exiting to prevent overwriting existing index"
    exit 2
else
    echo "The broad_resource_bundle directory was NOT found, creating and moving into it now"
    mkdir broad_resource_bundle
    cd broad_resource_bundle
fi

touch README
echo >> README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README


####################################
## Download broad resource bundle files
####################################

echo "####################################" >> README
echo "## Download broad resource bundle files " >> README
echo "####################################" >> README
echo >> README

echo "wget ${BROAD_BUNDLE_HAPMAP_DOWNLOAD_LINK}" >> README
echo "wget ${BROAD_BUNDLE_HAPMAP_DOWNLOAD_MD5SUM_LINK}" >> README
echo >> README
wget ${BROAD_BUNDLE_HAPMAP_DOWNLOAD_LINK}
wget ${BROAD_BUNDLE_HAPMAP_DOWNLOAD_MD5SUM_LINK}

HAPMAP=`basename ${BROAD_BUNDLE_HAPMAP_DOWNLOAD_LINK}`
if [ `md5sum ${HAPMAP} | cut -d " " -f 1` != `cut -d " " -f 1 ${HAPMAP}.md5` ]; then
    echo "The md5s do not match for ${HAPMAP}, please run the script to again."
    exit 1
else
    echo "md5sum ${HAPMAP} | cut -d " " -f 1" >> README
    md5sum ${HAPMAP} | cut -d " " -f 1 >> README
    md5sum ${HAPMAP} | cut -d " " -f 1
    echo "cut -d " " -f 1 ${HAPMAP}.md5" >> README
    cut -d " " -f 1 ${HAPMAP}.md5 >> README
    cut -d " " -f 1 ${HAPMAP}.md5
fi
echo >> README

echo "wget ${BROAD_BUNDLE_MILLS_DOWNLOAD_LINK}" >> README
echo "wget ${BROAD_BUNDLE_MILLS_MD5SUM_DOWNLOAD_LINK}" >> README
echo >> README
wget ${BROAD_BUNDLE_MILLS_DOWNLOAD_LINK}
wget ${BROAD_BUNDLE_MILLS_MD5SUM_DOWNLOAD_LINK}

MILLS=`basename ${BROAD_BUNDLE_MILLS_DOWNLOAD_LINK}`
if [ `md5sum ${MILLS} | cut -d " " -f 1` != `cut -d " " -f 1 ${MILLS}.md5` ]; then
    echo "The md5s do not match for ${MILLS}, please run the script to again."
    exit 1
else
    echo "md5sum ${MILLS} | cut -d " " -f 1" >> README
    md5sum ${MILLS} | cut -d " " -f 1 >> README
    md5sum ${MILLS} | cut -d " " -f 1
    echo "cut -d " " -f 1 ${MILLS}.md5" >> README
    cut -d " " -f 1 ${MILLS}.md5 >> README
    cut -d " " -f 1 ${MILLS}.md5
fi
echo >> README

echo "wget ${BROAD_BUNDLE_1000G_DOWNLOAD_LINK}" >> README
echo "wget ${BROAD_BUNDLE_1000G_DOWNLOAD_MD5SUM_LINK}" >> README
echo >> README
wget ${BROAD_BUNDLE_1000G_DOWNLOAD_LINK}
wget ${BROAD_BUNDLE_1000G_DOWNLOAD_MD5SUM_LINK}

G1000=`basename ${BROAD_BUNDLE_1000G_DOWNLOAD_LINK}`
if [ `md5sum ${G1000} | cut -d " " -f 1` != `cut -d " " -f 1 ${G1000}.md5` ]; then
    echo "The md5s do not match for ${G1000}, please run the script to again."
    exit 1
else
    echo "md5sum ${G1000} | cut -d " " -f 1" >> README
    md5sum ${G1000} | cut -d " " -f 1 >> README
    md5sum ${G1000} | cut -d " " -f 1
    echo "cut -d " " -f 1 ${G1000}.md5" >> README
    cut -d " " -f 1 ${G1000}.md5 >> README
    cut -d " " -f 1 ${G1000}.md5
fi
echo >> README

echo "wget ${BROAD_BUNDLE_1000Gphase3_DOWNLOAD_LINK}"
wget ${BROAD_BUNDLE_1000Gphase3_DOWNLOAD_LINK}
