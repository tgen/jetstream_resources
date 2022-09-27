#!/usr/bin/env bash

# Automated Script to download and build Topmed BCF annotation file for usage in Tempe workflow

# Usage: ./build_topmed.sh pipeline_resouces.ini

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
## Load Required Tools
###################################
if [ ${ENVIRONMENT} == "TGen" ]
then
  module load BCFtools/1.10.1-foss-2019a
elif [ ${ENVIRONMENT} == "LOCAL" ]
then
  echo
  echo "Assuming required tools are available in $PATH"
  echo
else
  echo "Unexpected Entry in ${WORKFLOW_NAME}_resources.ini Environment Variable"
  echo "Only TGen or LOCAL are supported"
  exit 1
fi

####################################
## Create Expected Folder Structure
###################################

# Make top level directory if not available
if [ -e ${PARENT_DIR} ]
then
    echo "Parent directory: ${PARENT_DIR} exists, moving into it"
    cd ${PARENT_DIR}
else
    echo "Parent directory NOT fount, creating and moving into it now"
    mkdir -p ${PARENT_DIR}
    cd ${PARENT_DIR}
fi

# Make public_databases folder if not available
if [ -e public_databases ]
then
    echo "Public Databases folder exists, moving into it"
    cd public_databases
else
    echo "Public Databases folder NOT fount, creating and moving into it now"
    mkdir -p public_databases
    cd public_databases
fi

# Make topmed folder if not available
if [ -e topmed ]
then
    echo "topmed folder exists, moving into it"
    cd topmed
else
    echo "topmed folder NOT fount, creating and moving into it now"
    mkdir -p topmed
    cd topmed
fi

# Make topmed release version folder if not available
if [ -e ${TOPMED_VERSION} ]
then
    echo "topmed ${TOPMED_VERSION} folder exists, exiting to prevent overwrite"
    exit 1
else
    echo "topmed ${TOPMED_VERSION} folder NOT fount, creating and moving into it now"
    mkdir -p ${TOPMED_VERSION}
    cd ${TOPMED_VERSION}
fi

####################################
## Parameterized Code
####################################

# Initialize a topmed README
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README

# Download the individual VCF files uploaded as part of freeze8
# STEP 1 - Download freeze 8 VCF files uploaded to dbSNP (https://bravo.sph.umich.edu/freeze8/hg38/downloads) for each chromosome

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/1' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr1.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/2' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr2.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/3' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr3.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/4' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr4.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/5' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr5.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/6' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr6.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/7' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr7.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/8' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr8.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/9' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr9.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/10' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr10.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/11' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr11.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/12' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr12.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/13' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr13.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/14' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr14.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/15' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr15.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/16' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr16.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/17' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr17.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/18' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr18.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/19' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr19.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/20' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr20.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/21' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr21.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/22' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chr22.BRAVO_TOPMed_Freeze_8.vcf.gz

curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/X' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: remember_token=jkeats@tgen.org|937c7fdbef17b651254e2f32cc7101ab058e8efe3d92a3cdfc08eb0edb8c1f426b9eca3d97642e4d213f6ff295cc4768112f907dfe1b3714cfc78aa5a3c61e64;session=.eJwdT8tqxDAQ-xefw2JPxvY4p0JLKKWkpYfuMTie8T7a7ELsPSyl_15TBNJFmpF-1Jw3KUc11O0mnZpPrAblenBWg42YyXmXE5GxOvbGQGhAJit-8aZvAosAOtShR0PgXPY2puzAU7BeAjAjp5B8bF6GYGnRicERLCYkyInR2QyEljNjNMk0Vp3aZJV1kW0ukq4XLmqg9kTvdKdKjVVay_dXj9P0bffHlZ7G8x1Pn3J_GcePt8f989Ru3ErL_w86f0ms5aEe5LK7bgf1-wemyEh3.YzJMMA.R_TEfWWO_97CbXnECyDEq-WFxsI;' --compressed > chrX.BRAVO_TOPMed_Freeze_8.vcf.gz


# To make a large annoation file the downloads need an index so that header contigs can be added by bcftools concat

for line in `ls *.BRAVO_TOPMed_Freeze_8.vcf.gz`
do

bcftools index --tbi --threads 4 ${line}

done

# Now make large annotation file

bcftools concat --threads 4 -O b -o All.BRAVO_TOPMed_Freeze_8.bcf chr1.BRAVO_TOPMed_Freeze_8.vcf.gz chr2.BRAVO_TOPMed_Freeze_8.vcf.gz chr3.BRAVO_TOPMed_Freeze_8.vcf.gz chr4.BRAVO_TOPMed_Freeze_8.vcf.gz chr5.BRAVO_TOPMed_Freeze_8.vcf.gz chr6.BRAVO_TOPMed_Freeze_8.vcf.gz chr7.BRAVO_TOPMed_Freeze_8.vcf.gz chr8.BRAVO_TOPMed_Freeze_8.vcf.gz chr9.BRAVO_TOPMed_Freeze_8.vcf.gz chr10.BRAVO_TOPMed_Freeze_8.vcf.gz chr11.BRAVO_TOPMed_Freeze_8.vcf.gz chr12.BRAVO_TOPMed_Freeze_8.vcf.gz chr13.BRAVO_TOPMed_Freeze_8.vcf.gz chr14.BRAVO_TOPMed_Freeze_8.vcf.gz chr15.BRAVO_TOPMed_Freeze_8.vcf.gz chr16.BRAVO_TOPMed_Freeze_8.vcf.gz chr17.BRAVO_TOPMed_Freeze_8.vcf.gz chr18.BRAVO_TOPMed_Freeze_8.vcf.gz chr19.BRAVO_TOPMed_Freeze_8.vcf.gz chr20.BRAVO_TOPMed_Freeze_8.vcf.gz chr21.BRAVO_TOPMed_Freeze_8.vcf.gz chr22.BRAVO_TOPMed_Freeze_8.vcf.gz chrX.BRAVO_TOPMed_Freeze_8.vcf.gz

# Make CSI index
bcftools index --threads 4 All.BRAVO_TOPMed_Freeze_8.bcf

# Generate a stats file for comparisons

bcftools stats --threads 4 All.BRAVO_TOPMed_Freeze_8.bcf > All.BRAVO_TOPMed_Freeze_8.stats

# Cleanup
rm *.BRAVO_TOPMed_Freeze_8.vcf.gz
rm *.BRAVO_TOPMed_Freeze_8.vcf.gz.tbi

echo "Process Complete"
