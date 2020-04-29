#!/usr/bin/env bash

# Usage: create_genome_reference.sh <Config.ini>

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

# Check that the reference genome was created successfully
if [ -e GENOME_FASTA_GENERATION_COMPLETE ]
then
    echo "Genome fasta exists, moving forward"
else
    echo "Genome fasta generation complete flag NOT found"
    echo "Try again later as this is required"
    exit 2
fi

if [ -e tool_resources ]
then
    echo "tool_resources directory exists, moving into it"
    cd tool_resources
else
    echo "tool_resources directory NOT found, creating and moving into it now"
    mkdir tool_resources
    cd tool_resources
fi

if [ -e "gatk_cnv" ]
then
    echo "The gatk_cnv directory exists, exiting to prevent overwriting existing interval lists."
    exit 2
else
    echo "The gatk_cnv directory was NOT found, creating and moving into it now"
    mkdir gatk_cnv
    cd gatk_cnv
fi

####################################
## Initialize a gatk_cnv README
####################################
touch README
echo >> README
echo "For details on file creation see the associated github repository:" >> README
echo "https://github.com/tgen/jetstream_resources/${WORKFLOW_NAME}" >> README
echo "Created and downloaded by ${CREATOR}" >> README
date >> README
echo >> README
echo "gatk cnv interval_list creation details:" >> README
echo "The interval lists are generated by essentially concatenating the reference.dict with the primary contig intervals" >> README
echo "These intervals are grabbed from the SN and LN of the reference.dict and then we apply the interval_list formatting" >> README
echo -e "e.g. SN\t1\tLN\t+\tSN" >> README
echo >> README

# Determine the full path to the reference genome dictionary file
GENOME_FASTA_BASENAME=`basename ${REFERENCE_DNA_GENOME_NAME} ".fa"`
GENOME_FASTA_DICT=${TOPLEVEL_DIR}/genome_reference/${GENOME_FASTA_BASENAME}.dict

# Place the reference dictionary as the header to the interval list
cat ${GENOME_FASTA_DICT} > ${GENOME_BUILD}.contigs.interval_list_header

# We only care about the sequence name and the length from the dictionary
# Take these values and apply interval_list formatting
cut -f2,3 ${GENOME_FASTA_DICT} \
  | \
  grep -v "VN" \
  | \
  sed 's/SN://g' \
  | \
  sed 's/LN://g' \
  | \
  awk '{ OFS="\t"; print $1,"1",$2,"+",$1 }' >> ${GENOME_BUILD}.contigs.interval_list

# Reducing the interval list to only include the primary contigs, to be generalizable, using a
# last primary contig variable to select the last contig that should be used. This contig will
# be removed in the female interval list
while read -r line; do
    contig_name=`echo ${line} | awk '{ print $1 }'`
    if [[ "${contig_name}" == "${LAST_PRIMARY_CONTIG}" ]]; then
        echo "$line" >> ${GENOME_BUILD}.primary.contigs.interval_list
        break
    else
        echo "$line" >> ${GENOME_BUILD}.primary.contigs.interval_list
    fi
done < ${GENOME_BUILD}.contigs.interval_list

# Creating the final interval lists
cat ${GENOME_BUILD}.contigs.interval_list_header ${GENOME_BUILD}.primary.contigs.interval_list \
    > ${GENOME_BUILD}.primary.contigs.male.interval_list

# Removing last line of male interval_list to create female
sed '$d' ${GENOME_BUILD}.primary.contigs.male.interval_list > ${GENOME_BUILD}.primary.contigs.female.interval_list

# Cleanup files that are no longer needed
rm ${GENOME_BUILD}.primary.contigs.interval_list ${GENOME_BUILD}.contigs.interval_list_header ${GENOME_BUILD}.contigs.interval_list