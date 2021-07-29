#!/bin/bash

set -o errexit -o pipefail -o noclobber

SEX=Female
REFERENCE=/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa
INTERVALS=/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/tool_resources/gatk_cnv/Homo_sapiens.GRCh38.primary.contigs.female.interval_list
ENCODE_DENY_LIST=/home/tgenref/homo_sapiens/grch38_hg38/public_databases/ncbi/Encode_deny_list_with_ncbi_centromere.merged.bed
GNOMAD=/home/tgenref/homo_sapiens/grch38_hg38/public_databases/gnomad/r3.0/gnomad.genomes.r3.0.sites.pass.AnnotationReference.bcf
OUTPUT_DIR=${PWD}
NAME=""
BIN_LENGTH=10000

# Reads each argument provided after calling the script
num_args_rem=$#
if [[ $num_args_rem == 0 ]] ; then
    printf  "
            Available arguments:
                    -i | --intervals        - Path to a interval list file, defaults to $INTERVALS
                    -d | --denylist         - Path to the denylist, defaults to $ENCODE_DENY_LIST
                    -g | --gnomad           - Path GnomAD Genome, defaults to ${GNOMAD}
                    -r | --reference        - Path to reference.fa, defaults to ${REFERENCE}
                    -l | --bin_length       - Bin Length, defaults to $BIN_LENGTH
                    -o | --out_dir          - Path to output directory, defaults to pwd
                    -n | --name             - Name/prefix of the output files
                    -s | --sex              - Sex of panel, default to $SEX
                    -b | --bams             - File contain bam file location(s) - They do not need to be in the same directory\n\n"
    exit
fi

# Grabbing the remaining arguments to pass on to operations
while [[ $num_args_rem -ne 0 ]] ; do
    case "$1" in
        -i|--intervals)
            INTERVALS="$2"
            shift 2 ;;
        -d|--denylist)
            ENCODE_DENY_LIST="$2"
            shift 2 ;;
        -g|--gnomad)
            GNOMAD="$2"
            shift 2 ;;
        -r|--reference)
            REFERENCE="$2"
            shift 2 ;;
        -l|--bin-length)
            BIN_LENGTH="$2"
            shift 2 ;;
        -t|--out_dir)
            OUTPUT_DIR="$2"
            shift 2 ;;
        -n|--name)
            NAME="$2"
            shift 2 ;;
        -s|--sex)
            SEX="$2"
            shift 2 ;;
        -b|--bams)
            BAMS="$2"
            shift 2 ;;
    esac
    num_args_rem=$((num_args_rem-1))
done

# Preparing input lists
hdf5_file_list=""
echo "Using the following bam files:"
while read -r bam; do
  if [[ -f "${bam}" ]]; then
    echo $bam
  else
    echo "Error - $bam not found"
    exit 1
  fi
  normal_name=$(basename ${bam} | cut -d. -f1 )
  hdf5_file_list+="--input ${OUTPUT_DIR}/${NAME}_hdf5_counts/${normal_name}.filter.counts.hdf5 "
done < ${BAMS}

echo "Starting analysis using the following arguments:"
echo "SEX=$SEX"
echo "REFERENCE=$REFERENCE"
echo "INTERVALS=$INTERVALS"
echo "ENCODE_DENY_LIST=$ENCODE_DENY_LIST"
echo "GNOMAD=$GNOMAD"
echo "OUTPUT_DIR=${OUTPUT_DIR}"
echo "NAME=${NAME}"
echo "BIN_LENGTH=$BIN_LENGTH"
echo

# TODO check if tools already exist
module load GATK/4.1.8.0-GCCcore-8.3.0-Java-1.8
module load BCFtools/1.10.2-GCC-8.2.0-2.31.1
module load R/3.6.1-phoenix
module load samtools/1.9

mkdir -p ${OUTPUT_DIR}/temp

gatk --java-options "-Xmx12g" PreprocessIntervals \
  --tmp-dir ${OUTPUT_DIR}/temp/ \
  --intervals ${INTERVALS} \
  --bin-length ${BIN_LENGTH} \
  --reference ${REFERENCE} \
  --output ${OUTPUT_DIR}/${NAME}_preprocessed.interval_list \
  --interval-merging-rule OVERLAPPING_ONLY

gatk --java-options "-Xmx12g" AnnotateIntervals \
    --tmp-dir ${OUTPUT_DIR}/temp/ \
    --mappability-track /home/tgenref/homo_sapiens/grch38_hg38/public_databases/bismap/k100.umap.no_header.bed \
    --intervals ${OUTPUT_DIR}/${NAME}_preprocessed.interval_list \
    --output ${OUTPUT_DIR}/${NAME}_PreFilter_anno_preprocessed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    --reference ${REFERENCE}

gatk --java-options "-Xmx12g" FilterIntervals \
    --tmp-dir ${OUTPUT_DIR}/temp/ \
    --intervals ${OUTPUT_DIR}/${NAME}_preprocessed.interval_list \
    --exclude-intervals ${ENCODE_DENY_LIST} \
    --interval-merging-rule OVERLAPPING_ONLY \
    --minimum-gc-content 0.1 \
    --maximum-gc-content 0.9 \
    --minimum-mappability 0.9 \
    --maximum-mappability 1.0 \
    --annotated-intervals ${OUTPUT_DIR}/${NAME}_PreFilter_anno_preprocessed.interval_list \
    --output ${OUTPUT_DIR}/${NAME}_preprocessed_filt_map.interval_list

gatk --java-options "-Xmx12g" AnnotateIntervals \
    --tmp-dir ${OUTPUT_DIR}/temp/ \
    --intervals ${OUTPUT_DIR}/${NAME}_preprocessed_filt_map.interval_list \
    --output ${OUTPUT_DIR}/${NAME}_anno_preprocessed_filt_map.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    --reference ${REFERENCE}

mkdir -p ${OUTPUT_DIR}/${NAME}_hdf5_counts

while read -r bam; do
  normal_name=$(basename ${bam} | cut -d. -f1 )
  mkdir ${OUTPUT_DIR}/${NAME}_temp_${normal_name}

  gatk --java-options "-Xmx12g" CollectReadCounts \
    --tmp-dir ${OUTPUT_DIR}/${NAME}_temp_${normal_name}/ \
    --input ${bam} \
    --intervals ${OUTPUT_DIR}/${NAME}_preprocessed_filt_map.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    --read-filter FirstOfPairReadFilter \
    --output ${OUTPUT_DIR}/${NAME}_hdf5_counts/${normal_name}.filter.counts.hdf5
done < ${BAMS}

# After getting all of the counts now we need to filter the intervals that have low or extream count
# This needs to be done by hand to get all of the inputs

gatk --java-options "-Xmx12g" FilterIntervals \
  --tmp-dir ${OUTPUT_DIR}/temp/ \
  --intervals ${OUTPUT_DIR}/${NAME}_preprocessed_filt_map.interval_list \
  --output ${OUTPUT_DIR}/${NAME}_preprocessed_filt_map_counts.interval_list \
  --interval-merging-rule OVERLAPPING_ONLY \
  $hdf5_file_list

gatk --java-options "-Xmx12g" AnnotateIntervals \
    --tmp-dir ${OUTPUT_DIR}/temp/ \
    --intervals ${OUTPUT_DIR}/${NAME}_preprocessed_filt_map_counts.interval_list \
    --output ${OUTPUT_DIR}/${NAME}_anno_preprocessed_filt_map_counts.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    --reference ${REFERENCE}

mv ${OUTPUT_DIR}/${NAME}_hdf5_counts ${OUTPUT_DIR}/${NAME}_hdf5_counts_preFilt
mkdir ${OUTPUT_DIR}/${NAME}_hdf5_counts

while read -r bam; do
  normal_name=$(basename ${bam} | cut -d. -f1 )
  mkdir -p ${OUTPUT_DIR}/${NAME}_temp_${normal_name} || true

  gatk --java-options "-Xmx12g" CollectReadCounts \
    --tmp-dir ${OUTPUT_DIR}/${NAME}_temp_${normal_name}/ \
    --input ${bam} \
    --intervals ${OUTPUT_DIR}/${NAME}_preprocessed_filt_map_counts.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    --read-filter FirstOfPairReadFilter \
    --output ${OUTPUT_DIR}/${NAME}_hdf5_counts/${normal_name}.filter.counts.hdf5
done < ${BAMS}

bcftools view \
    --exclude-types indels,mnps,ref,bnd,other \
    --genotype het \
    --targets ^chrX:10001-2781479,chrX:155701383-156030895,chrY:10001-2781479,chrY:56887903-57217415 \
    --exclude "AF<0.2 | AF>0.8" \
    --output-type v \
    --output-file ${OUTPUT_DIR}/${NAME}_GNOMAD.preFilt.snps.vcf \
    ${GNOMAD}

bcftools view \
    --targets-file ^${ENCODE_DENY_LIST} \
    --output-type v \
    --output-file ${OUTPUT_DIR}/${NAME}_filteredSnps.vcf \
    ${OUTPUT_DIR}/${NAME}_GNOMAD.preFilt.snps.vcf

rm ${OUTPUT_DIR}/${NAME}_GNOMAD.preFilt.snps.vcf

# create pon
mkdir -p "${OUTPUT_DIR}/${NAME}_temp3"

gatk --java-options "-Xmx12g" CreateReadCountPanelOfNormals \
    --tmp-dir ${OUTPUT_DIR}/${NAME}_temp3/ \
    --annotated-intervals ${OUTPUT_DIR}/${NAME}_anno_preprocessed_filt_map_counts.interval_list \
    --minimum-interval-median-percentile 10.0 \
    --output ${OUTPUT_DIR}/${NAME}_cnv.pon.hdf5 \
    $hdf5_file_list
