#!/bin/bash
#SBATCH -c 16

# Usage: make_index.sh <sjdbOverhang>

set -ux
module load STAR/2.6.1d

GTF=../../Homo_sapiens.GRCh38.95.ucsc.gtf.gz
FASTA=../../../../Homo_sapiens.GRCh38.noalt.fasta
SJDB_OVERHANG="${1}"
INDEX_DIR=./GRCh38.95_${1}bp/

if [ -d "${INDEX_DIR}" ]; then
  echo "Index already exists: ${INDEX_DIR}"
  exit 1
else
  mkdir -p "${INDEX_DIR}"
fi

GTF_TEMP=$(mktemp)
zcat "${GTF}" > "${GTF_TEMP}"

STAR \
  --runMode genomeGenerate \
  --outFileNamePrefix "${INDEX_DIR}" \
  --genomeDir "${INDEX_DIR}" \
  --runThreadN 16 \
  --sjdbOverhang "${SJDB_OVERHANG}" \
  --genomeFastaFiles "${FASTA}" \
  --sjdbGTFfile "${GTF_TEMP}"

EXIT_CODE=$?
rm "${GTF_TEMP}"
exit $EXIT_CODE