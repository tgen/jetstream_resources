#!/bin/bash
#SBATCH --job-name="cellranger_mkref"
#SBATCH --time=0-24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 20

##########################################
## Make Custom 10x Genomics Cellranger Reference Package
##########################################

REFERENCE_RNA_GENOME_FASTA=/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy.fa
INPUT_GTF=/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/Homo_sapiens.GRCh38.98.ucsc.gtf
CELLRANGER_VERSION="3.1.0"

# Load cellranger version
module load cellranger/${CELLRANGER_VERSION}

# Filter the input GTF as recommended by 10x Genomics
cellranger mkgtf ${INPUT_GTF} ${FILTERED_GTF} \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lincRNA \
    --attribute=gene_biotype:antisense \
    --attribute=gene_biotype:IG_LV_gene \
    --attribute=gene_biotype:IG_V_gene \
    --attribute=gene_biotype:IG_V_pseudogene \
    --attribute=gene_biotype:IG_D_gene \
    --attribute=gene_biotype:IG_J_gene \
    --attribute=gene_biotype:IG_J_pseudogene \
    --attribute=gene_biotype:IG_C_gene \
    --attribute=gene_biotype:IG_C_pseudogene \
    --attribute=gene_biotype:TR_V_gene \
    --attribute=gene_biotype:TR_V_pseudogene \
    --attribute=gene_biotype:TR_D_gene \
    --attribute=gene_biotype:TR_J_gene \
    --attribute=gene_biotype:TR_J_pseudogene \
    --attribute=gene_biotype:TR_C_gene

# Error Capture
if [ "$?" = "0" ]
then
    echo "PASSED_GTF_FILTERING" >> README
else
    touch FAILED_GTF_FILTERING
    echo "FAILED_GTF_FILTERING" >> README
    exit 1
fi

# Generate the index files
cellranger mkref \
    --nthreads 20 \
    --memgb 40 \
    --genome=GRCh38_hg38tgen.98 \
    --fasta=${REFERENCE_RNA_GENOME_FASTA} \
    --genes=${FILTERED_GTF}

# Error Capture
if [ "$?" = "0" ]
then
    echo "PASSED_REFPACKAGE_GENERATION" >> README
else
    touch FAILED_REFPACKAGE_GENERATION
    echo "FAILED_REFPACKAGE_GENERATION" >> README
    exit 1
fi
