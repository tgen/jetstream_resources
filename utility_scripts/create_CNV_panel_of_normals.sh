#!/bin/bash

# TODO check if tools already exist
module load GATK/4.1.8.0-GCCcore-8.3.0-Java-1.8
module load BCFtools/1.10.2-GCC-8.2.0-2.31.1
module load R/3.6.1-phoenix
module load samtools/1.9


# TODO user cli options

SEX=Female
INTERVALS=/home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/tool_resources/gatk_cnv/Homo_sapiens.GRCh38.primary.contigs.female.interval_list
ENCODE_DENY_LIST=/home/tgenref/homo_sapiens/grch38_hg38/public_databases/ncbi/Encode_deny_list_with_ncbi_centromere.merged.bed

temp_dir=test_001
# TODO Update for dynamic bin length
BIN_LENGTH=10000

mkdir -p ${temp_dir}/temp

gatk --java-options "-Xmx12g" PreprocessIntervals \
	--tmp-dir ${temp_dir}/temp/ \
	--intervals ${INTERVALS} \
	--bin-length ${BIN_LENGTH} \
	--reference /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa \
	--output ${temp_dir}/preprocessed.interval_list \
	--interval-merging-rule OVERLAPPING_ONLY

gatk --java-options "-Xmx12g" AnnotateIntervals \
    --tmp-dir ${temp_dir}/temp/ \
    --mappability-track /home/tgenref/homo_sapiens/grch38_hg38/public_databases/bismap/k100.umap.no_header.bed \
    --intervals ${temp_dir}/preprocessed.interval_list \
    --output ${temp_dir}/PreFilter_anno_preprocessed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    --reference /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa

gatk --java-options "-Xmx12g" FilterIntervals \
    --tmp-dir ${temp_dir}/temp/ \
    --intervals ${temp_dir}/preprocessed.interval_list \
	--exclude-intervals ${ENCODE_DENY_LIST} \
    --interval-merging-rule OVERLAPPING_ONLY \
    --minimum-gc-content 0.1 \
    --maximum-gc-content 0.9 \
    --minimum-mappability 0.9 \
    --maximum-mappability 1.0 \
    --annotated-intervals ${temp_dir}/PreFilter_anno_preprocessed.interval_list \
    --output ${temp_dir}/preprocessed_filt_map.interval_list

gatk --java-options "-Xmx12g" AnnotateIntervals \
    --tmp-dir ${temp_dir}/temp/ \
    --intervals ${temp_dir}/preprocessed_filt_map.interval_list \
    --output ${temp_dir}/anno_preprocessed_filt_map.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    --reference /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa


mkdir -p ${temp_dir}/hdf5_counts

for bam in $(ls bams/*.bam)
do
	normal_name=$(echo ${bam} | cut -d. -f1 | cut -d/ -f2 )

	mkdir ${temp_dir}/temp_${normal_name}

	gatk --java-options "-Xmx12g" CollectReadCounts \
		--tmp-dir ${temp_dir}/temp_${normal_name}/ \
		--input bams/${normal_name}.bwa.bam \
		--intervals ${temp_dir}/preprocessed_filt_map.interval_list \
		--interval-merging-rule OVERLAPPING_ONLY \
		--read-filter FirstOfPairReadFilter \
		--output ${temp_dir}/hdf5_counts/${normal_name}.filter.counts.hdf5
done

for bam in $(ls bams/*.bam)
do
	normal_name=$(echo ${bam} | cut -d. -f1 | cut -d/ -f2 )
	mkdir ${temp_dir}/temp_${normal_name}

	sbatch --cpus-per-task 10 /home/achristofferson/toolkits/toolkit_achristofferson/gatk_cna_pipeline/generate_CNV_panel_of_normals_CollectReadCounts.sh \
		${normal_name} \
		${temp_dir} \
		preprocessed_filt_map.interval_list
done

mkdir slurm_outs
mv slurm-* slurm_outs/

# After getting all of the counts now we need to filter the intervals that have low or extream count
# This needs to be done by hand to get all of the inputs

gatk --java-options "-Xmx12g" FilterIntervals \
	--tmp-dir ${temp_dir}/temp/ \
	--intervals ${temp_dir}/preprocessed_filt_map.interval_list \
	--output ${temp_dir}/preprocessed_filt_map_counts.interval_list \
	--interval-merging-rule OVERLAPPING_ONLY \
  # TODO dynamic number of inputs
	--input test_001/hdh5_counts/GM14632_CORIELL_p0_CL_Whole_C1_TWGSP.filter.counts.hdf5

gatk --java-options "-Xmx12g" AnnotateIntervals \
    --tmp-dir ${temp_dir}/temp/ \
    --intervals ${temp_dir}/preprocessed_filt_map_counts.interval_list \
    --output ${temp_dir}/anno_preprocessed_filt_map_counts.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    --reference /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa

mv test_001/hdh5_counts test_001/hdh5_counts_preFilt
mkdir test_001/hdh5_counts

for bam in $(ls bams/*.bam)
do
    normal_name=$(echo ${bam} | cut -d. -f1 | cut -d/ -f2 )
    mkdir ${temp_dir}/temp_${normal_name}

    sbatch --cpus-per-task 10 /home/achristofferson/toolkits/toolkit_achristofferson/gatk_cna_pipeline/generate_CNV_panel_of_normals_CollectReadCounts.sh \
        ${normal_name} \
        ${temp_dir} \
        preprocessed_filt_map_counts.interval_list
done

mv slurm-* slurm_outs/

# Filter the gnomad_genome_v3_0 vcf to get only high frequency snps

gnomad_genome_v3_0=/home/tgenref/homo_sapiens/grch38_hg38/public_databases/gnomad/r3.0/gnomad.genomes.r3.0.sites.pass.AnnotationReference.bcf

bcftools view \
    --exclude-types indels,mnps,ref,bnd,other \
    --genotype het \
    --targets ^chrX:10001-2781479,chrX:155701383-156030895,chrY:10001-2781479,chrY:56887903-57217415 \
    --exclude "AF<0.2 | AF>0.8" \
    --output-type v \
    --output-file ${temp_dir}/gnomad_genome_v3_0.preFilt.snps.vcf \
    ${gnomad_genome_v3_0}

bcftools view \
    --targets-file ^${ENCODE_DENY_LIST} \
    --output-type v \
    --output-file ${temp_dir}/filteredSnps.vcf \
    ${temp_dir}/gnomad_genome_v3_0.preFilt.snps.vcf

rm ${temp_dir}/gnomad_genome_v3_0.preFilt.snps.vcf

# create pon
mkdir -p "${temp_dir}/temp3"

gatk --java-options "-Xmx12g" CreateReadCountPanelOfNormals \
    --tmp-dir ${temp_dir}/temp3/ \
    # TODO dynamic number of inputs
    --input test_001/hdh5_counts/GM14632_CORIELL_p0_CL_Whole_C1_TWGSP.filter.counts.hdf5 \
    --annotated-intervals ${temp_dir}/anno_preprocessed_filt_map_counts.interval_list \
    --minimum-interval-median-percentile 10.0 \
    --output ${temp_dir}/cnv.pon.hdf5
