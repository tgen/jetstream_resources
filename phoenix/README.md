# Phoenix Analysis Pipeline (GRCh38)

These scripts are designed to create the full reference dataset needed for the Phoenix workflow. 
There is an expected order of script usage and each one checks that required steps are complete 
before they will start.

Because of difference between BWA for DNA alignment and STAR for RNA alignment this workflow leverages 
two different but related reference genomes for DNA, has ALT-contigs as BWA is ALT-aware and HLA contigs, 
while the reference used for RNA does not have ALT-contigs as STAR is not ALT-aware.

### Script Usage Order
* create_genome_reference.sh phoenix_resources.ini
  * create_bwa_genome_index.sh phoenix_resources.ini
  * create_gene_model.sh phoenix_resources.ini
    * create_exome_capture_resources.sh phoenix_resources.ini
    * create_salmon_index.sh
    * create_star_genome_index.sh
    * create_star-fusion_resource.sh
    * create_snpEff_db.sh
      * **NOTE**: make sure you update the snpEff.config in the utility_files directory BEFORE RUNNING
    * create_exome_capture_resources.sh phoenix_resources.ini
    * create_samtools_stats_non_N_region_file.sh phoenix_resources.ini
    * create_vep_database.sh phoenix_resources.ini
    * create_gatk_cnv_interval_list.sh phoenix_resources.ini
    * create_deepvariant_models.sh phoenix_resources.ini
* build_clinivar_20190715.sh
* build_cosmic_v90.sh
* build_dbSNP_b152.sh
  * create_genderCheck_SNP_list.sh
* build_delly_annotations_e97.sh
* build_gnomeAD_r3.0.sh
* build_encode_blacklist.sh
* build_lymphocyte_count_windows.sh
* build_UCSC2ensembl_crossmapping.sh
* build_centromere_and_heterochromatin_bed_files.sh
* build_Myeloma_FISH_Probe_Locations.sh

### Required binaries
* wget
* curl
* md5sum
* zcat
* python v3.7.2+
  * which modules are needed? We have lots in this module.
* bedtools v2.29.0
* GATK v4.1.4.0
* bwa v0.7.17+
* bowtie v2.3.5.1+
* star v2.7.3a+
* salmon v0.14.1+
* samtools v1.10.0+
* bcftools v1.10.1+
* snpEff v4_3t+
* R v3.6.1
  * library(tidyverse) # requires dplyr
* NCBI eUTILs
* JSON.awk
* gtfToGenePred ()
* cufflinks (cufflinks itself is not used but one of the provided tools, gffread, is used to create the transcriptome fasta)
  * gffread

NOTE: The "ENVIRONMENT" variable (TGen/LOCAL) in the phoenix_resourece.ini defines if tools are loaded into the path using a `module load` 
process, standard procedure at TGen, or if they are expected to be available in your $PATH at runtime. 

module load BEDTools/2.29.0-GCC-8.2.0-2.31.1
module load Python/3.7.2-foss-2019a