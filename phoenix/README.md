# Automated Scripts for Reference File Creation

These scripts are designed to create the full reference dataset needed for the Phoenix workflow. There is an expected order of script usage and each one checks that required steps are complete before they will start

### Script Usage Order
* create_genome_reference.sh (downloads fasta files for DNA and RNA alignments)
  * create_bwa_genome_index.sh
  * create_gene_model.sh
    * create_exome_capture_resources.sh
    * transcriptome_fasta_creation (automatic initiation by gene_model script)
      * create_salmon_index.sh
    * create_star_genome_index.sh
    * create_star-fusion_resource.sh
    * create_snpEff_db.sh
      * **NOTE**: make sure you update the snpEff.config in the utility_files directory BEFORE RUNNING
    * create_samtools_stats_non_N_region_file.sh
* create_vep_database.sh
* build_clinivar_20190715.sh
* build_cosmic_v90.sh
* build_dbSNP_b152.sh
  * create_genderCheck_SNP_list.sh
* build_delly_annotations_e97.sh (UPDATE, dependancies???)
* build_gnomeAD_r3.0.sh
* build_encode_blacklist.sh
* build_lymphocyte_count_windows.sh
* build_UCSC2ensembl_crossmapping.sh
* build_centromere_and_heterochromatin_bed_files.sh

### Required binaries
* wget
* md5sum
* bwa
* star
* salmon
* samtools
* bcftools

NOTE: Loaded using module load at TGen, but if available in your path things will work fine
