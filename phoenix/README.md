# Automated Scripts for Reference File Creation

These scripts are designed to create the full reference dataset needed for the Phoenix workflow. There is an expected order of script usage and each one checks that required steps are complete before they will start

### Script Usage Order
* create_genome_reference.sh (downloads fasta files for DNA and RNA alignments)
  * create_bwa_genome_index.sh
  * create_gene_model.sh
    * transcriptome_fasta_creation (automatic initiation by gene_model script)
      * create_salmon_index.sh
    * create_star_genome_index.sh 
    * create_snpEff_db.sh
* build_clinivar_20190715.sh
* build_cosmic_v90.sh
* build_dbSNP_b152.sh
  * create_genderCheck_SNP_list.sh
* build_delly_annotations_e97.sh (UPDATE, dependancies???)
* build_gnomeAD_r3.0.sh
* build lymphocyte_count_windows.sh
