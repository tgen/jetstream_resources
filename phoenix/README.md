# Automated Scripts for Reference File Creation

These scripts are designed to create the full reference dataset needed for the Phoenix workflow. There is an expected order of script usage and each one checks that required steps are complete before they will start

### Script Usage Order
* create_genome_reference.sh
  * create_bwa_genome_index.sh
  * create_gene_model.sh
    * --> transcriptome_fasta_creation
      * create_salmon_index.sh
    * create_star_genome_index.sh