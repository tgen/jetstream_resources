# Suncity Analysis Pipeline (GRCh37/b37/hg19 - hs37d5)

A JetStream workflow to support, old human genome things...


### This directory contains Automated scripts for building files need for Coyote pipeline

* ../shared_resource_creation_scripts/create_genome_reference.sh suncity_resources.ini
  * ../shared_resource_creation_scripts/create_bwa_genome_index.sh suncity_resources.ini
  * ../shared_resource_creation_scripts/create_gene_model.sh suncity_resources.ini
    * ../shared_resource_creation_scripts/create_star_genome_index.sh suncity_resources.ini star_index_lengths.csv 
    * ../shared_resource_creation_scripts/create_salmon_index.sh suncity_resources.ini 
    * ../shared_resource_creation_scripts/create_snpEff_db.sh suncity_resources.ini 
    * ../shared_resource_creation_scripts/create_vep_database.sh suncity_resources.ini 
    * create_sexCheck_SNP_list.sh suncity_resources.ini
     
    * create_exome_capture_resources.sh
    * create_star-fusion_resource.sh
    * create_samtools_stats_non_N_region_file.sh
    *   