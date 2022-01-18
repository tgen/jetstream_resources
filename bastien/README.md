# Bastien Analysis Pipeline (ROS Cfam 1.0)

### This directory contains Automated scripts for building files need for Bastien pipeline

* shared_resource_creation_script/create_genome_reference.sh bastien/bastien_resources.ini
  * shared_resource_creation_scripts/create_bwa_genome_index.sh bastien/bastien_resources.ini
  * shared_resource_creation_scripts/create_gene_model.sh bastien/bastien_resources.ini
    * shared_resource_creation_scripts/create_star_genome_index.sh bastien/bastien_resources.ini bastien/star_index_lengths.csv
    * shared_resource_creation_scripts/create_salmon_index.sh bastien/bastien_resources.ini
    * shared_resource_creation_scripts/create_snpEff_db.sh bastien/bastien_resources.ini
    * shared_resource_creation_scripts/create_vep_database.sh bastien/bastien_resources.ini
    * coyote/create_exome_capture_resources.sh bastien/bastien_resources.ini coyote/capture_kits.csv
    * coyote/create_star-fusion_resource.sh bastien/bastien_resources.ini
    * coyote/create_samtools_stats_non_N_region_file.sh bastien/bastien_resources.ini
