# Coyote Analysis Pipeline (canFAM3.1)

### This directory contains Automated scripts for building files need for Coyote pipeline

* ../shared_resource_creation_script/create_genome_reference.sh coyote_resources.ini
  * ../shared_resource_creation_scripts/create_bwa_genome_index.sh coyote_resources.ini
  * ../shared_resource_creation_scripts/create_gene_model.sh coyote_resources.ini
    * ../shared_resource_creation_scripts/create_star_genome_index.sh coyote_resources.ini star_index_lengths.csv
    * ../shared_resource_creation_scripts/create_salmon_index.sh coyote_resources.ini
    * ../shared_resource_creation_scripts/create_snpEff_db.sh coyote_resources.ini
        * THIS MIGHT HAVE FAILED, NO *.bin FILE CREATED
    * ../shared_resource_creation_scripts/create_vep_database.sh coyote_resources.ini
    * create_exome_capture_resources.sh
    * create_star-fusion_resource.sh
    * create_samtools_stats_non_N_region_file.sh
    *   
* next
