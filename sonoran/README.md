# Sonoran Analysis Pipeline (UU_Cfam_GSD_1.0)

### This directory contains Automated scripts for building files need for Sonoran pipeline

* sonoran/create_genome_reference.sh sonoran/sonoran_resources.ini
  * shared_resource_creation_scripts/create_bwa_genome_index.sh sonoran/sonoran_resources.ini
  * sonoran/create_gene_model.sh sonoran/sonoran_resources.ini
    * shared_resource_creation_scripts/create_star_genome_index.sh sonoran/sonoran_resources.ini sonoran/star_index_lengths.csv
    * shared_resource_creation_scripts/create_salmon_index.sh sonoran/sonoran_resources.ini
    * shared_resource_creation_scripts/create_snpEff_db.sh sonoran/sonoran_resources.ini
    * shared_resource_creation_scripts/create_vep_database.sh sonoran/sonoran_resources.ini
    * sonoran/create_exome_capture_resources.sh sonoran/sonoran_resources.ini sonoran/capture_kits.csv
    * sonoran/create_star-fusion_resource.sh sonoran/sonoran_resources.ini
    * sonoran/create_samtools_stats_non_N_region_file.sh sonoran/sonoran_resources.ini
