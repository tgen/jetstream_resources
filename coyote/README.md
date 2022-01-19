# Coyote Analysis Pipeline (Dog10K Boxer Tasha)*
\* Previously designed for CanFam3.1

### This directory contains Automated scripts for building files need for Coyote pipeline

* shared_resource_creation_script/create_genome_reference.sh coyote/coyote_resources.ini
  * shared_resource_creation_scripts/create_bwa_genome_index.sh coyote/coyote_resources.ini
  * shared_resource_creation_scripts/create_gene_model.sh coyote/coyote_resources.ini
    * shared_resource_creation_scripts/create_star_genome_index.sh coyote/coyote_resources.ini coyote/star_index_lengths.csv
    * shared_resource_creation_scripts/create_salmon_index.sh coyote/coyote_resources.ini
    * shared_resource_creation_scripts/create_snpEff_db.sh coyote/coyote_resources.ini
    * shared_resource_creation_scripts/create_vep_database.sh coyote/coyote_resources.ini
    * coyote/create_exome_capture_resources.sh coyote/coyote_resources.ini coyote/capture_kits.csv
    * coyote/create_star-fusion_resource.sh coyote/coyote_resources.ini
    * coyote/create_samtools_stats_non_N_region_file.sh coyote/coyote_resources.ini
