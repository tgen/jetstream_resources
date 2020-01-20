# Coyote Analysis Pipeline (canFAM3.1)

### This directory contains Automated scripts for building files need for Coyote pipeline

* ../shared_resource_creation_script/create_genome_reference.sh coyote_resources.ini
  * ../shared_resource_creation_scripts/create_bwa_genome_index.sh coyote_resources.ini
  * ../shared_resource_creation_scripts/create_gene_model.sh coyote_resources.ini
    * ../shared_resource_creation_scripts/create_star_genome_index.sh coyote_resources.ini star_index_lengths.csv
    * ../shared_resource_creation_scripts/create_salmon_index.sh coyote_resources.ini  
* next
