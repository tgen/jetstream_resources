# Mouse Pipeline supports mm10 mouse genome with ensembl 98 gene models


### This directory contains Automated scripts for building files need for cactusMouse pipeline

# ../shared_resource_creation_scripts/create_genome_reference.sh cactusMouse_resources.ini
# ../shared_resource_creation_scripts/create_bwa_genome_index_cactusMouse.sh cactusMouse_resources.ini 
# ../shared_resource_creation_scripts/create_gene_model.sh cactusMouse_resources.ini 
# ../shared_resource_creation_scripts/create_star_genome_index_catcusMouse.sh cactusMouse_resources.ini /home/kdrenner/jetstream_resources/cactusMouse/star_index_lengths.csv 
 ##### Had to update the create_star_genome_index.sh file. Received the error that not enough RAM was allocated.
 ####     Used --limitGenomeGenerateRAM 34000000000 and the code was able to run

# ../shared_resource_creation_scripts/create_salmon_index.sh cactusMouse_resources.ini
# ../shared_resource_creation_scripts/create_snpEff_db.sh cactusMouse_resources.ini

