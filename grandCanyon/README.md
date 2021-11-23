# GrandCanyon Analysis Pipeline (T2T Consortium - CHM13 v1.1 Human Reference)

A JetStream workflow to support the T2T human reference genome (WARNING: No chrY today)

### Script Usage Order
* shared_resource_creation_scripts/create_genome_reference.sh grandCanyon/grandcanyon_resources.ini 
  * shared_resource_creation_scripts/create_bwa_genome_index.sh grandCanyon/grandcanyon_resources.ini
  * shared_resource_creation_scripts/create_bwa_mem2_genome_index.sh grandCanyon/grandcanyon_resources.ini
  * shared_resource_creation_scripts/create_bowtie2_genome_index.sh grandCanyon/grandcanyon_resources.ini
  * shared_resource_creation_scripts/create_snpSniffer_references.sh suncity/grandcanyon_resources.ini
  * grandCanyon/create_gene_model.sh grandCanyon/grandcanyon_resources.ini
  

    * shared_resource_creation_scripts/create_star_genome_index.sh suncity/suncity_resources.ini suncity/star_index_lengths.csv
    * shared_resource_creation_scripts/create_salmon_index.sh suncity/suncity_resources.ini
    * shared_resource_creation_scripts/create_snpEff_db.sh suncity/suncity_resources.ini
    * shared_resource_creation_scripts/create_vep_database.sh suncity/suncity_resources.ini
    * shared_resource_creation_scripts/create_sexCheck_SNP_list.sh suncity/suncity_resources.ini
   
    * shared_resource_creation_scripts/create_exome_capture_resources.sh suncity/suncity_resources.ini suncity/capture_kits.csv
    * suncity/create_star-fusion_resource.sh suncity/suncity_resources.ini
    * suncity/create_samtools_stats_non_N_region_file.sh suncity/suncity_resources.ini
    * shared_resource_creation_scripts/create_gatk_cnv_interval_list.sh suncity/suncity_resources.ini
    * shared_resource_creation_scripts/create_deepvariant_models.sh suncity/suncity_resources.ini

    * suncity/build_annotation_files/build_dbSNP_b138_broadBundle.sh suncity/suncity_resources.ini
    * suncity/build_annotation_files/build_delly_annotations_e87.sh suncity/suncity_resources.ini
    * suncity/build_annotation_files/build_gnomAD_r2.1.1.sh suncity/suncity_resources.ini