# GrandCanyon Analysis Pipeline (T2T Consortium - CHM13 v2.0 Human Reference)

A JetStream workflow to support the T2T human reference genome

## Output location
`/home/tgenref/homo_sapiens/t2t_chm13/chm13_v2.0`

### Script Usage Order
* create_genome_reference.sh grandcanyon_resources.ini 
  * create_bwa_mem2_genome_index.sh grandcanyon_resources.ini
  * create_gene_model.sh grandcanyon_resources.ini
    * create_star_genome_index.sh grandcanyon_resources.ini star_index_lengths.csv
      * Calls sbatch for individual index creation using singularity container
    * create_vep_database.sh

#### NOTE: 
We don't create minimap2 indexes as a number of different configurations are used and if a script pulls the wrong 
pre-made index then the setting used in the index supersedes the one provided in th command. This ensures execution is 
performed as outlined in the command.  It does was some CPU cycles recreating the index each time
