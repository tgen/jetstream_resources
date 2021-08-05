# Phoenix Analysis Pipeline (GRCh38)

These scripts are designed to create the full reference dataset needed for the Phoenix workflow.
There is an expected order of script usage and each one checks that required steps are complete
before they will start.

Because of difference between BWA for DNA alignment and STAR for RNA alignment this workflow leverages
two different but related reference genomes for DNA, has ALT-contigs as BWA is ALT-aware and HLA contigs,
while the reference used for RNA does not have ALT-contigs as STAR is not ALT-aware.

### Script Usage Order
* create_genome_reference.sh phoenix_resources.ini  - RETEST-UPDATED/Pass
  * create_bwa_genome_index.sh phoenix_resources.ini  - Pass/Pass
  * ../shared_resource_creation_scripts/create_gatk_cnv_interval_list.sh phoenix_resources.ini  - Pass/Pass
  * create_samtools_stats_non_N_region_file.sh phoenix_resources.ini  - Pass/Pass
  * create_gene_model.sh phoenix_resources.ini phoenix_resources.ini - Pass/Pass
    * jetstream_resources/phoenix/create_disease_specific_resources.sh jetstream_resources/phoenix/phoenix_resources.ini - Testing/...
    * jetstream_resources/shared_resource_creation_scripts/create_salmon_index.sh jetstream_resources/phoenix/phoenix_resources.ini - Pass/Pass
    * jetstream_resources/shared_resource_creation_scripts/create_star_genome_index.sh jetstream_resources/phoenix/phoenix_resources.ini star_index_lengths.csv - Pass/Pass
    * jetstream_resources/phoenix/create_star-fusion_resource.sh jetstream_resources/phoenix/phoenix_resources.ini - Testing/... (hit walltime.. argh)
    * jetstream_resources/shared_resource_creation_scripts/create_snpEff_db.sh jetstream_resources/phoenix/phoenix_resources.ini - Pass/Pass
      * **NOTE**: make sure you update the snpEff.config in the utility_files directory BEFORE RUNNING
    * jetstream_resources/phoenix/create_vep_database.sh jetstream_resources/phoenix/phoenix_resources.ini  - Pass/Pass
    * jetstream_resources/shared_resource_creation_scripts/create_exome_capture_resources.sh jetstream_resources/phoenix/phoenix_resources.ini
      * NOTE: This is dependent on existing bed files, some of which are not freely available and must be downloaded in advance.
* ../shared_resource_creation_scripts/create_deepvariant_models.sh phoenix_resources.ini  - Pass/Pass
* ../shared_resource_creation_scripts/create_snpSniffer_references.sh phoenix_resources.ini
* build_annotation_files/build_delly_annotations_e98.sh   - Pass/Pass
* build_annotation_files/build_clinivar_20190715.sh  - Pass/Pass
* build_annotation_files/build_cosmic_v90.sh  - Pass/Pass
* build_annotation_files/build_dbSNP_b152.sh - Pass/Pass
  * create_genderCheck_SNP_list.sh - Pass/Pass
* build_annotation_files/build_gnomeAD_r2.1.1.sh  - Testing/...
* build_annotation_files/build_gnomeAD_r3.0.sh  - Testing/...
* build_annotation_files/build_encode_blacklist.sh  - Pass/Pass
  * build_centromere_and_heterochromatin_bed_files.sh  - Pass/Pass
* build_annotation_files/build_lymphocyte_count_windows.sh  - Pass/Pass
* build_annotation_files/build_UCSC2ensembl_crossmapping.sh  - Pass/Pass
* build_annotation_files/build_Myeloma_FISH_Probe_Locations.sh  - Pass/Pass


## Required Software
NOTE: The "ENVIRONMENT" variable (TGen/LOCAL) in the phoenix_resource.ini defines if tools are loaded into the path using a `module load`
process, standard procedure at TGen, or if they are expected to be available in your $PATH at runtime.

### Required Binaries that MUST be available in your $PATH
* module load SAMtools/1.10-GCC-8.2.0-2.31.1
* module load BWA/0.7.17-GCC-8.2.0-2.31.1
* module load cufflinks/2.2.1
  * NOTE: Provided "gffread" binary needed to build the transcriptome fasta
* module load Salmon/0.14.1-GCC-8.2.0-2.31.1
* module load STAR/2.7.3a
* module load snpEff/v4_3t
  * NOTE: requires JAVA v1.8
* module load STAR-Fusion/1.8.1-GCC-8.2.0-2.31.1-Perl-5.28.1-Python-3.7.2
  * NOTE: to build all starFusion resources the script downloads the required "dfamscan.pl" file
  * module load blast/2.7.1
  * module load hmmer/3.2.1
* module load BEDTools/2.29.0-GCC-8.2.0-2.31.1
* module load Python/3.7.2-foss-2019a
  * Needs: pandas

#### Standard UNIX Tools Used
* wget
* curl
* md5sum
* sum
* zcat
* python v3.7.2+
  * which modules are needed? We have lots in this module.
* bedtools v2.29.0
* GATK v4.1.4.0
* bwa v0.7.17+
* bowtie v2.3.5.1+
* star v2.7.3a+
* salmon v0.14.1+
* samtools v1.10.0+
* bcftools v1.10.1+
* snpEff v4_3t+
* R v3.6.1
  * library(tidyverse) # requires dplyr
* NCBI eUTILs
* JSON.awk
* gtfToGenePred ()
* cufflinks (cufflinks itself is not used but one of the provided tools, gffread, is used to create the transcriptome fasta)
  * gffread

### Required Binaries that MUST be defined as absolute paths in the phoenix.ini
NOTE: For local usage outside of TGen you must updates these variables before building resources
* snpEff.jar
    * https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip/download
* ncbi eUTILs
    * See (https://www.ncbi.nlm.nih.gov/books/NBK179288/) for details on installation
* JSON.awk
    * https://github.com/step-/JSON.awk/archive/1.3.tar.gz
* gtfToGenePred
    * UCSC (http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred)
* faToTwoBit
    * UCSC (http://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToTwoBit)
