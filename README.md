# Jetstream Resources
Collection of tools and scripts what can be used to summarize or build required annotation files for a JetStream workflow.

#### What is JetStream?
JetStream is the pipeline management system used at TGen. In a simple form it can be used to create analysis workflows
based on user defined tasks that are linked by a series of directives. In the advanced form, like the provided workflows,
project specific JSON files are integrated with template files to render a DAG that is then executed by a JetStream runner
on a local or distributed computing environment. For more details see the JetStream repository (https://github.com/tgen/jetstream) 

## JetStream Workflows
Multiple workflows exist to facilitate different analysis projects

### Phoenix (https://github.com/tgen/phoenix)
* Human analyis pipeline built around the GRCh38 reference genome and Ensembl version 98 gene models
  * Ensembl v98 is equivalent to Gencode v32
* Scripts used to build all required annotation files can be found in the "phoenix" directory

### Suncity (https://github.com/tgen/suncity)
* Developed for legacy project support
* Human analysis pipeline built around GRCh37, using the 1000 genomes hs37d5 reference, and ensembl GRCh37 archive gene models

### Coyote (https://github.com/tgen/coyote)
* Dog analysis pipeline built around CanFam3.1, using ensembl version 98 gene models

## Other Provided Tools and Resources
### Utility Files
* Files required for scripts in this repository

### Utility Scripts
* Fully parameterized scripts that can be used to create common outputs
  * annotation databases
  * alignment indexes
  * compile tools (used in testing only, production install by HPC)
  * transcriptome fasta