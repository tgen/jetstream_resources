#!/usr/bin/env Rscript

# Load required modules
library(tidyverse)
library(optparse)

# Define Options
option_list = list(
  make_option(c("-p", "--project"), 
              type="character", 
              default=NULL, 
              help="Project Name", 
              metavar="character"),
  make_option(c("-t", "--time_summary"), 
              type="character", 
              default=NULL, 
              help="Project time summary file (produced with report_cpu_usage.py", 
              metavar="character"),
  make_option(c("-a", "--task_summary"), 
              type="character", 
              default=NULL, 
              help="File with project task summary variables", 
              metavar="character"),
  make_option(c("-s", "--study_summary"), 
              type="character", 
              default=NULL, 
              help="File with study summary variables", 
              metavar="character")
  ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Validate arguments are provided
if (is.null(opt$project)){
  print_help(opt_parser)
  stop("You must provide a project string (ie. MMRF_1234) to -p/--project", call.=FALSE)
}

if (is.null(opt$time_summary)){
  print_help(opt_parser)
  stop("You must provide a time summary file to -t/--time_summary", call.=FALSE)
}

if (is.null(opt$task_summary)){
  print_help(opt_parser)
  stop("You must provide a study task file to -a/--task_summary", call.=FALSE)
}

if (is.null(opt$study_summary)){
  print_help(opt_parser)
  stop("You must provide a study project summary file to -s/--study_summary", call.=FALSE)
}

##############################################
# Plot out project run time for summary
##############################################

# Import data table
# Generated within project folder using <python3 /home/tgenjetstream/git_repositories/jetstream_resources/reporting_tools/report_cpu_usage.py>
data <- read_delim(opt$time_summary, delim = "\t")
project_name <- opt$project

## Update the tibble table

# Remove "<Task(complete): " and ">" from the TASK column
data <- data %>% mutate_if(is.character, str_replace_all, pattern = '<', replacement = "")
data <- data %>% mutate_if(is.character, str_replace_all, pattern = "[(]", replacement = "")
data <- data %>% mutate_if(is.character, str_replace_all, pattern = '[)]', replacement = "")
data <- data %>% mutate_if(is.character, str_replace_all, pattern = 'Taskcomplete: ', replacement = "")
data <- data %>% mutate_if(is.character, str_replace_all, pattern = '[>]', replacement = "")

# Parse the Elapsed time into hours, minutes, seconds
data <- data %>% separate(Elapsed, into = c("hours", "minutes", "seconds"), sep = ":", convert = TRUE, remove = FALSE)

# Add Summary Columns
data <- data %>% mutate(Group = case_when(str_detect(Task, "^copy_fastqs") ~ "Copy_Fastq",
                                          str_detect(Task, "^split_fastq") ~ "Split_Fastq",
                                          str_detect(Task, "^chunked_bwa_mem_samtools_fixmate") ~ "BWA_Align",
                                          str_detect(Task, "^chunked_samtools_merge_rg_bams") ~ "Samtools_Merge",
                                          str_detect(Task, "^samtools_markdup") ~ "Samtools_MarkDup",
                                          str_detect(Task, "^bam_to_cram") ~ "Samtools_BamCram",
                                          str_detect(Task, "^gatk_collectwgsmetrics") ~ "Picard_Metric",
                                          str_detect(Task, "^gatk_collectwgsmetricswithnonzerocoverage") ~ "Picard_Metric",
                                          str_detect(Task, "^gatk_collectrawwgsmetrics") ~ "Picard_Metric",
                                          str_detect(Task, "^gatk_collectmultiplemetrics") ~ "Picard_Metric",
                                          str_detect(Task, "^gatk_convertsequencingarrtifacttooxog") ~ "Picard_Metric",
                                          str_detect(Task, "^gatk_collecthsmetrics") ~ "Picard_Metric",
                                          str_detect(Task, "^gatk_collectrnaseqmetrics") ~ "Picard_Metric",
                                          str_detect(Task, "^samtools_stats") ~ "Samtools_Metric",
                                          str_detect(Task, "^samtools_flagstat") ~ "Samtools_Metric",
                                          str_detect(Task, "^samtools_idxstats") ~ "Samtools_Metric",
                                          str_detect(Task, "^verifybamid2") ~ "Random_Stat",
                                          str_detect(Task, "^freebayes_sex_check") ~ "Random_Stat",
                                          str_detect(Task, "^snpsniffer_geno") ~ "Random_Stat",
                                          str_detect(Task, "^hmmcopy_make_wig_bwa") ~ "iChor_CNA",
                                          str_detect(Task, "^ichor_cna_bwa") ~ "iChor_CNA",
                                          str_detect(Task, "^haplotypecaller_gvcf") ~ "HaplotypeCaller",
                                          str_detect(Task, "^haplotypecaller_gvcf_merge") ~ "HaplotypeCaller",
                                          str_detect(Task, "^manta") ~ "Manta_Strelka",
                                          str_detect(Task, "^strelka2_filter_variants") ~ "Variant_Filter",
                                          str_detect(Task, "^strelka2") ~ "Manta_Strelka",
                                          str_detect(Task, "^deepvariant_make_examples") ~ "Deepvariant",
                                          str_detect(Task, "^deepvariant_call_variants") ~ "Deepvariant",
                                          str_detect(Task, "^deepvariant_postprocess_variants") ~ "Deepvariant",
                                          str_detect(Task, "^deepvariant_filter_variants") ~ "Deepvariant",
                                          str_detect(Task, "^lancet_merge_chunks") ~ "Variant_Merge",
                                          str_detect(Task, "^lancet_filter_variants") ~ "Variant_Filter",
                                          str_detect(Task, "^lancet") ~ "Lancet",
                                          str_detect(Task, "^octopus_merge_chunks") ~ "Variant_Merge",
                                          str_detect(Task, "^octopus_filter_variants") ~ "Variant_Filter",
                                          str_detect(Task, "^octopus") ~ "Octopus", 
                                          str_detect(Task, "^vardict_merge_chunks") ~ "Variant_Merge",
                                          str_detect(Task, "^vardict_filter_variants") ~ "Variant_Filter",
                                          str_detect(Task, "^vardict") ~ "VarDictJava",
                                          str_detect(Task, "^mutect2_merge_chunks") ~ "Variant_Merge",
                                          str_detect(Task, "^mutect2_filter_variants") ~ "Variant_Filter",
                                          str_detect(Task, "^mutect2_filter_calls") ~ "Mutect",
                                          str_detect(Task, "^mutect2_calculate_contamination") ~ "Mutect",
                                          str_detect(Task, "^mutect2_merge_pileup_summaries") ~ "Mutect",
                                          str_detect(Task, "^mutect2_learn_readorientationmodel") ~ "Mutect",
                                          str_detect(Task, "^mutect2_merge_stats") ~ "Mutect",
                                          str_detect(Task, "^mutect2_merge_pileup_summaries") ~ "Mutect",
                                          str_detect(Task, "^mutect2_GetPileupSummaries") ~ "Mutect",
                                          str_detect(Task, "^mutect2") ~ "Mutect",
                                          str_detect(Task, "^vcfmerger2") ~ "VCFmerger",
                                          str_detect(Task, "^bcftools_annotate") ~ "Annotation",
                                          str_detect(Task, "^snpeff") ~ "Annotation",
                                          str_detect(Task, "^vep") ~ "Annotation",
                                          str_detect(Task, "^bcftools_annotate") ~ "Annotation",
                                          str_detect(Task, "^bcftools_annotate") ~ "Annotation",
                                          str_detect(Task, "^delly") ~ "Delly",
                                          str_detect(Task, "^gatk_call_cnv") ~ "GATK_CNV",
                                          str_detect(Task, "^add_matched_rna") ~ "RNA_Steps",
                                          str_detect(Task, "^add_rna_header_to_vcf") ~ "RNA_Steps",
                                          str_detect(Task, "^salmon_quant_cdna") ~ "RNA_Steps",
                                          str_detect(Task, "^star_quant") ~ "RNA_Steps",
                                          str_detect(Task, "^star_fusion") ~ "RNA_Steps",
                                          str_detect(Task, "^fixmate_sort_star") ~ "RNA_Steps",
                                          str_detect(Task, "^markduplicates_star_gatk") ~ "RNA_Steps",
                                          str_detect(Task, "^rna_getBTcellLociCounts") ~ "RNA_Steps",
                                          TRUE ~ "Misc"
                                          )
                        )


# Plot
ggplot(data, aes(x=Group, y=Elapsed, color=as.factor(CPUs))) +
  geom_jitter() + 
  scale_color_discrete() +
  coord_flip()
ggsave(file=paste(project_name, "_ElapsedTime_by_Task_per_Group.png", sep=""), dpi=150)

ggplot(data, aes(x=Group, y=Hours)) +
  geom_jitter() + 
  coord_flip()
ggsave(file=paste(project_name, "_CPUhours_by_Task_per_Group.png", sep=""), dpi=150)

# Group and summarize to get realtime and CPU hours by task Group
task_summary <- data %>% 
  group_by(Group) %>% 
  summarise(Tasks = n(),
            Total_CPU_Hours = sum(Hours), 
            Total_Elapsed_Hours = as.double(sum(Elapsed)/3600),
            ) %>% 
  mutate(PCT_CPU_Hours = Total_CPU_Hours/sum(Total_CPU_Hours)) %>% 
  mutate(PCT_Elapsed_Hours = Total_Elapsed_Hours/sum(Total_Elapsed_Hours))

# Add column with project
task_summary <- task_summary %>% 
  add_column(Project = project_name, .before = "Group")
  
# Plot Summary Data
ggplot(task_summary, aes(x=Group, y=Total_Elapsed_Hours)) +
  geom_bar(stat="identity") + 
  coord_flip()
ggsave(file=paste(project_name, "_ElapsedHours_by_TaskGroup.png", sep=""), dpi=150)

ggplot(task_summary, aes(x=Group, y=Total_CPU_Hours)) +
  geom_bar(stat="identity") + 
  coord_flip()
ggsave(file=paste(project_name, "_CPUhours_by_TaskGroup.png", sep=""), dpi=150)

# Generate Project Summary
project_summary <- task_summary %>% 
  group_by(Project) %>% 
  summarise(Tasks = sum(Tasks),
            Total_CPU_Hours = sum(Total_CPU_Hours), 
            Total_Elapsed_Hours = (sum(Total_Elapsed_Hours)),
  )

##############################################
# Group with other projects for Study Summary
##############################################
  
# Read in project summary files
study_task_summary <- read_delim(opt$task_summary, delim = "\t", col_types = "ccdddd")
study_project_summary <- read_delim(opt$study_summary, delim = "\t", col_types = "cdd")

# Concatonate into Project summary files
study_task_summary <- bind_rows(study_task_summary, task_summary)
study_project_summary <- bind_rows(study_project_summary, project_summary)

# Write out Project summary file
write_tsv(study_task_summary, "study_task_summary.txt")
write_tsv(study_project_summary, "study_project_summary.txt")
