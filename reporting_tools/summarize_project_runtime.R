# Plot out project run time for summary

# Load required modules
library(tidyverse)
library(lubridate)

# Import data table
# Generated within project folder using <python3 /home/tgenjetstream/git_repositories/jetstream_resources/reporting_tools/report_cpu_usage.py>
data <- read_delim("MMRF_1605_timing.txt", delim = "\t")
project_name <- "MMRF_1605"

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
  summarise(CPU_Hours = sum(Hours), 
            Total_Elapsed_Hours = sum(Elapsed)/3600,
            Total_Time = duration(hour = sum(hours), minute = sum(minutes), second = sum(seconds))
            )
  
# Plot Summary Data
ggplot(task_summary, aes(x=Group, y=Total_Elapsed_Hours)) +
  geom_bar(stat="identity") + 
  coord_flip()
ggsave(file=paste(project_name, "_ElapsedHours_by_TaskGroup.png", sep=""), dpi=150)

ggplot(task_summary, aes(x=Group, y=CPU_Hours)) +
  geom_bar(stat="identity") + 
  coord_flip()
ggsave(file=paste(project_name, "_CPUhours_by_TaskGroup.png", sep=""), dpi=150)

  
  
  
