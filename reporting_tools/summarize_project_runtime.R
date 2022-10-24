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

# Parse the Elapsed time into hours, minutes, seconds
data <- data %>% separate(Elapsed, into = c("hours", "minutes", "seconds"), sep = ":", convert = TRUE, remove = FALSE)

data <- data %>% select("Tags", everything())

# Plot
ggplot(data, aes(x=Tags, y=Elapsed, color=as.factor(CPUs))) +
  geom_jitter(size=0.5) + 
  scale_color_discrete() +
  coord_flip()
ggsave(file=paste(project_name, "_ElapsedTime_by_Task_per_Group.png", sep=""), dpi=150)

ggplot(data, aes(x=Tags, y=Hours)) +
  geom_jitter(size=0.5) + 
  coord_flip()
ggsave(file=paste(project_name, "_CPUhours_by_Task_per_Group.png", sep=""), dpi=150)

# Group and summarize to get realtime and CPU hours by task Group
task_summary <- data %>% 
  group_by(Tags) %>% 
  summarise(Tasks = n(),
            Total_CPU_Hours = sum(Hours), 
            Max_Task_CPU_Hours = max(Hours),
            Total_Elapsed_Hours = as.double(sum(Elapsed)/3600),
            Max_Task_Elapsed_Hours = max(Elapsed)
            ) %>% 
  mutate(PCT_CPU_Hours = Total_CPU_Hours/sum(Total_CPU_Hours)) %>% 
  mutate(PCT_Elapsed_Hours = Total_Elapsed_Hours/sum(Total_Elapsed_Hours))

# Add column with project
task_summary <- task_summary %>% 
  add_column(Project = project_name, .before = "Tags")
  
# Plot Summary Data
ggplot(task_summary, aes(x=Tags, y=Total_Elapsed_Hours)) +
  geom_bar(stat="identity") + 
  coord_flip()
ggsave(file=paste(project_name, "_ElapsedHours_by_TaskGroup.png", sep=""), dpi=150)

ggplot(task_summary, aes(x=Tags, y=Total_CPU_Hours)) +
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
study_task_summary <- read_delim(opt$task_summary, delim = "\t", col_types = "ccidddtdd")
study_project_summary <- read_delim(opt$study_summary, delim = "\t", col_types = "cidd")

# Concatonate into Project summary files
study_task_summary <- bind_rows(study_task_summary, task_summary)
study_project_summary <- bind_rows(study_project_summary, project_summary)

# Write out Project summary file
write_tsv(study_task_summary, "study_task_summary.txt")
write_tsv(study_project_summary, "study_project_summary.txt")
