# Graph Study Summary

library(tidyverse)
library(optparse)

# Define Options
option_list = list(
  make_option(c("-s", "--study_name"), 
              type="character", 
              default=NULL, 
              help="Study Name", 
              metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Validate arguments are provided
if (is.null(opt$study_name)){
  print_help(opt_parser)
  stop("You must provide a study name string (ie. MMRF) to -s/--study_name", call.=FALSE)
}

# Read in files

tasks <- read_tsv("study_task_summary.txt")
projects <- read_tsv("study_project_summary.txt")

ggplot(tasks, aes(x=Group, y=Total_CPU_Hours)) + 
  geom_jitter() + 
  coord_flip()
ggsave(file=paste(opt$study_name, "_TotalCPUhours_by_Group_per_Project.png", sep=""), dpi=150)

ggplot(tasks, aes(x=Group, y=Total_Elapsed_Hours)) + 
  geom_jitter() + 
  coord_flip()
ggsave(file=paste(opt$study_name, "_TotalElapsedHours_by_Group_per_Project.png", sep=""), dpi=150)

ggplot(tasks, aes(x=Group, y=Max_Task_Elapsed_Hours)) + 
  geom_jitter() + 
  coord_flip()
ggsave(file=paste(opt$study_name, "_MaxTaskElapsedHours_by_Group_per_Project.png", sep=""), dpi=150)

ggplot(tasks, aes(x=Group, y=Max_Task_CPU_Hours)) + 
  geom_jitter() + 
  coord_flip()
ggsave(file=paste(opt$study_name, "_MaxTaskCPUhours_by_Group_per_Project.png", sep=""), dpi=150)

ggplot(projects, aes(x=Total_CPU_Hours, y=Total_Elapsed_Hours, color=Tasks)) + 
  geom_point() + 
  scale_colour_gradientn(colours = rainbow(10))
ggsave(file=paste(opt$study_name, "_ElapsedTime_by_Task_per_Group.png", sep=""), dpi=150)
