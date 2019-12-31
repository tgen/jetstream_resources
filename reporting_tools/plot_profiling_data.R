#!/usr/bin/env Rscript --vanilla

# Load required modules
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(optparse))
#suppressPackageStartupMessages(require(cowplot))

# Define Options
option_list = list(
  make_option(c("-t", "--task_title"), 
              type="character", 
              default=NULL, 
              help="Summary task name, groups multiple inputs and titles graphs", 
              metavar="taskname"), 
  make_option(c("-i", "--monitor_interval"), 
              type="integer", 
              default=2, 
              help="Monitoring interval in seconds [2]"), 
  make_option(c("-p", "--max_cpu_pct"), 
              type="integer", 
              default=100, 
              help="Max CPU level for task, 1 core = 100 [100]"), 
  make_option(c("-e", "--max_mem_gb"), 
              type="integer", 
              default=2000, 
              help="Max memory in Gb expected, [2000]"),
  make_option(c("-n", "--number_of_inputs"), 
              type="integer", 
              default=1, 
              help="Number of input timing files to summarize for task [1]", 
              metavar="[1,2,3]"),
  make_option(c("-a", "--input1"), 
              type="character", 
              default=NULL, 
              help="Name of first input file", 
              metavar="filename"), 
  make_option(c("-x", "--input1_process"), 
              type="character", 
              default=NULL, 
              help="Name of first input process", 
              metavar="process_name (process1)"),
  make_option(c("-b", "--input2"), 
              type="character", 
              default=NULL, 
              help="Name of second input file", 
              metavar="filename"), 
  make_option(c("-y", "--input2_process"), 
              type="character", 
              default=NULL, 
              help="Name of second input process", 
              metavar="process_name (process2)"),
  make_option(c("-c", "--input3"), 
              type="character", 
              default=NULL, 
              help="Name of third input file", 
              metavar="filename"),
  make_option(c("-z", "--input3_process"), 
              type="character", 
              default=NULL, 
              help="Name of third input process", 
              metavar="process_name (process3)")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Validate arguments are provided
if (is.null(opt$task_title)){
  print_help(opt_parser)
  stop("You must provide a task title for summary (ie. bwa_fixmate_sort) to -t/--task_title", call.=FALSE)
}

if (is.null(opt$number_of_inputs)){
  print_help(opt_parser)
  stop("You must provide an integer indicating number of input files -n/--number_of_inputs", call.=FALSE)
}

if (is.null(opt$input1)){
  print_help(opt_parser)
  stop("You must provide a filename to -a/--input1", call.=FALSE)
}

if (is.null(opt$input1_process)){
  print_help(opt_parser)
  stop("You must provide a filename to -a/--input1", call.=FALSE)
}


##############################################
# Load data files add columns
##############################################

data1 <- read_table2(opt$input1, 
                    comment = "#", 
                    col_names = c("Time", "CPU", "Real_MEM", "Virtual_MEM"), 
                    col_types = "nnnn")
data1 <- data1 %>% 
  add_column(Process = opt$input1_process) %>% 
  add_column(Task = opt$task_title) %>% 
  rowid_to_column("ID") %>% 
  mutate(Timepoint = (ID-1)*opt$monitor_interval)

write_tsv(data1, "data1.tsv")

if(opt$number_of_inputs > 1){
  print("Two or three inputs provided")
  data2 <- read_delim(opt$input2, 
                      delim = " ", 
                      comment = "#", 
                      col_names = c("Time", "CPU", "Real_MEM", "Virtual_MEM"), 
                      col_types = "nnnn")
  data2 <- data2 %>% 
    add_column(Process = opt$input2_process) %>% 
    add_column(Task = opt$task_title) %>% 
    rowid_to_column("ID") %>% 
    mutate(Timepoint = (ID-1)*opt$monitor_interval)
  write_tsv(data2, "data2.tsv")
}

if(opt$number_of_inputs == 3){
  print("Three inputs provided")
  data3 <- read_delim(opt$input3, 
                      delim = " ", 
                      comment = "#", 
                      col_names = c("Time", "CPU", "Real_MEM", "Virtual_MEM"), 
                      col_types = "nnnn")
  data3 <- data3 %>% 
    add_column(Process = opt$input3_process) %>% 
    add_column(Task = opt$task_title) %>% 
    rowid_to_column("ID") %>% 
    mutate(Timepoint = (ID-1)*opt$monitor_interval)
  write_tsv(data3, "data3.tsv")
}

# Join tables
if(opt$number_of_inputs == 1){
  plot_data <- data1 
} else if(opt$number_of_inputs == 2){
  plot_data <- bind_rows(data1, data2)
  # calculate summary
  summary <- plot_data %>% 
    group_by(Timepoint) %>% 
    summarise(Time = mean(Time), 
              CPU = sum(CPU), 
              Real_MEM = sum(Real_MEM), 
              Virtual_MEM = sum(Virtual_MEM)) %>% 
    rowid_to_column("ID") %>% 
    add_column(Process = "Total") %>% 
    add_column(Task = opt$task_title) %>% 
    select(ID, Time, CPU, Real_MEM, Virtual_MEM, Process, Task, Timepoint)
  # combine with individual data
  plot_data <- bind_rows(plot_data, summary)
} else if(opt$number_of_inputs == 3){
  plot_data <- bind_rows(data1, data2, data3)
  # calculate summary
  summary <- plot_data %>% 
    group_by(Timepoint) %>% 
    summarise(Time = mean(Time), 
              CPU = sum(CPU), 
              Real_MEM = sum(Real_MEM), 
              Virtual_MEM = sum(Virtual_MEM)) %>% 
    rowid_to_column("ID") %>% 
    add_column(Process = "Total") %>% 
    add_column(Task = opt$task_title) %>% 
    select(ID, Time, CPU, Real_MEM, Virtual_MEM, Process, Task, Timepoint)
  # combine with individual data
  plot_data <- bind_rows(plot_data, summary)
}

##############################################
# Plot results
##############################################

filename_cpu <- paste(opt$input1_process, "cpu.png", sep = "_")

p_cpu <- ggplot(plot_data, aes(Time, CPU, color = Process)) + 
  geom_line() + 
  geom_hline(yintercept = opt$max_cpu_pct, color = "red", linetype = "dashed") + 
  xlim(c(0,250)) + 
  ggtitle(opt$task_title) + 
  theme(plot.title = element_text(hjust=0.5), legend.position="top")
p_cpu
ggsave(filename_cpu, width = 20)

filename_mem <- paste(opt$input1_process, "mem.png", sep = "_")

p_mem <- ggplot(plot_data, aes(Time, Real_MEM, color = Process)) + 
  geom_line() + 
  geom_hline(yintercept = opt$max_mem_gb, color = "red", linetype = "dashed") + 
  ggtitle(opt$task_title) + 
  theme(plot.title = element_text(hjust=0.5), legend.position="top")
p_mem
ggsave(filename_mem)

# plot as a single image to fit on our google slide template (8.5w, 4h)
#p_grid <- plot_grid(p_cpu, p_mem, ncol = 2, nrow = 1)
#save_plot(filename, p_grid, ncol = 2, base_asp = 1.2, base_width = 8.5)