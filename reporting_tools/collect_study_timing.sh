#!/usr/bin/env bash

# usage: collect_study_timing.sh <List_of_Projects>

# Load needed modules
module load R/3.6.1-phoenix
module load Python/3.7.2-foss-2019a

# make folder for summary
mkdir -p timing_summary

# enter summary folder and create initial study template files
cd timing_summary
echo -e Project"\t"Group"\t"Total_CPU_Hours"\t"Total_Elapsed_Hours"\t"PCT_CPU_Hours"\t"PCT_Elapsed_Hours > study_task_summary.txt
echo -e Project"\t"Total_CPU_Hours"\t"Total_Elapsed_Hours > study_project_summary.txt

echo
echo

for project in `cat $1`
do

    echo "-----------------------------------------"
    echo "Processing ${project}"
    # Enter project folder
    cd ../${project}

    # Generate timing result
    python3 /home/tgenjetstream/git_repositories/jetstream_resources/reporting_tools/report_cpu_usage.py > ../timing_summary/${project}_timing.txt

    # return to project summary to run R scripts
    cd ../timing_summary
    Rscript --vanilla \
        /home/tgenjetstream/git_repositories/jetstream_resources/reporting_tools/summarize_project_runtime.R \
        --project ${project} \
        --time_summary ${project}_timing.txt \
        --task_summary study_task_summary.txt \
        --study_summary study_project_summary.txt

done

echo
echo

