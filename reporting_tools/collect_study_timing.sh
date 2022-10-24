#!/usr/bin/env bash

# usage: collect_study_timing.sh <List_of_Projects> <Study>

# Requires: jetstream installed to whichever version of python is used
# Potential enhancements:
#   - singularity would benefit from dynamic arguments, especially for binding

set -e

trap "echo jetstream not found" ERR

jetstreams_python=$(head -n1 $(which jetstream) | grep -o "/.*")

# Unsetting the trap since it is no longer valid
trap - ERR

# Load needed modules
module load singularity/3.7.1-phoenix

# make folder for summary
mkdir -p timing_summary

# enter summary folder and create initial study template files
cd timing_summary
echo -e Project"\t"Tags"\t"Tasks"\t"Total_CPU_Hours"\t"Max_Task_CPU_Hours"\t"Total_Elapsed_Hours"\t"Max_Task_Elapsed_Hours"\t"PCT_CPU_Hours"\t"PCT_Elapsed_Hours > study_task_summary.txt
echo -e Project"\t"Tasks"\t"Total_CPU_Hours"\t"Total_Elapsed_Hours > study_project_summary.txt

echo
echo

for project in `cat $1`
do

    echo "-----------------------------------------"
    echo "Processing ${project}"
    # Enter project folder
    cd ../${project}

    # Generate timing result
    $jetstreams_python /home/tgenjetstream/git_repositories/jetstream_resources/reporting_tools/report_cpu_usage.py > ../timing_summary/${project}_timing.txt

    # return to project summary to run R scripts
    cd ../timing_summary
    # Singularity command needs to be more dynamic
    # Example for gemini PMED - singularity exec -e -B /coh_labs/PMED docker://ghcr.io/tgen/jetstream_containers/r-with_modules:3.6.1 Rscript --vanilla
    # This does ping ghcr.io, so a local image might be preferred.
    singularity exec /home/tgenref/containers/r-with_modules_3.6.1.sif Rscript --vanilla \
        /home/tgenjetstream/git_repositories/jetstream_resources/reporting_tools/summarize_project_runtime.R \
        --project ${project} \
        --time_summary ${project}_timing.txt \
        --task_summary study_task_summary.txt \
        --study_summary study_project_summary.txt

done

echo
echo

# Summarize the overall study results
echo "##################################################"
echo "Summarizing overall study timing results"

# See comment above for singularity usage caveat
singularity exec /home/tgenref/containers/r-with_modules_3.6.1.sif Rscript --vanilla \
    /home/tgenjetstream/git_repositories/jetstream_resources/reporting_tools/study_summary_Graphs.R \
    --study_name $2

