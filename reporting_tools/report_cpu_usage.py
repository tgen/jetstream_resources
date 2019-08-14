#!/usr/bin/env python3

###
# Copied from Ryan Richholt Toolkit
###

# Usage: report_cpu_usage.py <Path_to_Jetstream_Project>

# Generates a report of CPU Usage from a Jetstream project
# Future improvements : leverage tags for subsummary by subtype, RG_SM, Individual Tools

import argparse
import logging
import re
import sys

import jetstream
from datetime import datetime, timedelta

log = logging.getLogger(__name__)
ELAPSED_PAT = re.compile(r'((?P<days>\d*)-)?(?P<hours>\d*):(?P<minutes>\d*):(?P<seconds>\d*)$')


def arg_parser():
    parser = argparse.ArgumentParser(description='Report CPU hours for a project')
    parser.add_argument('project', nargs='*')
    return parser


def summarize_task_cpuhs(t):
    cpus = int(t.state.get('slurm_sacct', {}).get('AllocCPUS') or t.directives.get("cpus", 1))
    raw_elapsed = t.state.get('slurm_sacct', {}).get('Elapsed', '00:00:00')
    dt = {k: int(v or 0) for k, v in ELAPSED_PAT.match(raw_elapsed).groupdict().items()}
    elapsed = timedelta(**dt)
    seconds = elapsed.total_seconds()
    hours = seconds / 3600
    cpuh = cpus * hours
    return cpus, elapsed, cpuh


def report(project):
    total_cpu_hours = 0

    print('CPUs\tElapsed\tHours\tCumulative\tTask')
    workflow = project.load_workflow()
    log.critical(f'Reporting on: {workflow}')
    for name, t in project.load_workflow().tasks.items():
        cpus, task_elapsed, task_cpu_hours = summarize_task_cpuhs(t)
        total_cpu_hours += task_cpu_hours
        print(f'{cpus}\t{task_elapsed}\t{round(task_cpu_hours, 2)}\t{round(total_cpu_hours, 2)}\t{t}')


def main():
    parser = arg_parser()
    args = parser.parse_args()

    if args.project:
        for p in args.project:
            try:
                report(jetstream.Project(path=p))
            except Exception as e:
                log.exception(f'Failed to generate report for "{p}"')
    else:
        report(jetstream.Project())


if __name__ == '__main__':
    main()