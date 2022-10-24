#!/usr/bin/env python3

# Usage: report_cpu_usage.py --project <Path_to_Jetstream_Project> --ignore_tags <tag_to_ignore> <another_tag_to_ignore>
# Arguments are optional 

# Generates a report of CPU Usage from a Jetstream project

# Author: Ryan Richolt
# Updated by: Bryce Turner - bturner@tgen.org

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
    parser.add_argument('--project', nargs='*')
    parser.add_argument('--ignore_tags', nargs='+')
    return parser


def summarize_task_cpuhs(t):
    cpus = int(t.state.get('slurm_sacct', {}).get('AllocCPUS') or t.directives.get("cpus", 1))
    raw_elapsed = t.state.get('slurm_sacct', {}).get('Elapsed', '00:00:00')
    dt = {k: int(v or 0) for k, v in ELAPSED_PAT.match(raw_elapsed).groupdict().items()}
    elapsed = timedelta(**dt)
    total_seconds = int(elapsed.total_seconds())
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    fElapsed = f'{hours:02d}:{minutes:02d}:{seconds:02d}'
    cpuh = cpus * elapsed.total_seconds() / 3600
    return cpus, fElapsed, cpuh


def report(project, ignore_tags):
    total_cpu_hours = 0

    print('CPUs\tElapsed\tHours\tCumulative\tTask\tTags')
    workflow = project.load_workflow()
    log.critical(f'Reporting on: {workflow}')
    for name, t in workflow.tasks.items():
        cpus, task_elapsed, task_cpu_hours = summarize_task_cpuhs(t)
        total_cpu_hours += task_cpu_hours
        tags = t.directives.get('tags', ['Untagged'])
        ignore_tags_re_string = '|'.join(f'({tag})' for tag in ignore_tags).replace("-", "_")
        cleanre = re.compile(ignore_tags_re_string)
        # We perform a regex search because we expect the ignore tags to be partial
        # for example may want to ignore 'stats' tags, 'stats' may not be an actual tag
        # but it could be part of a tag, e.g. stats2json or stats2lims
        clean_tags = [tag for tag in tags if not re.search(cleanre, tag.replace("-", "_")) ]
        tags_string = '_'.join(clean_tags).replace(" ", "_")
        print(f'{cpus}\t{task_elapsed}\t{round(task_cpu_hours, 2)}\t{round(total_cpu_hours, 2)}\t{name}\t{tags_string}')


def main():
    parser = arg_parser()
    args = parser.parse_args()

    if args.ignore_tags:
        ignore_tags = args.ignore_tags
    else:
        ignore_tags = []

    if args.project:
        for p in args.project:
            try:
                report(jetstream.Project(path=p), ignore_tags)
            except Exception as e:
                log.exception(f'Failed to generate report for "{p}"')
    else:
        report(jetstream.Project(), ignore_tags)


if __name__ == '__main__':
    main()
