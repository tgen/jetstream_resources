#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 16:49:47 2016

Seg Stitch

@author: dpenaherrera
"""
import sys
import os

# ==============================================================================
# WHERE ALL THE MAGIC HAPPENS
# Grab all .seg files in path and convert them to .bed files


def seg_extend(file):
    with open(file, 'r') as f:
        seglist = f.readlines()

    # Parse begins
    all_lines = []

    for line in seglist:
        all_lines.append(line.rstrip().split('\t'))

    del seglist

    for idx in range(1, len(all_lines)):
        all_lines[idx][2:4] = [int(x) for x in all_lines[idx][2:4]]

# -------FINDING SEGMENTS THAT CROSS THE CENTROMERE-----------------------------
    for chrom in chromes:
        i = 1
        while i in range(1, len(all_lines)):
            line = all_lines[i]
            if line[1] == chrom:
                if (line[2] < centromeres[chrom][0]) and (line[3] > centromeres[chrom][1]):
                    newline = [line[0], line[1], centromeres[chrom][1], line[3], line[4], line[5]]
                    all_lines.insert(i+1, newline)
                    all_lines[i][3] = centromeres[chrom][0]
                    del newline
                    i += 1
                else:
                    i += 1
            else:
                i += 1

# ---STITCHING SEMGENTS---------------------------------------------------------
    for k in range(1, len(all_lines)-1):
        for chrom in chromes:
            if (all_lines[k][1] == chrom) and (all_lines[k+1][1] == chrom):
                if (all_lines[k][3] != centromeres[chrom][0]) and (all_lines[k+1][2] != centromeres[chrom][1]):
                    all_lines[k][3] = (all_lines[k][3]+all_lines[k+1][2])/2
                    all_lines[k+1][2] = all_lines[k][3]
                elif (all_lines[k][3] == centromeres[chrom][0]) and (all_lines[k+1][2] == centromeres[chrom][1]):
                    continue

# ---WRITE THE OUTPUTS TO FILE--------------------------------------------------
    for idx in range(1, len(all_lines)):
        all_lines[idx][2:4] = [str(x) for x in all_lines[idx][2:4]]

    for n in range(1, len(all_lines)):
        all_lines[n] = [all_lines[n][1], all_lines[n][2], all_lines[n][3], all_lines[n][4], all_lines[n][5]]

    del(all_lines[0])

    bedfile = os.path.splitext(file)[0]+".bed"

    with open(bedfile, 'w') as f:
        for line in all_lines:
            f.write("\t".join(elem for elem in line)+"\n")


# ==============================================================================
if __name__ == '__main__':
    centromeres = {}
    with open(sys.argv[1]) as centro:
        for interval in centro:
            (key, start, stop) = interval.rstrip().split('\t')
            centromeres[key] = [int(start), int(stop)]

    chromes = list(centromeres.keys())

    seg_extend(sys.argv[2])
