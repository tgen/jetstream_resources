#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 16:49:47 2016

Seg Stitch

@author: dpenaherrera
"""
import sys
import os

# PRELIMINARIES
# Chromosomes
chromes = ["chr1", "chr2", "chr3",
           "chr4", "chr5", "chr6",
           "chr7", "chr8", "chr9",
           "chr10", "chr11", "chr12",
           "chr13", "chr14", "chr15",
           "chr16", "chr17", "chr18",
           "chr19", "chr20", "chr21",
           "chr22", "chrX", "chrY"]

# Creating Centromere Dictionary
centromere = {
        'chr1': [121605200, 143562200],
        'chr2': [91969000, 94569500],
        'chr3': [90455400, 93749200],
        'chr4': [49077200, 51816100],
        'chr5': [46433900, 50165300],
        'chr6': [58432200, 60242300],
        'chr7': [58031800, 60997400],
        'chr8': [43937900, 45969600],
        'chr9': [43236168, 61518900],
        'chr10': [39686683, 41712900],
        'chr11': [50813600, 54425074],
        'chr12': [34715400, 37269100],
        'chr13': [16000001, 18171400],
        'chr14': [0, 18670900],
        'chr15': [0, 19790100],
        'chr16': [36225200, 46423000],
        'chr17': [22745200, 26987200],
        'chr18': [15367200, 20940300],
        'chr19': [24330100, 27274500],
        'chr20': [26364200, 30038348],
        'chr21': [10650900, 12965800],
        'chr22': [11976900, 15154400],
        'chrX': [58467900, 62522800],
        'chrY': [10246200, 11041200]
                }


# ==============================================================================
# WHERE ALL THE MAGIC HAPPENS
# Grab all .seg files in path and convert them to .bed files


def seg_extend(path):
    os.chdir(path)
    for file in os.listdir(path):
        if file.endswith(".seg"):
            with open(file, 'r') as f:
                seglist = f.readlines()

            # Parse begins
            all_lines = []

            for line in seglist:
                all_lines.append(line.split('\t'))

            del seglist

            for i in range(len(all_lines)):
                line = all_lines[i]
                line[5] = line[5].strip()
                all_lines[i] = line

            for idx in range(1, len(all_lines)):
                all_lines[idx][2:4] = [int(x) for x in all_lines[idx][2:4]]

# -------FINDING SEGMENTS THAT CROSS THE CENTROMERE-----------------------------
            for chrom in chromes:
                i = 1
                while i in range(1, len(all_lines)):
                    line = all_lines[i]
                    if line[1] == chrom:
                        if (line[2] < centromere[chrom][0]) and (line[3] > centromere[chrom][1]):
                            newline = [line[0], line[1], centromere[chrom][1], line[3], line[4], line[5]]
                            all_lines.insert(i+1, newline)
                            all_lines[i][3] = centromere[chrom][0]
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
                        if (all_lines[k][3] != centromere[chrom][0]) and (all_lines[k+1][2] != centromere[chrom][1]):
                            all_lines[k][3] = (all_lines[k][3]+all_lines[k+1][2])/2
                            all_lines[k+1][2] = all_lines[k][3]
                        elif (all_lines[k][3] == centromere[chrom][0]) and (all_lines[k+1][2] == centromere[chrom][1]):
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
    seg_extend(sys.argv[1])
