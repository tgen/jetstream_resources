#!/bin/bash

module load HMMCopy_Utils/1.0.0-5911bf6

gcCounter \
 --window 1000000 \
 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,X \
 /home/tgenref/canis_familiaris/canfam3.1/canfam3.1_tgen/genome_reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa > gc_canfam3.1_1000kb.wig##