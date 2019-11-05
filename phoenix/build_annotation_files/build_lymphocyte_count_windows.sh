#!/usr/bin/env bash

# Usage: ./build_lymphocyte_count_windows.sh

# Write this full document as a README
cat $0 > README

# Added user information and timestamp to README
USER=`whoami`
DATE=`date`
echo "Downloaded and Processed by:  ${USER}" >> README
echo ${DATE} >> README

## This will write out a bed file with the human b38 locations for the B-cell and T-cell receptor loci
## These are used to get a crude count to determine if a sample is B-cell or T-cell enriched

## Loci extracted from the ensembl v97 GTF by visual identification by Jonathan Keats

## IGK (IGKC - IGKV2-40
# chr2:88857161-89330429

## IGH (IGHA2 - IGHVIII-82)
# chr14:105583731-106879812

## IGL (IGLV4-69 - IGLC7)
# chr22:22030934-22923034

## TCR-alpha (TRAV1-1 - TRAC)
# chr14:21621838-22552156

## TCR-beta (TRBV1 - TRBC2) there is one V past the C2 but didn't include it
# chr7:142299177-142802748

### Print a bed file with the regions of interest (0-base corrected positions)
echo -e chr2"\t"88857160"\t"89330429"\t"IGK > lymphocyteReceptor_loci.bed
echo -e chr14"\t"105583730"\t"106879812"\t"IGH >> lymphocyteReceptor_loci.bed
echo -e chr22"\t"22030933"\t"22923034"\t"IGL >> lymphocyteReceptor_loci.bed
echo -e chr14"\t"21621837"\t"22552156"\t"TRA >> lymphocyteReceptor_loci.bed
echo -e chr7"\t"142299176"\t"142802748"\t"TRB >> lymphocyteReceptor_loci.bed
