#!/bin/bash
#SBATCH --job-name="MSIsensorPro_scan"
#SBATCH --time=0-1:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 1
#SBATCH --partition defq,overflow

## Script to create list of homopolymers and miscrosatelites from reference genome

## Usage:
## sbatch --export ALL,FASTA="genome.fa",MICROSATELLITES_LIST="GRCh38tgen_decoy_alts_hla.microsatellites.list",MSISENSOR_PRO_VERSION="1.1.a" create_msisensor-pro_ref_microsatellites_list.sh 

module load msisensor/${MSISENSOR_PRO_VERSION}

msisensor-pro scan -d "$FASTA" -o "$MICROSATELLITES_LIST"

# Error Capture
if [ "$?" = "0" ]
then
    echo "PASSED_MSISENSOR_PRO_SCAN" >> README
else
    touch FAILED_MSISENSOR_PRO_SCAN
    echo "FAILED_MSISENSOR_PRO_SCAN" >> README
    exit 1
fi
