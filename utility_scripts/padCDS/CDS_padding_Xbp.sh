#!/usr/bin/env bash


## NOTE: This file MUST be started from a DIR_WORK directory where temp and out files will be written;
## NOTE: minimum # of CPUs recommended 10
#set -eu -o pipefail


DIR_SCRIPTS=$(dirname $0)/scripts
## 'export' the script directory that comes with the tool installation
## expected scripts must be in the same directory as the current script.
export PATH=${DIR_SCRIPTS}:${PATH}

## Check requirements ; All MUST be in your PATH
type bedtools >/dev/null 2>&1 || { echo >&2 "require << bedtools >>  but it's not in PATH.  Aborting."; exit 1; }
type bgzip >/dev/null 2>&1 || { echo >&2 "require << bgzip >>  but it's not in PATH.  Aborting."; exit 1; }
type python3 >/dev/null 2>&1 || { echo >&2 "require << python3 >>  but it's not in PATH.  Aborting."; exit 1; }
type parse_gtf_to_bed_for_all_lines_and_features.py >/dev/null 2>&1 || { echo >&2 "require << parse_gtf_to_bed_for_all_lines_and_features.py >>  but it's not in PATH.  Aborting."; exit 1; }

source ${DIR_SCRIPTS}/functions.sh

init_some_vars
getOptions "$@"

##@@@@@@@@@@@@@@@@@@@@@@
## CHECKING INPUTS
##@@@@@@@@@@@@@@@@@@@@@@
if [[ ! -e ${GTF} || ! -f ${GTF} ]] ;
then
	echo -e "ERROR: FILE NOT FOUND: ${GTF}"
	usage ; requirements ;
	exit 1
fi

if [[ ${GTF##*.} != "gtf" ]] ;
then
	echo -e "ERROR: File extension expected is '.gtf'; Please check your input; Aborting" ;
	usage ; requirements ;
	exit 1
fi

if [[ ${CPUS} == "" || ! ${CPUS} =~ ^[0-9]*$ || ${CPUS} -lt 1 ]] ;
then
	echo -e "ERROR: CPUS must be an INTEGER and must be greater than 0; Aborting"
	usage ; requirements ;
	exit 1
fi

if [[ ( ${NBP} != "" && ! ${NBP} =~ ^[0-9]*$ ) || ${NBP} -lt 1 ]] ;
then
	echo -e "ERROR: NBP or Number of Padding must be an INTEGER and must be greater than 0; Aborting."
	usage ; requirements ;
	exit 1
fi

if [[ ${DIR_WORK} == "" ]]
then
	echo -e "ERROR: DIR WORK MUST be provided; Aborting" ;
	usage ; requirements ;
	exit 1
elif [[ ! -e ${DIR_WORK} || ! -d ${DIR_WORK} ]] ;
then
	echo -e "ERROR: DIR NOT FOUND: ${DIR_WORK} must be an EXISTING directory; Aborting."
	echo -e "you may run: mkdir -p ${DIR_WORK}"
	usage ; requirements ;
	exit 1
fi

## RECAP INPUTS after Checking
echo -e "recap input:"
echo "
DIR_WORK == ${DIR_WORK}
GTF == '${GTF}'
OUTFILE = '${OUTFILE}'
CLEAN_TEMP_FILES == ${CLEAN_TEMP_FILES}
CPUS == ${CPUS}
NBP == ${NBP}
"

##@@@@@@@@@@@@@@@@
## CONSTANTs
##@@@@@@@@@@@@@@@@
#(these can be transformed into dynamic variable or user given variable if requested)
GTF_BN=$(basename ${GTF})
BED=${GTF_BN/gtf/bed}
BED_BN=$(basename ${BED})
PADDED_BED=${BED_BN/.bed/.ss_padded.bed}
CDS_FILE=${BED_BN/.bed/.cds_only.7cols.bed}
MBED=isec__all_ensts_cds_start_stop.bed
MBEDCDS=${MBED}.padded_CDS.bed
MBEDSS=${MBED}.start_stop.bed
DEBUGGING="no"  ## HARDCODED ; developpers turn it on


echo -e "Entering the WORKING DIRECTORY ..."
cd ${DIR_WORK}
echo -e "curdir (or/and dir_work) is: ${PWD}"
echo -e "$(date)"


## split GTF to parallelize the GTF2BED conversion
echo -e "Removing Previous Split files ..."
rm sGTF__*  &> /dev/null
echo -e "split GTF ..."
split -n l/$( echo "scale=0 ; ($(cat ${GTF} | wc -l )/100000)+1" | bc -l ) ${GTF} sGTF__
check_ev $? "parallel GTF2BED"

echo -e "converting the input GTF file to BED file ..."
parallel --link -k --halt soon,fail=1 -j ${CPUS} parse_gtf_to_bed_for_all_lines_and_features.py --gtf {1} --out {2} ::: $(find -maxdepth 1 -type f -name "sGTF__*" | sort | tr '\n' ' ' ) ::: $(find -maxdepth 1 -type f -name "sGTF__*" | sort | sed 's/$/.bed/' | tr '\n' ' ' )
check_ev $? "parallel GTF2BED"

cat $(find -maxdepth 1 -type f -name "sGTF__*" | sort | tr '\n' ' ' ) | grep -vE "##seqname" | sort --parallel=${CPUS} > ${BED}
check_ev $? "cat sGTF___ files"


if [[ ! -e ${BED} || ! -f ${BED} ]] ;
then
	echo -e "ERROR: FILE NOT FOUND: ${BED} ; The file has not been created. Check if GTF was correct; Aborting"
	exit 1
fi


echo -e "get BED header line ..."
echo -e "$(head -n1 ${BED})\t$(head -n1 ${BED})\t$(head -n 1 ${BED} | cut -f1-6,21)" > header_line_50cols.txt

echo -e "extracting lines with CDS and making cds_only bed file ..."
sort --parallel=${CPUS} -k1,1V -k2,2n -k3,3n <( cat ${BED} | awk '$13=="CDS"' | sed 's/ \+/_/g' | cut -f1-6,11,13 | awk '{FS=OFS="\t" ; $5=$7 ; print}' | cut -f1-6,8) > ${CDS_FILE}
check_ev $? "subset CDS"

echo -e "extracting lines with start_codon and making 'start_codon.bed' file ..."
(head -n1 ${BED} ; sort --parallel=${CPUS} -k1,1V -k2,2n -k3,3n  ${BED}  | awk '$13=="start_codon"' )  | sed 's/ \+/_/g' | cut -f1-6,11,13 | awk '{FS=OFS="\t" ; $5=$7 ; print}' | cut -f1-6,8 > start_codon.bed

echo -e "extracting lines with stop_codon and making 'stop_codon.bed' file ..."
(head -n1 ${BED} ; sort --parallel=${CPUS} -k1,1V -k2,2n -k3,3n  ${BED}  | awk '$13=="stop_codon"' )  | sed 's/ \+/_/g' | cut -f1-6,11,13 | awk '{FS=OFS="\t" ; $5=$7 ; print}' | cut -f1-6,8 > stop_codon.bed

## need to remove or add three bp to STOP CODON, because stop codon does NOT overlap CDS since it is not a part of the CDS.
## removal or addition depends on the STRAND.  we do that on both positions to get a 100% 3bp overlap of the stop_codon to a CDS as the start_codon does since the start codon is the first codon in the first CDS.
echo -e "modifying stop codon coordinates to make them overlapping CDS region ... "
awk '{FS=OFS="\t" ; if($6=="+") {$2=$2-3 ; $3=$3-3} else {$2=$2+3 ; $3=$3+3} ; print }' stop_codon.bed > stop_codon.padded_3bp.bed
check_ev $? "add 3bp to stop codon"

## there is no overlapping of start and stop codon that belongs to the same transcript ;
## ALL lines are INDEPENDENT ;
echo -e "concatenating start and stop codons files ... -->  start_and_stop_codons.sorted.bed "
sort --parallel=${CPUS} -k1,1V -k2,2n -k3,3n <( cat start_codon.bed stop_codon.padded_3bp.bed ) | grep -v seqname > start_and_stop_codons.sorted.bed
check_ev $? "concat start and stop codons"


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ DEBUG @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Some transcripts do have the same start_codon or the same stop_codon
## So we should be able to group these BUT that depends on the CDS;
## Do these transcripts with the same start (or stop) use the same cds?
## Let's count how many start and stop codon per Transcript
## This debugging part is mostly for developper(s)
if [[ ${DEBUGGING} == "yes" ]]
	then
	echo -e "##@@@@@   entering DEBUG part   @@@@##"
	echo -e "capturing counts of CDS and start and Stop codons for DEBUGGING ..."
	for ENST in $(cat CDS_only.bed start_and_stop_codons.sorted.bed | cut -f5 | grep -v stable_id | sort -u | tr '\n' ' '  ) ;
	do
		echo ${ENST} ; fgrep ${ENST} CDS_only.bed start_and_stop_codons.sorted.bed | cut -f7 | sort -u | uniq -c ;
	done > counts_CDS_start_stop_codons.txt
	echo -e "debug counts files generated; 1) counts_canonical_CDS_start_stop_codons.txt 2) counts_CDS_start_stop_codons.txt "
fi
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ DEBUG @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## get list of ALL the transcripts
echo -e "making the list of unique ENSTS from the file ${BED}"
cat ${CDS_FILE} | cut -f5 | grep -v stable_id | sort -u > list_all_ENST_IDs.txt
check_ev $? "get unique ENTS ids list"
dos2unix list_all_ENST_IDs.txt


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Removing Previous directory if exist
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo -e "removing pre-existing directories if any ... "
for DIR in  enst_cds enst_start_stop isec_cds_start_stop
do
	if [[ -e ${DIR} ]]
	then
		echo -e "\tremoving directory ${DIR} ..."
		rm -r ${DIR}
	fi
done
mkdir -p enst_cds enst_start_stop isec_cds_start_stop

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## STRATEGY is HERE: to make bedtools intersect working as it should and
## without having false positive we MUST treat each ENST independently
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##  making files with only CDS information for each ENST [runtime ~ 5-10 mins ]
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo -e "$(date)"
echo -e "parallel #1 --> enst_cds "
parallel -k --halt soon,fail=1 -j ${CPUS} fgrep {1} ${BED_BN/.bed/.cds_only.7cols.bed} ">" enst_cds/{1}.CDS.bed "||" true ::: $(cat list_all_ENST_IDs.txt | tr '\n' ' ' )
ev=$?
echo -e "ev for parallel #1 ==> $ev"



##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## making files with only START/STOP codons information
## for each ENST [runtime ~ 5-10 mins ]
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo -e "$(date)"
echo -e "parallel #2 --> enst_start_stop "
parallel -k --halt soon,fail=1 -j ${CPUS} fgrep {1} start_and_stop_codons.sorted.bed ">" enst_start_stop/{1}.StartStopCodons.bed "||" true ::: $(cat list_all_ENST_IDs.txt | tr '\n' ' ' )
ev=$?
echo -e "ev = $ev"
check_ev $? "parallel #2 enst_start_stop_codon"


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##  BEDTOOLS INTERSECT  CDS and START and STOP codons of EACH ENST  [runtime ~ 5 mins ]
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo -e "$(date)"
echo -e "parallel #3 --> bedtools intersect "
parallel -k --halt soon,fail=1 -j ${CPUS} bedtools intersect -a enst_cds/{1}.CDS.bed -b enst_start_stop/{1}.StartStopCodons.bed -wao -s ">" isec_cds_start_stop/isec_{1}_cds_start_stop.bed "||" true ::: $(cat list_all_ENST_IDs.txt | tr '\n' ' ' )
ev=$?
echo -e "ev = $ev"
check_ev $ev "parallel #3 bedtools intersect by enst"
## Known issue: the command "parallel" above returns an exit value of "101";
## The explanation is still a mystery.
## and that causes to comment out the line `set -eu -o pipefile` at the beginning of that script
## otherwise the script fails here, because of that 101 ev value from parallel
##@@@@@@@@@@@@@@@-------------------------------------------@@@@@@@@@@@@@@@@@@@@
## NOTE: when a Transcript have one or more CDS listed but no start or stop codons
## described or listed for that transcripts, the 'bedtools intersect' outputs a 11-columns
## bedfile instead of a 15-columns bedfile when intersection with start and stop codons
##@@@@@@@@@@@@@@@-------------------------------------------@@@@@@@@@@@@@@@@@@@@

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Check if the 3 folders have exactly the same number of files
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
NCDS=$(ls enst_cds/ | wc -l)
NSS=$(ls enst_start_stop/ | wc -l)
NISEC=$(ls isec_cds_start_stop/ | wc -l)

if [[ ${NCDS} -ne ${NSS} ]]
then
	echo -e "ERROR: expected the same NUMBER of files between enst_cds (${NCDS}) and the folder enst_start_stop (${NSS}). FAILED; Aborting; Check what went wrong."
	exit 1
elif [[ ${NSS} -ne ${NISEC} ]]
then
	echo -e "ERROR: expected the same NUMBER of files between isec_cds_start_stop (${NISEC}) and the folder enst_start_stop (${NSS}). FAILED; Aborting; Check what went wrong."
	exit 1
fi
echo -e "${NCDS} =?= ${NSS} =?= ${NISEC}"



##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##  BEDTOOLS INTERSECT  CDS and START and STOP codons runtime ~ 5 mins
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo -e "$(date)"
echo -e "concatenating the isec files ..."

## why using the _for_ loop? because _if_ we use "cat  isec_cds_start_stop/isec_*_cds_start_stop.bed > ${MBED}"
## we get the ERROR MESSAGE:
## -bash: /usr/bin/cat: Argument list too long
## old method too slow: for F in $(find isec_cds_start_stop/ -name "isec*_cds_start_stop.bed" ) ; do cat ${F} ; done > ${MBED}  ;
>${MBED}
find isec_cds_start_stop/ -maxdepth 1 -type f -name 'isec*_cds_start_stop.bed' -print0 | xargs -0 cat -- >> ${MBED}
check_ev $? "concat isec_cds_start_stop files"
echo -e "$(date)"


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Processing the MBED file
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## first we check if unexpected difference exist; If so we abort.
echo -e "checking if input is correctly formatted ..."
awk 'NR>1 { if( $14!="" && $13 !="." && $13!=$6 ){ err=1 } } END { exit err}' ${MBED}
check_ev $? "Unexpected Difference in Strand between CDS and START/STOP codon within LINE ${LCOUNT}: ${LINE};\ncheck you input ; should not have happened ; Aborting"

echo -e "$(date)"
## first, let us process ALL the CDS without any start or stop codon associated
echo -e "awk to capture the CDS lines and pad them ... and ... subsets the start_stop for future padding ..."
echo -e "subsetting isec files ... to improve runtime"
echo -e "awk-ing to capture both the CDS and the start_stop and dump them into their respective file ..."
## for CDS only, we pad both side as it is a regular CDS (no start or stop associated)
## for start and stop, we put them in a new file we will process later in a while loop.
echo -e "dump CDS_without_SS  and  CDS_with_SS into 2 different files ..."
cat ${MBED} | awk -v NBP=$NBP '{ if( $14=="." || $14 =="-1" || $14=="" ){ FS=OFS="\t" ; $2=$2-NBP ; if($2<0){$2=0} ; $3=$3+NBP ; print > "'${MBEDCDS}'" } else { print >"'${MBEDSS}'" } }'
check_ev $? "cat ${MBED} and awk and CDS padding and StSt dump"
echo -e "$(date)"

## Merge SS to check if start and stop overlap
sort -k1,1V -k2,2n -k3,3n ${MBEDSS} | bedtools merge -s -c 4,5,6,7,8,9,10,11,12,13,14 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct > ${MBEDSS}.bt_merge.bed
#MBEDSS=${MBEDSS}.bt_merge.bed

## Splitting File for Parallel use
echo -e "remove previous split files if any remaining"
rm "ssforcds___*"  "padded_ssforcds___*" &>/dev/null

echo -e "split the '${MBEDSS}' file having only the start and stop codon information ..."
split -l 10000 ${MBEDSS} ssforcds___
check_ev $? "split ${MBEDSS}"

## Processing the Padding for CDS that have start or/and stop codons relationship
echo -e "parallel #4 ... parse_start_stop_for_padding.sh "
parallel --delay 10 --halt soon,fail=1 -k -j ${CPUS} unbuffer bash parse_start_stop_for_padding.sh {1} ${NBP} "||" true ::: $(find -type f -name "ssforcds___*" | tr '\n' ' ' )
check_ev $? "parallel parse_start_stop_for_padding.sh"


## Merge the files generated by parse_start_stop_for_padding.sh; files have prefix: padded_ssforcds___
>${PADDED_BED}
echo -e "merging padded_ssforcds___ files ..."
find . -maxdepth 1 -type f -name 'padded_ssforcds___*' -print0 | xargs -0 cat -- >> ${PADDED_BED}
check_ev $? "merging padded_ssforcds___"


##@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Finalizing Results File
##@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo -e "$(date)"

if [[ ${OUTFILE} == "" ]]
then
	OUTFILE=${PADDED_BED/.bed.cds_padded.bed/.cds_padded.final.bed}
elif [[ ${OUTFILE##*.} != "bed" ]]
then
	OUTFILE=${OUTFILE}.bed
fi
check_ev $? "set OUTFILE name"


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## MERGING OR NOT MERGING THAT IS THE CHOICE
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if [[ ${DO_NOT_MERGE} == "no" ]]
then
	echo -e "cat CDS and start_stops_codons, sorting, merging and finalizing output file that contains the ${NBP}bp-padded CDS lines for ALL the transcripts in the given input GTF..."
	cat ${MBEDCDS} ${PADDED_BED} \
	| cut -f1-6,14 \
	| sort --parallel=${CPUS} -k1,1V -k2,2n -k3,3n \
	| awk '{FS=OFS="\t" ; if(NF==6){$6=$6"\t."} ; print $0 }' \
	| bedtools merge -c 4,5,6,7 -o distinct,distinct,distinct,distinct \
	> ${OUTFILE}
	ev_final=$?
	check_ev $? "make OUTFILE with merging overllapping codons"
else
	echo -e "cat CDS and start_stops_codons, sorting and finalizing UNmerged output file that contains the ${NBP}bp-padded CDS lines for ALL the transcripts in the given input GTF..."
	cat ${MBEDCDS} ${PADDED_BED} \
	| cut -f1-6,14 \
	| sort --parallel=${CPUS} -k1,1V -k2,2n -k3,3n \
	| awk '{FS=OFS="\t" ; if(NF==6){$6=$6"\t."} ; print $0 }' \
	> ${OUTFILE}
	ev_final=$?
	check_ev $? "make OUTFILE without merging overllapping codons"
fi

echo -e "output file with all transcripts is: ${OUTFILE}"


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## CLEANING
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## default value for CLEAN_TEMP_FILES is "no"
if [[ $ev_final -eq 0 && ${ev_final_canonical} -eq 0 && ${CLEAN_TEMP_FILES} == "yes" ]]
then
	echo -e "$(date)"
	echo -e "cleaning intermediate or/and temporary files ..."
	rm ${BED} ${MBED} header_line_50cols.txt ${PADDED_BED}* ${CDS_FILE}
	rm start_codon.bed stop_codon.bed stop_codon.padded_3bp.bed
	rm isec__all_ensts_cds_start_stop.bed.start_stop.bed
	rm isec__all_ensts_cds_start_stop.bed.padded_CDS.bed
	rm start_and_stop_codons.sorted.bed
	rm ${MBEDSS}.bt_merge.bed
	rm list_all_ENST_IDs.txt
	rm -r enst_cds/ &
	rm -r enst_start_stop/ &
	rm -r isec_cds_start_stop/ &
	echo -e "removing temp SGTF___*"
	rm sGTF__* &>/dev/null
	echo -e "cleaning temp files padded_ssforcds___* and ssforcds___*  "
	rm padded_ssforcds___*  ssforcds___* &>/dev/null
	wait
elif [[ ${ev_final} -ne 0 ]]
then
	echo -e "ERROR: script $0 FAILED; Final expected BED file FAILED; Aborting cleanup step."
	exit 1
fi

echo -e "$(date)"
echo -e "$0 completed successfully"


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## RELATED ADD ON INFORMATION BELOW
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## HEADER in BED FILE generated from GTF to BED conversion
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##seqname
##start
##end
##gene_id__gene_name
##score
##strand
##frame
##gene_version
##gene_source
##gene_biotype
##transcript_id
##source
##feature
##transcript_version
##transcript_name
##transcript_source
##transcript_biotype
##tag
##transcript_support_level
##exon_number
##exon_id
##exon_version
##protein_id
##protein_version
##ccds_id
