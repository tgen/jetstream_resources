#!/usr/bin/env bash


##@@@@@@@@@@@@@@
## FUNCTION(S)
##@@@@@@@@@@@@@
function check_ev(){
	local ev=$1
	local msg=$2
	if [[ $ev -ne 0 ]]
	then
		echo -e "ERROR: ${msg} FAILED. Aborting";
		exit 1
	fi
}

## example of MBEDSS file content
## (15 or 11 columns after intersecting the CDS with their respective start and stop codons):
## chr19   35119376        35119416        ENSG00000089356___FXYD3 ENST00000603181 +       CDS     chr19   35119376        35119379        ENSG00000089356___FXYD3 ENST00000603181 +       start_codon     3
# chr19   35123440        35123454        ENSG00000089356___FXYD3 ENST00000603181 +       CDS     chr19   35123451        35123454        ENSG00000089356___FXYD3 ENST00000603181 +       stop_codon      3
# chr22   30299475        30299551        ENSG00000099992___TBC1D10A      ENST00000437122 -       CDS     chr22   30299475        30299478        ENSG00000099992___TBC1D10A      ENST00000437122 -       stop_codon      3
# chr22   30326672        30326881        ENSG00000099992___TBC1D10A      ENST00000437122 -       CDS     chr22   30326878        30326881        ENSG00000099992___TBC1D10A      ENST00000437122 -       start_codon     3

##@@@@@@@@@@@@@@@@@@@@
## MAIN INPUT FILE
##@@@@@@@@@@@@@@@@@@@@
MBEDSS=$1
NBP=$2

BN_MBEDSS=$(basename ${MBEDSS})
PADDED_BED=padded_${BN_MBEDSS}
## prepare data and variable for the while loop
>${PADDED_BED}  ; 	##init output file
TOT_LINES_MBEDSS=$(cat ${MBEDSS} | wc -l);	## total number of lines in start and stop file (MBEDSS)
LCOUNT=0  ; 	##init Line Counter
STEP=1000 ; 	## step counter ## HARDCODED

echo -e "parsing '${MBEDSS}' to add padding ... \nlooping over the ${TOT_LINES_MBEDSS} lines (slow process)"


while read LINE
do
	LCOUNT=$((LCOUNT+1))
	if [[ $((${LCOUNT} % ${STEP})) -eq 0 ]]
	then
		PCT=$( echo "${LCOUNT}/${TOT_LINES_MBEDSS}*100" | bc -l | awk '{printf "%0.2f",$1}' )
		echo -e "from ${BN_MBEDSS} : processed ... ${LCOUNT} lines ( ~${PCT}% ) " 1>&2
	fi

	S=$(echo -e "${LINE}" | cut -f14 ) ; 		## S for string of interest
	STRAND=$(echo -e "${LINE}" | cut -f6 ) ; 	## STRAND is the strand associated to the CDS
	STRAND2=$(echo -e "${LINE}" | cut -f13 ) ; 	## STRAND2 is the strand found for start and stop codons

	##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	## important NOTE: in the isec MBED file you will have lines with 15 columns
	## (meaning ENST had at least one start or one stop_codon)
	## if the isec line has only 11 columns this mean NO start or stop codon were
	## available for that ENST; therefore, we have a value "-1" instead of ".",
	## due to the way bedtools intersect manages the absence of feature in B file.
	##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	if [[ ${S} == "start_codon" ]]
	then
		#echo -e "start_codon found ... "
		if [[ ${STRAND} == "+" ]];
		then
			echo -e "${LINE}" | awk -v NBP=$NBP '{FS=OFS="\t" ; $3=$3+NBP ; print $0}'
		else
			echo -e "${LINE}" | awk -v NBP=$NBP '{FS=OFS="\t" ; $2=$2-NBP ; if($2<0){ $2=0 } ; print $0}'
		fi
	elif [[ ${S} == "stop_codon" ]]
	then
		# echo -e "stop_codon found ... "
		if [[ ${STRAND} == "+" ]];
		then
			echo -e "${LINE}" | awk -v NBP=$NBP '{FS=OFS="\t" ; $2=$2-NBP ; if($2<0){ $2=0 } ; print $0}'
		else
			echo -e "${LINE}" | awk -v NBP=$NBP '{FS=OFS="\t" ; $3=$3+NBP ; print $0}'
		fi
	elif [[ ${S} == "start_codon,stop_codon" || ${S} == "stop_codon,start_codon" ]]
	then
		## we just copy the line as is in the output file b/c the CDS has both start and stop codon
		echo -e "${LINE}"

	else
		echo -e "ERRORwithLINE: ${LINE}"  1>&2
		echo -e "ERROR: we should never have reached here; this means the value in expected column 13 is WRONG; check input; " 1>&2 ;
		echo -e "unexpected VALUE is: ${S} in  LINE: ${LINE}" 1>&2
		exit 1
	fi
done < ${MBEDSS} > ${PADDED_BED}
## we did not use the line ' done < ${MBEDSS} > ${PADDED_BED}' because
## it is actually slower than the current implementation
check_ev $? "while loop with ${MBEDSS}"

## cleaning the input file since we do not need it anymore
rm ${MBEDSS}
