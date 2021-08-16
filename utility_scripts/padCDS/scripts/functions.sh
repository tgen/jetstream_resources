##@@@@@@@@@@@@@@
## FUNCTION(S)
##@@@@@@@@@@@@@
function usage(){
	echo -e "USAGES"
	echo -e "USAGE1: $0 --dirwork \"\${DIR_WORK}\" --gtf '\${A_GTF_FILE}' --outfilename '\${OUTFILE}' --cpus '\$CPUS'  --nbp '\$NumberOfPaddingToAdd' --cleanup"
	echo -e "USAGE2: $0 -d \"\${DIR_WORK}\" -g '\${A_GTF_FILE}' -o '\${OUTFILE}' -t '\$CPUS'  -n '\$NumberOfPaddingToAdd' -c \n"
	echo -e "-d|--dirwork|--dir_work)   Working Directory where all the temp files are written and where the output file is written if no full/rel path for the outfilename is given"
	echo -e "-g|--gtf)   must be a GTF file"
	echo -e "-o|--outfilename)   is the name to the outfile name for the padded file with ALL the CDS from ALL the transcripts; extension recommended is '.bed'; default is: \${GTF/.gtf/.cds_padded.final.bed} "
	echo -e "-t|--cpus|--threads)   is the number of cpus allocated to the script to run faster the parallel steps; Default is 2 but recommended value is a minimum of 10; higher the better"
	echo -e "-n|--nbp)   is the number padding to be added to the CDS, start and stop codons: default is 2"
	echo -e "-c|--cleanup)   is The CLEAN_TEMP_FILES argument to specify if we want to remove ALL the temporary files and folders and only keep the padded files. Flag to say yes I want to keep the temp files and folders"
	echo -e "-m|--do_not_merge) will skip the last step of that tool which run 'bedtools merge' in order to merge the overlapping intervals; default is to merge overlapping intervals"
	exit
}

function requirements(){
	echo -e "bedtools, python3, bgzip, and parse_gtf_to_bed_for_all_lines_and_features.py (as executable) MUST be in your PATH."
	echo -e "make sure all tools are 'executables' with 'chmod a+x' "
}

function check_ev(){
	local ev=$1
	local msg=$2
	if [[ $ev -ne 0 ]]
	then
		echo -e "ERROR: ${msg} FAILED. Aborting";
		exit 1
	fi
}

function checkDir(){
	local D=$1
	if [[ ! -e ${D} ]] ; then echo -e "DIR NOT FOUND << ${D} >>; Aborting!" ; usage ; exit -1 ; fi ;
}
function checkFile(){
	local F=$1
	if [[ ! -e ${F} ]] ; then echo -e "FILE NOT FOUND << ${F} >>; Aborting!" ; usage ; exit -1 ; fi ;
}

function checkFileName(){

	for F in $@
	do
		if [[ ${F} == "" ]] ;
		then
			echo -e "ERROR: Missing filenames for one of the variable which is REQUIRED; Aborting."
			exit -1
		fi
	done
}

function init_some_vars(){
	LI="RECAP_INPUTS_USED:"
	FP_SNPEFF_JAR=""
	##@@@@@@@@@@@@@@@@@@@@@@@@@@@
	## DEFAULTS INPUTS VALUES
	##@@@@@@@@@@@@@@@@@@@@@@@@@@@
	DIR_WORK=""
	GTF=""
	OUTFILE=""
	CPUS=2
	NBP=2
	CLEAN_TEMP_FILES="no" 		## values can be 'yes' or 'no'
	DO_NOT_MERGE="no"		## Will not use bedtools merge to merge overlapping intervals; the user will have to do it manually if needed; Default is to merge the overlapping codons

}

function getOptions(){
# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -o hcd:g:l:t:o:n: -l help,dirwork:,dir_work:,gtf:,cpus:,threads:,nbp:,cleanup,do_not_merge,outfilename: -- "$@" )
	then
	# something went wrong, getopt will put out an error message for us
		echo "ERROR in Arguments" ; usage
		exit -1
	fi
	eval set -- "$options"
	while [[ $# -gt 0 ]]
	do
		# for options with required arguments, an additional shift is required
		case $1 in
		-d|--dirwork|--dir_work) DIR_WORK="$2" ; LI="${LI}\nDIR_WORK==\"${DIR_WORK}\"" ; shift ;;
		-g|--gtf)  GTF="$(echo ${2} | sed 's/,/ /g')" ;  LI="${LI}\nGTF==\"${GTF}\"";  shift ;;
		-t|--cpus|--threads)  CPUS="$(echo ${2} | sed 's/,/ /g')" ;  LI="${LI}\nCPUS==\"${CPUS}\"";  shift ;;
		-o|--outfilename) OUTFILE="$2" ; LI="${LI}\nOUTFILE==\"${OUTFILE}\"";  shift ;;
		-n|--nbp) NBP=$2 ; LI="${LI}\nNBP==\"${NBP}\"";  shift ;;
		-c|--cleanup)  CLEAN_TEMP_FILES="yes" ; LI="${LI}\nCLEAN_TEMP_FILES==\"${CLEAN_TEMP_FILES}\""  ;;
		--do_not_merge) DO_NOT_MERGE="yes" ; LI="${LI}\nDO_NOT_MERGE==\"${DO_NOT_MERGE}\""  ;;
		-h|--help) usage ; requirements ; exit ;;
		(--) shift ; break ; echo "--" ;;
		(-*) echo -e "$0: error - unrecognized option $1\n\n" 1>&2   ; usage; requirements ; exit 1  ;;
		(*) break ; echo "$0: error --- unrecognized option $1" 1>&2 ; usage; requirements ; exit 1  ;;
		esac
		shift
	done

	#input recap
	LI="${LI}\nCURDIR==\"${PWD}\""
	echo -e "\n\n+------------------------------------------------+\n${LI[@]}\n+------------------------------------------------+\n\n"

}
