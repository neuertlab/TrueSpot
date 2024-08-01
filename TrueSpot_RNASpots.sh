#!/bin/bash

#Determine directory script is running in...
#https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
get_script_dir(){
	#This version is supposed to resolve symlinks
	SRC=${BASH_SOURCE[0]}
	while [ -L "$SRC" ]; do
		MY_DIR=$( cd -P "$( dirname "$SRC" )" >/dev/null 2>&1 && pwd )
		SRC=$(readlink "$SRC")
		[[ $SRC != /* ]] && SRC=$MY_DIR/$SRC
	done
	MY_DIR=$( cd -P "$( dirname "$SRC" )" >/dev/null 2>&1 && pwd )
}

usage(){
	get_script_dir
	if [ -n "${MY_DIR}" ]; then
		cat ${MY_DIR}/doc/bash/bash_man_rnaspots.txt
	else
		echo "ERROR - Manual file could not be found!"
	fi
}

#Collect arguments and reformat into a string to pass to MATLAB
MAT_ARG_STR=""
while [[ $# -gt 0 ]]; do
	#Check to see if it is a special argument first.
	SKIP_ME=0
	case $1 in
		-h|-help|-man)
			usage
			exit 0
			;;
		-input)
			IN_PATH="$2"
			;;
		-log)
			LOG_PATH="$2"
			SKIP_ME=1
			shift
			;;
		-outdir)
			OUT_DIR="$2"
			;;
		-outstem)
			OUT_STEM="$2"
			;;
	esac
	
	if [[ $SKIP_ME -eq 0 ]]; then
		if [ -n "$NOT_FIRST" ]; then
			MAT_ARG_STR="${MAT_ARG_STR}, "
		fi
		MAT_ARG_STR=${MAT_ARG_STR}\'${1}\'
		NOT_FIRST="true"
	fi
	shift
done

if [ -z "$IN_PATH" ]; then
	echo "Input image path is required!"
	usage
	exit 1
fi

#Determine log file path
if [ -z "$LOG_PATH" ]; then
	echo "Log path not specified."
	if [ -n "$OUT_STEM" ]; then
		LOG_PATH="${OUT_STEM}_matlog.log"
	elif [ -n "$OUT_DIR" ]; then
		LOG_PATH="${OUT_DIR}/MAT_Main_RNASpots.log"
	else
		LOG_PATH="${IN_PATH}_MAT_rnaspots.log"
	fi
	echo "Log path has been set to: ${LOG_PATH}"
fi

#Call MATLAB
get_script_dir
SCRIPT_NAME="Main_RNASpots"
matlab -nodisplay -nosplash -logfile "${LOG_PATH}" -r "cd '${MY_DIR}/src'; ${SCRIPT_NAME}(${MAT_ARG_STR}); quit;"