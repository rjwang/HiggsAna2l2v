#!/bin/bash

###
# submitJobToBatch.sh
# wrapper for the submission of jobs for the lxbatch
###

USER=`whoami`
HOSTNAME=`hostname`
HOSTINFO=`host ${HOSTNAME}`
BASEDIR=`pwd`;

#
# PARSE CONFIGURATION PARAMETERS
#
for p in $*; 
    do
    if   [ "$p" != "${p/-q/}"  ]; then QUEUE="${p/-q/}"; 
    elif [ "$p" != "${p/-R/}"  ]; then REQUIREMENT="${p/-R/}"; 
    elif [ "$p" != "${p/-J/}"  ]; then JOBNAME="${p/-J/}";
    elif [ "$p" != "${p/.sh/}" ]; then SCRIPT="${p}"; 
    elif [ "$p" != "${p/-G/}"  ]; then LOG="${p/-G/}";
    else PARAMS="$PARAMS $p"
    fi
done

if [[ -z $SCRIPT ]]; then
    echo "Must provide a script to configure+run your job (~*.sh)"
    echo "SubmitJobToBatch.sh <-qQUEUE> <-Rrequirements> <-Jjobname> [SCRIPT TO RUN] <SCRIPT PARAMETERS> "
    echo ""
    exit -1;
else
    cd 
    if [[ -e $SCRIPT ]]; then
	#echo "$SCRIPT was found" 
	cd -
    else 
	SCRIPT="${BASEDIR}/${SCRIPT}"
	echo $SCRIPT
	if [[ -e $SCRIPT ]]; then
	    #echo "$SCRIPT was found"
	    cd -
	else
	    echo "Unable to find script to submit to batch"
	    cd -
	    exit -1 
	fi
    fi
fi

#
# CONFIGURE
#
if [[ -z $QUEUE ]]; then 
#    QUEUE="1nd"
    QUEUE="8nh"
#    QUEUE="2nd"
fi

#
# SUBMIT JOB
#
#echo "Submitting shell script: ${SCRIPT} with parameters: ${PARAMS} to queue $QUEUE"
#<<<<<<< submit2batch.sh
#bsub -q $QUEUE -R "type==SLC5_64 && pool>30000 && tmp>50000" `echo ${SCRIPT} ${PARAMS}`
echo -e "\e[0;35m Submitting shell script: bsub -q $QUEUE -o $LOG -R $REQUIREMENT -J $JOBNAME `echo ${SCRIPT} ${PARAMS}` \e[0m"
######
bsub -sp 99 -q $QUEUE -o $LOG -R "$REQUIREMENT" -J $JOBNAME `echo ${SCRIPT} ${PARAMS}`
######
#
# END
#
