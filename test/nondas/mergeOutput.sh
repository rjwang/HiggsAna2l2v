#!/bin/bash

totalOutputs=`ls ${1}/res/analysis_*.root | wc -l`
nsplit=$2
echo "****************************************"
echo "Integrating $1 PD"
echo "$totalOutputs found - will be merged into ${nsplit} file"
echo "****************************************"

let "step=${totalOutputs}/${nsplit}"
let "temp=${totalOutputs}%${nsplit}"
let "finalstep=(${totalOutputs}+$temp)/${nsplit}"

#echo "$finalstep"
#echo "total: ${totalOutputs}"
#echo "nsplit: ${nsplit}" 
#echo "temp: ${temp}" 
#echo "step: $step"


if [ $nsplit == 1 ]; then
    echo "nsplit: $nsplit"
    for i in `seq 0 0`; do
    	let "startOutput=${i}*$finalstep+1"
    	let "endOutput=(${i}+1)*$finalstep"

    	echo "${startOutput}-${endOutput}"

    	outputList=""
    	for j in `seq ${startOutput} ${endOutput}`; do
          newFile="${1}/res/analysis_${j}_*.root"
          if [ -e ${newFile} ]; then
            	outputList="${outputList} ${newFile}"
          fi
    	done
    	hadd -f /tmp/rewang/${1}.root ${outputList}
    done
fi


if [ $nsplit != 1 ]; then
   let "nn=${nsplit}-1"
   for i in `seq 0 ${nn}`; do
   	let "startOutput=${i}*$finalstep+1"
        let "endOutput=(${i}+1)*$finalstep"

    	echo "${startOutput}-${endOutput}"

    	outputList=""
    	for j in `seq ${startOutput} ${endOutput}`; do
	  newFile="${1}/res/analysis_${j}_*.root"
	  if [ -e ${newFile} ]; then
	    	outputList="${outputList} ${newFile}"
	  fi
    	done
    	hadd -f /tmp/rewang/${1}_${i}.root ${outputList}
   done
fi

