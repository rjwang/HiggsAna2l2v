#!bin/bash

for ifile in $(ls -d multicrab_*)
do 
	echo "***********************************"
	echo ${ifile}
	echo "***********************************"
	multicrab -status -c ${ifile}
done
