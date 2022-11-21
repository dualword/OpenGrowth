#!/bin/bash

#This scripts check if each jobs have succeeded, and writes
#in a new submission file (ending with New.sh) which one
#should be restarted.

for Protein in 1HPV
do
	for Iteration in `seq 1 50`
	do
		mv ${Protein}_${Iteration}.e ${Protein}_${Iteration}.o ${Protein}_${Iteration}/.
		Success=` grep "Finished" ${Protein}_${Iteration}/08-FullMD.log | wc -l `
		if [ ${Success} == 0 ]; then
		   echo "${Iteration}" >> Temp.txt
		fi
	done
	sed -i ':a;N;$!ba;s/\n/,/g' Temp.txt
	List=$(sed '1!d' Temp.txt)
	sed "s/--array=1-50/--array=${List}/g" Run_03-MD_${Protein}.sh > Run_03-MD_${Protein}_New.sh
	chmod u+x Run_03-MD_${Protein}_New.sh
	sleep 1s
	rm Temp.txt
done

