#!/bin/bash

#This scripts check if each jobs have succeeded, and writes
#in a new submission file (ending with New.sh) which one
#should be restarted.

for Protein in 1HPV
do
	for Iteration in `seq 1 50`
	do
		mv ${Protein}_${Iteration}_Analysis.e ${Protein}_${Iteration}_Analysis.o ${Protein}_${Iteration}/.
		Success=` ls ${Protein}_${Iteration}/MMPBSA/MMPBSA_${Protein}_${Iteration}.dat | wc -l `
		if [ ${Success} == 0 ]; then
		   echo "${Iteration}" >> Temp.txt
		fi
	done
	sed -i ':a;N;$!ba;s/\n/,/g' Temp.txt
	List=$(sed '1!d' Temp.txt)
	sed "s/--array=1-50/--array=${List}/g" Run_04-MMPBSA_${Protein}.sh > Run_04-MMPBSA_${Protein}_New.sh
	chmod u+x Run_04-MMPBSA_${Protein}_New.sh
	sleep 1s
	rm Temp.txt
done

