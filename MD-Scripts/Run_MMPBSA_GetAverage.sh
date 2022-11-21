#!/bin/bash

#This script computes the average energy values.

echo "Protein	#MMPBSA	MMPBSA	#S DeltaG"
for Protein in 1HPV
do
	SumMMPBSA=0
	SumEntropy=0
	xMMPBSA=0
	xEntropy=0
	for Iteration in `seq 1 50`
	do
		if [ -s ${Protein}_${Iteration}/MMPBSA/MMPBSA_${Protein}_${Iteration}.dat ]; then
			Value=` grep "DELTA TOTAL" ${Protein}_${Iteration}/MMPBSA/MMPBSA_${Protein}_${Iteration}.dat | awk '{ print $3 ; }' `
			SumMMPBSA=` echo "scale=4;(${SumMMPBSA} + ${Value})" | bc `
			xMMPBSA=`expr ${xMMPBSA} + 1`
		fi
		if [ -s ${Protein}_${Iteration}/MMPBSA_${Protein}_${Iteration}_Entropy.dat ]; then
			Value=` grep "DELTA S total" ${Protein}_${Iteration}/MMPBSA_${Protein}_${Iteration}_Entropy.dat | awk '{ print $4 ; }' `
			SumEntropy=` echo "scale=4;(${SumEntropy} + ${Value})" | bc `
			xEntropy=`expr ${xEntropy} + 1`
		fi
	done
	AverageMMPBSA=` echo "scale=6;(${SumMMPBSA} / ${xMMPBSA})" | bc `
	AverageEntropy=` echo "scale=6;(${SumEntropy} / ${xEntropy})" | bc `
	DeltaG=` echo "scale=2; ( ${AverageMMPBSA} - ${AverageEntropy} )" | bc `
	echo "${Protein}	${xMMPBSA}	${AverageMMPBSA}	${xEntropy} ${DeltaG}"
done

