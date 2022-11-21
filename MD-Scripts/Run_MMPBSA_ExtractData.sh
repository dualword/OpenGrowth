#!/bin/bash

#This script returns the different contributions for the MMPBSA and Entropy analysis.

echo "Protein	MMPBSAComplex	MMPBSAReceptor	MMPBSALigand	EntropyComplex	EntropyReceptor	EntropyLigand	DeltaG	#MMPBSA	#Entropy"
for Protein in 1HPV
do
	SumEntropyComplex=0
	SumEntropyReceptor=0
	SumEntropyLigand=0
	SumMMPBSAComplex=0
	SumMMPBSAReceptor=0
	SumMMPBSALigand=0
	x=0
	y=0
	for Iteration in `seq 1 50`
	do
		if [ -s ${Protein}_${Iteration}/MMPBSA/MMPBSA_${Protein}_${Iteration}.dat ]; then
			grep "TOTAL" ${Protein}_${Iteration}/MMPBSA/MMPBSA_${Protein}_${Iteration}.dat > MMPBSA.dat
			ArrayMMPBSA=(` awk '{ print $2 ; }' MMPBSA.dat `)
			SumMMPBSAComplex=` echo "scale=4;(${SumMMPBSAComplex} + ${ArrayMMPBSA[0]})" | bc `
			SumMMPBSAReceptor=` echo "scale=4;(${SumMMPBSAReceptor} + ${ArrayMMPBSA[1]})" | bc `
			SumMMPBSALigand=` echo "scale=4;(${SumMMPBSALigand} + ${ArrayMMPBSA[2]})" | bc `
			rm MMPBSA.dat
			x=`expr $x + 1`
		fi

		if [ -s ${Protein}_${Iteration}/MMPBSA/MMPBSA_${Protein}_${Iteration}_Entropy.dat ]; then
			grep "Total" ${Protein}_${Iteration}/MMPBSA/MMPBSA_${Protein}_${Iteration}_Entropy.dat > Entropy.dat
			ArrayEntropy=(` awk '{ print $2 ; }' Entropy.dat `)
			SumEntropyComplex=` echo "scale=4;(${SumEntropyComplex} + ${ArrayEntropy[0]})" | bc `
			SumEntropyReceptor=` echo "scale=4;(${SumEntropyReceptor} + ${ArrayEntropy[1]})" | bc `
			SumEntropyLigand=` echo "scale=4;(${SumEntropyLigand} + ${ArrayEntropy[2]})" | bc `
			rm Entropy.dat
			y=`expr $y + 1`
		fi
	done

	AverageEntropyComplex=` echo "scale=6;(${SumEntropyComplex} / ${y})" | bc `
	AverageEntropyReceptor=` echo "scale=6;(${SumEntropyReceptor} / ${y})" | bc `
	AverageEntropyLigand=` echo "scale=6;(${SumEntropyLigand} / ${y})" | bc `
	AverageMMPBSAComplex=` echo "scale=6;(${SumMMPBSAComplex} / ${x})" | bc `
	AverageMMPBSAReceptor=` echo "scale=6;(${SumMMPBSAReceptor} / ${x})" | bc `
	AverageMMPBSALigand=` echo "scale=6;(${SumMMPBSALigand} / ${x})" | bc `

	DeltaG=` echo "scale=2; ( (${AverageMMPBSAComplex} - ${AverageMMPBSAReceptor} - ${AverageMMPBSALigand}) - (${AverageEntropyComplex} - ${AverageEntropyReceptor} - ${AverageEntropyLigand}) )" | bc `
	echo "${Protein}	${AverageMMPBSAComplex}	${AverageMMPBSAReceptor}	${AverageMMPBSALigand}	${AverageEntropyComplex}	${AverageEntropyReceptor}	${AverageEntropyLigand}	${DeltaG}	${x}	${y}"
done

