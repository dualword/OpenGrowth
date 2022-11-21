#!/bin/bash
#SBATCH -n 1
#SBATCH -t 24:00:00 
#SBATCH -p general
#SBATCH --mem-per-cpu=2500MB
#SBATCH -e Analysis.e
#SBATCH -o Analysis.o
#SBATCH -J Analysis

#This scripts look in all the MMPBSA calculations, get all the 5000 data in a single file and then plot a distribution of the energies.
#It also plots the normal distribution with the same mean and standard deviation.

module load viz/gnuplot-4.6.0

mkdir 00-Analysis
cd 00-Analysis

for Protein in 1HPV
do
    for Iteration in `seq 1 50`
    do
      cp ../${Protein}_${Iteration}/MMPBSA/Files_${Protein}_${Iteration}.tar.gz .
      if [ -s Files_${Protein}_${Iteration}.tar.gz ]; then
	tar -zvxf Files_${Protein}_${Iteration}.tar.gz
	rm Files_${Protein}_${Iteration}.tar.gz

	for File in complex receptor ligand
	do
		#Assuming that the MMPBSA was made on 20 cores. Otherwise, the list of i must be adapted.
		for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
		do
			grep "BOND" Files_${Protein}_${Iteration}/_MMPBSA_${File}_pb.mdout.$i >> ${Protein}_${Iteration}_${File}_1.dat
			grep "VDWAALS" Files_${Protein}_${Iteration}/_MMPBSA_${File}_pb.mdout.$i >> ${Protein}_${Iteration}_${File}_2.dat
			grep "1-4 VDW" Files_${Protein}_${Iteration}/_MMPBSA_${File}_pb.mdout.$i >> ${Protein}_${Iteration}_${File}_3.dat
			grep "ECAVITY" Files_${Protein}_${Iteration}/_MMPBSA_${File}_pb.mdout.$i >> ${Protein}_${Iteration}_${File}_4.dat
		done
		ArrayBOND=(` awk '{ print $3 ; }' ${Protein}_${Iteration}_${File}_1.dat `)
		ArrayANGLE=(` awk '{ print $6 ; }' ${Protein}_${Iteration}_${File}_1.dat `)
		ArrayDIHED=(` awk '{ print $9 ; }' ${Protein}_${Iteration}_${File}_1.dat `)
		ArrayVDWAALS=(` awk '{ print $3 ; }' ${Protein}_${Iteration}_${File}_2.dat `)
		ArrayEEL=(` awk '{ print $6 ; }' ${Protein}_${Iteration}_${File}_2.dat `)
		ArrayEPB=(` awk '{ print $9 ; }' ${Protein}_${Iteration}_${File}_2.dat `)
		Array14VDW=(` awk '{ print $4 ; }' ${Protein}_${Iteration}_${File}_3.dat `)
		Array14EEL=(` awk '{ print $8 ; }' ${Protein}_${Iteration}_${File}_3.dat `)
		ArrayRESTRAINT=(` awk '{ print $11 ; }' ${Protein}_${Iteration}_${File}_3.dat `)
		ArrayECAVITY=(` awk '{ print $3 ; }' ${Protein}_${Iteration}_${File}_4.dat `)
		ArrayEDISPER=(` awk '{ print $6 ; }' ${Protein}_${Iteration}_${File}_4.dat `)
		###
		SizeFile=` wc -l ${Protein}_${Iteration}_${File}_1.dat | awk '{ print $1 ; }' `
		x=0
		echo "BOND		ANGLE		DIHED		VDWAALS		EEL		EPB		14VDW		14EEL		RESTRAINT	ECAVITY		EDISPER" > ${Protein}_${Iteration}_${File}_DATA.dat
		while [ $x -lt ${SizeFile} ]; do
			echo "${ArrayBOND[$x]}	${ArrayANGLE[$x]}	${ArrayDIHED[$x]}	${ArrayVDWAALS[$x]}	${ArrayEEL[$x]}	${ArrayEPB[$x]}	${Array14VDW[$x]}	${Array14EEL[$x]}	${ArrayRESTRAINT[$x]}		${ArrayECAVITY[$x]}	${ArrayEDISPER[$x]}" >> ${Protein}_${Iteration}_${File}_DATA.dat
			Sum=` echo "scale=4;(${ArrayBOND[$x]} + ${ArrayANGLE[$x]} + ${ArrayDIHED[$x]} + ${ArrayVDWAALS[$x]} + ${ArrayEEL[$x]} + ${ArrayEPB[$x]} + ${Array14VDW[$x]} + ${Array14EEL[$x]} + ${ArrayRESTRAINT[$x]} + ${ArrayECAVITY[$x]} + ${ArrayEDISPER[$x]})" | bc `
			echo "${Sum}" >> ${Protein}_${Iteration}_${File}_SUM.dat
			x=`expr $x + 1`
		done
		rm ${Protein}_${Iteration}_${File}_1.dat ${Protein}_${Iteration}_${File}_2.dat ${Protein}_${Iteration}_${File}_3.dat ${Protein}_${Iteration}_${File}_4.dat
	done

	########

	ArrayComplex=(` awk '{ print $1 ; }' ${Protein}_${Iteration}_complex_SUM.dat `)
	ArrayReceptor=(` awk '{ print $1 ; }' ${Protein}_${Iteration}_receptor_SUM.dat `)
	ArrayLigand=(` awk '{ print $1 ; }' ${Protein}_${Iteration}_ligand_SUM.dat `)
	SizeFile=` wc -l ${Protein}_${Iteration}_complex_SUM.dat | awk '{ print $1 ; }' `
	x=0
	while [ $x -lt ${SizeFile} ]; do
		Diff=` echo "scale=4;(${ArrayComplex[$x]} - ${ArrayReceptor[$x]} - ${ArrayLigand[$x]} )" | bc `
		echo "${Diff}" >> ${Protein}_${Iteration}_DeltaE.dat
		x=`expr $x + 1`
	done
	cat ${Protein}_${Iteration}_DeltaE.dat >> ${Protein}_All_DeltaE.dat
	rm ${Protein}_${Iteration}_DeltaE.dat

	########

	rm ${Protein}_${Iteration}_complex_SUM.dat ${Protein}_${Iteration}_receptor_SUM.dat ${Protein}_${Iteration}_ligand_SUM.dat
	rm ${Protein}_${Iteration}_complex_DATA.dat ${Protein}_${Iteration}_receptor_DATA.dat ${Protein}_${Iteration}_ligand_DATA.dat
	rm -rf Files_${Protein}_${Iteration}/
      fi
    done

    sed "s/LIGANDNAME/${Protein}/g" ../00-Files/TEMPLATE_Distribution.gplot > ${Protein}_All.gplot
    gnuplot ${Protein}_All.gplot
    sleep 1s
done

