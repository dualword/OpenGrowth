#!/bin/bash

# This script reads a file PrepareFragments.dat and first checks that the SMARTS patterns with the interesting atoms on the right or on the left give the same answer,
# and also that the SMARTS patterns with the interesting atoms on the left and the 3Mer patterns give the same answer. If it succeeds, it creates the needed files for
# the GUI. The PrepareFragments.dat must have been properly formatted, and the .xyz files for each fragment must have been prepared before and placed in the XYZ folder.
# An executable for OpenBabel must be available.

File=PrepareFragments
Library=ChEMBL_Subset-8018.smi
OPENBABEL_EXE=obabel

sed '1d' ${File}.dat > ${File}_Temp.dat
Size=` wc -l ${File}_Temp.dat | awk '{ print $1 ; }' `

Smiles=(` awk '{ print $1 ; }' ${File}_Temp.dat `)
Name=(` awk '{ print $2 ; }' ${File}_Temp.dat `)
MW=(` awk '{ print $3 ; }' ${File}_Temp.dat `)
SmartsRight=(` awk '{ print $4 ; }' ${File}_Temp.dat `)
SmartsLeft=(` awk '{ print $5 ; }' ${File}_Temp.dat `)

#We need 2 arrays for the 3Mers
SmartsLeft3Mer=(` awk '{ print $6 ; }' ${File}_Temp.dat `)
sed -i "s/(y)//g" ${File}_Temp.dat
SmartsLeft3MerNoY=(` awk '{ print $6 ; }' ${File}_Temp.dat `)
rm ${File}_Temp.dat

x=0
while [ $x -lt ${Size} ]; do
	./ProcessFragments.exe 2 ${SmartsRight[$x]} ${Library} >> ${File}_Right.txt
	./ProcessFragments.exe 2 ${SmartsLeft[$x]}  ${Library} >> ${File}_Left.txt
	./ProcessFragments.exe 2 ${SmartsLeft3MerNoY[$x]}  ${Library} >> ${File}_3Mer.txt
	x=`expr $x + 1`
done

NumberErrors1=` diff ${File}_Right.txt ${File}_Left.txt | wc -l `
NumberErrors2=` diff ${File}_Left.txt ${File}_3Mer.txt | wc -l `

if [[ ${NumberErrors1} == 0 && ${NumberErrors2} == 0 ]]; then
	mkdir Fragments-GUI
	x=0
	while [ $x -lt ${Size} ]; do
		#Get the short name of the fragment: for Ring2_54_1, we want to keep Ring2_54 (e.g.).
		EndOfName=${Name[$x]:6}
		y=`expr index "${EndOfName}" _`
		Position=`expr 6 + ${y} - 1`
		ShortName=${Name[$x]:0:${Position}}
		#Create the needed folder.
		mkdir Fragments-GUI/${Name[$x]}
		echo "${Smiles[$x]}" > Fragments-GUI/${Name[$x]}/${Name[$x]}.smi
		echo "${Name[$x]}	${MW[$x]}" > Fragments-GUI/${Name[$x]}/${Name[$x]}.dat
		echo "${SmartsRight[$x]}	${SmartsLeft[$x]}	${SmartsLeft3Mer[$x]}" > Fragments-GUI/${Name[$x]}/${Name[$x]}.smarts
		${OPENBABEL_EXE} Fragments-GUI/${Name[$x]}/${Name[$x]}.smi -O Fragments-GUI/${Name[$x]}/${ShortName}.svg
		sed -i '/rect /d' Fragments-GUI/${Name[$x]}/${ShortName}.svg
		cp XYZ/${Name[$x]}.xyz Fragments-GUI/${Name[$x]}/.
		echo "${Name[$x]}" >> Fragments-GUI/Fragments.dat
		x=`expr $x + 1`
	done
elif [ ${NumberErrors1} != 0 ]; then
	diff ${File}_Right.txt ${File}_Left.txt
elif [ ${NumberErrors2} != 0 ]; then
	diff ${File}_Left.txt ${File}_3Mer.txt
fi

