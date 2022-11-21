#!/bin/bash

#Define the command to launch OpenBabel
OPENBABEL_EXE=obabel
GAUSSIAN_EXE=g09

File=PrepareFragments
sed '1d' ${File}.dat > ${File}_Temp.dat

#Store the SMILES strings and the names
Size=` wc -l ${File}_Temp.dat | awk '{ print $1 ; }' `
Smiles=(` awk '{ print $1 ; }' ${File}_Temp.dat `)
Name=(` awk '{ print $2 ; }' ${File}_Temp.dat `)

#Prepare the files for the optimization of geometry
x=0
while [ $x -lt ${Size} ]; do
	### Convert to 3D with OpenBabel
	echo "${Smiles[$x]}	${Name[$x]}" > ${Name[$x]}.smi
	${OPENBABEL_EXE} ${Name[$x]}.smi -O ${Name[$x]}.xyz --gen3d -p 7
	### Get the charges. This doesn't always work so you should check both the structure and the charges
	${OPENBABEL_EXE} ${Name[$x]}.xyz -O ${Name[$x]}.mol2
        Array=(` awk '{ print $9 ; }' ${Name[$x]}.mol2 `)
        SizeMolecule=` wc -l ${Name[$x]}.mol2 | awk '{ print $1 ; }' `
        y=0
        Charge=0
        while [ $y -lt ${SizeMolecule} ]; do
                Charge=` echo "${Charge} ${Array[$y]}" | awk '{print ($1+$2)}'`
                y=`expr $y + 1`
        done
        Charge=` echo "${Charge}" | awk '{printf "%.0f", $1}'`
	### Create the Gaussian input files
	echo "%Nprocshared=16" >> ${Name[$x]}.com
	echo "%mem=10000MB" >> ${Name[$x]}.com
	echo "#P M062X/aug-cc-pVTZ Opt SCF=(MaxCycle=500)" >> ${Name[$x]}.com
	echo "# SCRF Test Units(Ang,Deg)" >> ${Name[$x]}.com
	echo " " >> ${Name[$x]}.com
	echo "${Name[$x]}" >> ${Name[$x]}.com
	echo " " >> ${Name[$x]}.com
	echo "${Charge}	1" >> ${Name[$x]}.com
	sed '1,2d' ${Name[$x]}.xyz >> ${Name[$x]}.com
	echo " " >> ${Name[$x]}.com
	### Clean 
	rm ${Name[$x]}.smi ${Name[$x]}.xyz ${Name[$x]}.mol2
	echo "${Name[$x]} done"
	x=`expr $x + 1`
done
rm ${File}_Temp.dat 

#We only need to optimize the fragments called _1.com. The others (_2.com for example) will have the same geometries
#So we keep only the *_1.com files
mkdir Temp3D
mv *_1.com Temp3D/.
rm *.com
mv Temp3D/*_1.com .
rmdir Temp3D/

#Optimize
for File in ` ls *.com `
do
	Name=${File/.com}
	${GAUSSIAN_EXE} ${Name}.com ${Name}.log
	sleep 5s
	${OPENBABEL_EXE} ${Name}.log -O ${Name}.xyz
done

