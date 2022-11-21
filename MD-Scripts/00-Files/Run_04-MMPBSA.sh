#!/bin/bash
#SBATCH -n 20
#SBATCH -t 01-00:00:00 
#SBATCH -p general
#SBATCH --mem-per-cpu=2500MB
#SBATCH --array=1-50
#SBATCH -J LIGANDNAME_Analysis
#SBATCH -e LIGANDNAME_%a_Analysis.e
#SBATCH -o LIGANDNAME_%a_Analysis.o

module load centos6/gcc-4.8.2
module load centos6/openmpi-1.6.5_gcc-4.8.0
source /n/home04/ncheron/Programs/Gromacs_5.0.2/MPI/bin/GMXRC
#This Gromacs version was installed with the following cmake command:
#cmake .. -DGMX_MPI=ON -DGMX_GPU=OFF -DGMX_THREAD_MPI=ON -DGMX_OPENMP=ON -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=/n/home04/ncheron/Programs/Gromacs_5.0.2/MPI/

Protein=PROTEINNAME
Ligand=LIGANDNAME
LigandCharge=CHARGEVALUE
Iteration=${SLURM_ARRAY_TASK_ID}
Extension=_mpi
File=${Ligand}
Cores=20

cd ${File}_${Iteration}
#Start by cleaning everything.
rm -rf MMPBSA Frames* Snapshot*

echo "#####################"
echo " "
echo "Get Frames for MMPBSA"
echo " "
echo "#####################"
#We extract 100 snapshots every 10ps.
TimeBegin=4010
TimeEnd=5000
gmx${Extension} trjconv -s 08-FullMD.tpr -n 04-index.ndx -f 08-FullMD.xtc -o Snapshot.pdb -dt 10 -sep -b ${TimeBegin} -e ${TimeEnd} << EOI
Protein_LIG
EOI
Time=${TimeBegin}
x=0
FileNumber=` ls Snapshot*.pdb | wc -l `
while [ $x -lt ${FileNumber} ]; do
        mv Snapshot$x.pdb Complex_${Time}ps.pdb
        Time=`expr ${Time} + 10`
        x=`expr $x + 1`
done
mkdir Frames_${File}_${Iteration}
mv Complex_*.pdb Frames_${File}_${Iteration}/.

module load centos6/amber12_openmpi-1.6.5_intel-13.0.079
AMBERHOME=/n/sw/centos6/amber12_openmpi-1.6.5_intel-13.0.079

echo "######"
echo " "
echo "MMPBSA"
echo " "
echo "######"
mkdir MMPBSA
cd MMPBSA
mkdir Mdcrd_${File}_${Iteration}
cp ../../00-Files/Structures/${Ligand}.mol2 LIG.mol2
cp ../../00-Files/Structures/${Ligand}_AC.frcmod LIG_AC.frcmod
sed -i "s/Cl1/CL1/g" LIG.mol2
sed -i "s/Br1/BR1/g" LIG.mol2
Time=${TimeBegin}
while [ ${Time} -le ${TimeEnd} ]; do
	cp ../Frames_${File}_${Iteration}/Complex_${Time}ps.pdb .
	#First all HIS are converted to HIE (since it is the protonation state most often found in proteins).
	sed -i "s/HIS/HIE/g" Complex_${Time}ps.pdb
	#If some HIS were not detected as HIE by Gromacs, it should be changed. Near the beginning of the .e outpout file, Gromacs will print information such as:
		#Will use HISE for residue 3
		#Will use HISE for residue 4
		#Will use HISE for residue 10
		#Will use HISE for residue 15
		#Will use HISD for residue 94
		#Will use HISD for residue 96
		#Will use HISH for residue 107
		#Will use HISE for residue 119
		#Will use HISD for residue 122
	#In such a case, you can use something like the following (please note that the above information are not the ones for 1HPV since this protein has only
	#two histidines and this is not enough for the example):
	#if [ ${Protein} == 1HPV ]; then
	#	sed -i -e "s/HIE A  94/HID A  94/g" -e "s/HIE A  96/HID A  96/g" -e "s/HIE A 107/HIP A 107/g" -e "s/HIE A 122/HID A 122/g" Complex_${Time}ps.pdb
	#fi

	#If some cysteines are involved in disulfide bonds, you can use:
	#sed -i "s/CYS A  60/CYX A  60/g" Complex_${Time}ps.pdb

	#This protein is made with two chains. The name of the oxygen of the terminal residue of each chain (PHE in both cases) must be changed.
	sed -i -e "s/OC1 PHE A  99/OXT PHE A  99/g" -e "s/OC2 PHE A  99/O   PHE A  99/g" -e "s/OC1 PHE B  99/OXT PHE B  99/g" -e "s/OC2 PHE B  99/O   PHE B  99/g" Complex_${Time}ps.pdb
	#Finish to prepare the complex file.
	sed -i -f ../../00-Files/Gromacs2Amber.sed Complex_${Time}ps.pdb
	sed -i -e '/REMARK/d' -e '/CRYST1/d' -e '/MODEL/d' -e '/ENDMDL/d' -e '/TITLE/d' Complex_${Time}ps.pdb

	#Prepare files for the ligand, the receptor, the complex.
	grep "LIG" Complex_${Time}ps.pdb > ${Ligand}.pdb
	sed -i -e "s/ATOM  /HETATM/g" -e "s/LIG B/LIG A/g" -e "s/LIG C/LIG A/g" ${Ligand}.pdb
	echo "TER" >> ${Ligand}.pdb
	sed -i "/LIG/d" Complex_${Time}ps.pdb
	cp Complex_${Time}ps.pdb Receptor.pdb
	cat ${Ligand}.pdb >> Complex_${Time}ps.pdb
	mv Complex_${Time}ps.pdb Complex.pdb
	mv ${Ligand}.pdb Ligand.pdb

	#Create coordinate files
	sed "s/NAMELIGAND/LIG/g" ../../00-Files/TEMPLATECOMPLEX_tleap.in > tleap.in
	tleap -f tleap.in

	#Clean
	mv LIG_Complex.inpcrd Mdcrd_${File}_${Iteration}/${Protein}_Complex_${Time}ps.mdcrd
	rm LIG_Receptor.inpcrd LIG_Ligand.inpcrd Complex.pdb Receptor.pdb Ligand.pdb tleap.in leap.log

	echo "${Time} Done"
	Time=`expr ${Time} + 10`
done

#Run MMPBSA
mpirun -np ${Cores} MMPBSA.py.MPI -O -i ../../00-Files/MMPBSA_Binding.in -o MMPBSA_${File}_${Iteration}.dat -cp Complex.prmtop -rp Receptor.prmtop -lp Ligand.prmtop -y Mdcrd_${File}_${Iteration}/*.mdcrd
sleep 10s

#Clean
mkdir Files_${File}_${Iteration}
mv _MMPBSA* Files_${File}_${Iteration}/.
#
tar czfv Files_${File}_${Iteration}.tar.gz Files_${File}_${Iteration}/
rm -rf Files_${File}_${Iteration}/ 
#
tar czfv Mdcrd_${File}_${Iteration}.tar.gz Mdcrd_${File}_${Iteration}/
rm -rf Mdcrd_${File}_${Iteration}/
#
cd ..
tar czfv Frames_${File}_${Iteration}.tar.gz Frames_${File}_${Iteration}/
rm -rf Frames_${File}_${Iteration}/

