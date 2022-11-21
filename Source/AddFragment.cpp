#include "OpenGrowth.h"

double dist(double const atom1[3], double const atom2[3]);
double Energy(Parameters const & parameters, Protein const & protein, Molecule & molecule);
float  distNB(char const type1[], char const type2[]);
int    StericClash(Parameters const & parameters, Protein const & protein, Molecule const & molecule);
void   OptimizationPosition(Parameters & parameters, Protein const & protein, Molecule & molecule);
float  OptimizationGeom(Parameters const & parameters, Protein const & protein, Molecule & molecule, int & errorFF, int const & optParameter);
int    Random(Parameters & parameters, double const probaFile[]);

// This function adds a new fragment to an existing molecule.
void  AddFragment(Parameters & parameters, Protein const protein[], Molecule molecule[], double & averageEnergy, string const typeOfSnapshots, Molecule const fragment[], double const probaFirstFrag[], double const probaTransition[][MAX_FRAGMENTS], int ligandDescription[], int & errorFF)
{
    if (parameters.verbose>=3 ) { cout << "\t\tStart adding a new fragment." << endl; }
    Molecule newFragment;
    unsigned int indexH_L=0, indexH_F=0, indexFragment=0;

    // Depending on which kind of snapshots we are using, there are not the same number of them.
    int snapshotNumber=0;
    if      ( typeOfSnapshots=="rotamers" )   { snapshotNumber = parameters.rotamersNumber;   }
    else if ( typeOfSnapshots=="conformers" ) { snapshotNumber = parameters.conformersNumber; }

    // Check if we can grow linearly or branchedly.
    int canGrowLinearly=0;
    int canGrowBranchedly=0;
    for (unsigned int j=0; j<molecule[0].growingSites.size() ; j++ ) {
        unsigned int indexAtom = molecule[0].growingSites[j];
        if (molecule[0].fragmentNeighbours[indexAtom-1]<=1) { canGrowLinearly++;   }
        if (molecule[0].fragmentNeighbours[indexAtom-1]>1)  { canGrowBranchedly++; }
        }

    // Pick a random hydrogen in the current ligand, or pick it from the regrow file. indexH_L will store the atom number in the molecule of this hydrogen.
    // If we are growing in conformers and we did a rotamer search before (or we are using REGROW), the hydrogen has already been chosen and we pick it from an array.
    if      ((parameters.rotamersNumber!=0 ) && (typeOfSnapshots=="conformers"))                                 { indexH_L = ligandDescription[3*(molecule[0].numberResidues-1)+1];  }
    else if (parameters.growthMode=="REGROW")                                                                    { indexH_L = ligandDescription[3*(molecule[0].numberResidues-1)+1];  }
    else if (parameters.growthMode=="RANDOM" || parameters.growthMode=="BIASED" || parameters.growthMode=="FOG") {
        // Decide which kind of growth we will perform and from which hydrogen.
        int badFragment=0;
        do {
            double randomNumber = (parameters.random() % 1000) / 1000.0 ;
            // Select indexH_L.
            if      ((randomNumber>=parameters.branchingProba) && (canGrowLinearly>=1))      {        // Grow linear
                do {
                    unsigned int indexSite = parameters.random() % molecule[0].growingSites.size();
                    indexH_L = molecule[0].growingSites[indexSite];
                    } while (molecule[0].fragmentNeighbours[indexH_L-1]>1);
                }
            else if (randomNumber>=parameters.branchingProba)                                {        // Can't grow linear, so we grow branched
                unsigned int indexSite = parameters.random() % molecule[0].growingSites.size();
                indexH_L = molecule[0].growingSites[indexSite];
                }
            else if ((randomNumber<parameters.branchingProba) && (canGrowBranchedly>=1))     {        // Grow branched
                do {
                    unsigned int indexSite = parameters.random() % molecule[0].growingSites.size();
                    indexH_L = molecule[0].growingSites[indexSite];
                    } while (molecule[0].fragmentNeighbours[indexH_L-1]<=1);
                }
            else if (randomNumber<parameters.branchingProba)                                 {        // Can't grow branched so grow linear
                unsigned int indexSite = parameters.random() % molecule[0].growingSites.size();
                indexH_L = molecule[0].growingSites[indexSite];
                }

            // Check that the fragment is allowed. This means that something can grow from here.
            badFragment=0;
            for (unsigned int j=0; j<parameters.forbiddenFragments.size() ; j++ ) {
                if (parameters.forbiddenFragments[j]==molecule[0].growth[indexH_L-1]) { badFragment++; }
                }
            } while (badFragment!=0);
        }
    
    // By reading in the growth array, we can find to which fragment this hydrogen belongs.
    int fragmentType=molecule[0].growth[indexH_L-1];
    if (parameters.verbose>=3) { cout << "\t\tindexH_L (" << indexH_L << ") was successfully chosen. It is from the fragment type " << fragmentType << "." << endl; }
    // Choose the new fragment, either from the previous growth, randomly, according to the FOG probabilities or according to a given REGROW_FILE.
    if      ((parameters.rotamersNumber!=0 ) && (typeOfSnapshots=="conformers"))  { indexFragment = ligandDescription[3*(molecule[0].numberResidues-1)+2]; }
    else if (parameters.growthMode=="RANDOM")                     { indexFragment = (parameters.random() % parameters.fragmentListSize)+1; }
    else if (parameters.growthMode=="BIASED")                     { indexFragment = Random(parameters, probaFirstFrag)+1;                  }
    else if (parameters.growthMode=="FOG")                        { indexFragment = Random(parameters, probaTransition[fragmentType-1])+1; }
    else if (parameters.growthMode=="REGROW")                     { indexFragment = ligandDescription[3*(molecule[0].numberResidues-1)+2]; }
    // Store the new fragment.
    newFragment=fragment[indexFragment-1];

    // Pick an hydrogen in the new fragment which describes the fragment properly i.e. which has a growth mode of the same value as the indexFragment, or pick it from the ligand description array.
    if      ((parameters.rotamersNumber!=0 ) && (typeOfSnapshots=="conformers"))                                 { indexH_F = ligandDescription[3*(molecule[0].numberResidues-1)+3]; }
    else if (parameters.growthMode=="REGROW")                                                                    { indexH_F = ligandDescription[3*(molecule[0].numberResidues-1)+3]; }
    else if (parameters.growthMode=="RANDOM" || parameters.growthMode=="BIASED" || parameters.growthMode=="FOG") {
        do {
            unsigned int indexSite = parameters.random() % newFragment.growingSites.size();
            indexH_F = newFragment.growingSites[indexSite];
            } while(newFragment.growth[indexH_F-1]!=indexFragment);
        if (parameters.verbose>=3) { cout << "\t\tindexH_F (" << indexH_F << ") was successfully chosen. It is from the fragment type " << indexFragment << " (" << fragment[indexFragment-1].obMol.GetTitle() << ")." << endl; }
        }

    // Get the atom from the selected hydrogen of the ligand and its neighbour.
    OpenBabel::OBAtom *atomH_L;
    OpenBabel::OBAtom *atomN_L;
    atomH_L = molecule[0].obMol.GetAtom(indexH_L);
    OpenBabel::OBAtomAtomIter nbrL(atomH_L);
    int indexN_L = nbrL->GetIdx();
    atomN_L = molecule[0].obMol.GetAtom(indexN_L);

    // Get the atom from the selected hydrogen of the fragment and its neighbour.
    OpenBabel::OBAtom *atomH_F;
    OpenBabel::OBAtom *atomN_F;
    atomH_F = newFragment.obMol.GetAtom(indexH_F);
    OpenBabel::OBAtomAtomIter nbrF(atomH_F);
    int indexN_F = nbrF->GetIdx();
    atomN_F = newFragment.obMol.GetAtom(indexN_F);

    // Store the coordinates from the new fragment.
    double fragment0[newFragment.obMol.NumAtoms()][3];
    for (unsigned int j=1; j<=newFragment.obMol.NumAtoms() ; j++ ) {
        OpenBabel::OBAtom *atom;
        atom = newFragment.obMol.GetAtom(j);
        fragment0[j-1][0]=atom->GetX();
        fragment0[j-1][1]=atom->GetY();
        fragment0[j-1][2]=atom->GetZ();
        }

    // Translation 1 - The selected Neighbour atom of the new Fragment is put in (0,0,0).
    double fragment1[newFragment.obMol.NumAtoms()][3];
    for (unsigned int j=1; j<=newFragment.obMol.NumAtoms() ; j++ ) {
        fragment1[j-1][0]=fragment0[j-1][0]-fragment0[indexN_F-1][0];
        fragment1[j-1][1]=fragment0[j-1][1]-fragment0[indexN_F-1][1];
        fragment1[j-1][2]=fragment0[j-1][2]-fragment0[indexN_F-1][2];
        }

    // We need to store all the new molecules to decide if we accept or not the new fragment.
    Molecule newMolecule[snapshotNumber];
    // Add the new fragment in every snapshot.
    for (int s=0 ; s < snapshotNumber ; s++) {
        // Store the coordinates from the ligand.
        double ligand[molecule[s].obMol.NumAtoms()][3];
        for (unsigned int j=1; j<=molecule[s].obMol.NumAtoms() ; j++ ) {
            OpenBabel::OBAtom *atom;
            atom = molecule[s].obMol.GetAtom(j);
            ligand[j-1][0]=atom->GetX();
            ligand[j-1][1]=atom->GetY();
            ligand[j-1][2]=atom->GetZ();
            }

        // Rotation 1.
        // We will do a first rotation to align the H_F/N_F axe with the N_L/H_L axe. Let's define two vectors U (fragment1[indexN_F-1]-fragment1[indexH_F-1]) and
        // V (ligand[indexH_L-1]-ligand[indexN_L-1]) and normalize them. After the rotation, U will be aligned with V. We define W=UxV/||UxV||.
        double uX = (fragment1[indexN_F-1][0]-fragment1[indexH_F-1][0])/dist(fragment1[indexN_F-1], fragment1[indexH_F-1]);
        double uY = (fragment1[indexN_F-1][1]-fragment1[indexH_F-1][1])/dist(fragment1[indexN_F-1], fragment1[indexH_F-1]);
        double uZ = (fragment1[indexN_F-1][2]-fragment1[indexH_F-1][2])/dist(fragment1[indexN_F-1], fragment1[indexH_F-1]);
        double vX = (ligand[indexH_L-1][0]-ligand[indexN_L-1][0])/dist(ligand[indexH_L-1], ligand[indexN_L-1]);
        double vY = (ligand[indexH_L-1][1]-ligand[indexN_L-1][1])/dist(ligand[indexH_L-1], ligand[indexN_L-1]);
        double vZ = (ligand[indexH_L-1][2]-ligand[indexN_L-1][2])/dist(ligand[indexH_L-1], ligand[indexN_L-1]);
        double wX = (uY*vZ - uZ*vY)/sqrt(pow(uY*vZ - uZ*vY,2) + pow(uZ*vX - uX*vZ,2) + pow(uX*vY - uY*vX,2));
        double wY = (uZ*vX - uX*vZ)/sqrt(pow(uY*vZ - uZ*vY,2) + pow(uZ*vX - uX*vZ,2) + pow(uX*vY - uY*vX,2));
        double wZ = (uX*vY - uY*vX)/sqrt(pow(uY*vZ - uZ*vY,2) + pow(uZ*vX - uX*vZ,2) + pow(uX*vY - uY*vX,2));
        // We need to know from which angle we want to rotate. We know with the dot and the cross products: U.V = U*V*cos(alpha) and UxV = U*V*sin(alpha). Since U and V are normalized to 1, the norm of
        // UxV gives us |sin(alpha)|. This is not enough, we need to know the sign of alpha. Thus, we use the following: sin(alpha) = det(W, U, V) where W=UxV/||UxV|| (the rule of Sarrus can be used).
        double cosAlpha = uX*vX + uY*vY + uZ*vZ;
        double sinAlpha = (wX*uY*vZ + uX*vY*wZ + vX*wY*uZ) - (vX*uY*wZ + wX*vY*uZ + uX*wY*vZ);
        // Let's write the rotation matrix around the vector W (http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle).
        double rotationMatrix1[3][3];
        rotationMatrix1[0][0] = cosAlpha + pow(wX,2)*(1-cosAlpha);
        rotationMatrix1[0][1] = wX*wY*(1-cosAlpha) + wZ*sinAlpha;
        rotationMatrix1[0][2] = wX*wZ*(1-cosAlpha) - wY*sinAlpha;     // rotationMatrix1 has the following form:
        rotationMatrix1[1][0] = wX*wY*(1-cosAlpha) - wZ*sinAlpha;     // [0][0] [1][0] [2][0]
        rotationMatrix1[1][1] = cosAlpha + pow(wY,2)*(1-cosAlpha);    // [0][1] [1][1] [2][1]
        rotationMatrix1[1][2] = wY*wZ*(1-cosAlpha) + wX*sinAlpha;     // [0][2] [1][2] [2][2]
        rotationMatrix1[2][0] = wX*wZ*(1-cosAlpha) + wY*sinAlpha;
        rotationMatrix1[2][1] = wY*wZ*(1-cosAlpha) - wX*sinAlpha;
        rotationMatrix1[2][2] = cosAlpha + pow(wZ,2)*(1-cosAlpha);

        // The alignement is made.
        double fragment2[newFragment.obMol.NumAtoms()][3];
        for (unsigned int j=1; j<=newFragment.obMol.NumAtoms() ; j++ ) {
            fragment2[j-1][0] = rotationMatrix1[0][0]*fragment1[j-1][0] + rotationMatrix1[1][0]*fragment1[j-1][1] + rotationMatrix1[2][0]*fragment1[j-1][2] ;
            fragment2[j-1][1] = rotationMatrix1[0][1]*fragment1[j-1][0] + rotationMatrix1[1][1]*fragment1[j-1][1] + rotationMatrix1[2][1]*fragment1[j-1][2] ;
            fragment2[j-1][2] = rotationMatrix1[0][2]*fragment1[j-1][0] + rotationMatrix1[1][2]*fragment1[j-1][1] + rotationMatrix1[2][2]*fragment1[j-1][2] ;
            }

        // Rotation 2 + Translation 2.
        // We will rotate around the axe fragment2[indexH_F-1]-fragment2[indexN_F-1] (vector R) every 2*PI/rotationPrecision.
        float rX = (fragment2[indexN_F-1][0]-fragment2[indexH_F-1][0])/dist(fragment2[indexN_F-1], fragment2[indexH_F-1]);
        float rY = (fragment2[indexN_F-1][1]-fragment2[indexH_F-1][1])/dist(fragment2[indexN_F-1], fragment2[indexH_F-1]);
        float rZ = (fragment2[indexN_F-1][2]-fragment2[indexH_F-1][2])/dist(fragment2[indexN_F-1], fragment2[indexH_F-1]);
        float Beta=2*PI/(parameters.rotationPrecision);
        float rotationMatrix2[3][3];
        // The selected atom of the fragment (indexN_F-1) is then put at the proper distance of the selected atom of the current ligand (indexN_L-1). The binding atom is aligned with the C-H
        // of the current ligand, but a little bit further (because the bound is longer than the C-H one). For example, for x : Binding_x = Neighbour_x + d(Neighbour-Binding)*(Hx-Nx)/d(H-N).
        // The distance Neighbour-Binding depends on the two atoms and is given by the function distNB which reads this information from databases.
        float distFragLig = 1.00*distNB(atomN_L->GetType(), atomN_F->GetType());
        float distFragLig_x = ligand[indexN_L-1][0] + distFragLig*(ligand[indexH_L-1][0]-ligand[indexN_L-1][0])/dist(ligand[indexH_L-1], ligand[indexN_L-1]) ;
        float distFragLig_y = ligand[indexN_L-1][1] + distFragLig*(ligand[indexH_L-1][1]-ligand[indexN_L-1][1])/dist(ligand[indexH_L-1], ligand[indexN_L-1]) ;
        float distFragLig_z = ligand[indexN_L-1][2] + distFragLig*(ligand[indexH_L-1][2]-ligand[indexN_L-1][2])/dist(ligand[indexH_L-1], ligand[indexN_L-1]) ;

        // Create new arrays for growth, fragmentIndex and fragmentNeighbours.
        newMolecule[s].growth.clear();
        newMolecule[s].fragmentIndex.clear();
        newMolecule[s].fragmentNeighbours.clear();
        for (unsigned int j=1; j<=molecule[s].obMol.NumAtoms() ; j++ ) {
            newMolecule[s].growth.push_back(molecule[s].growth[j-1]);
            newMolecule[s].fragmentIndex.push_back(molecule[s].fragmentIndex[j-1]);
            newMolecule[s].fragmentNeighbours.push_back(molecule[s].fragmentNeighbours[j-1]);
            }
        // The ligand neighbour atom get a new neighbour fragment. Do it only if we have at least 2 fragments.
        if (molecule[0].numberResidues>=2) {
            int startingFragment=molecule[0].fragmentIndex[indexN_L-1];
            for (unsigned int j=1; j<=molecule[s].obMol.NumAtoms() ; j++ ) {
                if(newMolecule[s].fragmentIndex[j-1]==startingFragment) { newMolecule[s].fragmentNeighbours[j-1]++; }
                }
            }
        // Add info from newFragment.
        for (unsigned int j=1; j<=newFragment.obMol.NumAtoms() ; j++ ) {
            newMolecule[s].growth.push_back(newFragment.growth[j-1]);
            newMolecule[s].fragmentIndex.push_back(molecule[s].numberResidues+1);
            newMolecule[s].fragmentNeighbours.push_back(1);
            }
        newMolecule[s].growth.erase(newMolecule[s].growth.begin()+indexH_L-1);
        newMolecule[s].growth.erase(newMolecule[s].growth.begin()+molecule[s].obMol.NumAtoms()-1+indexH_F-1);
        newMolecule[s].fragmentIndex.erase(newMolecule[s].fragmentIndex.begin()+indexH_L-1);
        newMolecule[s].fragmentIndex.erase(newMolecule[s].fragmentIndex.begin()+molecule[s].obMol.NumAtoms()-1+indexH_F-1);
        newMolecule[s].fragmentNeighbours.erase(newMolecule[s].fragmentNeighbours.begin()+indexH_L-1);
        newMolecule[s].fragmentNeighbours.erase(newMolecule[s].fragmentNeighbours.begin()+molecule[s].obMol.NumAtoms()-1+indexH_F-1);

        // Create new array for growingSites and update numberResidues.
        newMolecule[s].growingSites.clear();
        for (unsigned int j=0; j<molecule[s].growingSites.size() ; j++ ) {
            if       (molecule[s].growingSites[j]<indexH_L) { newMolecule[s].growingSites.push_back(molecule[s].growingSites[j]);   }
            else if  (molecule[s].growingSites[j]>indexH_L) { newMolecule[s].growingSites.push_back(molecule[s].growingSites[j]-1); }
            }
        for (unsigned int j=0; j<newFragment.growingSites.size(); j++ ) {
            if       (newFragment.growingSites[j]<indexH_F) { newMolecule[s].growingSites.push_back(molecule[s].obMol.NumAtoms()-1+newFragment.growingSites[j]);   }
            else if  (newFragment.growingSites[j]>indexH_F) { newMolecule[s].growingSites.push_back(molecule[s].obMol.NumAtoms()-1+newFragment.growingSites[j]-1); }
            }
        newMolecule[s].numberResidues = molecule[s].numberResidues+1;

        // We look for the best rotamer (with the lowest energy).
        for (int k=0 ; k<parameters.rotationPrecision ; k++) {
            if (parameters.verbose>=6 ) { cout << "\t\t\t--- Rotamer search: " << k << " ---" << endl; }
            // We will need a temporary molecule.
            Molecule tempMolecule;
            tempMolecule = newMolecule[s];
            tempMolecule.obMol.Clear();

            // Rotation 2 + Translation 2.
            rotationMatrix2[0][0] = cos(k*Beta) + pow(rX,2)*(1-cos(k*Beta));
            rotationMatrix2[0][1] = rX*rY*(1-cos(k*Beta)) + rZ*sin(k*Beta) ;
            rotationMatrix2[0][2] = rX*rZ*(1-cos(k*Beta)) - rY*sin(k*Beta) ;    // rotationMatrix2 has the following form:
            rotationMatrix2[1][0] = rX*rY*(1-cos(k*Beta)) - rZ*sin(k*Beta) ;    // [0][0] [1][0] [2][0]
            rotationMatrix2[1][1] = cos(k*Beta) + pow(rY,2)*(1-cos(k*Beta));    // [0][1] [1][1] [2][1]
            rotationMatrix2[1][2] = rY*rZ*(1-cos(k*Beta)) + rX*sin(k*Beta) ;    // [0][2] [1][2] [2][2]
            rotationMatrix2[2][0] = rX*rZ*(1-cos(k*Beta)) + rY*sin(k*Beta) ;
            rotationMatrix2[2][1] = rY*rZ*(1-cos(k*Beta)) - rX*sin(k*Beta) ;
            rotationMatrix2[2][2] = cos(k*Beta) + pow(rZ,2)*(1-cos(k*Beta));
            double fragment3[newFragment.obMol.NumAtoms()][3];
            for (unsigned int j=1; j<=newFragment.obMol.NumAtoms() ; j++ ) {
                fragment3[j-1][0] = rotationMatrix2[0][0]*fragment2[j-1][0] + rotationMatrix2[1][0]*fragment2[j-1][1] + rotationMatrix2[2][0]*fragment2[j-1][2] + distFragLig_x;
                fragment3[j-1][1] = rotationMatrix2[0][1]*fragment2[j-1][0] + rotationMatrix2[1][1]*fragment2[j-1][1] + rotationMatrix2[2][1]*fragment2[j-1][2] + distFragLig_y;
                fragment3[j-1][2] = rotationMatrix2[0][2]*fragment2[j-1][0] + rotationMatrix2[1][2]*fragment2[j-1][1] + rotationMatrix2[2][2]*fragment2[j-1][2] + distFragLig_z;
                }

            // Update newFragment.obMol
            for (unsigned int j=1; j<=newFragment.obMol.NumAtoms() ; j++ ) {
                OpenBabel::OBAtom *atom;
                atom = newFragment.obMol.GetAtom(j);
                atom->SetVector(fragment3[j-1][0], fragment3[j-1][1], fragment3[j-1][2]);
                }

            // Join the ligand and the fragment.
            tempMolecule.obMol = molecule[s].obMol;
            tempMolecule.obMol += newFragment.obMol;
            // Delete the atom from the ligand.
            OpenBabel::OBAtom *atomToDeleteLigand;
            atomToDeleteLigand = tempMolecule.obMol.GetAtom(indexH_L);
            tempMolecule.obMol.DeleteAtom(atomToDeleteLigand);
            // Delete the atom from the fragment.
            OpenBabel::OBAtom *atomToDeleteFragment;
            atomToDeleteFragment = tempMolecule.obMol.GetAtom(molecule[s].obMol.NumAtoms()-1+indexH_F);
            tempMolecule.obMol.DeleteAtom(atomToDeleteFragment);
            // Create a new bond.
            tempMolecule.obMol.AddBond(indexN_L, molecule[s].obMol.NumAtoms()-1+indexN_F, 1);

            // Optimize the molecule position for every rotamer.
            if ((parameters.optimizationMode/10)>=3) { OptimizationGeom(parameters, protein[s], tempMolecule, errorFF, 1); }
            if ((parameters.optimizationMode%10)>=3) { OptimizationPosition(parameters, protein[s], tempMolecule);         }

            // Calculate the energy and check for clashes.
            tempMolecule.stericClashes = StericClash(parameters, protein[s], tempMolecule);
            tempMolecule.interactionEnergy = Energy(parameters, protein[s], tempMolecule);

            // Accept a new orientation of the new fragment if there were steric clashes and there are not in the new orientation, or if the new energy is lower than the previous one (without clashes).
            if (k == 0)     { newMolecule[s] = tempMolecule; }
            else {
                if      ( !tempMolecule.stericClashes && newMolecule[s].stericClashes )                                   { newMolecule[s] = tempMolecule; }
                else if ( !tempMolecule.stericClashes && (tempMolecule.interactionEnergy<molecule[s].interactionEnergy) ) { newMolecule[s] = tempMolecule; }
                }
            }

        // Optimize the molecule position by allowing some translations or rotations in the active site, once the best rotamer is found.
        // It is done only if no intramolecular clashes are found (i.e. stericClashes==0 or 1).   //&& newMolecule[s].stericClashes<=1
        if ((parameters.optimizationMode/10)>=2) { OptimizationGeom(parameters, protein[s],  newMolecule[s], errorFF, 1); }
        if ((parameters.optimizationMode%10)>=2) { OptimizationPosition(parameters, protein[s],  newMolecule[s]);         }
        }

    // Check for steric clashes in all the molecules.
    int stericClashes=0, clashIntra=0, clashInter=0;
    for (int s=0 ; s < snapshotNumber ; s++) {
        if      ( newMolecule[s].stericClashes==1 ) { clashInter++;               }
        else if ( newMolecule[s].stericClashes==2 ) { clashIntra++;               }
        else if ( newMolecule[s].stericClashes==3 ) { clashInter++; clashIntra++; }
        }
    if      (clashIntra==0 && clashInter==0) { stericClashes=0; }
    else if (clashIntra==0 && clashInter>0)  { stericClashes=1; }
    else if (clashIntra>0  && clashInter==0) { stericClashes=2; }
    else if (clashIntra>0  && clashInter>0)  { stericClashes=3; }

    if (errorFF==0 && stericClashes==0) {
        long double newEnergy=0;
        if      (snapshotNumber==1)                   { newEnergy=newMolecule[0].interactionEnergy; }
        else if (parameters.averageType=="BOLTZMANN") {
            // Compute the Boltzmann average energy of the new molecules.
            long double partitionFunction=0;
            for (int s=0 ; s < snapshotNumber ; s++) {
                partitionFunction += exp(-BETA*newMolecule[s].interactionEnergy);
                newEnergy += newMolecule[s].interactionEnergy*exp(-BETA*newMolecule[s].interactionEnergy);
                }
            newEnergy = newEnergy/partitionFunction;
            }
        else if (parameters.averageType=="ARITHMETIC") {
            // Compute the arithmetic average energy of the new molecules.
            for (int s=0 ; s < snapshotNumber ; s++) {
                newEnergy += newMolecule[s].interactionEnergy;
                }
            newEnergy = newEnergy/snapshotNumber;
            }
        else if (parameters.averageType=="LOWESTSCORE") {
            // Keep only the lowest energy.
            newEnergy = molecule[0].interactionEnergy; 
            for (int s=1 ; s < snapshotNumber ; s++) {
                if(molecule[s].interactionEnergy<newEnergy) { newEnergy = molecule[s].interactionEnergy; }
                }
            }

        // Pick a random number for the growth algorithm.
        long double randomNumber = (parameters.random() % 100000000) / 100000000.0 ;
        double ratioMetropolis = exp( (averageEnergy-newEnergy)/parameters.MCTemp );

        // Accept the new fragment according to a Metropolis algorithm.
        if ((newEnergy<=averageEnergy) || (randomNumber<ratioMetropolis)) {
            // Save the information to describe the molecule afterwards.
            if (parameters.growthMode=="RANDOM" || parameters.growthMode=="BIASED" || parameters.growthMode=="FOG")  {
                ligandDescription[3*(molecule[0].numberResidues-1)+1]=indexH_L;
                ligandDescription[3*(molecule[0].numberResidues-1)+2]=indexFragment;
                ligandDescription[3*(molecule[0].numberResidues-1)+3]=indexH_F;
                }

            // If we are adding the second fragment, change the neighbours count for the first fragment.
            if (molecule[0].numberResidues==1) {
                for (int s=0 ; s < snapshotNumber ; s++) {
                    for (unsigned int j=1; j<=molecule[s].obMol.NumAtoms()-1 ; j++ ) {
                        newMolecule[s].fragmentNeighbours[j-1]++;
                        }
                    }
                }

            // The newMolecule becomes the current ligand.
            for (int s=0 ; s < snapshotNumber ; s++) {
                molecule[s] = newMolecule[s];
                }

            // Display results if we are dealing with rotamers and there will be no conformers search, or if we are dealing with conformers, or if verbose>=2.
            if ( (typeOfSnapshots=="rotamers" && parameters.conformersNumber==0) || (typeOfSnapshots=="conformers") ) {
                if (parameters.verbose>=1 ) { cout << "Fragment " << molecule[0].numberResidues << ": " << newFragment.obMol.GetTitle() << endl; }
                if (parameters.verbose>=3 ) {
                    cout << "\t\tAtom " << indexH_F << " of the new fragment will be used. It is bound to the atom " << indexN_F << "." << endl;
                    cout << "\t\tAtom " << indexH_L << " of the current ligand will be used. It is bound to the atom " << indexN_L << "." << endl;
                    }
                if (parameters.verbose>=2 ) {
                    if ( newEnergy>averageEnergy ) { cout << "\tNew Energy: " << newEnergy << " (Old energy: " << averageEnergy << ", ratioMetropolis: " << ratioMetropolis << ", random number: " << randomNumber << ")" << endl; }
                    else                           { cout << "\tNew Energy: " << newEnergy << endl;                                                                                                                                }
                    }
                }
            averageEnergy = newEnergy;
            }
        else {
            if (parameters.verbose>=3) { cout << "\t\tFragment " << newFragment.obMol.GetTitle() << " (atom " << indexH_F << ") on the atom " << indexH_L << " of the current ligand is rejected by Monte-Carlo." << endl; }
            }
        }
    else {
        if (parameters.verbose>=3) {
            if      (stericClashes==1) { cout << "\t\tAtom " << indexH_F << " of " << newFragment.obMol.GetTitle() << " on the atom " << indexH_L << " of the current ligand is rejected because of intermolecular clashes."           << endl; }
            else if (stericClashes==2) { cout << "\t\tAtom " << indexH_F << " of " << newFragment.obMol.GetTitle() << " on the atom " << indexH_L << " of the current ligand is rejected because of intramolecular clashes."           << endl; }
            else if (stericClashes==3) { cout << "\t\tAtom " << indexH_F << " of " << newFragment.obMol.GetTitle() << " on the atom " << indexH_L << " of the current ligand is rejected because of inter and intramolecular clashes." << endl; }
            }
        }
}

// This function returns the typical distance between two types of atoms.
// Check in http://publcif.iucr.org/cifmoldb/openbabel/share/openbabel/2.3.0/types.txt
// http://en.wikipedia.org/wiki/Chemical_bond and http://fr.wikipedia.org/wiki/Liaison_chimique
// Currently: Br C1 C2 C3 Cac Car Cl F HC HO I N1 N3 N3+ Nam Nar Nox Ntr O2 O3 O.co2 S2 S3 So2
float distNB(char const type1[], char const type2[])
{
    float Distance;

    // For first atom.
    char typeFirstAtom[6];
    if      (!strcmp(type1, "C1"))                                                     { strcpy(typeFirstAtom, "Csp");   }
    else if (!strcmp(type1, "C2") || !strcmp(type1, "Cac") || !strcmp(type1, "Car"))   { strcpy(typeFirstAtom, "Csp2");  }
    else if (!strcmp(type1, "C3"))                                                     { strcpy(typeFirstAtom, "Csp3");  }
    else if (!strcmp(type1, "N1") || !strcmp(type1, "N3")  || !strcmp(type1, "N3+") || !strcmp(type1, "Nam") || 
        !strcmp(type1, "Nar") || !strcmp(type1, "Nox") || !strcmp(type1, "Ntr"))       { strcpy(typeFirstAtom, "N");     }
    else if (!strcmp(type1 ,"O2") || !strcmp(type1, "O3")  || !strcmp(type1, "O.co2")) { strcpy(typeFirstAtom, "O");     }
    else if (!strcmp(type1, "S2") || !strcmp(type1, "S3")  || !strcmp(type1, "So2"))   { strcpy(typeFirstAtom, "S");     }
    else if (!strcmp(type1, "F"))                                                      { strcpy(typeFirstAtom, "F");     }
    else if (!strcmp(type1, "Cl"))                                                     { strcpy(typeFirstAtom, "Cl");    }
    else if (!strcmp(type1, "Br"))                                                     { strcpy(typeFirstAtom, "Br");    }
    else if (!strcmp(type1, "I"))                                                      { strcpy(typeFirstAtom, "I");     }

    // For second atom.
    char typeSecondAtom[6];
    if      (!strcmp(type2, "C1"))                                                     { strcpy(typeSecondAtom, "Csp");  }
    else if (!strcmp(type2, "C2") || !strcmp(type2, "Cac") || !strcmp(type2, "Car"))   { strcpy(typeSecondAtom, "Csp2"); }
    else if (!strcmp(type2, "C3"))                                                     { strcpy(typeSecondAtom, "Csp3"); }
    else if (!strcmp(type2, "N1") || !strcmp(type2, "N3")  || !strcmp(type2, "N3+") || !strcmp(type2,"Nam") || 
        !strcmp(type2, "Nar") || !strcmp(type2, "Nox") || !strcmp(type2, "Ntr"))       { strcpy(typeSecondAtom, "N");    }
    else if (!strcmp(type2, "O2") || !strcmp(type2, "O3")  || !strcmp(type2, "O.co2")) { strcpy(typeSecondAtom, "O");    }
    else if (!strcmp(type2, "S2") || !strcmp(type2, "S3")  || !strcmp(type2, "So2"))   { strcpy(typeSecondAtom, "S");    }
    else if (!strcmp(type2, "F"))                                                      { strcpy(typeSecondAtom, "F");    }
    else if (!strcmp(type2, "Cl"))                                                     { strcpy(typeSecondAtom, "Cl");   }
    else if (!strcmp(type2, "Br"))                                                     { strcpy(typeSecondAtom, "Br");   }
    else if (!strcmp(type2, "I"))                                                      { strcpy(typeSecondAtom, "I");    }

    // Assign values.
    if      (!strcmp(typeFirstAtom, "Csp3") && !strcmp(typeSecondAtom, "Csp3"))  { Distance=1.54; }
    else if (!strcmp(typeFirstAtom, "Csp2") && !strcmp(typeSecondAtom, "Csp2"))  { Distance=1.47; }
    else if (!strcmp(typeFirstAtom, "Csp")  && !strcmp(typeSecondAtom, "Csp"))   { Distance=1.37; }

    else if (!strcmp(typeFirstAtom, "Csp3") && !strcmp(typeSecondAtom, "Csp2"))  { Distance=1.50; }
    else if (!strcmp(typeFirstAtom, "Csp2") && !strcmp(typeSecondAtom, "Csp3"))  { Distance=1.50; }
    else if (!strcmp(typeFirstAtom, "Csp3") && !strcmp(typeSecondAtom, "Csp"))   { Distance=1.46; }
    else if (!strcmp(typeFirstAtom, "Csp")  && !strcmp(typeSecondAtom, "Csp3"))  { Distance=1.46; }
    else if (!strcmp(typeFirstAtom, "Csp2") && !strcmp(typeSecondAtom, "Csp"))   { Distance=1.43; }
    else if (!strcmp(typeFirstAtom, "Csp")  && !strcmp(typeSecondAtom, "Csp2"))  { Distance=1.43; }

    else if (!strcmp(typeFirstAtom, "Csp3") && !strcmp(typeSecondAtom, "N"))     { Distance=1.47; }
    else if (!strcmp(typeFirstAtom, "N")    && !strcmp(typeSecondAtom, "Csp3"))  { Distance=1.47; }
    else if (!strcmp(typeFirstAtom, "Csp2") && !strcmp(typeSecondAtom, "N"))     { Distance=1.36; }
    else if (!strcmp(typeFirstAtom, "N")    && !strcmp(typeSecondAtom, "Csp2"))  { Distance=1.36; }
    else if (!strcmp(typeFirstAtom, "Csp")  && !strcmp(typeSecondAtom, "N"))     { Distance=1.14; }
    else if (!strcmp(typeFirstAtom, "N")    && !strcmp(typeSecondAtom, "Csp"))   { Distance=1.14; }

    else if (!strcmp(typeFirstAtom, "Csp3") && !strcmp(typeSecondAtom, "O"))     { Distance=1.43; }
    else if (!strcmp(typeFirstAtom, "O")    && !strcmp(typeSecondAtom, "Csp3"))  { Distance=1.43; }
    else if (!strcmp(typeFirstAtom, "Csp2") && !strcmp(typeSecondAtom, "O"))     { Distance=1.43; }
    else if (!strcmp(typeFirstAtom, "O")    && !strcmp(typeSecondAtom, "Csp2"))  { Distance=1.43; }
    else if (!strcmp(typeFirstAtom, "Csp")  && !strcmp(typeSecondAtom, "O"))     { Distance=1.43; }
    else if (!strcmp(typeFirstAtom, "O")    && !strcmp(typeSecondAtom, "Csp"))   { Distance=1.43; }

    else if (!strcmp(typeFirstAtom, "Csp3") && !strcmp(typeSecondAtom, "S"))     { Distance=1.82; }
    else if (!strcmp(typeFirstAtom, "S")    && !strcmp(typeSecondAtom, "Csp3"))  { Distance=1.82; }
    else if (!strcmp(typeFirstAtom, "Csp2") && !strcmp(typeSecondAtom, "S"))     { Distance=1.82; }
    else if (!strcmp(typeFirstAtom, "S")    && !strcmp(typeSecondAtom, "Csp2"))  { Distance=1.82; }
    else if (!strcmp(typeFirstAtom, "Csp")  && !strcmp(typeSecondAtom, "S"))     { Distance=1.82; }
    else if (!strcmp(typeFirstAtom, "S")    && !strcmp(typeSecondAtom, "Csp"))   { Distance=1.82; }

    else if (!strcmp(typeFirstAtom, "Csp3") && !strcmp(typeSecondAtom, "F"))     { Distance=1.35; }
    else if (!strcmp(typeFirstAtom, "F")    && !strcmp(typeSecondAtom, "Csp3"))  { Distance=1.35; }
    else if (!strcmp(typeFirstAtom, "Csp2") && !strcmp(typeSecondAtom, "F"))     { Distance=1.35; }
    else if (!strcmp(typeFirstAtom, "F")    && !strcmp(typeSecondAtom, "Csp2"))  { Distance=1.35; }
    else if (!strcmp(typeFirstAtom, "Csp")  && !strcmp(typeSecondAtom, "F"))     { Distance=1.35; }
    else if (!strcmp(typeFirstAtom, "F")    && !strcmp(typeSecondAtom, "Csp"))   { Distance=1.35; }

    else if (!strcmp(typeFirstAtom, "Csp3") && !strcmp(typeSecondAtom, "Cl"))    { Distance=1.77; }
    else if (!strcmp(typeFirstAtom, "Cl")   && !strcmp(typeSecondAtom, "Csp3"))  { Distance=1.77; }
    else if (!strcmp(typeFirstAtom, "Csp2") && !strcmp(typeSecondAtom, "Cl"))    { Distance=1.77; }
    else if (!strcmp(typeFirstAtom, "Cl")   && !strcmp(typeSecondAtom, "Csp2"))  { Distance=1.77; }
    else if (!strcmp(typeFirstAtom, "Csp")  && !strcmp(typeSecondAtom, "Cl"))    { Distance=1.77; }
    else if (!strcmp(typeFirstAtom, "Cl")   && !strcmp(typeSecondAtom, "Csp"))   { Distance=1.77; }

    else if (!strcmp(typeFirstAtom, "Csp3") && !strcmp(typeSecondAtom, "Br"))    { Distance=1.94; }
    else if (!strcmp(typeFirstAtom, "Br")   && !strcmp(typeSecondAtom, "Csp3"))  { Distance=1.94; }
    else if (!strcmp(typeFirstAtom, "Csp2") && !strcmp(typeSecondAtom, "Br"))    { Distance=1.94; }
    else if (!strcmp(typeFirstAtom, "Br")   && !strcmp(typeSecondAtom, "Csp2"))  { Distance=1.94; }
    else if (!strcmp(typeFirstAtom, "Csp")  && !strcmp(typeSecondAtom, "Br"))    { Distance=1.94; }
    else if (!strcmp(typeFirstAtom, "Br")   && !strcmp(typeSecondAtom, "Csp"))   { Distance=1.94; }

    else if (!strcmp(typeFirstAtom, "Csp3") && !strcmp(typeSecondAtom, "I"))     { Distance=2.14; }
    else if (!strcmp(typeFirstAtom, "I")    && !strcmp(typeSecondAtom, "Csp3"))  { Distance=2.14; }
    else if (!strcmp(typeFirstAtom, "Csp2") && !strcmp(typeSecondAtom, "I"))     { Distance=2.14; }
    else if (!strcmp(typeFirstAtom, "I")    && !strcmp(typeSecondAtom, "Csp2"))  { Distance=2.14; }
    else if (!strcmp(typeFirstAtom, "Csp")  && !strcmp(typeSecondAtom, "I"))     { Distance=2.14; }
    else if (!strcmp(typeFirstAtom, "I")    && !strcmp(typeSecondAtom, "Csp"))   { Distance=2.14; }

    else if (!strcmp(typeFirstAtom,"N")     && !strcmp(typeSecondAtom,"N"))      { Distance=1.45; }
    else if (!strcmp(typeFirstAtom,"O")     && !strcmp(typeSecondAtom,"O"))      { Distance=1.48; }
    else                                                                         { Distance=1.50; }

    return Distance;
}

