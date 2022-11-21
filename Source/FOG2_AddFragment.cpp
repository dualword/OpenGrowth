#include "OpenGrowth.h"

double dist(double const atom1[3], double const atom2[3]);
float  distNB(char const type1[], char const type2[]);
float  Optimization(Parameters const & parameters, Molecule & molecule, int & errorFF, int const & optParameter=1);
int    Random(Parameters & parameters, double const probaFile[]);
int    StericClash(Parameters const & parameters, Molecule const & molecule);

// This function adds a new fragment to an existing molecule.
void AddFragment(Parameters & parameters, Molecule const fragment[], double const probaTransition[][MAX_FRAGMENTS], double const probaFirstFrag[], Molecule & molecule, int & errorFF)
{
    Molecule newFragment;
    unsigned int indexH_L=0, indexH_F=0, indexFragment=0;

    // Check if we can grow linearly or branchedly.
    int canGrowLinearly=0;
    int canGrowBranchedly=0;
    for (unsigned int j=0; j<molecule.growingSites.size() ; j++ ) {
        unsigned int indexAtom = molecule.growingSites[j];
        if (molecule.fragmentNeighbours[indexAtom-1]<=1) { canGrowLinearly++;   }
        if (molecule.fragmentNeighbours[indexAtom-1]>1)  { canGrowBranchedly++; }
        }

    // Decide which kind of growth we will perform and from which hydrogen. indexH_L will store the atom number in the molecule of this hydrogen.
    bool growLinearly=true;
    int badFragment;
    do {
        double randomNumber = (parameters.random() % 1000) / 1000.0 ;
        // Select indexH_L.
        if ( (randomNumber >= parameters.branchingProba) && (canGrowLinearly >= 1) )       {        // Grow linear
            do {
                unsigned int indexSite = parameters.random() % molecule.growingSites.size();
                indexH_L = molecule.growingSites[indexSite];
                } while (molecule.fragmentNeighbours[indexH_L-1] > 1 );
            growLinearly=true;
            }
        else if (randomNumber >= parameters.branchingProba)                                {        // Can't grow linear, so we grow branched
            unsigned int indexSite = parameters.random() % molecule.growingSites.size();
            indexH_L = molecule.growingSites[indexSite];
            growLinearly=false;
            }
        else if ( (randomNumber < parameters.branchingProba) && (canGrowBranchedly >= 1) ) {        // Grow branched
            do {
                unsigned int indexSite = parameters.random() % molecule.growingSites.size();
                indexH_L = molecule.growingSites[indexSite];
                } while (molecule.fragmentNeighbours[indexH_L-1] <= 1);
            growLinearly=false;
            }
        else if (randomNumber < parameters.branchingProba)                                 {        // Can't grow branched so grow linear
            unsigned int indexSite = parameters.random() % molecule.growingSites.size();
            indexH_L = molecule.growingSites[indexSite];
            growLinearly=true;
            }

        // Check that the fragment is allowed. This means that something can grow from here.
        badFragment=0;
        for (unsigned int j=0; j<parameters.forbiddenFragments.size() ; j++ ) {
            if (parameters.forbiddenFragments[j]==molecule.growth[indexH_L-1]) { badFragment++; }
            }
        } while (badFragment!=0);

    // By reading in the growth array, we can find to which fragment this hydrogen belongs.
    int fragmentType=molecule.growth[indexH_L-1];
    if (parameters.verbose>=3) { cout << "\t\tindexH_L (" << indexH_L << ") was successfully choosen. It is from the fragment type " << fragmentType << "." << endl; }
    // Choose the new fragment, either randomly or according to the FOG probabilities.
    if      (parameters.growthMode=="RANDOM") { indexFragment = (parameters.random() % parameters.fragmentListSize)+1;  }
    else if (parameters.growthMode=="BIASED") { indexFragment = Random(parameters, probaFirstFrag)+1;                   }
    else if (parameters.growthMode=="FOG")    { indexFragment = Random(parameters, probaTransition[fragmentType-1])+1;  }
    // Store the new fragment.
    newFragment=fragment[indexFragment-1];
    // Pick an hydrogen in the new fragment which describes the fragment properly i.e. which has a growth mode of the same value as the indexFragment.
    do {
        unsigned int indexSite = parameters.random() % newFragment.growingSites.size();
        indexH_F = newFragment.growingSites[indexSite];
        } while(newFragment.growth[indexH_F-1]!=indexFragment);
    if (parameters.verbose>=3) { cout << "\t\tindexH_F (" << indexH_F << ") was successfully choosen. It is from the fragment type " << indexFragment << "." << endl; }

    // Get the atom from the selected hydrogen of the ligand and its neighbour.
    OpenBabel::OBAtom *atomH_L;
    OpenBabel::OBAtom *atomN_L;
    atomH_L = molecule.obMol.GetAtom(indexH_L);
    OpenBabel::OBAtomAtomIter nbrL(atomH_L);
    int indexN_L = nbrL->GetIdx();
    atomN_L = molecule.obMol.GetAtom(indexN_L);

    // Get the atom from the selected hydrogen of the fragment and its neighbour.
    OpenBabel::OBAtom *atomH_F;
    OpenBabel::OBAtom *atomN_F;
    atomH_F = newFragment.obMol.GetAtom(indexH_F);
    OpenBabel::OBAtomAtomIter nbrF(atomH_F);
    int indexN_F = nbrF->GetIdx();
    atomN_F = newFragment.obMol.GetAtom(indexN_F);

    // Store the coordinates from the ligand.
    double ligand[molecule.obMol.NumAtoms()][3];
    for (unsigned int j=1; j<=molecule.obMol.NumAtoms() ; j++ ) {
        OpenBabel::OBAtom *atom;
        atom = molecule.obMol.GetAtom(j);
        ligand[j-1][0]=atom->GetX();
        ligand[j-1][1]=atom->GetY();
        ligand[j-1][2]=atom->GetZ();
        }

    // Store the coordinates from the fragment.
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

    // Rotation 1.
    // We will do a first rotation to align the H_F/N_F axe with the N_L/H_L axe. Let's define two vectors U (fragment1[indexN_F-1]-fragment1[indexH_F-1]) and
    // V (ligand[indexH_L-1]-ligand[indexN_L-1]) and normalize them. After the rotation, U will be aligned with V. We define W=UxV/||UxV||.
    float uX = (fragment1[indexN_F-1][0]-fragment1[indexH_F-1][0])/dist(fragment1[indexN_F-1], fragment1[indexH_F-1]);
    float uY = (fragment1[indexN_F-1][1]-fragment1[indexH_F-1][1])/dist(fragment1[indexN_F-1], fragment1[indexH_F-1]);
    float uZ = (fragment1[indexN_F-1][2]-fragment1[indexH_F-1][2])/dist(fragment1[indexN_F-1], fragment1[indexH_F-1]);
    float vX = (ligand[indexH_L-1][0]-ligand[indexN_L-1][0])/dist(ligand[indexH_L-1], ligand[indexN_L-1]);
    float vY = (ligand[indexH_L-1][1]-ligand[indexN_L-1][1])/dist(ligand[indexH_L-1], ligand[indexN_L-1]);
    float vZ = (ligand[indexH_L-1][2]-ligand[indexN_L-1][2])/dist(ligand[indexH_L-1], ligand[indexN_L-1]);
    float wX = (uY*vZ - uZ*vY)/sqrt(pow(uY*vZ - uZ*vY,2) + pow(uZ*vX - uX*vZ,2) + pow(uX*vY - uY*vX,2));
    float wY = (uZ*vX - uX*vZ)/sqrt(pow(uY*vZ - uZ*vY,2) + pow(uZ*vX - uX*vZ,2) + pow(uX*vY - uY*vX,2));
    float wZ = (uX*vY - uY*vX)/sqrt(pow(uY*vZ - uZ*vY,2) + pow(uZ*vX - uX*vZ,2) + pow(uX*vY - uY*vX,2));
    // We need to know from which angle we want to rotate. We know with the dot and the cross products: U.V = U*V*cos(alpha) and UxV = U*V*sin(alpha). Since U and V are normalized to 1, the norm of
    // UxV gives us |sin(alpha)|. This is not enough, we need to know the sign of alpha. Thus, we use the following: sin(alpha) = det(W, U, V) where W=UxV/||UxV|| (the rule of Sarrus can be used).
    float cosAlpha = uX*vX + uY*vY + uZ*vZ;
    float sinAlpha = (wX*uY*vZ + uX*vY*wZ + vX*wY*uZ) - (vX*uY*wZ + wX*vY*uZ + uX*wY*vZ);
    // Let's write the rotation matrix around the vector W (http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle).
    float rotationMatrix1[3][3];
    rotationMatrix1[0][0] = cosAlpha + pow(wX,2)*(1-cosAlpha);
    rotationMatrix1[0][1] = wX*wY*(1-cosAlpha) + wZ*sinAlpha;
    rotationMatrix1[0][2] = wX*wZ*(1-cosAlpha) - wY*sinAlpha;       // rotationMatrix1 has the following form:
    rotationMatrix1[1][0] = wX*wY*(1-cosAlpha) - wZ*sinAlpha;       // [0][0] [1][0] [2][0]
    rotationMatrix1[1][1] = cosAlpha + pow(wY,2)*(1-cosAlpha);      // [0][1] [1][1] [2][1]
    rotationMatrix1[1][2] = wY*wZ*(1-cosAlpha) + wX*sinAlpha;       // [0][2] [1][2] [2][2]
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

    // Define newMolecule.
    Molecule newMolecule;
    // Create new arrays for growth, fragmentIndex and fragmentNeighbours.
    newMolecule.growth.clear();
    newMolecule.fragmentIndex.clear();
    newMolecule.fragmentNeighbours.clear();
    for (unsigned int j=1; j<=molecule.obMol.NumAtoms() ; j++ ) {
        newMolecule.growth.push_back(molecule.growth[j-1]);
        newMolecule.fragmentIndex.push_back(molecule.fragmentIndex[j-1]);
        newMolecule.fragmentNeighbours.push_back(molecule.fragmentNeighbours[j-1]);
        }
    // Update the neighbour atoms.
    int startingFragment=molecule.fragmentIndex[indexN_L-1];
    for (unsigned int j=1; j<=newMolecule.obMol.NumAtoms() ; j++ ) {
        if(newMolecule.fragmentIndex[j-1]==startingFragment) { newMolecule.fragmentNeighbours[j-1]++; }
        }
    // Add info from newFragment.
    for (unsigned int j=1; j<=newFragment.obMol.NumAtoms() ; j++ ) {
        newMolecule.growth.push_back(newFragment.growth[j-1]);
        newMolecule.fragmentIndex.push_back(molecule.numberResidues+1);
        newMolecule.fragmentNeighbours.push_back(1);
        }
    newMolecule.growth.erase(newMolecule.growth.begin()+indexH_L-1);
    newMolecule.growth.erase(newMolecule.growth.begin()+molecule.obMol.NumAtoms()-1+indexH_F-1);
    newMolecule.fragmentIndex.erase(newMolecule.fragmentIndex.begin()+indexH_L-1);
    newMolecule.fragmentIndex.erase(newMolecule.fragmentIndex.begin()+molecule.obMol.NumAtoms()-1+indexH_F-1);
    newMolecule.fragmentNeighbours.erase(newMolecule.fragmentNeighbours.begin()+indexH_L-1);
    newMolecule.fragmentNeighbours.erase(newMolecule.fragmentNeighbours.begin()+molecule.obMol.NumAtoms()-1+indexH_F-1);

    // Create new array for growingSites and update numberResidues.
    newMolecule.growingSites.clear();
    for (unsigned int j=0; j<molecule.growingSites.size() ; j++ ) {
        if       (molecule.growingSites[j]<indexH_L) { newMolecule.growingSites.push_back(molecule.growingSites[j]);   }
        else if  (molecule.growingSites[j]>indexH_L) { newMolecule.growingSites.push_back(molecule.growingSites[j]-1); }
        }
    for (unsigned int j=0; j<newFragment.growingSites.size(); j++ ) {
        if       (newFragment.growingSites[j]<indexH_F) { newMolecule.growingSites.push_back(molecule.obMol.NumAtoms()-1+newFragment.growingSites[j]);   }
        else if  (newFragment.growingSites[j]>indexH_F) { newMolecule.growingSites.push_back(molecule.obMol.NumAtoms()-1+newFragment.growingSites[j]-1); }
        }
    newMolecule.numberResidues = molecule.numberResidues+1;

    // We look for the best rotamer (with the lowest energy).
    for (int k=0 ; k < parameters.rotationPrecision ; k++) {
        // We will need a temporary molecule.
        Molecule tempMolecule;
        tempMolecule = newMolecule;
        tempMolecule.obMol.Clear();

        // Rotation 2 + Translation 2.
        rotationMatrix2[0][0] = cos(k*Beta) + pow(rX,2)*(1-cos(k*Beta));
        rotationMatrix2[0][1] = rX*rY*(1-cos(k*Beta)) + rZ*sin(k*Beta) ;
        rotationMatrix2[0][2] = rX*rZ*(1-cos(k*Beta)) - rY*sin(k*Beta) ;       // rotationMatrix2 has the following form:
        rotationMatrix2[1][0] = rX*rY*(1-cos(k*Beta)) - rZ*sin(k*Beta) ;       // [0][0] [1][0] [2][0]
        rotationMatrix2[1][1] = cos(k*Beta) + pow(rY,2)*(1-cos(k*Beta));       // [0][1] [1][1] [2][1]
        rotationMatrix2[1][2] = rY*rZ*(1-cos(k*Beta)) + rX*sin(k*Beta) ;       // [0][2] [1][2] [2][2]
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
        tempMolecule.obMol = molecule.obMol;
        tempMolecule.obMol += newFragment.obMol;
        // Delete the atom from the ligand.
        OpenBabel::OBAtom *atomToDeleteLigand;
        atomToDeleteLigand = tempMolecule.obMol.GetAtom(indexH_L);
        tempMolecule.obMol.DeleteAtom(atomToDeleteLigand);
        // Delete the atom from the fragment.
        OpenBabel::OBAtom *atomToDeleteFragment;
        atomToDeleteFragment = tempMolecule.obMol.GetAtom(molecule.obMol.NumAtoms()-1+indexH_F);
        tempMolecule.obMol.DeleteAtom(atomToDeleteFragment);
        // Create a new bond.
        tempMolecule.obMol.AddBond(indexN_L, molecule.obMol.NumAtoms()-1+indexN_F, 1);
        // Check for clashes.
        tempMolecule.stericClashes = StericClash(parameters, tempMolecule);

        // Accept a new orientation of the new fragment if there were steric clashes and there are not in the new orientation, or if the new energy is lower than the previous one (without clashes).
        if (k == 0)     { newMolecule = tempMolecule; }
        else {
            if      ( !tempMolecule.stericClashes && newMolecule.stericClashes )                                                                              { newMolecule = tempMolecule; }
            else if ( !tempMolecule.stericClashes && (Optimization(parameters, tempMolecule, errorFF, 0)<Optimization(parameters, newMolecule, errorFF, 0)) ) { newMolecule = tempMolecule; }
            }

        tempMolecule.obMol.Clear();
        }

    if (errorFF==0 && !StericClash(parameters, newMolecule)) {
        // Optimize the molecule geometry, once the best rotamer is found.
        Optimization(parameters, newMolecule, errorFF, 1);

        // Save the new molecule.
        molecule.obMol.Clear();
        molecule = newMolecule;

        // Display information.
        if (parameters.verbose>=1) {
            cout << endl << "Fragment " << molecule.numberResidues << ": #" << indexFragment << ", " << newFragment.obMol.GetTitle() << endl;
            }
        if (parameters.verbose>=3) {
            cout << "\t\tAtom " << indexH_L << " of the current ligand will be used. It is bound to the atom " << indexN_L << " (fragment " << molecule.fragmentIndex[indexN_L-1] << ")." << endl;
            cout << "\t\tAtom " << indexH_F << " of the new fragment will be used. It is bound to the atom " << indexN_F << "." << endl;
            cout << "\t\tWe will grow " ; if(growLinearly==true) { cout << "linearly." << endl; } else { cout << "branchedly." << endl; }
            }
        }

    newFragment.obMol.Clear();
    newMolecule.obMol.Clear();
}

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
    else if (!strcmp(typeFirstAtom, "Cl")    && !strcmp(typeSecondAtom, "Csp3")) { Distance=1.77; }
    else if (!strcmp(typeFirstAtom, "Csp2") && !strcmp(typeSecondAtom, "Cl"))    { Distance=1.77; }
    else if (!strcmp(typeFirstAtom, "Cl")    && !strcmp(typeSecondAtom, "Csp2")) { Distance=1.77; }
    else if (!strcmp(typeFirstAtom, "Csp")  && !strcmp(typeSecondAtom, "Cl"))    { Distance=1.77; }
    else if (!strcmp(typeFirstAtom, "Cl")    && !strcmp(typeSecondAtom, "Csp"))  { Distance=1.77; }

    else if (!strcmp(typeFirstAtom, "Csp3") && !strcmp(typeSecondAtom, "Br"))    { Distance=1.94; }
    else if (!strcmp(typeFirstAtom, "Br")    && !strcmp(typeSecondAtom, "Csp3")) { Distance=1.94; }
    else if (!strcmp(typeFirstAtom, "Csp2") && !strcmp(typeSecondAtom, "Br"))    { Distance=1.94; }
    else if (!strcmp(typeFirstAtom, "Br")    && !strcmp(typeSecondAtom, "Csp2")) { Distance=1.94; }
    else if (!strcmp(typeFirstAtom, "Csp")  && !strcmp(typeSecondAtom, "Br"))    { Distance=1.94; }
    else if (!strcmp(typeFirstAtom, "Br")    && !strcmp(typeSecondAtom, "Csp"))  { Distance=1.94; }

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

