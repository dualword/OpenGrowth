#include "OpenGrowth.h"

double dist(double const atom1[3], double const atom2[3]);
int    Random(Parameters & parameters, double const probaFile[]);
double Energy(Parameters const & parameters, Protein const & protein, Molecule & molecule);
int    StericClash(Parameters const & parameters, Protein const & protein, Molecule const & molecule);
void   OptimizationPosition(Parameters & parameters, Protein const & protein, Molecule & molecule);
float  OptimizationGeom(Parameters const & parameters, Protein const & protein, Molecule & molecule, int & errorFF, int const & optParameter);

// This function adds the first fragment in the active site when DENOVO mode is selected.
void  FirstFragment(Parameters & parameters, Protein const protein[], Molecule molecule[], double & averageEnergy, string const typeOfSnapshots, Molecule const fragment[], double const probaFirstFrag[], double randomAxe[], double randomPoint[], int ligandDescription[], int & errorFF)
{
    // Depending on which kind of snapshots we are using, there are not the same number of them.
    int snapshotNumber=0;
    if      (typeOfSnapshots == "rotamers")   { snapshotNumber = parameters.rotamersNumber;   }
    else if (typeOfSnapshots == "conformers") { snapshotNumber = parameters.conformersNumber; }

    // Define some parameters.
    int indexFragment=0;
    int stericClashes=0;
    bool success=false;
    Molecule firstFragment;

    // There are two loops. We randomly choose the first fragment, and try to put it MAX_ITERATIONS times in the active site. If it fails, another fragment is chosen.
    // This should avoid the program to stay in an infinite loop without success.
    while(!success) {
        // Choose the first fragment, either from the previous growth, randomly, according to the FOG probabilities or according to a given REGROW_FILE.
        if      (parameters.rotamersNumber!=0 && typeOfSnapshots=="conformers") { indexFragment = ligandDescription[0];                                  }
        else if (parameters.growthMode=="RANDOM")                               { indexFragment = (parameters.random() % parameters.fragmentListSize)+1; }
        else if (parameters.growthMode=="BIASED")                               { indexFragment = Random(parameters, probaFirstFrag)+1;                  }
        else if (parameters.growthMode=="FOG")                                  { indexFragment = Random(parameters, probaFirstFrag)+1;                  }
        else if (parameters.growthMode=="REGROW")                               { indexFragment = ligandDescription[0];                                  }

        // Display information
        if (parameters.verbose>=3) { cout << "Fragment chosen for trial: #" << indexFragment << ", " << fragment[indexFragment-1].obMol.GetTitle() << endl; }

        // Store the first fragment and its coordinates.
        firstFragment=fragment[indexFragment-1];
        double fragment0[firstFragment.obMol.NumAtoms()][3];
        for (unsigned int j=1; j<=firstFragment.obMol.NumAtoms() ; j++ ) {
            OpenBabel::OBAtom *atom;
            atom = firstFragment.obMol.GetAtom(j);
            fragment0[j-1][0]=atom->GetX();
            fragment0[j-1][1]=atom->GetY();
            fragment0[j-1][2]=atom->GetZ();
            }

        // Translation 1 - The first atom of the first fragment is put in (0,0,0). This could be avoided since all the fragments have first atom in (0,0,0), but we keep it just in case.
        double fragment1[firstFragment.obMol.NumAtoms()][3];
        for (unsigned int j=1; j<=firstFragment.obMol.NumAtoms() ; j++ ) {
            fragment1[j-1][0]=fragment0[j-1][0]-fragment0[0][0];
            fragment1[j-1][1]=fragment0[j-1][1]-fragment0[0][1];
            fragment1[j-1][2]=fragment0[j-1][2]-fragment0[0][2];
            }

        int counterIterations=0;
        // Do the following block as long as steric clashes is found and we have not tried it too many times.
        while(!success && (counterIterations<parameters.maxIterations) ) {
            // We take random numbers between 0 and 99999, divide by 1000 to have them between 00.000 and 99.999, and center the distribution. It is done to have a random direction
            // for the second atom and we do it if it has not been done before.
            if ( (typeOfSnapshots == "rotamers") || ((typeOfSnapshots == "conformers") && (parameters.rotamersNumber == 0)) ) {
                 while (dist(randomAxe, fragment1[0])==0) {
                     randomAxe[0] = (parameters.random() % 100000) / 1000.0 - 50.0 ;
                     randomAxe[1] = (parameters.random() % 100000) / 1000.0 - 50.0 ;
                     randomAxe[2] = (parameters.random() % 100000) / 1000.0 - 50.0 ;
                     }
                 }

            // Rotation 1.
            // We will do a first rotation to align the 2 first atoms with the fragment1[0]-randomAxe axe. Let's define two vectors U (fragment1[1]-fragment1[0]) and
            // V (randomAxe-fragment1[0]) and normalize them. After the rotation, U will be aligned with V. We define W=UxV/||UxV||. 
            float uX = (fragment1[1][0]-fragment1[0][0])/dist(fragment1[1], fragment1[0]);
            float uY = (fragment1[1][1]-fragment1[0][1])/dist(fragment1[1], fragment1[0]);
            float uZ = (fragment1[1][2]-fragment1[0][2])/dist(fragment1[1], fragment1[0]);
            float vX = (randomAxe[0]-fragment1[0][0])/dist(randomAxe, fragment1[0]);
            float vY = (randomAxe[1]-fragment1[0][1])/dist(randomAxe, fragment1[0]);
            float vZ = (randomAxe[2]-fragment1[0][2])/dist(randomAxe, fragment1[0]);
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
            rotationMatrix1[1][0] = wX*wY*(1-cosAlpha) + wZ*sinAlpha;
            rotationMatrix1[2][0] = wX*wZ*(1-cosAlpha) - wY*sinAlpha;    // rotationMatrix1 has the following form:
            rotationMatrix1[0][1] = wX*wY*(1-cosAlpha) - wZ*sinAlpha;    // [0][0] [0][1] [0][2]
            rotationMatrix1[1][1] = cosAlpha + pow(wY,2)*(1-cosAlpha);   // [1][0] [1][1] [1][2]
            rotationMatrix1[2][1] = wY*wZ*(1-cosAlpha) + wX*sinAlpha;    // [2][0] [2][1] [2][2]
            rotationMatrix1[0][2] = wX*wZ*(1-cosAlpha) + wY*sinAlpha;
            rotationMatrix1[1][2] = wY*wZ*(1-cosAlpha) - wX*sinAlpha;
            rotationMatrix1[2][2] = cosAlpha + pow(wZ,2)*(1-cosAlpha);
            // The alignement is made.
            double fragment2[firstFragment.obMol.NumAtoms()][3];
            for (unsigned int j=1; j<=firstFragment.obMol.NumAtoms() ; j++ ) {
                fragment2[j-1][0] = rotationMatrix1[0][0]*fragment1[j-1][0] + rotationMatrix1[0][1]*fragment1[j-1][1] + rotationMatrix1[0][2]*fragment1[j-1][2] ;
                fragment2[j-1][1] = rotationMatrix1[1][0]*fragment1[j-1][0] + rotationMatrix1[1][1]*fragment1[j-1][1] + rotationMatrix1[1][2]*fragment1[j-1][2] ;
                fragment2[j-1][2] = rotationMatrix1[2][0]*fragment1[j-1][0] + rotationMatrix1[2][1]*fragment1[j-1][1] + rotationMatrix1[2][2]*fragment1[j-1][2] ;
                }

            // Rotation 2 + Translation 2.
            // We will rotate around the axe fragment2[0]-fragment2[1] (vector R) every 2*PI/rotationPrecision.
            float rX = (fragment2[0][0]-fragment2[1][0])/dist(fragment2[0], fragment2[1]);
            float rY = (fragment2[0][1]-fragment2[1][1])/dist(fragment2[0], fragment2[1]);
            float rZ = (fragment2[0][2]-fragment2[1][2])/dist(fragment2[0], fragment2[1]);
            float Beta=2*PI/(parameters.rotationPrecision);
            float rotationMatrix2[3][3];
            // The first atom will be put randomly in the active site. We pick a random point between 0 and 2*bindingSize, center around 0.00 and then around the binding site.
            // We do it if it has not been done before.
            if ( (typeOfSnapshots == "rotamers") || ((typeOfSnapshots == "conformers") && (parameters.rotamersNumber == 0)) ) {
                 int randomParameter = (int) 2*parameters.bindingSize*100;
                 randomPoint[0] = (parameters.random() % randomParameter) / 100.0 - parameters.bindingSize + parameters.bindingSite[0];
                 randomPoint[1] = (parameters.random() % randomParameter) / 100.0 - parameters.bindingSize + parameters.bindingSite[1];
                 randomPoint[2] = (parameters.random() % randomParameter) / 100.0 - parameters.bindingSize + parameters.bindingSite[2];
                 }

            // For storing the best orientation (with the lowest energy), we do 2 loops. We first do a loop on the rotationPrecision because it saves some
            // computational time (fragment3 is computed only once and not for every conformers).
            for (int k=0 ; k < parameters.rotationPrecision ; k++) {
                // Rotation 2 + Translation 2.
                rotationMatrix2[0][0] = cos(k*Beta) + pow(rX,2)*(1-cos(k*Beta));
                rotationMatrix2[1][0] = rX*rY*(1-cos(k*Beta)) + rZ*sin(k*Beta) ;
                rotationMatrix2[2][0] = rX*rZ*(1-cos(k*Beta)) - rY*sin(k*Beta) ;    // rotationMatrix2 has the following form:
                rotationMatrix2[0][1] = rX*rY*(1-cos(k*Beta)) - rZ*sin(k*Beta) ;    // [0][0] [0][1] [0][2]
                rotationMatrix2[1][1] = cos(k*Beta) + pow(rY,2)*(1-cos(k*Beta));    // [1][0] [1][1] [1][2]
                rotationMatrix2[2][1] = rY*rZ*(1-cos(k*Beta)) + rX*sin(k*Beta) ;    // [2][0] [2][1] [2][2]
                rotationMatrix2[0][2] = rX*rZ*(1-cos(k*Beta)) + rY*sin(k*Beta) ;
                rotationMatrix2[1][2] = rY*rZ*(1-cos(k*Beta)) - rX*sin(k*Beta) ;
                rotationMatrix2[2][2] = cos(k*Beta) + pow(rZ,2)*(1-cos(k*Beta));
                double fragment3[firstFragment.obMol.NumAtoms()][3];
                for (unsigned int j=1; j<=firstFragment.obMol.NumAtoms() ; j++ ) {
                    fragment3[j-1][0] = rotationMatrix2[0][0]*fragment2[j-1][0] + rotationMatrix2[0][1]*fragment2[j-1][1] + rotationMatrix2[0][2]*fragment2[j-1][2] + randomPoint[0];
                    fragment3[j-1][1] = rotationMatrix2[1][0]*fragment2[j-1][0] + rotationMatrix2[1][1]*fragment2[j-1][1] + rotationMatrix2[1][2]*fragment2[j-1][2] + randomPoint[1];
                    fragment3[j-1][2] = rotationMatrix2[2][0]*fragment2[j-1][0] + rotationMatrix2[2][1]*fragment2[j-1][1] + rotationMatrix2[2][2]*fragment2[j-1][2] + randomPoint[2];
                    }

                // For each snapshot, look for the best rotamer.
                for (int s=0 ; s < snapshotNumber ; s++) {
                    // We will need a temporary molecule, and we update its coordinates.
                    Molecule tempMolecule;
                    tempMolecule=firstFragment;
                    for (unsigned int j=1; j<=tempMolecule.obMol.NumAtoms() ; j++ ) {
                        OpenBabel::OBAtom *atom;
                        atom = tempMolecule.obMol.GetAtom(j);
                        atom->SetVector(fragment3[j-1][0], fragment3[j-1][1], fragment3[j-1][2]);
                        }

                    // Optimize the molecule position for every rotamer.
                    if ((parameters.optimizationMode/10)>=3) { OptimizationGeom(parameters, protein[s], tempMolecule, errorFF, 1);}
                    if ((parameters.optimizationMode%10)>=3) { OptimizationPosition(parameters, protein[s], tempMolecule);        }

                    // Compute some information.
                    tempMolecule.stericClashes = StericClash(parameters, protein[s], tempMolecule);
                    tempMolecule.interactionEnergy = Energy(parameters, protein[s], tempMolecule);

                    // We save the tempMolecule if it is "better" than the old one.
                    if (k == 0)     { molecule[s]=tempMolecule; }
                    else {
                        if      ( !tempMolecule.stericClashes && molecule[s].stericClashes   )                                      { molecule[s] = tempMolecule; }
                        else if ( !tempMolecule.stericClashes && (tempMolecule.interactionEnergy < molecule[s].interactionEnergy) ) { molecule[s] = tempMolecule; }
                        }
                    }
                }

            // Check for steric clashes in all the molecules.
            stericClashes=0;
            for (int s=0 ; s < snapshotNumber ; s++) {
                if ( molecule[s].stericClashes ) { stericClashes++; }
                }
            if (stericClashes==0) { success=true; }
            counterIterations++;
            }
        }
    if (parameters.verbose>=3) { cout << "\t\tThe first atom will be put in " << randomPoint[0] << " " << randomPoint[1] << " " << randomPoint[2] << "." << endl; }

    // Optimize the molecule position, once the best rotamer is found.
    for (int s=0 ; s < snapshotNumber ; s++) {
        if ((parameters.optimizationMode/10)>=2) { OptimizationGeom(parameters, protein[s], molecule[s], errorFF, 1); }
        if ((parameters.optimizationMode%10)>=2) { OptimizationPosition(parameters, protein[s], molecule[s]);         }
        }

    // Compute the average energy of the molecules.
    averageEnergy=0;
    long double sumEnergy=0;
    if      (snapshotNumber==1)                   { averageEnergy=molecule[0].interactionEnergy; }
    else if (parameters.averageType=="BOLTZMANN") {
        // Compute the Boltzmann average energy.
        long double partitionFunction=0;
        for (int s=0 ; s < snapshotNumber ; s++) {
            partitionFunction += exp(-BETA*molecule[s].interactionEnergy);
            sumEnergy += molecule[s].interactionEnergy*exp(-BETA*molecule[s].interactionEnergy);
            }
        averageEnergy = sumEnergy/partitionFunction;
        }
    else if (parameters.averageType=="ARITHMETIC") {
        // Compute the arithmetic average energy.
        for (int s=0 ; s < snapshotNumber ; s++) {
            averageEnergy += molecule[s].interactionEnergy;
            }
        averageEnergy = averageEnergy/snapshotNumber;
        }
    else if (parameters.averageType=="LOWESTSCORE") {
        // Keep only the lowest energy.
        averageEnergy = molecule[0].interactionEnergy; 
        for (int s=1 ; s < snapshotNumber ; s++) {
            if(molecule[s].interactionEnergy<averageEnergy) { averageEnergy = molecule[s].interactionEnergy; }
            }
        }

    // Save the information to describe the molecule afterwards.
    if (parameters.growthMode=="RANDOM" || parameters.growthMode=="BIASED" || parameters.growthMode=="FOG") {
        ligandDescription[0]=indexFragment;
        }

    // Display results if we are dealing with rotamers and there will be no conformers search, or if we are dealing with conformers.
    if ( ((typeOfSnapshots == "rotamers") && (parameters.conformersNumber == 0)) || (typeOfSnapshots == "conformers") ) {
        if (parameters.verbose>=1) { cout << "Fragment " << molecule[0].numberResidues << ": " << fragment[indexFragment-1].obMol.GetTitle() << endl; }
        if (parameters.verbose>=2) { cout << "\tEnergy: " << averageEnergy << endl; }
        }
}

