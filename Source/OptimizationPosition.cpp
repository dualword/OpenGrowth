#include "OpenGrowth.h"

double dist(double const atom1[3], double const atom2[3]);
double Energy(Parameters const & parameters, Protein const & protein, Molecule & molecule);
int    StericClash(Parameters const & parameters, Protein const & protein, Molecule const & molecule);
void   Translation(Parameters & parameters, Protein const & protein, Molecule & molecule);
void   Rotation(Parameters & parameters, Protein const & protein, Molecule & molecule);
void   CenterOfMolecule(Molecule const & molecule, double moleculeCenter[]);

// This functions optimizes the position of a molecule as a whole to see if the energy decreases. It randomly translates or rotates the molecule along/around
// a random axe chosen OPTIMIZATION_NUMBER times, and OPTIMIZATION_ITERATIONS moves are made (either of OPTIMIZATION_DISTANCE or OPTIMIZATION_ANGLE).
void OptimizationPosition(Parameters & parameters, Protein const & protein, Molecule & molecule)
{
    for (int i=0 ; i < parameters.optimizationNumber ; i++) {
        int choice = parameters.random() % 2;
        if      (choice==0) { Translation(parameters, protein, molecule); }
        else if (choice==1) { Rotation(parameters, protein, molecule);    }
        }
}

// This function makes a translation of a molecule as a whole along a random direction.
void Translation(Parameters & parameters, Protein const & protein, Molecule & molecule)
{
    // Find the center of the molecule.
    double moleculeCenter[3];
    CenterOfMolecule(molecule, moleculeCenter);

    // Define a random direction toward which the molecule will be translated.
    double randomDirection[3];
    do {
        randomDirection[0] = (parameters.random() % 100000) / 1000.0 - 50.0 + parameters.bindingSite[0];
        randomDirection[1] = (parameters.random() % 100000) / 1000.0 - 50.0 + parameters.bindingSite[1];
        randomDirection[2] = (parameters.random() % 100000) / 1000.0 - 50.0 + parameters.bindingSite[2];
        } while (dist(moleculeCenter, randomDirection)==0);

    // Define a new molecule.
    Molecule translatedMolecule;
    translatedMolecule=molecule;

    // We do one loop between 1 and parameters.optimizationIterations, but every time we check for two translations.
    for (int j=1 ; j <= parameters.optimizationIterations ; j++) {
        // Translation vector.
        double translation[3];
        translation[0] = j*parameters.optimizationDistance*(randomDirection[0] - moleculeCenter[0])/dist(randomDirection, moleculeCenter);
        translation[1] = j*parameters.optimizationDistance*(randomDirection[1] - moleculeCenter[1])/dist(randomDirection, moleculeCenter);
        translation[2] = j*parameters.optimizationDistance*(randomDirection[2] - moleculeCenter[2])/dist(randomDirection, moleculeCenter);

        // Translate molecule in one direction.
        for (unsigned int i=1 ; i <= translatedMolecule.obMol.NumAtoms() ; i++) {
            OpenBabel::OBAtom *atom;
            atom = translatedMolecule.obMol.GetAtom(i);
            atom->SetVector(atom->GetX()-translation[0], atom->GetY()-translation[1], atom->GetZ()-translation[2]);
            }

        translatedMolecule.interactionEnergy = Energy(parameters, protein, translatedMolecule);
        translatedMolecule.stericClashes = StericClash(parameters, protein, translatedMolecule);
        // Keep the translated molecule: 1) if there were steric clashes, and there are no more. 2) if the energy is lower than the previous best one and there are no clashes.
        if      ( !translatedMolecule.stericClashes && molecule.stericClashes )                                              { molecule = translatedMolecule; }
        else if ( !translatedMolecule.stericClashes && (translatedMolecule.interactionEnergy < molecule.interactionEnergy) ) { molecule = translatedMolecule; }

        // Translate molecule in the other direction.
        for (unsigned int i=1 ; i <= translatedMolecule.obMol.NumAtoms() ; i++) {
            OpenBabel::OBAtom *atom;
            atom = translatedMolecule.obMol.GetAtom(i);
            atom->SetVector(atom->GetX()+translation[0], atom->GetY()+translation[1], atom->GetZ()+translation[2]);
            }

        translatedMolecule.interactionEnergy = Energy(parameters, protein, translatedMolecule);
        translatedMolecule.stericClashes = StericClash(parameters, protein, translatedMolecule);
        // Keep the translated molecule: 1) if there were steric clashes, and there are no more. 2) if the energy is lower than the previous best one and there are no clashes.
        if      ( !translatedMolecule.stericClashes && molecule.stericClashes )                                              { molecule = translatedMolecule; }
        else if ( !translatedMolecule.stericClashes && (translatedMolecule.interactionEnergy < molecule.interactionEnergy) ) { molecule = translatedMolecule; }
        }
}

// This function makes a rotation of a molecule as a whole around a random direction.
void Rotation(Parameters & parameters, Protein const & protein, Molecule & molecule)
{
    // Find the center of the molecule.
    double moleculeCenter[3];
    CenterOfMolecule(molecule, moleculeCenter);

    // Define a random direction around which the molecule will be rotated.
    double origin[3] = {0.0, 0.0, 0.0};
    double randomDirection[3];
    do {
        randomDirection[0] = (parameters.random() % 100000) / 1000.0 - 50.0;
        randomDirection[1] = (parameters.random() % 100000) / 1000.0 - 50.0;
        randomDirection[2] = (parameters.random() % 100000) / 1000.0 - 50.0;
        } while (dist(randomDirection, origin)==0);

    // Normalized axe of rotation.
    double rotation[3];
    rotation[0] = (randomDirection[0])/dist(randomDirection, origin);
    rotation[1] = (randomDirection[1])/dist(randomDirection, origin);
    rotation[2] = (randomDirection[2])/dist(randomDirection, origin);

    // Store the coordinates from the fragment.
    double molecule0[molecule.obMol.NumAtoms()][3];
    for (unsigned int j=1; j<=molecule.obMol.NumAtoms() ; j++ ) {
        OpenBabel::OBAtom *atom;
        atom = molecule.obMol.GetAtom(j);
        molecule0[j-1][0]=atom->GetX();
        molecule0[j-1][1]=atom->GetY();
        molecule0[j-1][2]=atom->GetZ();
        }

    // We need to know the borders of the molecule for the next step. Get coordinates of the first atom, then compare values.
    double Xmin=molecule0[0][0], Xmax=molecule0[0][0];
    double Ymin=molecule0[0][1], Ymax=molecule0[0][1];
    double Zmin=molecule0[0][2], Zmax=molecule0[0][2];
    for (unsigned int j=1; j<=molecule.obMol.NumAtoms() ; j++ ) {
        if (molecule0[j-1][0]<Xmin) { Xmin = molecule0[j-1][0]; }
        if (molecule0[j-1][0]>Xmax) { Xmax = molecule0[j-1][0]; }
        if (molecule0[j-1][1]<Ymin) { Ymin = molecule0[j-1][1]; }
        if (molecule0[j-1][1]>Ymax) { Ymax = molecule0[j-1][1]; }
        if (molecule0[j-1][2]<Zmin) { Zmin = molecule0[j-1][2]; }
        if (molecule0[j-1][2]>Zmax) { Zmax = molecule0[j-1][2]; }
        }

    // We will then define a random point within the molecule cartesian space.
    double randomPoint[3];
    randomPoint[0] = (parameters.random() % 100000) * (Xmax-Xmin) / 100000.0 + Xmin;
    randomPoint[1] = (parameters.random() % 100000) * (Ymax-Ymin) / 100000.0 + Ymin;
    randomPoint[2] = (parameters.random() % 100000) * (Zmax-Zmin) / 100000.0 + Zmin;

    // Translate molecule in (0,0,0).
    double translatedMolecule[molecule.obMol.NumAtoms()][3];
    for (unsigned int j=1; j<=molecule.obMol.NumAtoms() ; j++ ) {
        translatedMolecule[j-1][0] = molecule0[j-1][0] - randomPoint[0];
        translatedMolecule[j-1][1] = molecule0[j-1][1] - randomPoint[1];
        translatedMolecule[j-1][2] = molecule0[j-1][2] - randomPoint[2];
        }

    // Define a new molecule.
    Molecule rotatedMolecule;
    rotatedMolecule=molecule;

    float Beta=parameters.optimizationAngle*PI/180;
    // We do one loop between 1 and parameters.optimizationIterations, but every time we check for two rotations.
    for (int j=1 ; j <= parameters.optimizationIterations ; j++) {
        // Rotation in one direction.
        float rotationMatrix1[3][3];
        rotationMatrix1[0][0] = cos(j*Beta) + pow(rotation[0],2)*(1-cos(j*Beta));
        rotationMatrix1[1][0] = rotation[0]*rotation[1]*(1-cos(j*Beta)) + rotation[2]*sin(j*Beta) ;
        rotationMatrix1[2][0] = rotation[0]*rotation[2]*(1-cos(j*Beta)) - rotation[1]*sin(j*Beta) ;          // rotationMatrix1 has the following form:
        rotationMatrix1[0][1] = rotation[0]*rotation[1]*(1-cos(j*Beta)) - rotation[2]*sin(j*Beta) ;          // [0][0] [0][1] [0][2]
        rotationMatrix1[1][1] = cos(j*Beta) + pow(rotation[1],2)*(1-cos(j*Beta));                            // [1][0] [1][1] [1][2]
        rotationMatrix1[2][1] = rotation[1]*rotation[2]*(1-cos(j*Beta)) + rotation[0]*sin(j*Beta) ;          // [2][0] [2][1] [2][2]
        rotationMatrix1[0][2] = rotation[0]*rotation[2]*(1-cos(j*Beta)) + rotation[1]*sin(j*Beta) ;
        rotationMatrix1[1][2] = rotation[1]*rotation[2]*(1-cos(j*Beta)) - rotation[0]*sin(j*Beta) ;
        rotationMatrix1[2][2] = cos(j*Beta) + pow(rotation[2],2)*(1-cos(j*Beta));
        double molecule1[molecule.obMol.NumAtoms()][3];
        for (unsigned int i=1; i<=molecule.obMol.NumAtoms() ; i++ ) {
            molecule1[i-1][0] = rotationMatrix1[0][0]*translatedMolecule[i-1][0] + rotationMatrix1[0][1]*translatedMolecule[i-1][1] + rotationMatrix1[0][2]*translatedMolecule[i-1][2] + randomPoint[0];
            molecule1[i-1][1] = rotationMatrix1[1][0]*translatedMolecule[i-1][0] + rotationMatrix1[1][1]*translatedMolecule[i-1][1] + rotationMatrix1[1][2]*translatedMolecule[i-1][2] + randomPoint[1];
            molecule1[i-1][2] = rotationMatrix1[2][0]*translatedMolecule[i-1][0] + rotationMatrix1[2][1]*translatedMolecule[i-1][1] + rotationMatrix1[2][2]*translatedMolecule[i-1][2] + randomPoint[2];
            }
        // Update molecule.
        for (unsigned int i=1 ; i <= rotatedMolecule.obMol.NumAtoms() ; i++) {
            OpenBabel::OBAtom *atom;
            atom = rotatedMolecule.obMol.GetAtom(i);
            atom->SetVector(molecule1[i-1][0], molecule1[i-1][1], molecule1[i-1][2]);
            }
        rotatedMolecule.interactionEnergy = Energy(parameters, protein, rotatedMolecule);
        rotatedMolecule.stericClashes = StericClash(parameters, protein, rotatedMolecule);
        // Keep the translated molecule: 1) if there were steric clashes, and there are no more. 2) if the energy is lower than the previous best one and there are no clashes.
        if      ( !rotatedMolecule.stericClashes && molecule.stericClashes )                                           { molecule = rotatedMolecule; }
        else if ( !rotatedMolecule.stericClashes && (rotatedMolecule.interactionEnergy < molecule.interactionEnergy) ) { molecule = rotatedMolecule; }

        // Rotation in the other direction.
        float rotationMatrix2[3][3];
        rotationMatrix2[0][0] = cos(-j*Beta) + pow(rotation[0],2)*(1-cos(-j*Beta));
        rotationMatrix2[1][0] = rotation[0]*rotation[1]*(1-cos(-j*Beta)) + rotation[2]*sin(-j*Beta) ;
        rotationMatrix2[2][0] = rotation[0]*rotation[2]*(1-cos(-j*Beta)) - rotation[1]*sin(-j*Beta) ;          // rotationMatrix2 has the following form:
        rotationMatrix2[0][1] = rotation[0]*rotation[1]*(1-cos(-j*Beta)) - rotation[2]*sin(-j*Beta) ;          // [0][0] [0][1] [0][2]
        rotationMatrix2[1][1] = cos(-j*Beta) + pow(rotation[1],2)*(1-cos(-j*Beta));                            // [1][0] [1][1] [1][2]
        rotationMatrix2[2][1] = rotation[1]*rotation[2]*(1-cos(-j*Beta)) + rotation[0]*sin(-j*Beta) ;          // [2][0] [2][1] [2][2]
        rotationMatrix2[0][2] = rotation[0]*rotation[2]*(1-cos(-j*Beta)) + rotation[1]*sin(-j*Beta) ;
        rotationMatrix2[1][2] = rotation[1]*rotation[2]*(1-cos(-j*Beta)) - rotation[0]*sin(-j*Beta) ;
        rotationMatrix2[2][2] = cos(-j*Beta) + pow(rotation[2],2)*(1-cos(-j*Beta));
        double molecule2[molecule.obMol.NumAtoms()][3];
        for (unsigned int i=1; i<=molecule.obMol.NumAtoms() ; i++ ) {
            molecule2[i-1][0] = rotationMatrix2[0][0]*translatedMolecule[i-1][0] + rotationMatrix2[0][1]*translatedMolecule[i-1][1] + rotationMatrix2[0][2]*translatedMolecule[i-1][2] + randomPoint[0];
            molecule2[i-1][1] = rotationMatrix2[1][0]*translatedMolecule[i-1][0] + rotationMatrix2[1][1]*translatedMolecule[i-1][1] + rotationMatrix2[1][2]*translatedMolecule[i-1][2] + randomPoint[1];
            molecule2[i-1][2] = rotationMatrix2[2][0]*translatedMolecule[i-1][0] + rotationMatrix2[2][1]*translatedMolecule[i-1][1] + rotationMatrix2[2][2]*translatedMolecule[i-1][2] + randomPoint[2];
            }
        // Update molecule.
        for (unsigned int i=1 ; i <= rotatedMolecule.obMol.NumAtoms() ; i++) {
            OpenBabel::OBAtom *atom;
            atom = rotatedMolecule.obMol.GetAtom(i);
            atom->SetVector(molecule2[i-1][0], molecule2[i-1][1], molecule2[i-1][2]);
            }
        rotatedMolecule.interactionEnergy = Energy(parameters, protein, rotatedMolecule);
        rotatedMolecule.stericClashes = StericClash(parameters, protein, rotatedMolecule);
        // Keep the translated molecule: 1) if there were steric clashes, and there are no more. 2) if the energy is lower than the previous best one and there are no clashes.
        if      ( !rotatedMolecule.stericClashes && molecule.stericClashes )                                           { molecule = rotatedMolecule; }
        else if ( !rotatedMolecule.stericClashes && (rotatedMolecule.interactionEnergy < molecule.interactionEnergy) ) { molecule = rotatedMolecule; }
        }
}

// This function finds the cartesian center of a molecule.
void CenterOfMolecule(Molecule const & molecule, double moleculeCenter[])
{
    // Get coordinates of the first atom.
    OpenBabel::OBAtom *atom;
    atom = molecule.obMol.GetAtom(1);
    double Xmin=atom->GetX(), Xmax=atom->GetX();
    double Ymin=atom->GetY(), Ymax=atom->GetY();
    double Zmin=atom->GetZ(), Zmax=atom->GetZ();

    // Compare values.
    for (unsigned int j=1; j<=molecule.obMol.NumAtoms() ; j++ ) {
        atom = molecule.obMol.GetAtom(j);
        if ( (atom->GetX()) < Xmin ) { Xmin = atom->GetX(); }
        if ( (atom->GetX()) > Xmax ) { Xmax = atom->GetX(); }
        if ( (atom->GetY()) < Ymin ) { Ymin = atom->GetY(); }
        if ( (atom->GetY()) > Ymax ) { Ymax = atom->GetY(); }
        if ( (atom->GetZ()) < Zmin ) { Zmin = atom->GetZ(); }
        if ( (atom->GetZ()) > Zmax ) { Zmax = atom->GetZ(); }
        }

    moleculeCenter[0]=(Xmin+Xmax)/2;
    moleculeCenter[1]=(Ymin+Ymax)/2;
    moleculeCenter[2]=(Zmin+Zmax)/2;
}

