#include "OpenGrowth.h"

double Energy(Parameters const & parameters, Protein const & protein, Molecule & molecule);
int    StericClash(Parameters const & parameters, Protein const & protein, Molecule const & molecule);
double dist2(double const atom1[3], double const atom2[3]);

// This function optimizes the geometry of a ligand in the active site of the protein.
float OptimizationGeom(Parameters const & parameters, Protein const & protein, Molecule & molecule, int & errorFF, int const & optParameter)
{
    // Define a temporary molecule and OBMol for the complex.
    Molecule tempMolecule;
    tempMolecule=molecule;
    OpenBabel::OBMol obMolComplex;
    obMolComplex = tempMolecule.obMol;
    obMolComplex += protein.obMol;

    // Create tables for ligand + pocket
    OpenBabel::OBBitVec ligand;
    OpenBabel::OBBitVec pocket;
    ligand.SetRangeOn(1, tempMolecule.obMol.NumAtoms());
    pocket.SetRangeOn(tempMolecule.obMol.NumAtoms()+1, tempMolecule.obMol.NumAtoms()+protein.obMol.NumAtoms());

    // Define the forcefield
    OpenBabel::OBForceField *forceField = OpenBabel::OBForceField::FindForceField(parameters.optimizationForceField);
    if (!forceField) {
        cerr << "ERROR: Could not find forcefield " << parameters.optimizationForceField << "." << endl;
        errorFF=1;
        }
    forceField->SetLogFile(&cout);
    if (parameters.verbose>=6) { forceField->SetLogLevel(OBFF_LOGLVL_MEDIUM); } // To set the log level.
    else                       { forceField->SetLogLevel(OBFF_LOGLVL_NONE);   } // NONE / LOW / MEDIUM / HIGH

    // Fix the atom which is the closest from the binding site. To save computational time, we compare the square of distances.
    OpenBabel::OBFFConstraints constraints;
    constraints.Clear();
    int atomToFix=1;
    double minimalDistance2=0;
    for (unsigned int i=1; i<=tempMolecule.obMol.NumAtoms() ; i++ ) {
        OpenBabel::OBAtom *atom;
        atom = tempMolecule.obMol.GetAtom(i);
        double coordinatesAtom[3]={atom->GetX(), atom->GetY(), atom->GetZ()};
        if (i==1) { minimalDistance2=dist2(coordinatesAtom, parameters.bindingSite); }
        if (dist2(coordinatesAtom, parameters.bindingSite)<=minimalDistance2) { minimalDistance2=dist2(coordinatesAtom, parameters.bindingSite); atomToFix=i; }
        }
    if (parameters.verbose>=4) { cout << "\t\t\tWe will fix atom " << atomToFix << " which is at " << sqrt(minimalDistance2) << " Angstrom of the binding site." << endl; }
    constraints.AddAtomConstraint(atomToFix);

    // Fix protein atoms.
    for (unsigned int i=tempMolecule.obMol.NumAtoms()+1; i<=tempMolecule.obMol.NumAtoms()+protein.obMol.NumAtoms() ; i++ ) {
        if ( pocket.BitIsOn(i) ) { constraints.AddAtomConstraint(i); }
        }

    // Specify the interacting groups. The pocket atoms are fixed, so there is no need to calculate intra- and inter-molecular interactions for the binding pocket.
    forceField->ClearGroups();
    forceField->AddIntraGroup(ligand);            // bonded interactions in the ligand
    forceField->AddInterGroup(ligand);            // non-bonded between ligand-ligand atoms
    forceField->AddInterGroups(ligand, pocket);   // non-bonded between ligand and pocket atoms

    // We pass the constraints as argument for Setup().
    if (!forceField->Setup(obMolComplex, constraints)) {
        cerr << "ERROR: Could not setup force field." << endl;
        errorFF=1;
        }

    forceField->EnableCutOff(true);
    forceField->SetVDWCutOff(parameters.optimizationVdwCutOff);
    forceField->SetElectrostaticCutOff(parameters.optimizationElecCutOff);

    // Print the energy and perform the minimization.
    if (parameters.verbose>=4) { cout << "\t\t\tBefore geometry optimization: " << forceField->Energy() << " " << forceField->GetUnit() << ". " ; }
    forceField->SteepestDescent(optParameter*parameters.optimizationSteepestDescent);
    forceField->ConjugateGradients(optParameter*parameters.optimizationConjugateGradient);
    if (parameters.verbose>=4) { cout << "After: " << forceField->Energy() << " " << forceField->GetUnit() << "." << endl; }

    // Write the output.
    forceField->GetCoordinates(obMolComplex);
    for (unsigned int i=1; i<=tempMolecule.obMol.NumAtoms() ; i++ ) {
        OpenBabel::OBAtom *atomComplex;
        atomComplex = obMolComplex.GetAtom(i);
        OpenBabel::OBAtom *atomTempMolecule;
        atomTempMolecule = tempMolecule.obMol.GetAtom(i);
        atomTempMolecule->SetVector(atomComplex->GetX(), atomComplex->GetY(), atomComplex->GetZ());
        }

    tempMolecule.interactionEnergy = Energy(parameters, protein, tempMolecule);
    tempMolecule.stericClashes = StericClash(parameters, protein, tempMolecule);
    // Keep the optimized molecule if there are no steric clashes (we don't consider the interaction energy criteria here).
    if (!tempMolecule.stericClashes) { molecule = tempMolecule; }
    obMolComplex.Clear();

    float energy=forceField->Energy();
    return energy;
}

