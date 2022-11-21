#include "OpenGrowth.h"

// This function optimizes the geometry of the molecule and compute the difference of energy.
float LigandConstraints(Parameters const & parameters, Molecule molecule)
{
    // Define the forcefield
    OpenBabel::OBForceField *forceField = OpenBabel::OBForceField::FindForceField(parameters.optimizationForceField);
    if (!forceField) {
        cerr << "ERROR: Could not find forcefield " << parameters.optimizationForceField << "." << endl;
        exit (-1);
        }
    forceField->SetLogFile(&cout);
    forceField->SetLogLevel(OBFF_LOGLVL_NONE);    // NONE / LOW / MEDIUM / HIGH

    // Setup forcefield.
    if (!forceField->Setup(molecule.obMol)) {
        cerr << "ERROR: Could not setup force field " << parameters.optimizationForceField << "." << endl;
        exit (-1);
        }

    float energyBefore = forceField->Energy();
    forceField->SteepestDescent(10*parameters.optimizationSteepestDescent);
    forceField->ConjugateGradients(10*parameters.optimizationConjugateGradient);
    float energyAfter = forceField->Energy();
    return energyBefore-energyAfter;
}

