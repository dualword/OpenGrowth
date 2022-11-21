#include "OpenGrowth.h"

float Optimization(Parameters const & parameters, Molecule & molecule, int & errorFF, int const & optParameter)
{
    float energy=0;

    // Define the forcefield
    OpenBabel::OBForceField *forceField = OpenBabel::OBForceField::FindForceField(parameters.optimizationForceField);
    if (!forceField) {
        cerr << "ERROR: Could not find forcefield." << endl;
        exit (-1);
        }
    forceField->SetLogFile(&cout);
    if (parameters.verbose>=6) { forceField->SetLogLevel(OBFF_LOGLVL_MEDIUM); } // To set the log level.
    else                       { forceField->SetLogLevel(OBFF_LOGLVL_NONE);   } // NONE / LOW / MEDIUM / HIGH

    // Set the constraints, to minimize the structure with fixing the position of atom with index 1.
    OpenBabel::OBFFConstraints constraints;
    constraints.AddAtomConstraint(1);
    
    // We pass the constraints as argument for Setup()
    if (!forceField->Setup(molecule.obMol, constraints)) {
        cerr << "ERROR: Could not setup force field." << endl;
        errorFF=1;
        }

    forceField->EnableCutOff(true);
    forceField->SetVDWCutOff(parameters.optimizationVdwCutOff);
    forceField->SetElectrostaticCutOff(parameters.optimizationElecCutOff);

    if (errorFF==0) {
        // Print the energy and unit and perform the minimization
        if (parameters.verbose>=4) { cout << "\t\t\tEnergy before: " << forceField->Energy() << " " << forceField->GetUnit() << ". " ; }
        forceField->SteepestDescent(optParameter*parameters.optimizationSteepestDescent);
        forceField->ConjugateGradients(optParameter*parameters.optimizationConjugateGradient);
        if (parameters.verbose>=4) { cout << "Energy after: " << forceField->Energy() << " " << forceField->GetUnit() << "." << endl; }

        // Save the output coordinates.
        forceField->GetCoordinates(molecule.obMol);

        energy=forceField->Energy();
        }

    return energy;
}

