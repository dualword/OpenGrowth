#include "OpenGrowth.h"

void   LigandTypeSMOG2001(Molecule & molecule);
float  SMOG2001(Protein const & protein, Molecule const & molecule);
void   LigandTypeSMOG2016(Molecule & molecule);
double SMOG2016(Parameters const & parameters, Protein const & protein, Molecule & molecule);

// This function computes the binding score. Depending on the choice made in the input file, the energy is computed in different ways.
double Energy(Parameters const & parameters, Protein const & protein, Molecule & molecule)
{
    double energy=0;
    if (parameters.scoringFunction == "SMOG2001") {
        LigandTypeSMOG2001(molecule);
        energy=SMOG2001(protein, molecule);
        }
    if (parameters.scoringFunction == "SMOG2016") {
        LigandTypeSMOG2016(molecule);
        energy=SMOG2016(parameters, protein, molecule);
        }
    return energy;
}

