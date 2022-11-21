#include "OpenGrowth.h"

double dist2(double const atom1[3], double const atom2[3]);
bool   IsThereABond(Molecule const & molecule, unsigned int const i, unsigned int const j);
float  rVDW(OpenBabel::OBAtom *atom);

// This function looks for steric clash between the ligand and the protein and inside the ligand. Depending if no/only inter/only intra/both inter and intra clashes are found, the returned value is different.
int StericClash(Parameters const & parameters, Protein const & protein, Molecule const & molecule)
{
    int clashIntra=0;
    int clashInter=0;
    for (unsigned int i=1; i<=molecule.obMol.NumAtoms() ; i++ ) {
        // Get info for first atom.
        OpenBabel::OBAtom *firstAtom;
        firstAtom = molecule.obMol.GetAtom(i);
        double firstAtomCoord[3];
        firstAtomCoord[0]=firstAtom->GetX();
        firstAtomCoord[1]=firstAtom->GetY();
        firstAtomCoord[2]=firstAtom->GetZ();

        float vdWScaleIntra2=pow(parameters.vdWScaleIntra,2);
        float vdWScaleInter2=pow(parameters.vdWScaleInter,2);

        // INTRA
        // If the distance between two atoms of the ligand (which are not from the same fragment and which are not bounded) is lower than VDW_SCALE_INTRA*(sum of the vdW radii), there is an
        // intramolecular clash. To save computational time, we compare the square of the distances.
        // If we want to only check it for the newly added fragment, add : (molecule.fragmentIndex[i-1]==molecule.numberResidues).
        for (unsigned int j=1; j<i && !clashIntra && vdWScaleIntra2; j++ ) {
            // Get info for second atom.
            OpenBabel::OBAtom *secondAtom;
            secondAtom = molecule.obMol.GetAtom(j);
            double secondAtomCoord[3];
            secondAtomCoord[0]=secondAtom->GetX();
            secondAtomCoord[1]=secondAtom->GetY();
            secondAtomCoord[2]=secondAtom->GetZ();
            float sumVDW2=pow(rVDW(firstAtom)+rVDW(secondAtom),2);
            // Check clash.
            if ( (molecule.fragmentIndex[i-1]!=molecule.fragmentIndex[j-1]) && !IsThereABond(molecule, i, j) &&
                 (dist2(firstAtomCoord, secondAtomCoord)<(vdWScaleIntra2*sumVDW2)) ) { clashIntra++ ; }
            }

        // INTER
        // If the distance between an heavy atom of the ligand and an heavy atom of the protein is lower than VDW_SCALE_INTER*(sum of the vdW radii), there is a clash. To save computational
        // time, we compare the square of the distances. We don't do it if VDW_SCALE_INTER=0, i.e. for functions where it has been disabled because there is repulsion included.
        for (unsigned int j=1; j<=protein.atom.size() && !clashInter && vdWScaleInter2; j++ ) {
            // Get info for second atom.
            OpenBabel::OBAtom *secondAtom;
            secondAtom = protein.obMol.GetAtom(protein.atom[j-1].index);
            double secondAtomCoord[3];
            secondAtomCoord[0]=secondAtom->GetX();
            secondAtomCoord[1]=secondAtom->GetY();
            secondAtomCoord[2]=secondAtom->GetZ();
            float sumVDW2=pow(rVDW(firstAtom)+rVDW(secondAtom),2);
            // Check clash.
            if ( dist2(firstAtomCoord, secondAtomCoord)<(vdWScaleInter2*sumVDW2) ) { clashInter++ ; }
            }
        }

    // Return the clash value, depending on the status for inter and intramolecular clashes.
    int clash=0;
    if      (clashIntra==0 && clashInter>0)  { clash=1; }
    else if (clashIntra>0  && clashInter==0) { clash=2; }
    else if (clashIntra>0  && clashInter>0)  { clash=3; }
    return clash;
}

// This function checks if there is a bond between atoms i and j.
bool IsThereABond(Molecule const & molecule, unsigned int const i, unsigned int const j)
{
    OpenBabel::OBAtom *firstAtom;
    firstAtom = molecule.obMol.GetAtom(i);
    OpenBabel::OBAtom *secondAtom;
    secondAtom = molecule.obMol.GetAtom(j);

    OpenBabel::OBBond *bond;
    bond=molecule.obMol.GetBond(firstAtom, secondAtom);

    if ( bond != NULL ) { return true;  }
    else                { return false; }
}

// This function returns the van der Waals radii of an atom (from http://periodictable.com/Properties/A/VanDerWaalsRadius.v.html).
float rVDW(OpenBabel::OBAtom *atom)
{
    float radii;

    // Raddi by Bondi (http://pubs.acs.org/doi/abs/10.1021/j100785a001), except B by Truhlar (http://pubs.acs.org/doi/full/10.1021/jp8111556)
    if      (atom->GetAtomicNum()==1)  { radii=1.10; }    // H
    else if (atom->GetAtomicNum()==5)  { radii=1.92; }    // B
    else if (atom->GetAtomicNum()==6)  { radii=1.70; }    // C
    else if (atom->GetAtomicNum()==7)  { radii=1.55; }    // N
    else if (atom->GetAtomicNum()==8)  { radii=1.52; }    // O
    else if (atom->GetAtomicNum()==9)  { radii=1.47; }    // F
    else if (atom->GetAtomicNum()==15) { radii=1.80; }    // P
    else if (atom->GetAtomicNum()==16) { radii=1.80; }    // S
    else if (atom->GetAtomicNum()==17) { radii=1.75; }    // Cl
    else if (atom->GetAtomicNum()==35) { radii=1.83; }    // Br
    else if (atom->GetAtomicNum()==53) { radii=1.98; }    // I
    else                               { radii=1.67; }    // Average of the previous values, just in case it is needed and even if it doesn't mean anything.

    return radii;
}

