#include "OpenGrowth.h"
#include <iomanip>            // for setprecision

string atomName(OpenBabel::OBAtom *atom);

// This function is mainly used for debugging. It displays information for a given molecule.
void CheckGrowth(Molecule const molecule)
{
    cout << fixed << setprecision(5) ;
    cout << "\t\t**********" << endl;
    cout << "\t\tGrowingSites:                           "; for (unsigned int j=0; j<molecule.growingSites.size() ; j++ ) { cout << molecule.growingSites[j] << "  " ; } cout << endl;
    cout << "\t\tSize (all numbers should be the same):  " << molecule.obMol.NumAtoms() << "  " << molecule.growth.size() << "  " << molecule.fragmentIndex.size() << "  " << molecule.fragmentNeighbours.size() << endl;
    cout << "\t\tAtom	X            Y            Z     	index    growth   fragmentIndex   fragmentNeighbours" << endl;
    for (unsigned int j=1; j<=molecule.obMol.NumAtoms() ; j++ ) {
        OpenBabel::OBAtom *atom;
        atom = molecule.obMol.GetAtom(j);
        cout << "\t\t" << atomName(atom) << "\t" << atom->GetX() << "     " << atom->GetY() << "     " << atom->GetZ() << "\t" << j << "        " << molecule.growth[j-1] << "        " << molecule.fragmentIndex[j-1] << "               " << molecule.fragmentNeighbours[j-1] << endl;
        }
    cout << "\t\t**********" << endl;
}

string atomName(OpenBabel::OBAtom *atom)
{
    string name;

    if      (atom->GetAtomicNum()==1)  { name="H";  }
    else if (atom->GetAtomicNum()==5)  { name="B";  }
    else if (atom->GetAtomicNum()==6)  { name="C";  }
    else if (atom->GetAtomicNum()==7)  { name="N";  }
    else if (atom->GetAtomicNum()==8)  { name="O";  }
    else if (atom->GetAtomicNum()==9)  { name="F";  }
    else if (atom->GetAtomicNum()==15) { name="P";  }
    else if (atom->GetAtomicNum()==16) { name="S";  }
    else if (atom->GetAtomicNum()==17) { name="Cl"; }
    else if (atom->GetAtomicNum()==35) { name="Br"; }
    else if (atom->GetAtomicNum()==53) { name="I";  }
    else                               { name="X";  }

    return name;
}

