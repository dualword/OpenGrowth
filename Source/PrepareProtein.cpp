#include "OpenGrowth.h"
#include <openbabel/residue.h>

double dist2(double const atom1[3], double const atom2[3]);
void   ProteinTypeSMOG2001(string const residue, string const atomID, string & atomType);
void   ProteinTypeSMOG2016(string const residue, string const atomID, int & atomTypeNumber, float & LJr, float & LJe);

// This function prepares the protein files by storing the atom coordinates and the atom types.
void PrepareProtein(Parameters const & parameters, Protein protein[], string const snapshotNames[], string const typeOfSnapshots)
{
    cout << endl << "**********************************************************************" << endl;

    // Depending on which kind of snapshots we are using, there are not the same number of them.
    int snapshotNumber=0;
    if      (typeOfSnapshots == "rotamers")   { snapshotNumber = parameters.rotamersNumber;   }
    else if (typeOfSnapshots == "conformers") { snapshotNumber = parameters.conformersNumber; }

    for (int s=0 ; s < snapshotNumber ; s++) {
        // Open the .pdb file.
        OpenBabel::OBConversion obConversion;
        OpenBabel::OBFormat *format = obConversion.FormatFromExt(snapshotNames[s].c_str());
        obConversion.SetInFormat(format);
        obConversion.ReadFile(&protein[s].obMol, snapshotNames[s].c_str());

        // This part assign the potential atomtypes for the protein atoms. We do it only for atoms within the PROTEIN_RANGE from the center of the BINDING_SITE.
        for (unsigned int i=1; i<=protein[s].obMol.NumAtoms(); i++) {
            OpenBabel::OBAtom *atomProt;
            atomProt = protein[s].obMol.GetAtom(i);
            OpenBabel::OBResidue *obResidue;
            obResidue = atomProt->GetResidue();
            string residue=obResidue->GetName();
            string atomID=obResidue->GetAtomID(atomProt);
            residue.erase(remove(residue.begin(), residue.end(), ' '), residue.end());
            atomID.erase(remove(atomID.begin(), atomID.end(), ' '), atomID.end());
            double coordinatesAtom[3] = {atomProt->GetX(), atomProt->GetY(), atomProt->GetZ()};
            if ( dist2(coordinatesAtom, parameters.bindingSite) <= pow(parameters.proteinRange,2) ) {
                Atom tempAtom;
                string atomType="";
                int atomTypeNumber=0;
                float LJe=0.0;
                float LJr=0.0;
                if      (parameters.scoringFunction == "SMOG2001") { ProteinTypeSMOG2001(residue, atomID, atomType);                 }
                else if (parameters.scoringFunction == "SMOG2016") { ProteinTypeSMOG2016(residue, atomID, atomTypeNumber, LJr, LJe); }
                tempAtom.coordinates[0] = atomProt->GetX();
                tempAtom.coordinates[1] = atomProt->GetY();
                tempAtom.coordinates[2] = atomProt->GetZ();
                tempAtom.Type = atomType;
                tempAtom.TypeNumber = atomTypeNumber;
                tempAtom.LJr = LJr;
                tempAtom.LJe = LJe;
                tempAtom.index = i;
                protein[s].atom.push_back(tempAtom);
                }
            }

        cout << "Preparation of " << protein[s].obMol.GetTitle() << " (" << protein[s].obMol.NumAtoms() << " atoms in total, " << protein[s].atom.size() << " kept) successful" << endl;
        }

    cout << "**********************************************************************" << endl;
}

