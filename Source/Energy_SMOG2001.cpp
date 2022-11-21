#include "OpenGrowth.h"

double dist2(double const atom1[3], double const atom2[3]);

// This file contains three functions: LigandTypeSMOG2001, ProteinTypeSMOG2001, SMOG2001.
// Everything related to SMOG2001 can be found here.

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function assigns the atom-types of the ligand atoms. It also stores in the dynamic array molecule.atom the coordinates and the atom-types to facilitate the calculation of energy.
void LigandTypeSMOG2001(Molecule & molecule)
{
    molecule.atom.clear();
    for (unsigned int i=1; i<=molecule.obMol.NumAtoms() ; i++ ) {
        OpenBabel::OBAtom *atom;
        atom = molecule.obMol.GetAtom(i);
        // Check neighbor atoms to see if they are different from C or H. These 6 lines are coming from typer.cpp from OpenBabel.
        int polarity = 0;
        OpenBabel::OBAtom *nbr;
        vector<OpenBabel::OBBond*>::iterator k;
        for (nbr = atom->BeginNbrAtom(k); nbr; nbr = atom->NextNbrAtom(k)) {
            if (nbr->IsNotCorH()) { polarity++; }
            }

        // Begin atom-typing
        string atomType;
        if      (atom->IsHydrogen())                                                             { atomType = "H";  }
        else if (atom->IsSulfur() || atom->IsPhosphorus())                                       { atomType = "SL"; }
        else if (atom->IsAmideNitrogen())                                                        { atomType = "NM"; }
        else if (atom->IsNitrogen() && (atom->GetValence() > atom->GetHyb()) )                   { atomType = "NC"; }  // A nitrogen is charged if its valence is higher than its hybridation
        else if (atom->IsNitrogen() && atom->IsHbondDonor())                                     { atomType = "ND"; }
        else if (atom->IsNitrogen() && atom->IsHbondAcceptor())                                  { atomType = "NA"; }
        else if (atom->IsCarboxylOxygen() || atom->IsPhosphateOxygen() || atom->IsNitroOxygen()) { atomType = "OC"; }
        else if (atom->IsOxygen() && atom->MatchesSMARTS("[O;$(O=*)]"))                          { atomType = "OB"; }
        else if (atom->IsOxygen() && atom->IsHbondAcceptor())                                    { atomType = "OA"; }
        else if (atom->IsOxygen() && atom->IsHbondDonor())                                       { atomType = "OD"; }
        else if (atom->IsCarbon() && (atom->GetHyb()==3) && polarity==0)                         { atomType = "C3"; }
        else if (atom->IsCarbon() && (atom->GetHyb()==2) && polarity==0)                         { atomType = "C2"; }
        else if (atom->IsCarbon() && (atom->GetHyb()==1) && polarity==0)                         { atomType = "C2"; }
        else if (atom->IsCarbon() && atom->MatchesSMARTS("[#6;$(C=O)]"))                         { atomType = "CC"; }
        else if (atom->IsCarbon() && atom->MatchesSMARTS("[#6;$(C(=N)(N)(N))]"))                 { atomType = "CC"; }
        else if (atom->IsCarbon() && polarity!=0)                                                { atomType = "CP"; }
        else                                                                                     { atomType = "H";  }

        // Create a new Atom in the array.
        Atom tempAtom;
        tempAtom.coordinates[0] = atom->GetX();
        tempAtom.coordinates[1] = atom->GetY();
        tempAtom.coordinates[2] = atom->GetZ();
        tempAtom.Type = atomType;
        tempAtom.index = i;
        molecule.atom.push_back(tempAtom);
        }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function assigns the protein atom types for the SMOG2001 scoring function.
void ProteinTypeSMOG2001(string const residue, string const atomID, string & atomType)
{
    // ALA : N, CA, C, O, CB, H
    if      ((residue == "ALA") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "ALA") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "ALA") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "ALA") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "ALA") && (atomID == "CB"))  { atomType = "C3"; }
    // ARG : N, CA, C, O, CB, CG, CD, NE, CZ, NH1, NH2, H, HE, HH11, HH12, HH21, HH22
    else if ((residue == "ARG") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "ARG") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "ARG") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "ARG") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "ARG") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "ARG") && (atomID == "CG"))  { atomType = "C3"; }
    else if ((residue == "ARG") && (atomID == "CD"))  { atomType = "CP"; }
    else if ((residue == "ARG") && (atomID == "NE"))  { atomType = "NC"; }
    else if ((residue == "ARG") && (atomID == "CZ"))  { atomType = "CC"; }
    else if ((residue == "ARG") && (atomID == "NH1")) { atomType = "NC"; }
    else if ((residue == "ARG") && (atomID == "NH2")) { atomType = "NC"; }
    // ARN = ARG-H
    else if ((residue == "ARN") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "ARN") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "ARN") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "ARN") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "ARN") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "ARN") && (atomID == "CG"))  { atomType = "C3"; }
    else if ((residue == "ARN") && (atomID == "CD"))  { atomType = "CP"; }
    else if ((residue == "ARN") && (atomID == "NE"))  { atomType = "NC"; }
    else if ((residue == "ARN") && (atomID == "CZ"))  { atomType = "CC"; }
    else if ((residue == "ARN") && (atomID == "NH1")) { atomType = "NC"; }
    else if ((residue == "ARN") && (atomID == "NH2")) { atomType = "NC"; }
    // ASN : N, CA, C, O, CB, CG, OD1, ND2, H, HD21, HD22
    else if ((residue == "ASN") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "ASN") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "ASN") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "ASN") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "ASN") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "ASN") && (atomID == "CG"))  { atomType = "CC"; }
    else if ((residue == "ASN") && (atomID == "OD1")) { atomType = "OB"; }
    else if ((residue == "ASN") && (atomID == "ND2")) { atomType = "ND"; }
    // ASP : N, CA, C, O, CB, CG, OD1, OD2, H
    else if ((residue == "ASP") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "ASP") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "ASP") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "ASP") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "ASP") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "ASP") && (atomID == "CG"))  { atomType = "CP"; }
    else if ((residue == "ASP") && (atomID == "OD1")) { atomType = "OC"; }
    else if ((residue == "ASP") && (atomID == "OD2")) { atomType = "OC"; }
    // ASH = ASP+H
    else if ((residue == "ASH") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "ASH") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "ASH") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "ASH") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "ASH") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "ASH") && (atomID == "CG"))  { atomType = "CP"; }
    else if ((residue == "ASH") && (atomID == "OD1")) { atomType = "OC"; }
    else if ((residue == "ASH") && (atomID == "OD2")) { atomType = "OC"; }
    // CYS : N, CA, C, O, CB, SG, H
    else if ((residue == "CYS") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "CYS") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "CYS") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "CYS") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "CYS") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "CYS") && (atomID == "SG"))  { atomType = "SP"; }
    // CYX = CYS for S-S bridge
    else if ((residue == "CYX") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "CYX") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "CYX") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "CYX") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "CYX") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "CYX") && (atomID == "SG"))  { atomType = "SP"; }
    // GLU : N, CA, C, O, CB, CG, CD, OE1, OE2, H
    else if ((residue == "GLU") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "GLU") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "GLU") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "GLU") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "GLU") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "GLU") && (atomID == "CG"))  { atomType = "C3"; }
    else if ((residue == "GLU") && (atomID == "CD"))  { atomType = "CP"; }
    else if ((residue == "GLU") && (atomID == "OE1")) { atomType = "OC"; }
    else if ((residue == "GLU") && (atomID == "OE2")) { atomType = "OC"; }
    // GLH = GLU+H
    else if ((residue == "GLH") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "GLH") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "GLH") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "GLH") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "GLH") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "GLH") && (atomID == "CG"))  { atomType = "C3"; }
    else if ((residue == "GLH") && (atomID == "CD"))  { atomType = "CP"; }
    else if ((residue == "GLH") && (atomID == "OE1")) { atomType = "OC"; }
    else if ((residue == "GLH") && (atomID == "OE2")) { atomType = "OC"; }
    // GLN : N, CA, C, O, CB, CG, CD, OE1, NE2, H, HE21, HE22
    else if ((residue == "GLN") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "GLN") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "GLN") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "GLN") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "GLN") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "GLN") && (atomID == "CG"))  { atomType = "C3"; }
    else if ((residue == "GLN") && (atomID == "CD"))  { atomType = "CC"; }
    else if ((residue == "GLN") && (atomID == "OE1")) { atomType = "OB"; }
    else if ((residue == "GLN") && (atomID == "NE2")) { atomType = "ND"; }
    // GLY : N, CA, C, O, H
    else if ((residue == "GLY") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "GLY") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "GLY") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "GLY") && (atomID == "O"))   { atomType = "OB"; }
    // HIS : N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2, H, HD1, HE2
    else if ((residue == "HIS") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "HIS") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "HIS") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "HIS") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "HIS") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "HIS") && (atomID == "CG"))  { atomType = "CP"; }
    else if ((residue == "HIS") && (atomID == "ND1")) { atomType = "NC"; }
    else if ((residue == "HIS") && (atomID == "CD2")) { atomType = "CP"; }
    else if ((residue == "HIS") && (atomID == "CE1")) { atomType = "CP"; }
    else if ((residue == "HIS") && (atomID == "NE2")) { atomType = "NC"; }
    // HID = HIS protonated on delta
    else if ((residue == "HID") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "HID") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "HID") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "HID") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "HID") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "HID") && (atomID == "CG"))  { atomType = "CP"; }
    else if ((residue == "HID") && (atomID == "ND1")) { atomType = "NC"; }
    else if ((residue == "HID") && (atomID == "CD2")) { atomType = "CP"; }
    else if ((residue == "HID") && (atomID == "CE1")) { atomType = "CP"; }
    else if ((residue == "HID") && (atomID == "NE2")) { atomType = "NC"; }
    // HIE = HIS protonated on epsilon
    else if ((residue == "HIE") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "HIE") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "HIE") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "HIE") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "HIE") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "HIE") && (atomID == "CG"))  { atomType = "CP"; }
    else if ((residue == "HIE") && (atomID == "ND1")) { atomType = "NC"; }
    else if ((residue == "HIE") && (atomID == "CD2")) { atomType = "CP"; }
    else if ((residue == "HIE") && (atomID == "CE1")) { atomType = "CP"; }
    else if ((residue == "HIE") && (atomID == "NE2")) { atomType = "NC"; }
    // HIP = HIS doubly protonated
    else if ((residue == "HIP") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "HIP") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "HIP") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "HIP") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "HIP") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "HIP") && (atomID == "CG"))  { atomType = "CP"; }
    else if ((residue == "HIP") && (atomID == "ND1")) { atomType = "NC"; }
    else if ((residue == "HIP") && (atomID == "CD2")) { atomType = "CP"; }
    else if ((residue == "HIP") && (atomID == "CE1")) { atomType = "CP"; }
    else if ((residue == "HIP") && (atomID == "NE2")) { atomType = "NC"; }
    // ILE : N, CA, C, O, CB, CG1, CG2, CD1, H
    else if ((residue == "ILE") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "ILE") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "ILE") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "ILE") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "ILE") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "ILE") && (atomID == "CG1")) { atomType = "C3"; }
    else if ((residue == "ILE") && (atomID == "CG2")) { atomType = "C3"; }
    else if ((residue == "ILE") && (atomID == "CD1")) { atomType = "C3"; }
    else if ((residue == "ILE") && (atomID == "CD"))  { atomType = "C3"; }
    // LEU : N, CA, C, O, CB, CG, CD1, CD2, H
    else if ((residue == "LEU") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "LEU") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "LEU") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "LEU") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "LEU") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "LEU") && (atomID == "CG"))  { atomType = "C3"; }
    else if ((residue == "LEU") && (atomID == "CD1")) { atomType = "C3"; }
    else if ((residue == "LEU") && (atomID == "CD2")) { atomType = "C3"; }
    // LYS : N, CA, C, O, CB, CG, CD, CE, NZ, H, HZ1, HZ2, HZ3
    else if ((residue == "LYS") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "LYS") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "LYS") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "LYS") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "LYS") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "LYS") && (atomID == "CG"))  { atomType = "C3"; }
    else if ((residue == "LYS") && (atomID == "CD"))  { atomType = "C3"; }
    else if ((residue == "LYS") && (atomID == "CE"))  { atomType = "CP"; }
    else if ((residue == "LYS") && (atomID == "NZ"))  { atomType = "NC"; }
    // LYN = LYS-H
    else if ((residue == "LYN") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "LYN") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "LYN") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "LYN") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "LYN") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "LYN") && (atomID == "CG"))  { atomType = "C3"; }
    else if ((residue == "LYN") && (atomID == "CD"))  { atomType = "C3"; }
    else if ((residue == "LYN") && (atomID == "CE"))  { atomType = "CP"; }
    else if ((residue == "LYN") && (atomID == "NZ"))  { atomType = "NC"; }
    // MET : N, CA, C, O, CB, CG, SD, CE, H
    else if ((residue == "MET") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "MET") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "MET") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "MET") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "MET") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "MET") && (atomID == "CG"))  { atomType = "C3"; }
    else if ((residue == "MET") && (atomID == "SD"))  { atomType = "SP"; }
    else if ((residue == "MET") && (atomID == "CE"))  { atomType = "C3"; }
    // PHE : N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, H
    else if ((residue == "PHE") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "PHE") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "PHE") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "PHE") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "PHE") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "PHE") && (atomID == "CG"))  { atomType = "C2"; }
    else if ((residue == "PHE") && (atomID == "CD1")) { atomType = "C2"; }
    else if ((residue == "PHE") && (atomID == "CD2")) { atomType = "C2"; }
    else if ((residue == "PHE") && (atomID == "CE1")) { atomType = "C2"; }
    else if ((residue == "PHE") && (atomID == "CE2")) { atomType = "C2"; }
    else if ((residue == "PHE") && (atomID == "CZ"))  { atomType = "C2"; }
    // PRO : N, CA, C, O, CB, CG, CD, H2, H3
    else if ((residue == "PRO") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "PRO") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "PRO") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "PRO") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "PRO") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "PRO") && (atomID == "CG"))  { atomType = "C3"; }
    else if ((residue == "PRO") && (atomID == "CD"))  { atomType = "CP"; }
    // SER : N, CA, C, O, CB, OG
    else if ((residue == "SER") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "SER") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "SER") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "SER") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "SER") && (atomID == "CB"))  { atomType = "CP"; }
    else if ((residue == "SER") && (atomID == "OG"))  { atomType = "OD"; }
    // THR : N, CA, C, O, CB, OG1, CG2, H, HG1
    else if ((residue == "THR") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "THR") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "THR") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "THR") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "THR") && (atomID == "CB"))  { atomType = "CP"; }
    else if ((residue == "THR") && (atomID == "OG1")) { atomType = "OD"; }
    else if ((residue == "THR") && (atomID == "CG2")) { atomType = "C3"; }
    // TRP : N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2, H, HE1
    else if ((residue == "TRP") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "TRP") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "TRP") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "TRP") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "TRP") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "TRP") && (atomID == "CG"))  { atomType = "C2"; }
    else if ((residue == "TRP") && (atomID == "CD1")) { atomType = "CP"; }
    else if ((residue == "TRP") && (atomID == "CD2")) { atomType = "C2"; }
    else if ((residue == "TRP") && (atomID == "NE1")) { atomType = "ND"; }
    else if ((residue == "TRP") && (atomID == "CE2")) { atomType = "CP"; }
    else if ((residue == "TRP") && (atomID == "CE3")) { atomType = "C2"; }
    else if ((residue == "TRP") && (atomID == "CZ2")) { atomType = "C2"; }
    else if ((residue == "TRP") && (atomID == "CZ3")) { atomType = "C2"; }
    else if ((residue == "TRP") && (atomID == "CH2")) { atomType = "C2"; }
    // TYR : N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH, H, HH
    else if ((residue == "TYR") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "TYR") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "TYR") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "TYR") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "TYR") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "TYR") && (atomID == "CG"))  { atomType = "C2"; }
    else if ((residue == "TYR") && (atomID == "CD1")) { atomType = "C2"; }
    else if ((residue == "TYR") && (atomID == "CD2")) { atomType = "C2"; }
    else if ((residue == "TYR") && (atomID == "CE1")) { atomType = "C2"; }
    else if ((residue == "TYR") && (atomID == "CE2")) { atomType = "C2"; }
    else if ((residue == "TYR") && (atomID == "CZ"))  { atomType = "CP"; }
    else if ((residue == "TYR") && (atomID == "OH"))  { atomType = "OD"; }
    // VAL : N, CA, C, O, CB, CG1, CG2, H
    else if ((residue == "VAL") && (atomID == "N"))   { atomType = "NM"; }
    else if ((residue == "VAL") && (atomID == "CA"))  { atomType = "CA"; }
    else if ((residue == "VAL") && (atomID == "C"))   { atomType = "CC"; }
    else if ((residue == "VAL") && (atomID == "O"))   { atomType = "OB"; }
    else if ((residue == "VAL") && (atomID == "CB"))  { atomType = "C3"; }
    else if ((residue == "VAL") && (atomID == "CG1")) { atomType = "C3"; }
    else if ((residue == "VAL") && (atomID == "CG2")) { atomType = "C3"; }
    // Metals: CA, ZN, HG, MG, CD, NI, MN, NA
    else if ((residue == "CA")  && (atomID == "CA"))  { atomType = "M";  }
    else if ((residue == "ZN")  && (atomID == "ZN"))  { atomType = "M";  }
    else if ((residue == "HG")  && (atomID == "HG"))  { atomType = "M";  }
    else if ((residue == "MG")  && (atomID == "MG"))  { atomType = "M";  }
    else if ((residue == "CD")  && (atomID == "CD"))  { atomType = "M";  }
    else if ((residue == "NI")  && (atomID == "NI"))  { atomType = "M";  }
    else if ((residue == "MN")  && (atomID == "MN"))  { atomType = "M";  }
    else if ((residue == "NA")  && (atomID == "NA"))  { atomType = "M";  }
    // H
    else                                              { atomType = "H";  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function computes the interaction energy using SMoG2001. The parameter file was originally called params_3.5_4.5_0.9.
// If the distance between an atom of the ligand and an atom of the protein is lower than 3.5 Angstrom, the interaction energy
// has a given value. If it is between 3.5 and 4.5 Angstrom, it has a different value. Above 4.5 Angstrom, they are not in contact
// and the energy is 0. To save computational time, we compare the square of the distance with 3.5^2=12.25 and 4.5^2=20.25.
// There are 13 atom types for molecule atoms, and 13 for protein atoms.
// Ref.: http://pubs.acs.org/doi/abs/10.1021/jm0105833
float SMOG2001(Protein const & protein, Molecule const & molecule)
{
    float energy=0.0;

    //   L I G A N D
    // P
    // R
    // O
    // T
    // E
    // I
    // N

    //From 0 to 3.5:     C3      C2    CP     CC     OC     OB     OD     OA     NC     NM     ND     NA     SL 
    float shortC3[13]={ -0.80, -0.55,  0.40,  0.48, -0.33, -0.53, -0.29, -0.16,  0.52,  1.29, -0.28, -0.55,  0.57 };
    float shortC2[13]={ -0.92, -0.88, -0.14,  0.21,  0.13, -0.48, -0.27,  0.33,  0.63,  0.63, -0.14, -0.12,  1.44 };
    float shortCP[13]={ -0.24,  0.74,  0.60,  0.23, -1.16, -0.46, -0.57, -0.22,  1.02,  1.06,  0.67,  0.19,  0.65 };
    float shortCC[13]={ -0.69,  0.96,  0.33,  0.78, -0.24, -0.30, -1.10, -0.34, -0.58,  0.78, -0.17, -0.13,  1.74 };
    float shortCA[13]={  0.59,  0.49,  1.31,  0.86, -0.79, -1.07, -0.56, -0.75,  0.07,  2.87,  0.10, -0.12,  0.40 };
    float shortOC[13]={ -1.02,  0.44, -0.69,  0.22,  0.57,  0.21, -1.61,  0.32, -0.90, -0.53, -0.89, -0.39,  0.79 };
    float shortOB[13]={ -0.84, -0.33, -0.77, -0.12,  0.49, -0.19, -0.32,  0.31, -1.44, -1.22, -0.91, -0.64,  0.82 };
    float shortOD[13]={ -0.59,  0.13, -0.48, -0.75, -1.19, -0.36, -0.66, -0.32, -0.84, -0.48, -0.25, -0.18, -0.70 };
    float shortNC[13]={ -0.01,  0.12,  0.15, -0.56, -1.57, -0.38, -0.87, -0.70,  1.40,  0.86,  0.48, -0.17, -0.85 };
    float shortNM[13]={  0.64,  0.79,  1.02,  0.38, -1.17, -1.05, -0.58, -0.38,  0.46,  1.69, -0.06, -0.52,  0.09 };
    float shortND[13]={ -0.05,  0.62, -0.62, -0.59, -0.91, -1.10, -0.68, -0.86,  0.41, -0.39, -0.16, -0.15,  0.14 };
    float shortSP[13]={ -1.25, -0.68, -0.40, -0.08,  1.79,  1.06,  0.01,  0.63,  1.15,  0.10,  0.87, -0.40, -0.61 };
    float shortM[13]={  -1.65,  2.69,  1.37,  0.66, -0.22,  1.80,  0.86, -0.37,  0.91,  1.25,  2.01,  2.87, -1.06 };
    //From 3.5 to 4.5:   C3      C2    CP     CC     OC     OB     OD     OA     NC     NM     ND     NA     SL 
    float longC3[13]={  -0.45, -0.61, -0.08,  0.13, -0.05, -0.10, -0.13,  0.06, -0.05,  0.47, -0.59, -0.57,  0.34 };
    float longC2[13]={  -0.66, -0.67, -0.53, -0.08,  0.31, -0.29, -0.02,  0.42,  0.15, -0.03, -0.50, -0.20,  0.96 };
    float longCP[13]={  -0.10,  0.26, -0.06, -0.20, -0.61, -0.00, -0.23, -0.05,  0.18,  0.11, -0.03,  0.07, -0.57 };
    float longCC[13]={  -0.03,  0.25,  0.05,  0.31, -0.59, -0.24, -0.40,  0.05, -0.10,  0.04, -0.42,  0.11, -0.06 };
    float longCA[13]={   0.07,  0.08,  0.11,  0.12, -0.58, -0.13, -0.19, -0.17,  0.15,  0.47, -0.43, -0.08, -0.43 };
    float longOC[13]={  -0.36,  0.57, -0.57,  0.10,  0.50,  0.10, -0.36,  0.11,  0.42,  0.28, -0.11,  0.18,  0.69 };
    float longOB[13]={  -0.36, -0.23, -0.30, -0.42, -0.12, -0.18, -0.07, -0.04, -0.01,  0.34, -0.46, -0.25,  0.67 };
    float longOD[13]={  -0.27, -0.05, -0.23, -0.40, -0.33, -0.06,  0.11, -0.17,  0.09,  0.03,  0.06,  0.05, -0.88 };
    float longNC[13]={   0.07,  0.17, -0.20, -0.38, -0.83,  0.48, -0.08, -0.06,  1.05,  0.45,  0.25,  0.02, -1.23 };
    float longNM[13]={   0.13,  0.22,  0.19,  0.05, -0.30,  0.37, -0.08, -0.32, -0.18,  0.52, -0.25,  0.17, -0.84 };
    float longND[13]={  -0.20,  0.09, -0.39, -0.52, -0.26,  0.08, -0.18,  0.29, -0.04,  0.14,  0.00,  0.33, -0.32 };
    float longSP[13]={  -0.06, -0.60, -0.36,  0.41,  1.84,  0.06, -0.18,  0.58, -0.09,  0.26, -0.51, -0.50,  2.04 };
    float longM[13]={    0.60,  1.59,  0.65,  0.68, -0.49,  0.99,  0.29, -0.64,  0.18,  0.86,  0.89,  2.04, -1.03 };

    // We do a double loop in the protein atom, then in the ligand atom. Depending on the distance, one of the value is read.
    for (unsigned int i=0 ; i < protein.atom.size() ; i++) {
        for (unsigned int j=0 ; j < molecule.atom.size() ; j++) {
            if (dist2(protein.atom[i].coordinates, molecule.atom[j].coordinates) <= 12.25) {
                if (protein.atom[i].Type == "C3")      {
                    if      (molecule.atom[j].Type == "C3") { energy += shortC3[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortC3[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortC3[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortC3[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortC3[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortC3[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortC3[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortC3[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortC3[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortC3[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortC3[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortC3[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortC3[12]; }
                    }
                else if (protein.atom[i].Type == "C2") {
                    if      (molecule.atom[j].Type == "C3") { energy += shortC2[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortC2[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortC2[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortC2[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortC2[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortC2[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortC2[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortC2[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortC2[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortC2[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortC2[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortC2[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortC2[12]; }
                    }
                else if (protein.atom[i].Type == "CP") {
                    if      (molecule.atom[j].Type == "C3") { energy += shortCP[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortCP[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortCP[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortCP[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortCP[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortCP[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortCP[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortCP[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortCP[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortCP[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortCP[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortCP[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortCP[12]; }
                    }
                else if (protein.atom[i].Type == "CC") {
                    if      (molecule.atom[j].Type == "C3") { energy += shortCC[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortCC[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortCC[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortCC[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortCC[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortCC[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortCC[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortCC[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortCC[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortCC[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortCC[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortCC[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortCC[12]; }
                    }
                else if (protein.atom[i].Type == "CA") {
                    if      (molecule.atom[j].Type == "C3") { energy += shortCA[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortCA[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortCA[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortCA[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortCA[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortCA[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortCA[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortCA[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortCA[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortCA[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortCA[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortCA[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortCA[12]; }
                    }
                else if (protein.atom[i].Type == "OC") {
                    if      (molecule.atom[j].Type == "C3") { energy += shortOC[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortOC[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortOC[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortOC[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortOC[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortOC[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortOC[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortOC[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortOC[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortOC[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortOC[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortOC[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortOC[12]; }
                    }
                else if (protein.atom[i].Type == "OB") {
                    if      (molecule.atom[j].Type == "C3") { energy += shortOB[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortOB[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortOB[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortOB[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortOB[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortOB[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortOB[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortOB[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortOB[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortOB[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortOB[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortOB[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortOB[12]; }
                    }
                else if (protein.atom[i].Type == "OD") {
                    if      (molecule.atom[j].Type == "C3") { energy += shortOD[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortOD[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortOD[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortOD[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortOD[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortOD[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortOD[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortOD[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortOD[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortOD[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortOD[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortOD[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortOD[12]; }
                    }
                else if (protein.atom[i].Type == "NC") {
                    if      (molecule.atom[j].Type == "C3") { energy += shortNC[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortNC[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortNC[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortNC[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortNC[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortNC[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortNC[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortNC[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortNC[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortNC[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortNC[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortNC[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortNC[12]; }
                    }
                else if (protein.atom[i].Type == "NM") {
                    if      (molecule.atom[j].Type == "C3") { energy += shortNM[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortNM[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortNM[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortNM[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortNM[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortNM[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortNM[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortNM[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortNM[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortNM[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortNM[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortNM[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortNM[12]; }
                    }
                else if (protein.atom[i].Type == "ND") {
                    if      (molecule.atom[j].Type == "C3") { energy += shortND[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortND[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortND[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortND[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortND[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortND[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortND[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortND[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortND[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortND[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortND[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortND[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortND[12]; }
                    }
                else if (protein.atom[i].Type == "SP") {
                    if      (molecule.atom[j].Type == "C3") { energy += shortSP[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortSP[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortSP[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortSP[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortSP[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortSP[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortSP[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortSP[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortSP[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortSP[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortSP[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortSP[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortSP[12]; }
                    }
                else if (protein.atom[i].Type == "M")  {
                    if      (molecule.atom[j].Type == "C3") { energy += shortM[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += shortM[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += shortM[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += shortM[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += shortM[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += shortM[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += shortM[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += shortM[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += shortM[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += shortM[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += shortM[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += shortM[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += shortM[12]; }
                    }
                }
            else if (dist2(protein.atom[i].coordinates, molecule.atom[j].coordinates) <= 20.25) {
                if (protein.atom[i].Type == "C3")      {
                    if      (molecule.atom[j].Type == "C3") { energy += longC3[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longC3[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longC3[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longC3[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longC3[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longC3[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longC3[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longC3[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longC3[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longC3[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longC3[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longC3[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longC3[12]; }
                    }
                else if (protein.atom[i].Type == "C2") {
                    if      (molecule.atom[j].Type == "C3") { energy += longC2[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longC2[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longC2[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longC2[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longC2[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longC2[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longC2[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longC2[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longC2[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longC2[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longC2[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longC2[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longC2[12]; }
                    }
                else if (protein.atom[i].Type == "CP") {
                    if      (molecule.atom[j].Type == "C3") { energy += longCP[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longCP[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longCP[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longCP[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longCP[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longCP[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longCP[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longCP[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longCP[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longCP[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longCP[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longCP[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longCP[12]; }
                    }
                else if (protein.atom[i].Type == "CC") {
                    if      (molecule.atom[j].Type == "C3") { energy += longCC[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longCC[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longCC[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longCC[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longCC[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longCC[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longCC[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longCC[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longCC[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longCC[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longCC[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longCC[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longCC[12]; }
                    }
                else if (protein.atom[i].Type == "CA") {
                    if      (molecule.atom[j].Type == "C3") { energy += longCA[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longCA[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longCA[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longCA[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longCA[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longCA[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longCA[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longCA[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longCA[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longCA[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longCA[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longCA[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longCA[12]; }
                    }
                else if (protein.atom[i].Type == "OC") {
                    if      (molecule.atom[j].Type == "C3") { energy += longOC[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longOC[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longOC[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longOC[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longOC[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longOC[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longOC[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longOC[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longOC[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longOC[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longOC[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longOC[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longOC[12]; }
                    }
                else if (protein.atom[i].Type == "OB") {
                    if      (molecule.atom[j].Type == "C3") { energy += longOB[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longOB[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longOB[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longOB[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longOB[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longOB[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longOB[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longOB[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longOB[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longOB[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longOB[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longOB[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longOB[12]; }
                    }
                else if (protein.atom[i].Type == "OD") {
                    if      (molecule.atom[j].Type == "C3") { energy += longOD[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longOD[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longOD[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longOD[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longOD[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longOD[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longOD[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longOD[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longOD[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longOD[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longOD[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longOD[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longOD[12]; }
                    }
                else if (protein.atom[i].Type == "NC") {
                    if      (molecule.atom[j].Type == "C3") { energy += longNC[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longNC[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longNC[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longNC[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longNC[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longNC[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longNC[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longNC[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longNC[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longNC[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longNC[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longNC[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longNC[12]; }
                    }
                else if (protein.atom[i].Type == "NM") {
                    if      (molecule.atom[j].Type == "C3") { energy += longNM[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longNM[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longNM[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longNM[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longNM[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longNM[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longNM[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longNM[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longNM[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longNM[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longNM[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longNM[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longNM[12]; }
                    }
                else if (protein.atom[i].Type == "ND") {
                    if      (molecule.atom[j].Type == "C3") { energy += longND[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longND[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longND[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longND[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longND[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longND[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longND[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longND[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longND[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longND[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longND[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longND[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longND[12]; }
                    }
                else if (protein.atom[i].Type == "SP") {
                    if      (molecule.atom[j].Type == "C3") { energy += longSP[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longSP[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longSP[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longSP[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longSP[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longSP[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longSP[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longSP[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longSP[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longSP[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longSP[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longSP[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longSP[12]; }
                    }
                else if (protein.atom[i].Type == "M") {
                    if      (molecule.atom[j].Type == "C3") { energy += longM[0];  }
                    else if (molecule.atom[j].Type == "C2") { energy += longM[1];  }
                    else if (molecule.atom[j].Type == "CP") { energy += longM[2];  }
                    else if (molecule.atom[j].Type == "CC") { energy += longM[3];  }
                    else if (molecule.atom[j].Type == "OC") { energy += longM[4];  }
                    else if (molecule.atom[j].Type == "OB") { energy += longM[5];  }
                    else if (molecule.atom[j].Type == "OD") { energy += longM[6];  }
                    else if (molecule.atom[j].Type == "OA") { energy += longM[7];  }
                    else if (molecule.atom[j].Type == "NC") { energy += longM[8];  }
                    else if (molecule.atom[j].Type == "NM") { energy += longM[9];  }
                    else if (molecule.atom[j].Type == "ND") { energy += longM[10]; }
                    else if (molecule.atom[j].Type == "NA") { energy += longM[11]; }
                    else if (molecule.atom[j].Type == "SL") { energy += longM[12]; }
                    }
                }
            }
        }

    return energy;
}

