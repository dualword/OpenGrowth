#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/atom.h>
#include <iomanip>

using namespace std;
struct Parameters{
    string       energyFile;       // To store the file with the knowledge-based potential.
};
struct Atom{
    double       coordinates[3];   // Coordinates of the atom.
    string       Type;             // Atom type for the scoring function when it is assigned as a string.
    int          TypeNumber;       // Atom type for the scoring function when it is assigned as a number.
    float        LJr;              // r0 value for the Lennard-Jones potential.
    float        LJe;              // Epsilon value for the Lennard-Jones potential.
    int          index;            // Atom number in the protein (we need it for StericClash).
};
struct Molecule{
    OpenBabel::OBMol     obMol;    // Molecule from OpenBabel.
    vector<Atom>         atom;     // To store the potential types of the atoms and their coordinates.
};
struct Protein{
    OpenBabel::OBMol obMol;        // Molecule from OpenBabel.
    vector<Atom>     atom;         // To store the potential types of the protein atoms and their coordinates.
};

double dist(double const atom1[3], double const atom2[3]);
bool   IsThereABond(Molecule const & molecule, unsigned int const i, unsigned int const j);
double KBP2016(Parameters const & parameters, Protein const & protein, Molecule & molecule);
double LJP(Parameters const & parameters, Protein const & protein, Molecule & molecule);
double RotorXScore(OpenBabel::OBMol obMol);
void   LigandTypeSMOG2016(Molecule & molecule);
void   ProteinTypeSMOG2016(string const residue, string const atomID, int & atomTypeNumber, float & LJr, float & LJe);

// This programs computes single point scores with the SMoG2016 function.
int main(int nbarg, char * argv[])
{
    if(nbarg!=4) {
        cout << "Use with: './SMoG2016.exe Protein.pdb Protein.sdf DeltaG' where DeltaG is the experimental value. The purpose of this is to display it in the output. If you don't have it, write anything instead ('DG' e.g.)." << endl;
        exit(0);
        }

    // Define the energy file.
    Parameters parameters;
    parameters.energyFile="KBP-3.0-5.0-8.5.dat";

    // Define a Protein structure.
    Protein protein;
    // Create an OpenBabel object for the protein and open it.
    OpenBabel::OBConversion obConversion;
    OpenBabel::OBFormat *format = obConversion.FormatFromExt(argv[1]);
    obConversion.SetInFormat(format);
    obConversion.ReadFile(&protein.obMol, argv[1]);
    // This part assign the potential atomtypes for the protein atoms.
    for (unsigned int i=1; i<=protein.obMol.NumAtoms(); i++) {
        OpenBabel::OBAtom *atomProt;
        atomProt = protein.obMol.GetAtom(i);
        OpenBabel::OBResidue *obResidue;
        obResidue = atomProt->GetResidue();
        string residue=obResidue->GetName();
        string atomID=obResidue->GetAtomID(atomProt);
        residue.erase(remove(residue.begin(), residue.end(), ' '), residue.end());
        atomID.erase(remove(atomID.begin(), atomID.end(), ' '), atomID.end());
        // Get the parameters.
        Atom tempAtom;
        string atomType="";
        int atomTypeNumber=0;
        float LJe=0.0;
        float LJr=0.0;
        ProteinTypeSMOG2016(residue, atomID, atomTypeNumber, LJr, LJe);
        tempAtom.coordinates[0] = atomProt->GetX();
        tempAtom.coordinates[1] = atomProt->GetY();
        tempAtom.coordinates[2] = atomProt->GetZ();
        tempAtom.Type = atomType;
        tempAtom.TypeNumber = atomTypeNumber;
        tempAtom.LJr = LJr;
        tempAtom.LJe = LJe;
        tempAtom.index = i;
        protein.atom.push_back(tempAtom);
        }
    // Get the PDB ID.
    string PDBID=argv[1];
    unsigned int positionToErase = PDBID.find_last_of("_");
    PDBID.erase(positionToErase);

    // Prepare the ligand.
    Molecule moleculeRef;
    format = obConversion.FormatFromExt(argv[2]);
    obConversion.SetInFormat(format);
    obConversion.ReadFile(&moleculeRef.obMol, argv[2]);
    LigandTypeSMOG2016(moleculeRef);
    Molecule molecule=moleculeRef;

    //double ratioMass=(moleculeRef.obMol.GetExactMass()*protein.obMol.GetExactMass())/(moleculeRef.obMol.GetExactMass()+protein.obMol.GetExactMass());
    // Compute the energy.
    cout << fixed << setprecision(3) ;
    double energy = 0.032*(KBP2016(parameters, protein, molecule)+0.535487*LJP(parameters, protein, molecule)+1.91331*molecule.obMol.NumRotors()-21.9741*log(molecule.obMol.GetExactMass()));
    cout << PDBID << "\t" << argv[3] << "\t" << energy <<endl ;

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////

// We define a function which compute the distance between two atoms.
double dist(double const atom1[3], double const atom2[3])
{
    double Distance=sqrt( pow(atom1[0]-atom2[0],2) + pow(atom1[1]-atom2[1],2) + pow(atom1[2]-atom2[2],2) );
    return Distance;
}

/////////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function computes the energy between a protein and a ligand. It calls other functions defined below.
double KBP2016(Parameters const & parameters, Protein const & protein, Molecule & molecule)
{
    // Initialize parameters. In SMoG2016 we use a continuous description of the contacts. For example, if one shell starts at 4.5A,
    // if d<4.5-range there is a contact of 0, if d>4.5+range there is a contact of 1 (until the shell stops) and between the contact
    // value varies linearly: contact=(4.5+range-distance)/(2*range).
    float range = 0.1;
    int NbrLigandAtomType  = 14;    // Number of ligand atom types
    int NbrProteinAtomType = 30;    // Number of protein atom types

    // Define the arrays for the KBP.
    float KBP_1st_Shell[NbrProteinAtomType][NbrLigandAtomType];
    float KBP_2nd_Shell[NbrProteinAtomType][NbrLigandAtomType];
    float KBP_3rd_Shell[NbrProteinAtomType][NbrLigandAtomType];
    for(int i=0; i<NbrProteinAtomType; i++)  {
        for(int j=0; j<NbrLigandAtomType; j++) { KBP_1st_Shell[i][j]=0; }
        }
    for(int i=0; i<NbrProteinAtomType; i++)  {
        for(int j=0; j<NbrLigandAtomType; j++) { KBP_2nd_Shell[i][j]=0; }
        }
    for(int i=0; i<NbrProteinAtomType; i++)  {
        for(int j=0; j<NbrLigandAtomType; j++) { KBP_3rd_Shell[i][j]=0; }
        }

    // Open and read the KBP.dat file.
    //   L I G A N D
    // P
    // R
    // O
    // T
    // E
    // I
    // N

    // Open the energy file.
    ifstream KBP(parameters.energyFile.c_str(), ios::in) ;
    if(!KBP) {
        cerr << "ERROR1: The file " << parameters.energyFile << " is missing, or I can't read it." << endl;
        exit (-1);
        }
    // Size of the shells. It is read from the first line of the energyFile. This code works for 3 shells.
    float sizeShell_1=0.0;
    float sizeShell_2=0.0;
    float sizeShell_3=0.0;
    KBP >> sizeShell_1 >> sizeShell_2 >> sizeShell_3;

    // First shell
    for (int i=0; i<NbrProteinAtomType; i++) {
        for(int j=0; j<NbrLigandAtomType; j++) { KBP >> KBP_1st_Shell[i][j]; }
        }
    // Second shell
    for (int i=0; i<NbrProteinAtomType; i++) {
        for(int j=0; j<NbrLigandAtomType; j++) { KBP >> KBP_2nd_Shell[i][j]; }
        }
    // Third shell
    for (int i=0; i<NbrProteinAtomType; i++) {
        for(int j=0; j<NbrLigandAtomType; j++) { KBP >> KBP_3rd_Shell[i][j]; }
        }
    KBP.close();

    // Compute the KBP energy
    double ScoreShell1=0;
    double ScoreShell2=0;
    double ScoreShell3=0;
    for (unsigned int i=0 ; i < protein.atom.size() ; i++) {
        for (unsigned int j=0 ; j < molecule.atom.size() ; j++) {
            if   ( (protein.atom[i].TypeNumber == 0) || (molecule.atom[j].TypeNumber == 0) ) { }
            else {
                double distance = dist(protein.atom[i].coordinates, molecule.atom[j].coordinates);
                // Shell1
                if      (distance <= sizeShell_1-range) {
                    ScoreShell1 += KBP_1st_Shell[protein.atom[i].TypeNumber-1][molecule.atom[j].TypeNumber-1];
                    }
                // Shell1-Shell2
                else if (distance <= sizeShell_1+range) {
                    double c1 = (sizeShell_1+range-distance)/(2*range); 
                    double c2 = (distance-sizeShell_1+range)/(2*range); 
                    ScoreShell1 += c1*KBP_1st_Shell[protein.atom[i].TypeNumber-1][molecule.atom[j].TypeNumber-1] ;
                    ScoreShell2 += c2*KBP_2nd_Shell[protein.atom[i].TypeNumber-1][molecule.atom[j].TypeNumber-1];
                    }
                // Shell2
                else if (distance <= sizeShell_2-range) {
                    ScoreShell2 += KBP_2nd_Shell[protein.atom[i].TypeNumber-1][molecule.atom[j].TypeNumber-1];
                    }
                // Shell2-Shell3
                else if (distance <= sizeShell_2+range) {
                    double c1 = (sizeShell_2+range-distance)/(2*range); 
                    double c2 = (distance-sizeShell_2+range)/(2*range); 
                    ScoreShell2 += c1*KBP_2nd_Shell[protein.atom[i].TypeNumber-1][molecule.atom[j].TypeNumber-1];
                    ScoreShell3 += c2*KBP_3rd_Shell[protein.atom[i].TypeNumber-1][molecule.atom[j].TypeNumber-1];
                    }
                // Shell3
                else if (distance <= sizeShell_3-range) {
                    ScoreShell3 += KBP_3rd_Shell[protein.atom[i].TypeNumber-1][molecule.atom[j].TypeNumber-1];
                    }
                // Shell3-EndOfShell3
                else if (distance <= sizeShell_3+range) {
                    double c1 = (sizeShell_3+range-distance)/(2*range); 
                    ScoreShell3 += c1*KBP_3rd_Shell[protein.atom[i].TypeNumber-1][molecule.atom[j].TypeNumber-1];
                    }
                }
            }
        }
    // We sum the contributions of the three shells. Each score is divided by the average distance in the shells.
    double KBPScore = ScoreShell1/((0.0+sizeShell_1)/2) + ScoreShell2/((sizeShell_1+sizeShell_2)/2) + ScoreShell3/((sizeShell_2+sizeShell_3)/2);

    return KBPScore;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function computes the repulsion between a protein and a ligand.
double LJP(Parameters const & parameters, Protein const & protein, Molecule & molecule)
{
    // Get the repulsion part (1/r12) of the Lennard-Jones potential
    double LJPinter=0.0;
    for(unsigned int i=0; i < protein.atom.size(); i++)  {
        for(unsigned int j=0; j < molecule.atom.size(); j++)  {
            double distance = dist(protein.atom[i].coordinates, molecule.atom[j].coordinates);
            if ( (distance!=0) && (distance<5.0) && (protein.atom[i].LJe!=0) && (molecule.atom[j].LJe!=0) && (protein.atom[i].LJr!=0) && (molecule.atom[j].LJr!=0) ) { 
                LJPinter += sqrt(protein.atom[i].LJe*molecule.atom[j].LJe) *pow((protein.atom[i].LJr+molecule.atom[j].LJr)/(distance), 12);
                }
            }
        }
    double LJPintra=0.0;
    for (unsigned int i=1; i<=molecule.obMol.NumAtoms() ; i++ ) {
        //Compute repulsion with other atoms.
        for (unsigned int j=i+1; j<=molecule.obMol.NumAtoms() ; j++ ) {
            int goFurther=1;
            //Check if atom j is part of 1-2 or 1-3 interactions.
            if (IsThereABond(molecule, i, j)) { goFurther=0; }
            else {
                for (unsigned int k=1; k<=molecule.obMol.NumAtoms(); k++ ) {
                    if (k!=i && k!=j && IsThereABond(molecule, i, k) && IsThereABond(molecule, k, j)) { goFurther=0; }
                    }
                }
            //If everything is OK, compute repulsion.
            if ( (goFurther==1) && (molecule.atom[i-1].LJe!=0) && (molecule.atom[j-1].LJe!=0) && (molecule.atom[i-1].LJr!=0) && (molecule.atom[j-1].LJr!=0) ) { 
                double distance = dist(molecule.atom[i-1].coordinates, molecule.atom[j-1].coordinates);
                LJPintra += sqrt(molecule.atom[i-1].LJe*molecule.atom[j-1].LJe)*pow((molecule.atom[i-1].LJr+molecule.atom[j-1].LJr)/(distance), 12); 
                }
            }
        }

    double LJP = (LJPinter+LJPintra);
    return LJP;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function calculates a ligand entropy as defined in XScore.
double RotorXScore(OpenBabel::OBMol obMol) {
    //double RotorScore=obMol.NumRotors();
    double RotorScore=0.0;

    // Count in how many rotors each atom is involved.
    int RotorPerAtom[obMol.NumAtoms()]={};
    for (unsigned int k=0; k<obMol.NumBonds(); k++) {
        OpenBabel::OBBond *obBond;
        obBond=obMol.GetBond(k);
        if (obBond->IsRotor()) {
            RotorPerAtom[obBond->GetBeginAtomIdx()-1]++;
            RotorPerAtom[obBond->GetEndAtomIdx()-1]++;
            }
        }

    // Assign an energy to each atom depending on how many rotors it is involved in.
    for (unsigned int k=1; k<=obMol.NumAtoms(); k++) {
        double ScorePerAtom=0.0;
        if      (RotorPerAtom[k] == 1) { ScorePerAtom=0.5; }
        else if (RotorPerAtom[k] == 2) { ScorePerAtom=1;   }
        else if (RotorPerAtom[k] >  2) { ScorePerAtom=0.5; }
        RotorScore += ScorePerAtom;
        }

    return RotorScore;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function assigns the ligand atom type for SMOG2016. It also stores in the dynamic array molecule.atom the coordinates and the atom-types to facilitate the calculation of energy. We want here to point out that the order in which atom types are assigned is important: for example, a hydrogen bond donor oxygen (hydroxyl) can also be a hydrogen bond acceptor oxygen and it is important to use the same order as the one we have used to have consistent results.
void LigandTypeSMOG2016(Molecule & molecule)
{
    // Define the forcefield to get atom types.
    string forceFieldName="GAFF";
    OpenBabel::OBForceField *forceField=OpenBabel::OBForceField::FindForceField(forceFieldName);
    if (!forceField) {
        cerr << "ERROR2: Could not find forcefield " << forceFieldName << "." << endl;
        }
    forceField->SetLogFile(&cout);
    forceField->SetLogLevel(OBFF_LOGLVL_NONE);      // NONE / LOW / MEDIUM / HIGH
    // Setup the forcefield
    if (!forceField->Setup(molecule.obMol)) {
        cout << "ERROR3: Could not setup force field." << endl;
        }
    forceField->GetAtomTypes(molecule.obMol); 

    // Create an array of atoms.
    molecule.atom.clear();
    for (unsigned int i=1; i<=molecule.obMol.NumAtoms() ; i++ ) {
        OpenBabel::OBAtom *atom;
        atom = molecule.obMol.GetAtom(i);
        OpenBabel::OBPairData *type = (OpenBabel::OBPairData*) atom->GetData("FFAtomType");

        // Check neighbor atoms to see if they are different from C or H. These 6 lines are coming from typer.cpp from OpenBabel.
        int polarity = 0;
        OpenBabel::OBAtom *nbr;
        vector<OpenBabel::OBBond*>::iterator k;
        for (nbr = atom->BeginNbrAtom(k); nbr; nbr = atom->NextNbrAtom(k)) {
            if (nbr->IsNotCorH()) { polarity++; }
            }

        // Begin atom-typing for the scoring function.
        string atomType=type->GetValue();
        int atomTypeNumber=0;
        if      (atom->IsCarbon() && (atom->GetHyb()==3) && polarity==0)         { atomTypeNumber = 1;  }
        else if (atom->IsCarbon() && (atom->GetHyb()==2) && polarity==0)         { atomTypeNumber = 2;  }
        else if (atom->IsCarbon() && (atom->GetHyb()==1) && polarity==0)         { atomTypeNumber = 2;  }
        else if (atom->IsCarbon() && atom->MatchesSMARTS("[#6;$(C=O)]"))         { atomTypeNumber = 3;  }
        else if (atom->IsCarbon() && atom->MatchesSMARTS("[#6;$(C(=N)(N)(N))]")) { atomTypeNumber = 3;  }
        else if (atom->IsCarbon() && polarity!=0)                                { atomTypeNumber = 4;  }

        else if (atomType == "c")                                                { atomTypeNumber = 3;  }
        else if (atomType == "c1")                                               { atomTypeNumber = 2;  }
        else if (atomType == "c2")                                               { atomTypeNumber = 2;  }
        else if (atomType == "c3" && polarity!=0)                                { atomTypeNumber = 4;  }
        else if (atomType == "c3")                                               { atomTypeNumber = 1;  }
        else if (atomType == "ca")                                               { atomTypeNumber = 2;  }

        else if (atomType == "cc")                                               { atomTypeNumber = 2;  }
        else if (atomType == "cd")                                               { atomTypeNumber = 2;  }
        else if (atomType == "ce")                                               { atomTypeNumber = 2;  }
        else if (atomType == "cf")                                               { atomTypeNumber = 2;  }
        else if (atomType == "cu")                                               { atomTypeNumber = 2;  }
        else if (atomType == "cv")                                               { atomTypeNumber = 2;  }
        else if (atomType == "cx")                                               { atomTypeNumber = 1;  }
        else if (atomType == "cy")                                               { atomTypeNumber = 1;  }
        else if (atomType == "cz")                                               { atomTypeNumber = 2;  }
        else if (atomType == "cg")                                               { atomTypeNumber = 2;  } // sp1 bound to a C

        else if (atomType == "n")                                                { atomTypeNumber = 7;  }
        else if (atom->IsAmideNitrogen())                                        { atomTypeNumber = 7;  }
        else if (atom->IsNitrogen() && (atom->GetValence() > atom->GetHyb()))    { atomTypeNumber = 5;  }
        else if (atom->IsNitrogen() && atom->IsHbondDonor())                     { atomTypeNumber = 5;  }
        else if (atom->IsNitrogen() && atom->IsHbondAcceptor())                  { atomTypeNumber = 6;  }

        else if (atom->IsCarboxylOxygen() || atom->IsPhosphateOxygen() || atom->IsNitroOxygen()) { atomTypeNumber = 11; }
        else if (atom->IsOxygen() && atom->MatchesSMARTS("[O;$(O=*)]"))          { atomTypeNumber = 8;  }
        else if (atom->IsOxygen() && atom->IsHbondDonor())                       { atomTypeNumber = 9;  }
        else if (atom->IsOxygen() && atom->IsHbondAcceptor())                    { atomTypeNumber = 10; }
        else if (atomType == "o")                                                { atomTypeNumber = 10; }
        else if (atomType == "oh")                                               { atomTypeNumber = 9;  }
        else if (atomType == "os")                                               { atomTypeNumber = 10; }

        else if (atomType == "p3")                                               { atomTypeNumber = 12; }
        else if (atomType == "p5")                                               { atomTypeNumber = 12; }
        else if (atomType == "py")                                               { atomTypeNumber = 12; }

        else if (atomType == "s")                                                { atomTypeNumber = 13; }
        else if (atomType == "s2")                                               { atomTypeNumber = 13; }
        else if (atomType == "s4")                                               { atomTypeNumber = 13; }
        else if (atomType == "s6")                                               { atomTypeNumber = 13; }
        else if (atomType == "sh")                                               { atomTypeNumber = 13; }
        else if (atomType == "ss")                                               { atomTypeNumber = 13; }
        else if (atomType == "Sy")                                               { atomTypeNumber = 13; }

        else if (atomType == "br")                                               { atomTypeNumber = 14; }
        else if (atomType == "f")                                                { atomTypeNumber = 14; }
        else if (atomType == "i")                                                { atomTypeNumber = 14; }
        else if (atomType == "cl")                                               { atomTypeNumber = 14; }

        else if (atomType == "Si")                                               { atomTypeNumber = 0;  }
        else if (atomType == "B")                                                { atomTypeNumber = 0;  }
        else if (atomType == "Al")                                               { atomTypeNumber = 0;  }
        else if (atomType == "Pt")                                               { atomTypeNumber = 0;  }
        else if (atomType == "As")                                               { atomTypeNumber = 0;  }
        else if (atomType == "Ru")                                               { atomTypeNumber = 0;  }
        else if (atomType == "V")                                                { atomTypeNumber = 0;  }
        else if (atomType == "Se")                                               { atomTypeNumber = 0;  }
        else if (atomType == "Cu")                                               { atomTypeNumber = 0;  }
        else if (atomType == "Fe")                                               { atomTypeNumber = 0;  }
        else if (atomType == "Hg")                                               { atomTypeNumber = 0;  }

        else if (atomType == "h1")                                               { atomTypeNumber = 0;  }
        else if (atomType == "h2")                                               { atomTypeNumber = 0;  }
        else if (atomType == "h3")                                               { atomTypeNumber = 0;  }
        else if (atomType == "h4")                                               { atomTypeNumber = 0;  }
        else if (atomType == "h5")                                               { atomTypeNumber = 0;  }
        else if (atomType == "ha")                                               { atomTypeNumber = 0;  }
        else if (atomType == "hc")                                               { atomTypeNumber = 0;  }
        else if (atomType == "hn")                                               { atomTypeNumber = 0;  }
        else if (atomType == "ho")                                               { atomTypeNumber = 0;  }
        else if (atomType == "hp")                                               { atomTypeNumber = 0;  }
        else if (atomType == "hs")                                               { atomTypeNumber = 0;  }
        else if (atomType == "X")                                                { atomTypeNumber = 0;  }
        else                                                                     { atomTypeNumber = 0;  }

        // Begin atom-typing for the Lennard-Jones potential.
        float LJe=0.0;
        float LJr=0.0;
        // Open the VDW parameters file.
        ifstream VDWParameters("VDWParameters.txt");
        if(!VDWParameters) {
            cerr << "ERROR4: The file VDWParameters.txt is missing, or I can't read it." << endl;
            exit (-1);
            }
        // Read the file and look for the atom type to get epsilon and r0 values.
        string line;
        int parametersHasBeenSet=0;
        while(getline(VDWParameters, line)) {
            if (parametersHasBeenSet==0) {
                istringstream input(line);
                int atomTypeNumberLJ=0;
                string atomTypeLJ="";
                float rValueLJ=0.0;
                float eValueLJ=0.0;
                input >> atomTypeNumberLJ >> atomTypeLJ >> rValueLJ >> eValueLJ;
                if (atomType==atomTypeLJ) { LJr=rValueLJ; LJe=eValueLJ; parametersHasBeenSet++; }
                line="";
                }
            }
        VDWParameters.close();
        // In case no atom type was found by OpenBabel.
        if (atomType=="X") {
            if      (atom->IsHydrogen())       { LJe=0.0;    LJr=0.0;    }
            else if (atom->IsSulfur())         { LJe=0.25;   LJr=2.0;    }
            else if (atom->IsPhosphorus())     { LJe=0.2;    LJr=2.1;    }
            else if (atom->IsNitrogen())       { LJe=0.17;   LJr=1.85;   }
            else if (atom->IsCarboxylOxygen()) { LJe=0.21;   LJr=1.6612; }
            else if (atom->IsOxygen())         { LJe=0.1968; LJr=1.6886; }
            else if (atom->IsCarbon())         { LJe=0.0977; LJr=1.908;  }
            else                               { LJe=0.0;    LJr=0.0;    }
            }

        // Create a new Atom in the array.
        Atom tempAtom;
        tempAtom.coordinates[0] = atom->GetX();
        tempAtom.coordinates[1] = atom->GetY();
        tempAtom.coordinates[2] = atom->GetZ();
        tempAtom.Type = atomType;
        tempAtom.TypeNumber = atomTypeNumber;
        tempAtom.LJr = LJr;
        tempAtom.LJe = LJe;
        tempAtom.index = i;
        molecule.atom.push_back(tempAtom);
        }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function assigns the protein atom type for SMOG2016. These atom types are inspired by Chen2005: http://dx.doi.org/10.1110/ps.051440705
void ProteinTypeSMOG2016(string const residue, string const atomID, int & atomTypeNumber, float & LJr, float & LJe)
{
    // ALA : N, CA, C, O, CB, H
    if      ((residue == "ALA") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "ALA") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "ALA") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "ALA") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "ALA") && (atomID == "CB"))  { atomTypeNumber = 1;  }
    // ARG : N, CA, C, O, CB, CG, CD, NE, CZ, NH1, NH2, H, HE, HH11, HH12, HH21, HH22
    else if ((residue == "ARG") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "ARG") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "ARG") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "ARG") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "ARG") && (atomID == "CB"))  { atomTypeNumber = 15; }
    else if ((residue == "ARG") && (atomID == "CG"))  { atomTypeNumber = 16; }
    else if ((residue == "ARG") && (atomID == "CD"))  { atomTypeNumber = 16; }
    else if ((residue == "ARG") && (atomID == "NE"))  { atomTypeNumber = 19; }
    else if ((residue == "ARG") && (atomID == "CZ"))  { atomTypeNumber = 18; }
    else if ((residue == "ARG") && (atomID == "NH1")) { atomTypeNumber = 30; }
    else if ((residue == "ARG") && (atomID == "NH2")) { atomTypeNumber = 30; }
    // ARN = ARG-H
    else if ((residue == "ARN") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "ARN") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "ARN") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "ARN") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "ARN") && (atomID == "CB"))  { atomTypeNumber = 15; }
    else if ((residue == "ARN") && (atomID == "CG"))  { atomTypeNumber = 16; }
    else if ((residue == "ARN") && (atomID == "CD"))  { atomTypeNumber = 16; }
    else if ((residue == "ARN") && (atomID == "NE"))  { atomTypeNumber = 19; }
    else if ((residue == "ARN") && (atomID == "CZ"))  { atomTypeNumber = 18; }
    else if ((residue == "ARN") && (atomID == "NH1")) { atomTypeNumber = 30; }
    else if ((residue == "ARN") && (atomID == "NH2")) { atomTypeNumber = 30; }
    // ASN : N, CA, C, O, CB, CG, OD1, ND2, H, HD21, HD22
    else if ((residue == "ASN") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "ASN") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "ASN") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "ASN") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "ASN") && (atomID == "CB"))  { atomTypeNumber = 10; }
    else if ((residue == "ASN") && (atomID == "CG"))  { atomTypeNumber = 13; }
    else if ((residue == "ASN") && (atomID == "OD1")) { atomTypeNumber = 29; }
    else if ((residue == "ASN") && (atomID == "ND2")) { atomTypeNumber = 14; }
    // ASP : N, CA, C, O, CB, CG, OD1, OD2, H
    else if ((residue == "ASP") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "ASP") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "ASP") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "ASP") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "ASP") && (atomID == "CB"))  { atomTypeNumber = 20; }
    else if ((residue == "ASP") && (atomID == "CG"))  { atomTypeNumber = 21; }
    else if ((residue == "ASP") && (atomID == "OD1")) { atomTypeNumber = 22; }
    else if ((residue == "ASP") && (atomID == "OD2")) { atomTypeNumber = 22; }
    // ASH = ASP+H
    else if ((residue == "ASH") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "ASH") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "ASH") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "ASH") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "ASH") && (atomID == "CB"))  { atomTypeNumber = 20; }
    else if ((residue == "ASH") && (atomID == "CG"))  { atomTypeNumber = 21; }
    else if ((residue == "ASH") && (atomID == "OD1")) { atomTypeNumber = 22; }
    else if ((residue == "ASH") && (atomID == "OD2")) { atomTypeNumber = 22; }
    // CYS : N, CA, C, O, CB, SG, H
    else if ((residue == "CYS") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "CYS") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "CYS") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "CYS") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "CYS") && (atomID == "CB"))  { atomTypeNumber = 1;  }
    else if ((residue == "CYS") && (atomID == "SG"))  { atomTypeNumber = 4;  }
    // CYX = CYS for S-S bridge
    else if ((residue == "CYX") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "CYX") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "CYX") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "CYX") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "CYX") && (atomID == "CB"))  { atomTypeNumber = 1;  }
    else if ((residue == "CYX") && (atomID == "SG"))  { atomTypeNumber = 4;  }
    // GLU : N, CA, C, O, CB, CG, CD, OE1, OE2, H
    else if ((residue == "GLU") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "GLU") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "GLU") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "GLU") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "GLU") && (atomID == "CB"))  { atomTypeNumber = 20; }
    else if ((residue == "GLU") && (atomID == "CG"))  { atomTypeNumber = 20; }
    else if ((residue == "GLU") && (atomID == "CD"))  { atomTypeNumber = 21; }
    else if ((residue == "GLU") && (atomID == "OE1")) { atomTypeNumber = 22; }
    else if ((residue == "GLU") && (atomID == "OE2")) { atomTypeNumber = 22; }
    // GLH = GLU+H
    else if ((residue == "GLH") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "GLH") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "GLH") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "GLH") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "GLH") && (atomID == "CB"))  { atomTypeNumber = 20; }
    else if ((residue == "GLH") && (atomID == "CG"))  { atomTypeNumber = 20; }
    else if ((residue == "GLH") && (atomID == "CD"))  { atomTypeNumber = 21; }
    else if ((residue == "GLH") && (atomID == "OE1")) { atomTypeNumber = 22; }
    else if ((residue == "GLH") && (atomID == "OE2")) { atomTypeNumber = 22; }
    // GLN : N, CA, C, O, CB, CG, CD, OE1, NE2, H, HE21, HE22
    else if ((residue == "GLN") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "GLN") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "GLN") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "GLN") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "GLN") && (atomID == "CB"))  { atomTypeNumber = 10; }
    else if ((residue == "GLN") && (atomID == "CG"))  { atomTypeNumber = 12; }
    else if ((residue == "GLN") && (atomID == "CD"))  { atomTypeNumber = 13; }
    else if ((residue == "GLN") && (atomID == "OE1")) { atomTypeNumber = 29; }
    else if ((residue == "GLN") && (atomID == "NE2")) { atomTypeNumber = 14; }
    // GLY : N, CA, C, O, H
    else if ((residue == "GLY") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "GLY") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "GLY") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "GLY") && (atomID == "O"))   { atomTypeNumber = 28; }
    // HIS : N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2, H, HD1, HE2
    else if ((residue == "HIS") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "HIS") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "HIS") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "HIS") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "HIS") && (atomID == "CB"))  { atomTypeNumber = 5;  }
    else if ((residue == "HIS") && (atomID == "CG"))  { atomTypeNumber = 6;  }
    else if ((residue == "HIS") && (atomID == "ND1")) { atomTypeNumber = 9;  }
    else if ((residue == "HIS") && (atomID == "CD2")) { atomTypeNumber = 6;  }
    else if ((residue == "HIS") && (atomID == "CE1")) { atomTypeNumber = 6;  }
    else if ((residue == "HIS") && (atomID == "NE2")) { atomTypeNumber = 9;  }
    // HID = HIS protonated on delta
    else if ((residue == "HID") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "HID") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "HID") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "HID") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "HID") && (atomID == "CB"))  { atomTypeNumber = 5;  }
    else if ((residue == "HID") && (atomID == "CG"))  { atomTypeNumber = 6;  }
    else if ((residue == "HID") && (atomID == "ND1")) { atomTypeNumber = 9;  }
    else if ((residue == "HID") && (atomID == "CD2")) { atomTypeNumber = 6;  }
    else if ((residue == "HID") && (atomID == "CE1")) { atomTypeNumber = 6;  }
    else if ((residue == "HID") && (atomID == "NE2")) { atomTypeNumber = 9;  }
    // HIE = HIS protonated on epsilon
    else if ((residue == "HIE") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "HIE") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "HIE") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "HIE") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "HIE") && (atomID == "CB"))  { atomTypeNumber = 5;  }
    else if ((residue == "HIE") && (atomID == "CG"))  { atomTypeNumber = 6;  }
    else if ((residue == "HIE") && (atomID == "ND1")) { atomTypeNumber = 9;  }
    else if ((residue == "HIE") && (atomID == "CD2")) { atomTypeNumber = 6;  }
    else if ((residue == "HIE") && (atomID == "CE1")) { atomTypeNumber = 6;  }
    else if ((residue == "HIE") && (atomID == "NE2")) { atomTypeNumber = 9;  }
    // HIP = HIS doubly protonated 
    else if ((residue == "HIP") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "HIP") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "HIP") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "HIP") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "HIP") && (atomID == "CB"))  { atomTypeNumber = 5;  }
    else if ((residue == "HIP") && (atomID == "CG"))  { atomTypeNumber = 6;  }
    else if ((residue == "HIP") && (atomID == "ND1")) { atomTypeNumber = 9;  }
    else if ((residue == "HIP") && (atomID == "CD2")) { atomTypeNumber = 9;  }
    else if ((residue == "HIP") && (atomID == "CE1")) { atomTypeNumber = 6;  }
    else if ((residue == "HIP") && (atomID == "NE2")) { atomTypeNumber = 9;  }
    // ILE : N, CA, C, O, CB, CG1, CG2, CD1, H
    else if ((residue == "ILE") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "ILE") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "ILE") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "ILE") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "ILE") && (atomID == "CB"))  { atomTypeNumber = 1;  }
    else if ((residue == "ILE") && (atomID == "CG1")) { atomTypeNumber = 2;  }
    else if ((residue == "ILE") && (atomID == "CG2")) { atomTypeNumber = 2;  }
    else if ((residue == "ILE") && (atomID == "CD1")) { atomTypeNumber = 2;  }
    else if ((residue == "ILE") && (atomID == "CD"))  { atomTypeNumber = 2;  }
    // LEU : N, CA, C, O, CB, CG, CD1, CD2, H
    else if ((residue == "LEU") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "LEU") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "LEU") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "LEU") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "LEU") && (atomID == "CB"))  { atomTypeNumber = 1;  }
    else if ((residue == "LEU") && (atomID == "CG"))  { atomTypeNumber = 2;  }
    else if ((residue == "LEU") && (atomID == "CD1")) { atomTypeNumber = 2;  }
    else if ((residue == "LEU") && (atomID == "CD2")) { atomTypeNumber = 2;  }
    // LYS : N, CA, C, O, CB, CG, CD, CE, NZ, H, HZ1, HZ2, HZ3
    else if ((residue == "LYS") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "LYS") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "LYS") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "LYS") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "LYS") && (atomID == "CB"))  { atomTypeNumber = 15; }
    else if ((residue == "LYS") && (atomID == "CG"))  { atomTypeNumber = 16; }
    else if ((residue == "LYS") && (atomID == "CD"))  { atomTypeNumber = 16; }
    else if ((residue == "LYS") && (atomID == "CE"))  { atomTypeNumber = 16; }
    else if ((residue == "LYS") && (atomID == "NZ"))  { atomTypeNumber = 17; }
    // LYN = LYS-H
    else if ((residue == "LYN") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "LYN") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "LYN") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "LYN") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "LYN") && (atomID == "CB"))  { atomTypeNumber = 15; }
    else if ((residue == "LYN") && (atomID == "CG"))  { atomTypeNumber = 16; }
    else if ((residue == "LYN") && (atomID == "CD"))  { atomTypeNumber = 16; }
    else if ((residue == "LYN") && (atomID == "CE"))  { atomTypeNumber = 16; }
    else if ((residue == "LYN") && (atomID == "NZ"))  { atomTypeNumber = 17; }
    // MET : N, CA, C, O, CB, CG, SD, CE, H
    else if ((residue == "MET") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "MET") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "MET") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "MET") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "MET") && (atomID == "CB"))  { atomTypeNumber = 1;  }
    else if ((residue == "MET") && (atomID == "CG"))  { atomTypeNumber = 2;  }
    else if ((residue == "MET") && (atomID == "SD"))  { atomTypeNumber = 3;  }
    else if ((residue == "MET") && (atomID == "CE"))  { atomTypeNumber = 2;  }
    // PHE : N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, H
    else if ((residue == "PHE") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "PHE") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "PHE") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "PHE") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "PHE") && (atomID == "CB"))  { atomTypeNumber = 5;  }
    else if ((residue == "PHE") && (atomID == "CG"))  { atomTypeNumber = 6;  }
    else if ((residue == "PHE") && (atomID == "CD1")) { atomTypeNumber = 6;  }
    else if ((residue == "PHE") && (atomID == "CD2")) { atomTypeNumber = 6;  }
    else if ((residue == "PHE") && (atomID == "CE1")) { atomTypeNumber = 6;  }
    else if ((residue == "PHE") && (atomID == "CE2")) { atomTypeNumber = 6;  }
    else if ((residue == "PHE") && (atomID == "CZ"))  { atomTypeNumber = 6;  }
    // PRO : N, CA, C, O, CB, CG, CD, H2, H3
    else if ((residue == "PRO") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "PRO") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "PRO") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "PRO") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "PRO") && (atomID == "CB"))  { atomTypeNumber = 23; }
    else if ((residue == "PRO") && (atomID == "CG"))  { atomTypeNumber = 23; }
    else if ((residue == "PRO") && (atomID == "CD"))  { atomTypeNumber = 23; }
    // SER : N, CA, C, O, CB, OG
    else if ((residue == "SER") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "SER") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "SER") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "SER") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "SER") && (atomID == "CB"))  { atomTypeNumber = 10; }
    else if ((residue == "SER") && (atomID == "OG"))  { atomTypeNumber = 11; }
    // THR : N, CA, C, O, CB, OG1, CG2, H, HG1
    else if ((residue == "THR") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "THR") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "THR") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "THR") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "THR") && (atomID == "CB"))  { atomTypeNumber = 10; }
    else if ((residue == "THR") && (atomID == "OG1")) { atomTypeNumber = 11; }
    else if ((residue == "THR") && (atomID == "CG2")) { atomTypeNumber = 12; }
    // TRP : N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2, H, HE1
    else if ((residue == "TRP") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "TRP") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "TRP") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "TRP") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "TRP") && (atomID == "CB"))  { atomTypeNumber = 5;  }
    else if ((residue == "TRP") && (atomID == "CG"))  { atomTypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "CD1")) { atomTypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "CD2")) { atomTypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "NE1")) { atomTypeNumber = 7;  }
    else if ((residue == "TRP") && (atomID == "CE2")) { atomTypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "CE3")) { atomTypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "CZ2")) { atomTypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "CZ3")) { atomTypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "CH2")) { atomTypeNumber = 6;  }
    // TYR : N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH, H, HH
    else if ((residue == "TYR") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "TYR") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "TYR") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "TYR") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "TYR") && (atomID == "CB"))  { atomTypeNumber = 5;  }
    else if ((residue == "TYR") && (atomID == "CG"))  { atomTypeNumber = 6;  }
    else if ((residue == "TYR") && (atomID == "CD1")) { atomTypeNumber = 6;  }
    else if ((residue == "TYR") && (atomID == "CD2")) { atomTypeNumber = 6;  }
    else if ((residue == "TYR") && (atomID == "CE1")) { atomTypeNumber = 6;  }
    else if ((residue == "TYR") && (atomID == "CE2")) { atomTypeNumber = 6;  }
    else if ((residue == "TYR") && (atomID == "CZ"))  { atomTypeNumber = 6;  }
    else if ((residue == "TYR") && (atomID == "OH"))  { atomTypeNumber = 8;  }
    // VAL : N, CA, C, O, CB, CG1, CG2, H
    else if ((residue == "VAL") && (atomID == "N"))   { atomTypeNumber = 25; }
    else if ((residue == "VAL") && (atomID == "CA"))  { atomTypeNumber = 26; }
    else if ((residue == "VAL") && (atomID == "C"))   { atomTypeNumber = 27; }
    else if ((residue == "VAL") && (atomID == "O"))   { atomTypeNumber = 28; }
    else if ((residue == "VAL") && (atomID == "CB"))  { atomTypeNumber = 1;  }
    else if ((residue == "VAL") && (atomID == "CG1")) { atomTypeNumber = 2;  }
    else if ((residue == "VAL") && (atomID == "CG2")) { atomTypeNumber = 2;  }
    // Metals: CA, ZN, HG, MG, CD, NI, MN, NA
    else if ((residue == "CA")  && (atomID == "CA")) { atomTypeNumber = 24; }
    else if ((residue == "ZN")  && (atomID == "ZN")) { atomTypeNumber = 24; }
    else if ((residue == "HG")  && (atomID == "HG")) { atomTypeNumber = 24; }
    else if ((residue == "MG")  && (atomID == "MG")) { atomTypeNumber = 24; }
    else if ((residue == "CD")  && (atomID == "CD")) { atomTypeNumber = 24; }
    else if ((residue == "NI")  && (atomID == "NI")) { atomTypeNumber = 24; }
    else if ((residue == "MN")  && (atomID == "MN")) { atomTypeNumber = 24; }
    else if ((residue == "NA")  && (atomID == "NA")) { atomTypeNumber = 24; }
    // H
    else                                             { atomTypeNumber = 0;  }

    // Begin atom-typing for the Lennard-Jones potential.
    // Open the parameter file for the Amber atom types.
    ifstream AmberTypes("AmberAtomTypes.txt");
    if(!AmberTypes) {
        cerr << "ERROR5: The file AmberAtomTypes.txt is missing, or I can't read it." << endl;
        exit (-1);
        }
    // Parse the parameterFile by looking for specific words.
    string lineAmberTypes;
    string atomType="";
    int typeHasBeenSet=0;
    while(getline(AmberTypes, lineAmberTypes)) {
        if (typeHasBeenSet==0) {
            istringstream input(lineAmberTypes);
            int number=0;
            int typeValue=0;
            string atomTypeTemp="";
            string residueName="";
            string atomIDName="";
            input >> number >> typeValue >> atomTypeTemp >> residueName >> atomIDName;
            if ((residueName==residue) && (atomIDName==atomID)) { atomType=atomTypeTemp; typeHasBeenSet++; }
            lineAmberTypes="";
            }
        }
    AmberTypes.close();

    // Open the VDW parameters file.
    ifstream VDWParameters("VDWParameters.txt");
    if(!VDWParameters) {
        cerr << "ERROR6: The file VDWParameters.txt is missing, or I can't read it." << endl;
        exit (-1);
        }
    // Read the file and look for the atom type to get epsilon and r0 values.
    string lineVDWParameters="";
    int parametersHasBeenSet=0;
    while(getline(VDWParameters, lineVDWParameters)) {
        if (parametersHasBeenSet==0) {
        istringstream input(lineVDWParameters);
        int atomTypeNumberLJ=0;
        string atomTypeLJ="";
        float rValueLJ=0.0;
        float eValueLJ=0.0;
        input >> atomTypeNumberLJ >> atomTypeLJ >> rValueLJ >> eValueLJ;
        if (atomType==atomTypeLJ) { LJr=rValueLJ; LJe=eValueLJ; parametersHasBeenSet++; }
        lineVDWParameters="";
        }
        }
    VDWParameters.close();
}
