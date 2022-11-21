#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/atom.h>
#include <iomanip>

using namespace std;
struct Atom{
    double       coordinates[3];   // Coordinates of the atom.
    string       Type;             // Atom type for the scoring function when it is assigned as a string.
    int          TypeNumber;       // Atom type for the scoring function when it is assigned as a number.
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
void ProteinTypeSMOG2016(string const residue, string const atomID, int & TypeNumber);
void LigandTypeSMOG2016(Molecule & molecule);

int main(int nbarg, char * argv[])
{
    if(nbarg!=6) {
           cout << "Use with: './KBP2016-Training.exe 3.0 5.0 8.5 TrainingList.txt Library/'." << endl;
           exit(0);
           }

    double range = 0.1 ;           // range is a parameter for the smoothing of the shells. Between two shells the contact is smoothed by a linear function between SizeShell-range and SizeShell+range.
    double alpha = 0.9 ;           // alpha and beta are two empirical parameters used in th Knowlege Based Potential Training function.  
    double beta  = 0.9 ;           // They were Optimized in the SMoG2001 function (ref above), the same values are used in this function.
    int NbrLigandAtomType  = 14 ;  // How many atom types for the ligand.
    int NbrProteinAtomType = 30 ;  // How many atom types for the protein.

    //First shell
    double TotalShell_1 = 0;                                            //This variable is used to store the number of protein-ligand atom contacts in the first shell.
    double NbrProteinAtomOfTypeShell_1[NbrProteinAtomType+1];           //This table is used to store the number of contacts of each protein atom type in the first shell.
    double NbrLigandAtomOfTypeShell_1[NbrLigandAtomType+1];             //This table is udes to store the number of contacts of each ligand atom type in the first shell.
    double sizeShell_1 = atof(argv[1]);                                 //End of the first shell.
    double contactsShell_1[NbrProteinAtomType+1][NbrLigandAtomType+1];  // This table is used to store the number of contact of each ligand-protein pairs.

    //Second shell
    double TotalShell_2 = 0;
    double NbrProteinAtomOfTypeShell_2[NbrProteinAtomType+1];
    double NbrLigandAtomOfTypeShell_2[NbrLigandAtomType+1];
    double sizeShell_2 = atof(argv[2]);
    double contactsShell_2[NbrProteinAtomType+1][NbrLigandAtomType+1];

    //Third shell
    double TotalShell_3 = 0;
    double NbrProteinAtomOfTypeShell_3[NbrProteinAtomType+1];
    double NbrLigandAtomOfTypeShell_3[NbrLigandAtomType+1];
    double sizeShell_3 = atof(argv[3]);
    double contactsShell_3[NbrProteinAtomType+1][NbrLigandAtomType+1];

    //Assign 0 to all tables
    for (int i =0 ; i < NbrProteinAtomType+1 ; i++) {
         NbrProteinAtomOfTypeShell_1[i] = 0; 
         NbrProteinAtomOfTypeShell_2[i] = 0; 
         NbrProteinAtomOfTypeShell_3[i] = 0; 
         }
    for (int i =0 ; i < NbrLigandAtomType+1 ; i++) {
        NbrLigandAtomOfTypeShell_1[i] = 0;
        NbrLigandAtomOfTypeShell_2[i] = 0;
        NbrLigandAtomOfTypeShell_3[i] = 0;
        }
    for(int i =0 ; i < NbrProteinAtomType+1 ; i++)  {
        for(int j =0 ; j < NbrLigandAtomType+1 ; j++)  {
            contactsShell_1[i][j] = 0; 
            contactsShell_2[i][j] = 0; 
            contactsShell_3[i][j] = 0; 
            }
        }

    //Open the training list file.
    ifstream TrainingFlux(argv[4]);
    if (!TrainingFlux) {
        cerr << "ERROR: The file " << argv[4] << " is missing, or I can't read it.";
        exit(-1);
        }     

    //Loop over all the pdb names listed in the training list.
    string TrainingFile;
    while(getline(TrainingFlux,TrainingFile)) {
        string LibraryPath=argv[5];
        string proteinFile = LibraryPath + "/" + TrainingFile + "/" + TrainingFile + "_protein.pdb";
        string ligandFile = LibraryPath + "/" + TrainingFile  + "/" + TrainingFile + "_ligand.sdf";

        //Preparation of the protein. Start by creating an OpenBabel object for the protein.
        Protein protein;
        OpenBabel::OBConversion obConversion;
        OpenBabel::OBFormat *format = obConversion.FormatFromExt(proteinFile.c_str());
        obConversion.SetInFormat(format);
        obConversion.ReadFile(&protein.obMol, proteinFile.c_str());
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
            int TypeNumber=0;
            ProteinTypeSMOG2016(residue, atomID, TypeNumber);
            tempAtom.coordinates[0] = atomProt->GetX();
            tempAtom.coordinates[1] = atomProt->GetY();
            tempAtom.coordinates[2] = atomProt->GetZ();
            tempAtom.Type = atomType;
            tempAtom.TypeNumber = TypeNumber;
            protein.atom.push_back(tempAtom);
            }
        
        //Preparation of the ligand.
        Molecule ligand;
        format = obConversion.FormatFromExt(ligandFile);
        obConversion.SetInFormat(format);
        obConversion.ReadFile(&ligand.obMol, ligandFile);
        LigandTypeSMOG2016(ligand);

        //Compute the number of contacts in each shell.
        for(unsigned int i=0 ; i < ligand.obMol.NumAtoms(); i++) {
           for(unsigned int j=0; j< protein.obMol.NumAtoms(); j++) {
              if((ligand.atom[i].TypeNumber == 0) || (protein.atom[j].TypeNumber == 0)){;}        //To exclude all atoms that do not have atom types (hydrogens, metals in ligands..)
                 else  {    
                    double distance = dist(protein.atom[j].coordinates, ligand.atom[i].coordinates);            
                    if (distance <= sizeShell_1-range)    {                       
                       contactsShell_1[protein.atom[j].TypeNumber][ligand.atom[i].TypeNumber]++;   //Number of contacts between each atom types pairs in shell 1.
                       TotalShell_1++ ;                                                            //Total number of contacts in shell 1.
                       NbrProteinAtomOfTypeShell_1[protein.atom[j].TypeNumber]++;                  //Number of protein atoms of type TypeNumber involved in a contact in shell 1.
                       NbrLigandAtomOfTypeShell_1[ligand.atom[i].TypeNumber]++;                    //Number of ligand atoms of type TypeNumber involved in a contact in shell 1.                  
                       }

                    else if (distance <= sizeShell_1+range) {                             //From sizeShell1-range to sizeShell1+range the function is smoothed.
                       double c1 = (sizeShell_1+range-distance)/(2*range);                //When a contact is found there, c1 is added instead of 1. c1 decreases from 1 to 0.
                       double c2 = (distance-sizeShell_1+range)/(2*range);                //For Shell2, c2 implements the counting. Its behavior is the opposite of c1: c2 linearly increases from 0 to 1.
                       contactsShell_1[protein.atom[j].TypeNumber][ligand.atom[i].TypeNumber] += c1;     
                       TotalShell_1 += c1 ;                                                            
                       NbrProteinAtomOfTypeShell_1[protein.atom[j].TypeNumber] += c1;                      
                       NbrLigandAtomOfTypeShell_1[ligand.atom[i].TypeNumber] += c1;
                       contactsShell_2[protein.atom[j].TypeNumber][ligand.atom[i].TypeNumber] += c2;
                       TotalShell_2 += c2 ;
                       NbrProteinAtomOfTypeShell_2[protein.atom[j].TypeNumber] += c2; 
                       NbrLigandAtomOfTypeShell_2[ligand.atom[i].TypeNumber] += c2;
                       }

                    else if (distance <= sizeShell_2-range)  {
                       contactsShell_2[protein.atom[j].TypeNumber][ligand.atom[i].TypeNumber]++;
                       TotalShell_2++ ;
                       NbrProteinAtomOfTypeShell_2[protein.atom[j].TypeNumber]++; 
                       NbrLigandAtomOfTypeShell_2[ligand.atom[i].TypeNumber]++;
                       }

                    else if (distance <= sizeShell_2+range) {
                       double c1 = (sizeShell_2+range-distance)/(2*range);
                       double c2 = (distance-sizeShell_2+range)/(2*range);
                       contactsShell_2[protein.atom[j].TypeNumber][ligand.atom[i].TypeNumber] += c1;     
                       TotalShell_2 += c1 ;                                                            
                       NbrProteinAtomOfTypeShell_2[protein.atom[j].TypeNumber] += c1;                      
                       NbrLigandAtomOfTypeShell_2[ligand.atom[i].TypeNumber] += c1;
                       contactsShell_3[protein.atom[j].TypeNumber][ligand.atom[i].TypeNumber] += c2;
                       TotalShell_3 += c2 ;
                       NbrProteinAtomOfTypeShell_3[protein.atom[j].TypeNumber] += c2; 
                       NbrLigandAtomOfTypeShell_3[ligand.atom[i].TypeNumber] += c2;
                       }

                    else if (distance <= sizeShell_3-range)  {
                       contactsShell_3[protein.atom[j].TypeNumber][ligand.atom[i].TypeNumber]++;
                       TotalShell_3++ ;
                       NbrProteinAtomOfTypeShell_3[protein.atom[j].TypeNumber]++; 
                       NbrLigandAtomOfTypeShell_3[ligand.atom[i].TypeNumber]++;
                       }

                    else if (distance <= sizeShell_3+range) {
                       double c1 = (sizeShell_3+range-distance)/(2*range);
                       contactsShell_3[protein.atom[j].TypeNumber][ligand.atom[i].TypeNumber] += c1;     
                       TotalShell_3 += c1 ;                                                            
                       NbrProteinAtomOfTypeShell_3[protein.atom[j].TypeNumber] += c1;                      
                       NbrLigandAtomOfTypeShell_3[ligand.atom[i].TypeNumber] += c1;
                       }
                    } //End of loop on protein.
                } //End of loop on ligand.
            } //End of the loop over all protein-ligand atoms pairs.
        } //End of the loop over all the pdb names.

    //To display the number of contacts.
    /*
    for(int i=1 ; i<=NbrProteinAtomType ;i++) {
        for(int j=1 ; j<=NbrLigandAtomType ;j++) {
            cout << fixed << setprecision(2) << contactsShell_1[i][j]<< "\t" << "\t" ;            
            }
            cout << endl ;
        }
    cout << endl;
    */

    //Free energy calculation for each protein-ligand atom pair. From the contacts counting, the energy associated of 14x30 types of contact (NumberOfLigandAtomType x NumberOfProteinAtomType)
    //is computed in each shell. The form of the energy function is the same as the SMoG2001 function.

    //Calculation of the normalization constants for each shells.
    double NormConst_1 = 0.; 
    double NormConst_2 = 0.;
    double NormConst_3 = 0.;
    for(int i=1; i<=NbrProteinAtomType; i++)  {
        for(int j=1; j<=NbrLigandAtomType; j++)  {
            if(contactsShell_1[i][j] == 0)   {contactsShell_1[i][j] = 1;}    // This is to avoid an error (log[0]) if the number of contact is 0.
            NormConst_1 += 2.302585*log(contactsShell_1[i][j]/(pow(NbrProteinAtomOfTypeShell_1[i],alpha)*pow(NbrLigandAtomOfTypeShell_1[j],beta)*TotalShell_1)) ;
            if(contactsShell_2[i][j] == 0)   {contactsShell_2[i][j] = 1;}
            NormConst_2 += 2.302585*log(contactsShell_2[i][j]/(pow(NbrProteinAtomOfTypeShell_2[i],alpha)*pow(NbrLigandAtomOfTypeShell_2[j],beta)*TotalShell_2)) ;
            if(contactsShell_3[i][j] == 0)   {contactsShell_3[i][j] = 1;}
            NormConst_3 += 2.302585*log(contactsShell_3[i][j]/(pow(NbrProteinAtomOfTypeShell_3[i],alpha)*pow(NbrLigandAtomOfTypeShell_3[j],beta)*TotalShell_3)) ;
            }
        }
    NormConst_1 = NormConst_1/(NbrProteinAtomType*NbrLigandAtomType);
    NormConst_2 = NormConst_2/(NbrProteinAtomType*NbrLigandAtomType);
    NormConst_3 = NormConst_3/(NbrProteinAtomType*NbrLigandAtomType);

    //Calculation of the free energy for each types of contact.
    double F_1[NbrProteinAtomType+1][NbrLigandAtomType+1] ;
    double F_2[NbrProteinAtomType+1][NbrLigandAtomType+1] ;
    double F_3[NbrProteinAtomType+1][NbrLigandAtomType+1] ;
    for(int i=1 ;i<NbrProteinAtomType+1 ;i++)  {
        for(int j=1 ;j<NbrLigandAtomType+1 ;j++)  {
            F_1[i][j] = NormConst_1 - 2.302585*log(contactsShell_1[i][j]/(pow(NbrProteinAtomOfTypeShell_1[i],alpha)*pow(NbrLigandAtomOfTypeShell_1[j],beta)*TotalShell_1)) ; 
            F_2[i][j] = NormConst_2 - 2.302585*log(contactsShell_2[i][j]/(pow(NbrProteinAtomOfTypeShell_2[i],alpha)*pow(NbrLigandAtomOfTypeShell_2[j],beta)*TotalShell_2)) ; 
            F_3[i][j] = NormConst_3 - 2.302585*log(contactsShell_3[i][j]/(pow(NbrProteinAtomOfTypeShell_3[i],alpha)*pow(NbrLigandAtomOfTypeShell_3[j],beta)*TotalShell_3)) ; 
            }
        }

    //Calculation of the F0 term.
    double F0_1 = 0;
    double F0_2 = 0;
    double F0_3 = 0;
    for(int i=1 ;i<NbrProteinAtomType+1 ;i++)  {
        for(int j=1 ;j<NbrLigandAtomType+1 ;j++)  {
            F0_1 += F_1[i][j] ;
            F0_2 += F_2[i][j] ;
            F0_3 += F_3[i][j] ;
            }
        }
    F0_1 /= NbrProteinAtomType*NbrLigandAtomType ;
    F0_2 /= NbrProteinAtomType*NbrLigandAtomType ;
    F0_3 /= NbrProteinAtomType*NbrLigandAtomType ;

    //Correct the free energy terms.
    for(int i=1 ;i<NbrProteinAtomType+1 ;i++)  {
        for(int j=1 ;j<NbrLigandAtomType+1 ;j++)  {
            F_1[i][j] += F0_1 ;
            F_2[i][j] += F0_2 ;
            F_3[i][j] += F0_3 ;
            }
        }    

    //Write output file.
    string firstShell=argv[1]; 
    string secondShell=argv[2]; 
    string thirdShell=argv[3]; 
    string KBPOutFile = "KBP-" + firstShell + "-" + secondShell + "-" + thirdShell + ".dat";
    ofstream KBPOutFlux(KBPOutFile.c_str());
    //Shell 1
    KBPOutFlux << argv[1] << " " << argv[2] << " " << argv[3] << endl;                //First line is the shell configuration.
    for(int i=1 ; i<=NbrProteinAtomType ;i++) {
        for(int j=1 ; j<=NbrLigandAtomType ;j++)  {
            KBPOutFlux << fixed << setprecision (4) << F_1[i][j] << "\t" << "\t";    //Display the free energies for each type of contact. 30 lines for the protein atom types, 14 columns for the ligand atom types.
            }
            KBPOutFlux << endl ;
        }
    KBPOutFlux << endl;
    //Shell 2
    for(int i=1 ; i<=NbrProteinAtomType ;i++) {
        for(int j=1 ; j<=NbrLigandAtomType ;j++)  {
            KBPOutFlux  << F_2[i][j] << "\t" << "\t";            
            }
            KBPOutFlux << endl ;
        }
    KBPOutFlux << endl ;
    //Shell 3
    for(int i=1 ; i<=NbrProteinAtomType ;i++) {
        for(int j=1 ; j<=NbrLigandAtomType ;j++)  {
            KBPOutFlux  << F_3[i][j] << "\t" << "\t";            
            }
            KBPOutFlux << endl ;
        }
    KBPOutFlux << endl ;

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////

// We define a function which compute the distance between two atoms.
double dist(double const atom1[3], double const atom2[3])
{
    double Distance=sqrt( pow(atom1[0]-atom2[0],2) + pow(atom1[1]-atom2[1],2) + pow(atom1[2]-atom2[2],2) );
    return Distance;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function assigns the protein atom type for SMOG2016. These atom types are inspired by Chen2005: http://dx.doi.org/10.1110/ps.051440705
void ProteinTypeSMOG2016(string const residue, string const atomID, int & TypeNumber)
{
    // ALA : N, CA, C, O, CB, H
    if      ((residue == "ALA") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "ALA") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "ALA") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "ALA") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "ALA") && (atomID == "CB"))  { TypeNumber = 1;  }
    // ARG : N, CA, C, O, CB, CG, CD, NE, CZ, NH1, NH2, H, HE, HH11, HH12, HH21, HH22
    else if ((residue == "ARG") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "ARG") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "ARG") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "ARG") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "ARG") && (atomID == "CB"))  { TypeNumber = 15; }
    else if ((residue == "ARG") && (atomID == "CG"))  { TypeNumber = 16; }
    else if ((residue == "ARG") && (atomID == "CD"))  { TypeNumber = 16; }
    else if ((residue == "ARG") && (atomID == "NE"))  { TypeNumber = 19; }
    else if ((residue == "ARG") && (atomID == "CZ"))  { TypeNumber = 18; }
    else if ((residue == "ARG") && (atomID == "NH1")) { TypeNumber = 30; }
    else if ((residue == "ARG") && (atomID == "NH2")) { TypeNumber = 30; }
    // ARN = ARG-H
    else if ((residue == "ARN") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "ARN") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "ARN") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "ARN") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "ARN") && (atomID == "CB"))  { TypeNumber = 15; }
    else if ((residue == "ARN") && (atomID == "CG"))  { TypeNumber = 16; }
    else if ((residue == "ARN") && (atomID == "CD"))  { TypeNumber = 16; }
    else if ((residue == "ARN") && (atomID == "NE"))  { TypeNumber = 19; }
    else if ((residue == "ARN") && (atomID == "CZ"))  { TypeNumber = 18; }
    else if ((residue == "ARN") && (atomID == "NH1")) { TypeNumber = 30; }
    else if ((residue == "ARN") && (atomID == "NH2")) { TypeNumber = 30; }
    // ASN : N, CA, C, O, CB, CG, OD1, ND2, H, HD21, HD22
    else if ((residue == "ASN") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "ASN") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "ASN") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "ASN") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "ASN") && (atomID == "CB"))  { TypeNumber = 10; }
    else if ((residue == "ASN") && (atomID == "CG"))  { TypeNumber = 13; }
    else if ((residue == "ASN") && (atomID == "OD1")) { TypeNumber = 29; }
    else if ((residue == "ASN") && (atomID == "ND2")) { TypeNumber = 14; }
    // ASP : N, CA, C, O, CB, CG, OD1, OD2, H
    else if ((residue == "ASP") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "ASP") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "ASP") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "ASP") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "ASP") && (atomID == "CB"))  { TypeNumber = 20; }
    else if ((residue == "ASP") && (atomID == "CG"))  { TypeNumber = 21; }
    else if ((residue == "ASP") && (atomID == "OD1")) { TypeNumber = 22; }
    else if ((residue == "ASP") && (atomID == "OD2")) { TypeNumber = 22; }
    // ASH = ASP+H
    else if ((residue == "ASH") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "ASH") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "ASH") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "ASH") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "ASH") && (atomID == "CB"))  { TypeNumber = 20; }
    else if ((residue == "ASH") && (atomID == "CG"))  { TypeNumber = 21; }
    else if ((residue == "ASH") && (atomID == "OD1")) { TypeNumber = 22; }
    else if ((residue == "ASH") && (atomID == "OD2")) { TypeNumber = 22; }
    // CYS : N, CA, C, O, CB, SG, H
    else if ((residue == "CYS") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "CYS") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "CYS") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "CYS") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "CYS") && (atomID == "CB"))  { TypeNumber = 1;  }
    else if ((residue == "CYS") && (atomID == "SG"))  { TypeNumber = 4;  }
    // CYX = CYS for S-S bridge
    else if ((residue == "CYX") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "CYX") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "CYX") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "CYX") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "CYX") && (atomID == "CB"))  { TypeNumber = 1;  }
    else if ((residue == "CYX") && (atomID == "SG"))  { TypeNumber = 4;  }
    // GLU : N, CA, C, O, CB, CG, CD, OE1, OE2, H
    else if ((residue == "GLU") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "GLU") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "GLU") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "GLU") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "GLU") && (atomID == "CB"))  { TypeNumber = 20; }
    else if ((residue == "GLU") && (atomID == "CG"))  { TypeNumber = 20; }
    else if ((residue == "GLU") && (atomID == "CD"))  { TypeNumber = 21; }
    else if ((residue == "GLU") && (atomID == "OE1")) { TypeNumber = 22; }
    else if ((residue == "GLU") && (atomID == "OE2")) { TypeNumber = 22; }
    // GLH = GLU+H
    else if ((residue == "GLH") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "GLH") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "GLH") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "GLH") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "GLH") && (atomID == "CB"))  { TypeNumber = 20; }
    else if ((residue == "GLH") && (atomID == "CG"))  { TypeNumber = 20; }
    else if ((residue == "GLH") && (atomID == "CD"))  { TypeNumber = 21; }
    else if ((residue == "GLH") && (atomID == "OE1")) { TypeNumber = 22; }
    else if ((residue == "GLH") && (atomID == "OE2")) { TypeNumber = 22; }
    // GLN : N, CA, C, O, CB, CG, CD, OE1, NE2, H, HE21, HE22
    else if ((residue == "GLN") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "GLN") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "GLN") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "GLN") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "GLN") && (atomID == "CB"))  { TypeNumber = 10; }
    else if ((residue == "GLN") && (atomID == "CG"))  { TypeNumber = 12; }
    else if ((residue == "GLN") && (atomID == "CD"))  { TypeNumber = 13; }
    else if ((residue == "GLN") && (atomID == "OE1")) { TypeNumber = 29; }
    else if ((residue == "GLN") && (atomID == "NE2")) { TypeNumber = 14; }
    // GLY : N, CA, C, O, H
    else if ((residue == "GLY") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "GLY") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "GLY") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "GLY") && (atomID == "O"))   { TypeNumber = 28; }
    // HIS : N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2, H, HD1, HE2
    else if ((residue == "HIS") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "HIS") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "HIS") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "HIS") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "HIS") && (atomID == "CB"))  { TypeNumber = 5;  }
    else if ((residue == "HIS") && (atomID == "CG"))  { TypeNumber = 6;  }
    else if ((residue == "HIS") && (atomID == "ND1")) { TypeNumber = 9;  }
    else if ((residue == "HIS") && (atomID == "CD2")) { TypeNumber = 6;  }
    else if ((residue == "HIS") && (atomID == "CE1")) { TypeNumber = 6;  }
    else if ((residue == "HIS") && (atomID == "NE2")) { TypeNumber = 9;  }
    // HID = HIS protonated on delta
    else if ((residue == "HID") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "HID") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "HID") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "HID") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "HID") && (atomID == "CB"))  { TypeNumber = 5;  }
    else if ((residue == "HID") && (atomID == "CG"))  { TypeNumber = 6;  }
    else if ((residue == "HID") && (atomID == "ND1")) { TypeNumber = 9;  }
    else if ((residue == "HID") && (atomID == "CD2")) { TypeNumber = 6;  }
    else if ((residue == "HID") && (atomID == "CE1")) { TypeNumber = 6;  }
    else if ((residue == "HID") && (atomID == "NE2")) { TypeNumber = 9;  }
    // HIE = HIS protonated on epsilon
    else if ((residue == "HIE") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "HIE") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "HIE") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "HIE") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "HIE") && (atomID == "CB"))  { TypeNumber = 5;  }
    else if ((residue == "HIE") && (atomID == "CG"))  { TypeNumber = 6;  }
    else if ((residue == "HIE") && (atomID == "ND1")) { TypeNumber = 9;  }
    else if ((residue == "HIE") && (atomID == "CD2")) { TypeNumber = 6;  }
    else if ((residue == "HIE") && (atomID == "CE1")) { TypeNumber = 6;  }
    else if ((residue == "HIE") && (atomID == "NE2")) { TypeNumber = 9;  }
    // HIP = HIS doubly protonated 
    else if ((residue == "HIP") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "HIP") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "HIP") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "HIP") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "HIP") && (atomID == "CB"))  { TypeNumber = 5;  }
    else if ((residue == "HIP") && (atomID == "CG"))  { TypeNumber = 6;  }
    else if ((residue == "HIP") && (atomID == "ND1")) { TypeNumber = 9;  }
    else if ((residue == "HIP") && (atomID == "CD2")) { TypeNumber = 9;  }
    else if ((residue == "HIP") && (atomID == "CE1")) { TypeNumber = 6;  }
    else if ((residue == "HIP") && (atomID == "NE2")) { TypeNumber = 9;  }
    // ILE : N, CA, C, O, CB, CG1, CG2, CD1, H
    else if ((residue == "ILE") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "ILE") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "ILE") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "ILE") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "ILE") && (atomID == "CB"))  { TypeNumber = 1;  }
    else if ((residue == "ILE") && (atomID == "CG1")) { TypeNumber = 2;  }
    else if ((residue == "ILE") && (atomID == "CG2")) { TypeNumber = 2;  }
    else if ((residue == "ILE") && (atomID == "CD1")) { TypeNumber = 2;  }
    else if ((residue == "ILE") && (atomID == "CD"))  { TypeNumber = 2;  }
    // LEU : N, CA, C, O, CB, CG, CD1, CD2, H
    else if ((residue == "LEU") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "LEU") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "LEU") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "LEU") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "LEU") && (atomID == "CB"))  { TypeNumber = 1;  }
    else if ((residue == "LEU") && (atomID == "CG"))  { TypeNumber = 2;  }
    else if ((residue == "LEU") && (atomID == "CD1")) { TypeNumber = 2;  }
    else if ((residue == "LEU") && (atomID == "CD2")) { TypeNumber = 2;  }
    // LYS : N, CA, C, O, CB, CG, CD, CE, NZ, H, HZ1, HZ2, HZ3
    else if ((residue == "LYS") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "LYS") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "LYS") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "LYS") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "LYS") && (atomID == "CB"))  { TypeNumber = 15; }
    else if ((residue == "LYS") && (atomID == "CG"))  { TypeNumber = 16; }
    else if ((residue == "LYS") && (atomID == "CD"))  { TypeNumber = 16; }
    else if ((residue == "LYS") && (atomID == "CE"))  { TypeNumber = 16; }
    else if ((residue == "LYS") && (atomID == "NZ"))  { TypeNumber = 17; }
    // LYN = LYS-H
    else if ((residue == "LYN") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "LYN") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "LYN") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "LYN") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "LYN") && (atomID == "CB"))  { TypeNumber = 15; }
    else if ((residue == "LYN") && (atomID == "CG"))  { TypeNumber = 16; }
    else if ((residue == "LYN") && (atomID == "CD"))  { TypeNumber = 16; }
    else if ((residue == "LYN") && (atomID == "CE"))  { TypeNumber = 16; }
    else if ((residue == "LYN") && (atomID == "NZ"))  { TypeNumber = 17; }
    // MET : N, CA, C, O, CB, CG, SD, CE, H
    else if ((residue == "MET") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "MET") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "MET") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "MET") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "MET") && (atomID == "CB"))  { TypeNumber = 1;  }
    else if ((residue == "MET") && (atomID == "CG"))  { TypeNumber = 2;  }
    else if ((residue == "MET") && (atomID == "SD"))  { TypeNumber = 3;  }
    else if ((residue == "MET") && (atomID == "CE"))  { TypeNumber = 2;  }
    // PHE : N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, H
    else if ((residue == "PHE") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "PHE") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "PHE") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "PHE") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "PHE") && (atomID == "CB"))  { TypeNumber = 5;  }
    else if ((residue == "PHE") && (atomID == "CG"))  { TypeNumber = 6;  }
    else if ((residue == "PHE") && (atomID == "CD1")) { TypeNumber = 6;  }
    else if ((residue == "PHE") && (atomID == "CD2")) { TypeNumber = 6;  }
    else if ((residue == "PHE") && (atomID == "CE1")) { TypeNumber = 6;  }
    else if ((residue == "PHE") && (atomID == "CE2")) { TypeNumber = 6;  }
    else if ((residue == "PHE") && (atomID == "CZ"))  { TypeNumber = 6;  }
    // PRO : N, CA, C, O, CB, CG, CD, H2, H3
    else if ((residue == "PRO") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "PRO") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "PRO") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "PRO") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "PRO") && (atomID == "CB"))  { TypeNumber = 23; }
    else if ((residue == "PRO") && (atomID == "CG"))  { TypeNumber = 23; }
    else if ((residue == "PRO") && (atomID == "CD"))  { TypeNumber = 23; }
    // SER : N, CA, C, O, CB, OG
    else if ((residue == "SER") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "SER") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "SER") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "SER") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "SER") && (atomID == "CB"))  { TypeNumber = 10; }
    else if ((residue == "SER") && (atomID == "OG"))  { TypeNumber = 11; }
    // THR : N, CA, C, O, CB, OG1, CG2, H, HG1
    else if ((residue == "THR") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "THR") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "THR") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "THR") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "THR") && (atomID == "CB"))  { TypeNumber = 10; }
    else if ((residue == "THR") && (atomID == "OG1")) { TypeNumber = 11; }
    else if ((residue == "THR") && (atomID == "CG2")) { TypeNumber = 12; }
    // TRP : N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2, H, HE1
    else if ((residue == "TRP") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "TRP") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "TRP") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "TRP") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "TRP") && (atomID == "CB"))  { TypeNumber = 5;  }
    else if ((residue == "TRP") && (atomID == "CG"))  { TypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "CD1")) { TypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "CD2")) { TypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "NE1")) { TypeNumber = 7;  }
    else if ((residue == "TRP") && (atomID == "CE2")) { TypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "CE3")) { TypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "CZ2")) { TypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "CZ3")) { TypeNumber = 6;  }
    else if ((residue == "TRP") && (atomID == "CH2")) { TypeNumber = 6;  }
    // TYR : N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH, H, HH
    else if ((residue == "TYR") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "TYR") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "TYR") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "TYR") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "TYR") && (atomID == "CB"))  { TypeNumber = 5;  }
    else if ((residue == "TYR") && (atomID == "CG"))  { TypeNumber = 6;  }
    else if ((residue == "TYR") && (atomID == "CD1")) { TypeNumber = 6;  }
    else if ((residue == "TYR") && (atomID == "CD2")) { TypeNumber = 6;  }
    else if ((residue == "TYR") && (atomID == "CE1")) { TypeNumber = 6;  }
    else if ((residue == "TYR") && (atomID == "CE2")) { TypeNumber = 6;  }
    else if ((residue == "TYR") && (atomID == "CZ"))  { TypeNumber = 6;  }
    else if ((residue == "TYR") && (atomID == "OH"))  { TypeNumber = 8;  }
    // VAL : N, CA, C, O, CB, CG1, CG2, H
    else if ((residue == "VAL") && (atomID == "N"))   { TypeNumber = 25; }
    else if ((residue == "VAL") && (atomID == "CA"))  { TypeNumber = 26; }
    else if ((residue == "VAL") && (atomID == "C"))   { TypeNumber = 27; }
    else if ((residue == "VAL") && (atomID == "O"))   { TypeNumber = 28; }
    else if ((residue == "VAL") && (atomID == "CB"))  { TypeNumber = 1;  }
    else if ((residue == "VAL") && (atomID == "CG1")) { TypeNumber = 2;  }
    else if ((residue == "VAL") && (atomID == "CG2")) { TypeNumber = 2;  }
    // Metals: CA, ZN, HG, MG, CD, NI, MN, NA
    else if ((residue == "CA")  && (atomID == "CA")) { TypeNumber = 24; }
    else if ((residue == "ZN")  && (atomID == "ZN")) { TypeNumber = 24; }
    else if ((residue == "HG")  && (atomID == "HG")) { TypeNumber = 24; }
    else if ((residue == "MG")  && (atomID == "MG")) { TypeNumber = 24; }
    else if ((residue == "CD")  && (atomID == "CD")) { TypeNumber = 24; }
    else if ((residue == "NI")  && (atomID == "NI")) { TypeNumber = 24; }
    else if ((residue == "MN")  && (atomID == "MN")) { TypeNumber = 24; }
    else if ((residue == "NA")  && (atomID == "NA")) { TypeNumber = 24; }
    // H
    else                                             { TypeNumber = 0;  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function assigns the ligand atom type for SMOG2016. It also stores in the dynamic array molecule.atom the coordinates and the atom-types to facilitate the calculation of energy. We want here to point out that the order in which atom types are assigned is important: for example, a hydrogen bond donor oxygen (hydroxyl) can also be a hydrogen bond acceptor oxygen and it is important to use the same order as the one we have used to have consistent results.
void LigandTypeSMOG2016(Molecule & molecule)
{
    // Define the forcefield to get atom types.
    string forceFieldName="GAFF";
    OpenBabel::OBForceField *forceField=OpenBabel::OBForceField::FindForceField(forceFieldName);
    if (!forceField) {
        cerr << "ERROR: Could not find forcefield " << forceFieldName << "." << endl;
        }
    forceField->SetLogFile(&cout);
    forceField->SetLogLevel(OBFF_LOGLVL_NONE);      // NONE / LOW / MEDIUM / HIGH
    // Setup the forcefield
    if (!forceField->Setup(molecule.obMol)) {
        cout << "ERROR: Could not setup force field." << endl;
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
        int TypeNumber=0;
        if      (atom->IsCarbon() && (atom->GetHyb()==3) && polarity==0)         { TypeNumber = 1;  }
        else if (atom->IsCarbon() && (atom->GetHyb()==2) && polarity==0)         { TypeNumber = 2;  }
        else if (atom->IsCarbon() && (atom->GetHyb()==1) && polarity==0)         { TypeNumber = 2;  }
        else if (atom->IsCarbon() && atom->MatchesSMARTS("[#6;$(C=O)]"))         { TypeNumber = 3;  }
        else if (atom->IsCarbon() && atom->MatchesSMARTS("[#6;$(C(=N)(N)(N))]")) { TypeNumber = 3;  }
        else if (atom->IsCarbon() && polarity!=0)                                { TypeNumber = 4;  }

        else if (atomType == "c")                                                { TypeNumber = 3;  }
        else if (atomType == "c1")                                               { TypeNumber = 2;  }
        else if (atomType == "c2")                                               { TypeNumber = 2;  }
        else if (atomType == "c3" && polarity!=0)                                { TypeNumber = 4;  }
        else if (atomType == "c3")                                               { TypeNumber = 1;  }
        else if (atomType == "ca")                                               { TypeNumber = 2;  }

        else if (atomType == "cc")                                               { TypeNumber = 2;  }
        else if (atomType == "cd")                                               { TypeNumber = 2;  }
        else if (atomType == "ce")                                               { TypeNumber = 2;  }
        else if (atomType == "cf")                                               { TypeNumber = 2;  }
        else if (atomType == "cu")                                               { TypeNumber = 2;  }
        else if (atomType == "cv")                                               { TypeNumber = 2;  }
        else if (atomType == "cx")                                               { TypeNumber = 1;  }
        else if (atomType == "cy")                                               { TypeNumber = 1;  }
        else if (atomType == "cz")                                               { TypeNumber = 2;  }
        else if (atomType == "cg")                                               { TypeNumber = 2;  } // sp1 bound to a C

        else if (atomType == "n")                                                { TypeNumber = 7;  }
        else if (atom->IsAmideNitrogen())                                        { TypeNumber = 7;  }
        else if (atom->IsNitrogen() && (atom->GetValence() > atom->GetHyb()))    { TypeNumber = 5;  }
        else if (atom->IsNitrogen() && atom->IsHbondDonor())                     { TypeNumber = 5;  }
        else if (atom->IsNitrogen() && atom->IsHbondAcceptor())                  { TypeNumber = 6;  }

        else if (atom->IsCarboxylOxygen() || atom->IsPhosphateOxygen() || atom->IsNitroOxygen()) { TypeNumber = 11; }
        else if (atom->IsOxygen() && atom->MatchesSMARTS("[O;$(O=*)]"))          { TypeNumber = 8;  }
        else if (atom->IsOxygen() && atom->IsHbondDonor())                       { TypeNumber = 9;  }
        else if (atom->IsOxygen() && atom->IsHbondAcceptor())                    { TypeNumber = 10; }
        else if (atomType == "o")                                                { TypeNumber = 10; }
        else if (atomType == "oh")                                               { TypeNumber = 9;  }
        else if (atomType == "os")                                               { TypeNumber = 10; }

        else if (atomType == "p3")                                               { TypeNumber = 12; }
        else if (atomType == "p5")                                               { TypeNumber = 12; }
        else if (atomType == "py")                                               { TypeNumber = 12; }

        else if (atomType == "s")                                                { TypeNumber = 13; }
        else if (atomType == "s2")                                               { TypeNumber = 13; }
        else if (atomType == "s4")                                               { TypeNumber = 13; }
        else if (atomType == "s6")                                               { TypeNumber = 13; }
        else if (atomType == "sh")                                               { TypeNumber = 13; }
        else if (atomType == "ss")                                               { TypeNumber = 13; }
        else if (atomType == "Sy")                                               { TypeNumber = 13; }

        else if (atomType == "br")                                               { TypeNumber = 14; }
        else if (atomType == "f")                                                { TypeNumber = 14; }
        else if (atomType == "i")                                                { TypeNumber = 14; }
        else if (atomType == "cl")                                               { TypeNumber = 14; }

        else if (atomType == "Si")                                               { TypeNumber = 0;  }
        else if (atomType == "B")                                                { TypeNumber = 0;  }
        else if (atomType == "Al")                                               { TypeNumber = 0;  }
        else if (atomType == "Pt")                                               { TypeNumber = 0;  }
        else if (atomType == "As")                                               { TypeNumber = 0;  }
        else if (atomType == "Ru")                                               { TypeNumber = 0;  }
        else if (atomType == "V")                                                { TypeNumber = 0;  }
        else if (atomType == "Se")                                               { TypeNumber = 0;  }
        else if (atomType == "Cu")                                               { TypeNumber = 0;  }
        else if (atomType == "Fe")                                               { TypeNumber = 0;  }
        else if (atomType == "Hg")                                               { TypeNumber = 0;  }

        else if (atomType == "h1")                                               { TypeNumber = 0;  }
        else if (atomType == "h2")                                               { TypeNumber = 0;  }
        else if (atomType == "h3")                                               { TypeNumber = 0;  }
        else if (atomType == "h4")                                               { TypeNumber = 0;  }
        else if (atomType == "h5")                                               { TypeNumber = 0;  }
        else if (atomType == "ha")                                               { TypeNumber = 0;  }
        else if (atomType == "hc")                                               { TypeNumber = 0;  }
        else if (atomType == "hn")                                               { TypeNumber = 0;  }
        else if (atomType == "ho")                                               { TypeNumber = 0;  }
        else if (atomType == "hp")                                               { TypeNumber = 0;  }
        else if (atomType == "hs")                                               { TypeNumber = 0;  }
        else if (atomType == "X")                                                { TypeNumber = 0;  }
        else                                                                     { TypeNumber = 0;  }

        // Create a new Atom in the array.
        Atom tempAtom;
        tempAtom.coordinates[0] = atom->GetX();
        tempAtom.coordinates[1] = atom->GetY();
        tempAtom.coordinates[2] = atom->GetZ();
        tempAtom.Type = atomType;
        tempAtom.TypeNumber = TypeNumber;
        molecule.atom.push_back(tempAtom);
        }
}
