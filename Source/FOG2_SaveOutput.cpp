#include "OpenGrowth.h"
#include <iomanip>                // for setprecision
#include <openbabel/parsmart.h>

int SizeFile(string const fileName);

// This function saves a molecule to an .xyz output file. The name is the OUTPUT parameter, then "_x.xyz" where x is a counter.
void SaveOutput(Parameters const & parameters, Molecule & molecule, int const & counterOutput, string const threeMerList[])
{
    // Open the summary file.
    string outputSummary = parameters.outputName + "_Summary.txt";
    ofstream outputSum(outputSummary.c_str(), ios::out|ios::app);
    if(!outputSum) {
        cerr << "ERROR: I cannot write in the output file " << outputSummary << "." << endl;
        exit (-1);
        }
    // Open the summary SMILES file.
    string outputSummarySMILES = parameters.outputName + "_Summary.smi";
    ofstream outputSumSMILES(outputSummarySMILES.c_str(), ios::out|ios::app);
    if(!outputSumSMILES) {
        cerr << "ERROR: I cannot write in the output file " << outputSummarySMILES << "." << endl;
        exit (-1);
        }

    // Counters for the output file.
    ostringstream oss_counterOutput;
    oss_counterOutput << counterOutput;
    string outputNumber = oss_counterOutput.str();

    string outputFile;
    string outputFileName;
    outputFile = parameters.outputName + "_" + outputNumber;
    outputFileName = parameters.outputName + "_" + outputNumber + ".xyz";

    OpenBabel::OBConversion obConversion;
    if(parameters.smilesOnly==0) {
        obConversion.SetOutFormat("xyz");
        molecule.obMol.SetTitle(outputFile);
        obConversion.WriteFile(&molecule.obMol, outputFileName);
        }

    // Get the SMILES string.
    obConversion.SetOutFormat("smi");
    obConversion.AddOption("n",OpenBabel::OBConversion::OUTOPTIONS); 
    string smilesString = obConversion.WriteString(&molecule.obMol);
    smilesString[smilesString.length()-1]=' ';
    outputSumSMILES << smilesString << "\t" << parameters.outputName + "_" + outputNumber << endl;

    // 3Mer-Screen
    int Counting=0;
    if (parameters.threeMerFile!="") {
        for (int z=0 ; z < parameters.threeMerSize && !Counting ; z++) {
            OpenBabel::OBSmartsPattern obSmarts;
            obSmarts.Init(threeMerList[z]);
            if (obSmarts.HasMatch(molecule.obMol)) { Counting++; }
            }
        }

    // Write in the summary file (the first line is for the setprecision below). It is not very beautiful, but I am doing it this way because I am trying to have a
    // formatted output file. I am not using \t because depending on the text editor it is not of the same size.
    outputSum << fixed;
    // Name
    outputSum << outputFile ;
    if      ( counterOutput < 10 )       { outputSum << "         " ; }
    else if ( counterOutput < 100 )      { outputSum << "        " ;  }
    else if ( counterOutput < 1000 )     { outputSum << "       " ;   }
    else if ( counterOutput < 10000 )    { outputSum << "      " ;    }
    else if ( counterOutput < 100000 )   { outputSum << "     " ;     }
    else if ( counterOutput < 1000000 )  { outputSum << "    " ;      }
    else                                 { outputSum << "   " ;       }
    outputSum << setprecision(2) ;
    // #Residues
    outputSum << molecule.numberResidues ;
    if      ( molecule.numberResidues <= 9 )  { outputSum << "             " ; }
    else if ( molecule.numberResidues <= 99 ) { outputSum << "            " ;  }
    else                                      { outputSum << "           " ;   }
    // #HeavyAtoms
    outputSum << molecule.obMol.NumHvyAtoms() ;
    if      ( molecule.obMol.NumHvyAtoms() <= 9 )  { outputSum << "               " ; }
    else if ( molecule.obMol.NumHvyAtoms() <= 99 ) { outputSum << "              " ;  }
    else                                           { outputSum << "             " ;   }
    // Molecular weight
    outputSum << molecule.obMol.GetMolWt() ;
    if      ( molecule.obMol.GetMolWt() < 100 )  { outputSum << "                " ; }
    else if ( molecule.obMol.GetMolWt() < 1000 ) { outputSum << "               " ;  }
    else                                         { outputSum << "              " ;   }
    // 3Mer
    if      (parameters.threeMerFile!="" && !Counting) { outputSum << "PASS     " ; }
    else if (parameters.threeMerFile!="" && Counting)  { outputSum << "FAIL     " ; }
    // SMILES
    outputSum << "\t" << smilesString << endl; 

    // Close the summary file and the input file.
    outputSum.close();
    outputSumSMILES.close();
}

