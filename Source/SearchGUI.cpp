#include <iomanip>       // For setprecision
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>

using namespace std;

int  SizeFile(char const fileName[]);
int  findFragment(string databasepath, string fragment, char runType);
void computeProbabilities(string databasepath, string resourcepath, string outputpath, char runType);
void computeProbabilitiesSplitted(string databasepath, string resourcepath, string outputpath, string counter);
void findThreemers(string databasepath, string resourcepath,  string outputpath, char runType);

int main(int argc, char * argv[])
{
    if (argc < 3) {
        cerr << "Missing arguments" << endl;
        return 1;
        }
    
    char process = argv[1][0];
    string databasepath = argv[2];
    
    // Standard run type if not defined. If runType == 'q' program communicates with the GUI.
    char runType = 's';

    // Find Fragments
    if (process == 'f') {
        string fragment = argv[3];
        if (argc > 4) {
            runType = argv[4][0];
            }
        int number = findFragment(databasepath, fragment, runType);
        cout << number << endl;
        }

    // Threemer Search
    else if (process == 't') {
        string resourcepath = argv[3];
        string outputpath;
        if (argc > 5) {
            runType = argv[5][0];
            outputpath = argv[4];
            }
        else if (argc > 4) {
            outputpath = argv[4];
            }
        findThreemers(databasepath, resourcepath, outputpath, runType);
        }

    // FirstFrag/Transition Probabilities
    else if (process == 'p') {
        string resourcepath = argv[3];
        string outputpath;
        if (argc > 5) {
            runType = argv[5][0];
            outputpath = argv[4];
            }
        else if (argc > 4) {
            outputpath = argv[4];
            }
        computeProbabilities(databasepath, resourcepath, outputpath, runType);
        }

    // FirstFrag/Transition Probabilities IF we work in a splitted mode (i.e. we want to compute probabilites for one fragment, and we will then gather everything).
    // We don't use runType here because this is not made to be used with the GUI.
    else if (process == 's') {
        string resourcepath = argv[3];
        string outputpath;
        outputpath = argv[4];
        string counter(argv[5]);
        computeProbabilitiesSplitted(databasepath, resourcepath, outputpath, counter);
        }

    else {
        cerr << "Process code not present: 'f' to find fragment, 'p' to compute probabilities, 's' to compute probabilities in a splitted mode, and 't' to find unfound 3mers" << endl;
        }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////

// This function calculates the size of a file.
int SizeFile(char const fileName[])
{
    // Open the file and check it can be opened.
    ifstream listFile(fileName);
    if(!listFile) 
        {
        cerr << "ERROR: The file " << fileName << " is missing, or I can't read it." << endl;
        exit (-1);
        }
    // If a read line is not empty, increase the counter
    int fragmentListSize=0;
    string line;
    while(getline(listFile, line)) 
        {
        if(strcmp(line.c_str(),"")) { fragmentListSize++; }
        }
    
    return fragmentListSize;
}

/////////////////////////////////////////////////////////////////////////////////////

// Find fragment takes a path to a database and a SMARTS fragment query, and computes the number of times the fragment is found in that database
int findFragment(string databasepath, string fragment, char runType)
{
    // OpenBabel setup
    OpenBabel::OBConversion obconversion;
    OpenBabel::OBMol mol;
    obconversion.SetInFormat("smi");
    OpenBabel::OBSmartsPattern smarts;

    int total = 0;

    // Loads the SMARTS query
    smarts.Init(fragment);
    std::vector<std::vector<int> > maplist;
    bool notatend = obconversion.ReadFile(&mol, databasepath.c_str());

    // Loops through database, counting number of matches found
    while (notatend) {
    smarts.Match(mol);
    maplist = smarts.GetUMapList();
    total += maplist.size();
    mol.Clear();
    notatend = obconversion.Read(&mol);
    }
    return total;
}

/////////////////////////////////////////////////////////////////////////////////////

// Computes probability files Proba_FirstFrag.dat and Proba_Transition.dat
void computeProbabilities(string databasepath, string resourcepath, string outputpath, char runType)
{
    // runState accepts communication from GUI
    string runState;

    string libraryFileName = databasepath;
    if (*resourcepath.rbegin() != '/') {
        resourcepath = resourcepath + "/";
        }
    //if (outputpath == "") {
    //    outputpath = resourcepath + "Output/";
    //    }
    if (*outputpath.rbegin() != '/') {
        outputpath = outputpath + "/";
        }

    string fragmentFileName = outputpath + "selected_fragments.dat";
    string probaFirstFragFileName = outputpath + "Proba_FirstFrag.dat";
    string probaTransitionFileName = outputpath + "Proba_Transition.dat";

    // Check whether probability computation has been started, stores number of completed fragments
    int firstFragCompleted = 0;
    int transitionCompleted = 0;
    if (FILE *FirstFragFile = fopen(probaFirstFragFileName.c_str(), "r")) {
        fclose(FirstFragFile);
        firstFragCompleted = SizeFile(probaFirstFragFileName.c_str());
        }
    
    if (FILE *TransitionFile = fopen(probaTransitionFileName.c_str(), "r")) {
        fclose(TransitionFile);
        transitionCompleted = SizeFile(probaTransitionFileName.c_str());
        }

    // Define some parameters.
    int numberFragments = SizeFile(fragmentFileName.c_str());
    vector<string> smilesString (numberFragments);
    vector<string> smartsAtomOnRight (numberFragments);
    vector<string> smartsAtomOnLeft (numberFragments);
    vector<string> nameFragment (numberFragments);
    vector<string> smartsmiddlePattern (numberFragments);
    vector<double> occurence (numberFragments);
    vector<vector<int> > connection (numberFragments);

    // Dynamically set size of connection matrix 
    for (int i = 0; i < numberFragments; i++) {
        connection[i].resize(numberFragments);
        }

    // pos sets the position in the probaTransition file that we start writing from. 
    // If we already have some lines completed, it will point to the end of those lines.
    // flag monitors whether Proba_FirstFrag.dat has been computed already
    long pos=0;
    int flag=0;

    if (firstFragCompleted != numberFragments) {
        //start from beginning
        flag = 0;
        pos = 0;
        }
    else if (transitionCompleted == 0) {
        //start from transitions
        flag = 1;
        pos = 0;
        }
    else if (transitionCompleted < numberFragments) {
        //start from transitionCompleted
        flag = 1;
        string line;
        ifstream probaTransitionFile(probaTransitionFileName.c_str());
        for (int i = 0; i < transitionCompleted; i++) {
            getline(probaTransitionFile, line);
            }
        pos = probaTransitionFile.tellg();
        }
    else if (transitionCompleted == numberFragments) {
        cout << "Files exist already!" << endl;
        return;
        }

    // Sets GUI progress bar value
    if (transitionCompleted) {
        cout << numberFragments*(transitionCompleted + 2) << endl;
        }

    // Open the fragment file, and check if it has been opened.
    ifstream fragmentFilestream(fragmentFileName.c_str());
    if(fragmentFilestream.is_open()) {
        for (int i=0 ; i < numberFragments ; i++) {
            fragmentFilestream >> nameFragment[i];
            }
        }
    else {
        cerr << "Can't find selected_fragments file" << endl;
        }

    // Read each fragment file
    for (int i = 0; i < numberFragments; i++) {
        // Smarts
        string fragFileName = resourcepath + "/" + nameFragment[i] + "/" + nameFragment[i] + ".smarts";
        ifstream fragFileStream(fragFileName.c_str());
        if (fragFileStream.is_open()) {
            fragFileStream >> smartsAtomOnRight[i];
            fragFileStream >> smartsAtomOnLeft[i];
            fragFileStream >> smartsmiddlePattern[i];
            }
        else {
            cerr << "Can't find smarts patterns for fragment " << nameFragment[i] << "!" << endl;
            }

        // Smiles
        string smiFileName = resourcepath + "/" + nameFragment[i] + "/" + nameFragment[i] + ".smi";
        ifstream smiFileStream(smiFileName.c_str());
        if (smiFileStream.is_open()) {
            smiFileStream >> smilesString[i];
            }
        else {
            cerr << "Can't find smiles pattern for fragment " << nameFragment[i] << "!" << endl;
            }
        }

    // Process with OpenBabel.
    OpenBabel::OBConversion obconversion;
    OpenBabel::OBMol mol;
    obconversion.SetInFormat("smi");

    if (flag == 0) {
        // Open the proba file for the first fragment, and check if it has been opened.
        ofstream probaFirstFragFile(probaFirstFragFileName.c_str(), ios::out);
        if(!probaFirstFragFile) {
            cerr << "ERROR: The file " << probaFirstFragFileName.c_str() << " is missing, or I can't read it." << endl;
            exit (-1);   
            }

        for (int i=0 ; i < numberFragments ; i++) {
            // Count how many times the pattern is found
            OpenBabel::OBSmartsPattern smarts;
            string rightPattern(smartsAtomOnRight[i]);
            smarts.Init(rightPattern.c_str());
            vector<vector<int> > maplist;
            int total = 0;
            bool notatend = obconversion.ReadFile(&mol,libraryFileName.c_str());
            while (notatend) {
                smarts.Match(mol);
                maplist = smarts.GetUMapList();
                total += maplist.size();
                mol.Clear();
                notatend = obconversion.Read(&mol);
                }

            // Check if the next pattern is a duplicate of the current one
            OpenBabel::OBSmartsPattern smartsDuplicate;
            smartsDuplicate.Init(rightPattern.c_str());
            int NumberDuplicates=1;
            for (int j=1 ; j < numberFragments-i ; j++) {
                string duplicate(smilesString[i+j]);
                OpenBabel::OBConversion obconversionDuplicate;
                obconversionDuplicate.SetInFormat("smi");
                OpenBabel::OBMol molDuplicate;
                obconversionDuplicate.ReadString(&molDuplicate, duplicate.c_str());
                // Compare the two patterns
                if (!smartsDuplicate.HasMatch(molDuplicate)) { j += numberFragments;  }
                else                                         { NumberDuplicates += 1; }
                molDuplicate.Clear();
                }
            for (int k=0 ; k < NumberDuplicates ; k++) { occurence[i+k] = 1.0*total/NumberDuplicates; }
            i += (NumberDuplicates-1);
            cout << i+1 << endl;
            }

        // Compute the probabilities.
        float SumOccurence=0;
        for (int i=0 ; i < numberFragments ; i++) { SumOccurence += occurence[i]; }
        probaFirstFragFile << fixed << setprecision(8) ;
        for (int i=0 ; i < numberFragments ; i++) { 
            probaFirstFragFile << occurence[i]/SumOccurence << endl;
            }
        }

    // Open the proba transition file, and check if it has been opened.
    ofstream probaTransitionFile;
    if (pos == 0) {                                                                                         // If the transition file needs to be created
        probaTransitionFile.open(probaTransitionFileName.c_str(), ios::out);
        }
    else {
        probaTransitionFile.open(probaTransitionFileName.c_str(), ios::in | ios::out | ios::binary);        // If it's already made
        }
    
    if(!probaTransitionFile) {
        cerr << "ERROR: The file " << probaTransitionFileName.c_str() << " is missing, or I can't read it." << endl;
        exit (-1);
        }
    
    probaTransitionFile << fixed << setprecision(8) ;
    if (pos != 0) {
        probaTransitionFile.seekp(pos);
        }

    for (int i=transitionCompleted; i < numberFragments ; i++) {
        // If we called the program using the Qt GUI, we need to check whether user wants to pause
        if (runType == 'q') {
            cout << "ping" << endl;
            cin >> runState;
            int runStateInt = atoi(runState.c_str());
            if (runStateInt == 1) {
                probaTransitionFile.close();
                cout << "Paused!" << endl;
                return;
                }
            }

        OpenBabel::OBSmartsPattern smarts;
        string rightPattern(smartsAtomOnRight[i]);
        int Sum=0;
        vector<double> probaTransition (numberFragments);
        mol.Clear();

        // Calculate the number of connections.
        for (int j=0 ; j < numberFragments ; j++) {
            string leftPattern(smartsAtomOnLeft[j]);
            string smartsPattern = rightPattern + leftPattern;
            smarts.Init(smartsPattern.c_str());
            vector<vector<int> > maplist;
            int total = 0;
            bool notatend = obconversion.ReadFile(&mol,libraryFileName.c_str());
            while (notatend) {
                smarts.Match(mol);
                maplist = smarts.GetUMapList();
                total += maplist.size();
                mol.Clear();
                notatend = obconversion.Read(&mol);
                }
            connection[i][j] = total;
            Sum += connection[i][j]; 
            }

        // Update progress bar
        int progress = numberFragments*(i+2);
        cout << progress << endl;

        // Compute the probabilities.
        for (int j=0 ; j < numberFragments ; j++) { 
            if(Sum != 0) { probaTransition[j] = 1.0*connection[i][j] / Sum; }
            else         { probaTransition[j] = 0.0;                        }
            probaTransitionFile << probaTransition[j] << "\t";
            }

        probaTransitionFile << endl;
        }

    cout << "Done" << endl;
}

/////////////////////////////////////////////////////////////////////////////////////

void computeProbabilitiesSplitted(string databasepath, string resourcepath, string outputpath, string counter)
{
    string libraryFileName = databasepath;
    if (*resourcepath.rbegin() != '/') {
        resourcepath = resourcepath + "/";
        }
    if (*outputpath.rbegin() != '/') {
        outputpath = outputpath + "/";
        }
    
    string fragmentFileName = outputpath + "selected_fragments.dat";
    string probaTransitionFileName = outputpath + "Proba_Transition" + "_" + counter + ".dat";;

    // Define some parameters.
    int numberFragments = SizeFile(fragmentFileName.c_str());
    vector<string> smilesString (numberFragments);
    vector<string> smartsAtomOnRight (numberFragments);
    vector<string> smartsAtomOnLeft (numberFragments);
    vector<string> nameFragment (numberFragments);
    vector<string> smartsmiddlePattern (numberFragments);
    vector<double> occurence (numberFragments);
    vector<vector<int> > connection (numberFragments);

    // Dynamically set size of connection matrix 
    for (int i = 0; i < numberFragments; i++) {
        connection[i].resize(numberFragments);
        }

    // Open the fragment file, and check if it has been opened.
    ifstream fragmentFilestream(fragmentFileName.c_str());
    if(fragmentFilestream.is_open()) {
        for (int i=0 ; i < numberFragments ; i++) {
            fragmentFilestream >> nameFragment[i];
            }
        }
    else {
        cerr << "Can't find selected_fragments file" << endl;
        }

    // Read each fragment file
    for (int i = 0; i < numberFragments; i++) {
        // Smarts
        string fragFileName = resourcepath + "/" + nameFragment[i] + "/" + nameFragment[i] + ".smarts";
        ifstream fragFileStream(fragFileName.c_str());
        if (fragFileStream.is_open()) {
            fragFileStream >> smartsAtomOnRight[i];
            fragFileStream >> smartsAtomOnLeft[i];
            fragFileStream >> smartsmiddlePattern[i];
            }
        else {
            cerr << "Can't find smarts patterns for fragment " << nameFragment[i] << "!" << endl;
            }

        // Smiles
        string smiFileName = resourcepath + "/" + nameFragment[i] + "/" + nameFragment[i] + ".smi";
        ifstream smiFileStream(smiFileName.c_str());
        if (smiFileStream.is_open()) {
            smiFileStream >> smilesString[i];
            }
        else {
            cerr << "Can't find smiles pattern for fragment " << nameFragment[i] << "!" << endl;
            }
        }

    // Process with OpenBabel.
    OpenBabel::OBConversion obconversion;
    OpenBabel::OBMol mol;
    obconversion.SetInFormat("smi");

    // Open the proba transition file, and check if it has been opened.
    ofstream probaTransitionFile;
    probaTransitionFile.open(probaTransitionFileName.c_str(), ios::out);
    if(!probaTransitionFile) {
        cerr << "ERROR: The file " << probaTransitionFileName.c_str() << " is missing, or I can't read it." << endl;
        exit (-1);
        }

    probaTransitionFile << fixed << setprecision(8) ;
    int inputValue = atoi(counter.c_str());
    for (int i=inputValue ; i<inputValue+1 ; i++) {
        OpenBabel::OBSmartsPattern smarts;
        string rightPattern(smartsAtomOnRight[i]);
        int Sum=0;
        vector<double> probaTransition (numberFragments);
        mol.Clear();

        // Calculate the number of connections.
        cout << "Line " << i << endl ;
        for (int j=0 ; j < numberFragments ; j++) {
            string leftPattern(smartsAtomOnLeft[j]);
            string smartsPattern = rightPattern + "-" + leftPattern;
            smarts.Init(smartsPattern.c_str());
            vector<vector<int> > maplist;
            int total = 0;
            bool notatend = obconversion.ReadFile(&mol,libraryFileName.c_str());
            while (notatend) {
                smarts.Match(mol);
                maplist = smarts.GetUMapList();
                total += maplist.size();
                mol.Clear();
                notatend = obconversion.Read(&mol);
                }
            connection[i][j] = total;
            Sum += connection[i][j];
            cout << connection[i][j] << endl ;
            }

        // Compute the probabilities.
        for (int j=0 ; j < numberFragments ; j++) { 
            if(Sum != 0) { probaTransition[j] = 1.0*connection[i][j] / Sum; }
            else     { probaTransition[j] = 0.0;                            }
            probaTransitionFile << probaTransition[j] << "\t";
            }

        probaTransitionFile << endl;
        }

    cout << "Done" << endl;
}

/////////////////////////////////////////////////////////////////////////////////////

// Searches for 3mers that are never found in the given database
void findThreemers(string databasepath, string resourcepath, string outputpath, char runType)
{
    if (*resourcepath.rbegin() != '/') {
        resourcepath = resourcepath + "/";
        }
    if (*outputpath.rbegin() != '/') {
        outputpath = outputpath + "/";
        }

    // OpenBabel Setup
    OpenBabel::OBConversion obconversion;
    OpenBabel::OBMol mol;
    OpenBabel::OBSmartsPattern smarts;
    obconversion.SetInFormat("smi");

    int pauseflag = 0;

    string libraryFileName = databasepath;
    string fragmentFileName = outputpath + "selected_fragments.dat";
    string probaTransitionFileName = outputpath + "Proba_Transition.dat";
    string threemersleftname = outputpath + "3mers_left.dat";
    string threemerfile = outputpath + "unfound_3mers.dat";
    string threemersToSearch = outputpath + "3mers_to_search.dat";
    
    if (FILE *threemersleft = fopen(threemersleftname.c_str(), "r")) {
        pauseflag = 1;
        fclose(threemersleft);
        }

    // Stores list of possible 3mers
    list<string> threemers;

    // If we need to start from the beginning, we must compute the list of 3mers to be searched
    if (pauseflag == 0) {
        ifstream threemerstream(threemersToSearch.c_str());
        if (threemerstream) {
            string smartsPattern;
            while (threemerstream >> smartsPattern) {
                threemers.push_back(smartsPattern);
                }
            }
        else {
            // Define some parameters.
            int numberFragments = SizeFile(fragmentFileName.c_str());
            vector<string> smartsAtomOnRight (numberFragments);
            vector<string> smartsAtomOnLeft (numberFragments);
            vector<string> nameFragment (numberFragments);
            vector<string> smartsMiddleFragment (numberFragments);

            // Open the fragment file, and check if it has been opened.
            ifstream fragmentFilestream(fragmentFileName.c_str());
            if(fragmentFilestream.is_open()) {
                for (int i=0 ; i < numberFragments ; i++) {
                    fragmentFilestream >> nameFragment[i];
                    }
                }

            // Read each fragment file
            for (int i = 0; i < numberFragments; i++) {
                // Smarts
                string fragFileName = resourcepath + "/" + nameFragment[i] + "/" + nameFragment[i] + ".smarts";
                ifstream fragFileStream(fragFileName.c_str());
                if (fragFileStream.is_open()) {
                    fragFileStream >> smartsAtomOnRight[i];
                    fragFileStream >> smartsAtomOnLeft[i];
                    fragFileStream >> smartsMiddleFragment[i];
                    }
                }
   
            // Read from transition probability file
            vector<vector<float> > data (numberFragments);
            for (int i = 0; i < numberFragments; i++) {
                data[i].resize(numberFragments);
                }
            ifstream file;
            file.open(probaTransitionFileName.c_str());
            for (int i = 0; i < numberFragments; ++i) {
                for (int j = 0; j < numberFragments; ++j) {
                    file >> data[i][j];
                    }
                }

            file.close();

            float AB = 0;
            std::size_t found = 0;

            ofstream toSearch(threemersToSearch.c_str(), ios::out); 
            if (toSearch) {
                for (int i=0 ; i < numberFragments ; i++) {
                    string rightPattern(smartsAtomOnRight[i]);

                    // Calculate the number of connections.
                    for (int j=0 ; j < numberFragments ; j++) {

                        string middlePattern(smartsMiddleFragment[j]);
                        int numberofSpots = std::count(middlePattern.begin(), middlePattern.end(), 'y');

                            for (int k=0; k < numberFragments; k++) {
                                string leftPattern(smartsAtomOnLeft[k]);
                                for (int j2=0; j2 < numberofSpots; j2++) {
                                    string newmiddlePattern = middlePattern;

                                    for (int j3 = 0; j3 < numberofSpots; j3++) {
                                        found = newmiddlePattern.find('y');
                                        if (j3 == j2) {
                                            newmiddlePattern.replace(found, 1, leftPattern);
                                            }
                                        else {
                                            newmiddlePattern.replace(found-1, 3, "");
                                            }
                                        }

                                    AB = data[i][j];
                                    string smartsPattern = rightPattern + "-" + newmiddlePattern;

                                    // Checks whether transition A-B exists; if not, we don't need to search this set of 3mers.
                                    if (AB) {
                                        threemers.push_back(smartsPattern);
                                        toSearch << smartsPattern << endl;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
    // If we are resuming from a pause, we can load the fragments left to search
    else {
        ifstream leftstream;
        leftstream.open(threemersleftname.c_str());
        string smartsPattern;

        while (leftstream >> smartsPattern) {
            threemers.push_back(smartsPattern);
            }

        if(remove(threemersleftname.c_str()) != 0) {
            cout << "Can't delete 3mers_left file" << endl;
            }
        }

    // Tell GUI total number of threemers, to set progress bar maximum
    int numberThreemers = threemers.size();
    cout << "#" + to_string(static_cast<long long>(numberThreemers)) << endl;

    ofstream out(threemerfile.c_str(), ios::out | ios::app);
    if (out) {
        //Process threemers with openbabel
        list<string>::iterator it = threemers.begin();

        int total = 0;
        vector<vector<int> > maplist;
        string runState;

        while (it != threemers.end()) {
            if (runType == 'q') {
                cout << "ping" << endl;    
                cin >> runState;
                int runStateInt = atoi(runState.c_str());
                if (runStateInt == 1) {
                    ofstream leftstream;
                    leftstream.open(threemersleftname.c_str());
                    while (it != threemers.end()) {
                        string smartsPattern = *it;
                        leftstream << smartsPattern << endl;
                        it++;
                        }
                    leftstream.close();
                    cout << "Paused!" << endl;
                    return;
                    }
                }

            string smartsPattern = *it;
            smarts.Init(smartsPattern.c_str());
            bool notatend = obconversion.ReadFile(&mol,libraryFileName.c_str());
                while (notatend) {
                    smarts.Match(mol);
                    maplist = smarts.GetUMapList();
                    total += maplist.size();
                    mol.Clear();
                    notatend = obconversion.Read(&mol);
                    }

            if (total != 0) {
                }
            else {
                out << smartsPattern << endl;
                }

            list<string>::iterator start = threemers.begin();
            cout << distance(start, it) << endl;

            total = 0;
            it++;
            }
        }

    cout << "Done" << endl;
    ifstream file(threemersleftname.c_str()); 
    if (file.is_open()) {
        file.close();
        remove(threemersleftname.c_str());
        }
}

