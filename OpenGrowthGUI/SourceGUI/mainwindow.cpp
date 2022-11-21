#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QToolButton>
#include <qapplication.h>
#include <QProgressBar>
#include <QProcess>
#include <QMultiMap>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <QFile>
#include <QtGui>
#include <ctime>

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this); // Sets up UI from designer

    numNoRings = 0;
    numOneRings = 0;
    numTwoRings = 0;

    ui->progressBar->setFixedWidth(200);
    ui->runButton->setFixedSize(90,50);
    ui->createFragmentFileButton->setFixedSize(100,50);
    ui->threemerlistpushButton->setFixedSize(80,50);
    ui->threemerpushButton->setFixedSize(70,50);
    ui->killButton->setFixedSize(60,50);
}

MainWindow::~MainWindow()
{
    delete ui;
}

// Select "Resource Directory" button, allows user to select directory containing fragments, icons, etc.
void MainWindow::on_resourcedirpushButton_clicked()
{
    QString directory = QFileDialog::getExistingDirectory(this,                                                   // Prompt for file entry
                                                          tr("Open Directory"),                                   // Prompt Text
                                                          "C://Users/",                                           // Default Directory Shown
                                                          QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks
                                                          );
    // If user cancels out of database selection, database in textedit is not changed
    if (directory != "")
    {
        ui->resourcedirlineEdit->setText(directory);
        if (ui->databaselineEdit->text() == "")
            ui->progressLabel->setText("Please choose output directory or database");
        else ui->progressLabel->setText("Press update to load fragments");
    }
}

// Select "Output directory" button, allows user to select output directory
void MainWindow::on_outputdirpushButton_clicked()
{
    QString directory = QFileDialog::getExistingDirectory(this,                                                   // Prompt for file entry
                                                          tr("Open Directory"),                                        // Prompt Text
                                                          "C://Users/",                                                // Default Directory Shown
                                                          QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks
                                                          );
    if (directory != "")
    {
        ui->outputdirlineEdit->setText(directory);
        if (ui->databaselineEdit->text() == "")
            ui->progressLabel->setText("Please choose database");
    }
}

// Select "Database" button, allows user to select database file
void MainWindow::on_selectDatabase_clicked()
{
    QString filename = QFileDialog::getOpenFileName(this,                                                   // Prompt for file entry
                                                    tr("Open File"),                                        // Prompt Text
                                                    "C://Users/",                                           // Default Directory Shown            ***************CHANGE WHEN SWITCHING COMPUTER*****************
                                                    "All files (*.*);;Smiles File (*.smi)"                  // Filter View Options
                                                    );
    if (filename != "")                                                                                     // If user cancels out of database selection, database in textedit is not changed
    {
        ui->databaselineEdit->setText(filename);                                                            // Sets database name on UI
        if (ui->resourcedirlineEdit->text() == "")
            ui->progressLabel->setText("Please set resource directory");
        else ui->progressLabel->setText("Press update to load fragments");
    }
}

// Select "Update" for the database
void MainWindow::on_updateDatabase_clicked()
{
    if (ui->resourcedirlineEdit->text() != "")
    {
        if (ui->databaselineEdit->text() != "")
        {
            // Deletes old buttons from GUI
            deleteButtons();
            // Reset number of rings and ring map
            numOneRings = 0;
            ring1map.clear();
            numTwoRings = 0;
            ring2map.clear();

            // Reads fragments from fragment folder in resource directory
            loadFragments();
            // Sets progress bar min/max, based on number of fragments (global numOneRings and numTwoRings are set in loadFragments)
            ui->progressBar->setMinimum(1);
            ui->progressBar->setMaximum(numOneRings+numTwoRings);

            fragmentNo = 1;

            // Searches for fragments in database
            findFragments();
        }else ui->progressLabel->setText("No database selected");
    }else ui->progressLabel->setText("No resource directory selected");
}

void MainWindow::deleteButtons()
{
    for(int i = 0; i < numOneRings; i++){
        // Find each button by its name (assigned when making the button)
        QToolButton *button = ui->oneringScrollArea->findChild<QToolButton *>(ring1[i].second);
        delete button;
    }

    for(int i = 0; i < numTwoRings; i++){
        QToolButton *button = ui->tworingScrollArea->findChild<QToolButton *>(ring2[i].second);
        delete button;
    }
}

void MainWindow::loadFragments()
{
    QString qpath = ui->resourcedirlineEdit->text() + "/";
    QFile fragments(qpath + "Fragments.dat");

    if (fragments.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QTextStream stream(&fragments);
        QString fragment;
        while (!stream.atEnd())
        {
            fragment = stream.readLine();
            int   lastindex = fragment.lastIndexOf("_");
            int   patternNumber = fragment.at(lastindex+1).digitValue();
            QChar testEndOfLine = fragment.at(lastindex+2);

            // Only load the smarts pattern for the fragment X_X_1 and check if the name stops (to avoid X_X_10)
            if (patternNumber==1 && testEndOfLine=='\0')
            {
                int ringNumber = fragment.at(4).digitValue();
                QString ring = fragment.mid(0,lastindex);

                QFile smarts(qpath + "/" + fragment + "/" + fragment + ".smarts");
                if (smarts.open(QIODevice::ReadOnly | QIODevice::Text))
                {
                    QTextStream smartsStream(&smarts);
                    QString pattern;
                    smartsStream >> pattern;

                    switch (ringNumber)
                    {
                    case 0:
                        ring0 << QPair<QString, QString>(pattern, ring);
                        numNoRings++;
                        break;
                    case 1:
                        ring1 << QPair<QString, QString>(pattern, ring);
                        numOneRings++;
                        break;
                    case 2:
                        ring2 << QPair<QString, QString>(pattern, ring);
                        numTwoRings++;
                        break;
                    default:    break;
                    }
                }else ui->progressLabel->setText("Can't find smarts pattern for fragment: " + fragment);
            }
        }

    }else{
        ui->progressLabel->setText("Can't load fragment file, check resource directory");
    }
}


// findFragments starts the search for fragment matches in the database
// searchForFragment is the function that starts the lookup process. On completion of the first lookup, the loadprocesslog slot starts the lookup for the next fragment.
// This way, the lookup QProcesses can occur sequentially without a waitForFinished() or a waiting QEventLoop
void MainWindow::findFragments()
{
    ringType = 1;
    string firstFragment = ring1[0].first.toUtf8().constData();
    searchForFragment(firstFragment);
}

// Counts number of times that each fragment is found in the selected database. Uses external program SearchGUI.exe.
void MainWindow::searchForFragment(std::string smartsString)
{
    // Sets progress label to display the molecule that is being calculated
    ui->progressLabel->setText("Scanning fragment: " + QString::number(fragmentNo) + " out of " + QString::number(numOneRings+numTwoRings));

    // Sets the progress bar to a value corresponding to molecule#
    ui->progressBar->setValue(fragmentNo);

    QString program = "./SearchGUI.exe";
    QString database = ui->databaselineEdit->text();
    QStringList arguments;
    QString fragment = QString::fromStdString(smartsString);

    arguments << "f" << database << fragment << "q";

    myProcess = new QProcess(this);
    connect(myProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(loadprocesslog()));
    connect(myProcess, SIGNAL(readyReadStandardError()), this, SLOT(loadprocesserr()));
    connect(myProcess, SIGNAL(finished(int)), this, SLOT(killProcess()));
    myProcess->start(program, arguments);
    if (!myProcess->pid()){
        ui->progressLabel->setText("Can't find and/or start SearchGUI.exe. It must be located in the same folder as OpenGrowthGUI.exe.");
    }
}

// Slot for searchForFragment stdout
void MainWindow::loadprocesslog()
{
    QString matches = myProcess->readAllStandardOutput();
    int nummatches = matches.toInt();

    if (ringType == 1)
    {
        // Add number of matches, fragment name to ring1map
        ring1map.insert(nummatches, ring1[fragmentNo-1].second);
    }else if (ringType == 2){
        ring2map.insert(nummatches, ring2[fragmentNo-numOneRings-1].second);
    }
}

// Slot for searchForFragment stderr
void MainWindow::loadprocesserr()
{
    QString string = myProcess->readAllStandardError();
    ui->progressLabel->setText(string);
}

// Slot to kill process when finished
void MainWindow::killProcess()
{
    string nextFragment;

    if (ringType == 1)
    {

        // If there are still fragments to be searched
        if (fragmentNo < (numOneRings+numTwoRings))
        {
            // As long as the next fragment is still a one ring fragment
            if (fragmentNo != numOneRings)
            {
                nextFragment = ring1[fragmentNo].first.toUtf8().constData();
            }else{
                ringType = 2;
                nextFragment = ring2[fragmentNo-numOneRings].first.toUtf8().constData();
            }

            fragmentNo++;

            // Search for the next fragment
            searchForFragment(nextFragment);
        }
    }else if (ringType == 2){
        if (fragmentNo < (numOneRings+numTwoRings))
        {
            nextFragment = ring2[fragmentNo-numOneRings].first.toUtf8().constData();
            fragmentNo++;
            searchForFragment(nextFragment);
        }else{
            // Adds fragment buttons to GUI
            showButtons();
            ui->progressLabel->setText("Ready");
        }
    }
}

void MainWindow::showButtons()                                                      // After database is updated, this function displays the fragment buttons that can be selected
{                                                                                   // Fragment buttons ordered by # of times they appear in database
    // Show buttons for onering
    int columnCounter = 1;
    int rowCounter = 1;
    QString iconPath = ui->resourcedirlineEdit->text() + "/";

    // ring1map<int matches, string ring>
    QMultiMap<int, QString>::const_iterator rit = ring1map.constEnd();
    while (rit != ring1map.constBegin()){

        --rit;

        // If we have 5 in a row, switch to the next row
        if (columnCounter > 5)
        {
            rowCounter++;
            columnCounter = 1;
        }

        QToolButton *button=new QToolButton();                                                      // Creates new button

        button->setCheckable(true);                                                                 // On/off states
        button->setObjectName(rit.value());                                                         // Gives the button a name, eg Ring1_12
        button->setMinimumWidth(115);
        button->setMinimumHeight(125);
        button->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);                      // Sets button size policy, to make sure it expands with a window expand
        QFile iconfile(iconPath + rit.value() + "_1/" + rit.value() + ".svg");
        if(iconfile.exists())
        {
            button->setIcon(QIcon(iconPath + rit.value() + "_1/" + rit.value() + ".svg"));
            button->setIconSize(QSize(115,100));                                                    // Set fixed icon size
            button->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);                                // Sets tool button text underneath icon. Text is number of times found in database.
        }else{
            QPixmap pm(100,100);
            pm.fill(Qt::transparent);
            QPainter *painter = new QPainter(&pm);
            painter->drawText(QPoint(20,40), rit.value());
            delete painter;
            button->setIcon(QIcon(pm));
            button->setIconSize(QSize(100,100));
            button->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
        }
        button->setText(QString::number(rit.key()));                                               // first value of ring1map is number of times the molecule was found in database
        ui->oneringButtonGrid->addWidget(button, rowCounter, columnCounter);                       // Adds button to grid layout, addWidget(widget, row, column)
        columnCounter++;
    }

    // Show buttons for tworing
    columnCounter = 1;
    rowCounter = 1;

    QMultiMap<int, QString>::const_iterator rit2 = ring2map.constEnd();

    while (rit2 != ring2map.constBegin()){

        --rit2;

        // If we have 3 in a row, switch to the next row
        if (columnCounter > 3)
        {
            rowCounter++;
            columnCounter = 1;
        }

        QToolButton *button=new QToolButton();
        ui->tworingButtonGrid->addWidget(button, rowCounter, columnCounter);
        button->setCheckable(true);
        button->setObjectName(rit2.value());
        button->setMinimumWidth(150);
        button->setMinimumHeight(200);
        button->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        QFile iconfile(iconPath + rit2.value() + "_1/" + rit2.value() + ".svg");
        if(iconfile.exists())
        {
            button->setIcon(QIcon(iconPath + rit2.value() +"_1/" + rit2.value() + ".svg"));
            button->setIconSize(QSize(150,150));
            button->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
        }else{
            QPixmap pm(180,180);
            pm.fill(Qt::transparent);
            QPainter *painter = new QPainter(&pm);
            painter->drawText(QPoint(40,90), rit2.value());
            delete painter;
            button->setIcon(QIcon(pm));
            button->setIconSize(QSize(180,180));
            button->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
        }
        button->setText(QString::number(rit2.key()));
        columnCounter++;
    }
}

// Creates file selected_fragments in output folder
void MainWindow::on_createFragmentFileButton_clicked()
{
    if (ui->resourcedirlineEdit->text() != "")
    {
        if (ui->databaselineEdit->text() != "")
        {
            if (numNoRings || numOneRings || numTwoRings)
            {
                QString qpath;
                if (ui->outputdirlineEdit->text() != "")
                    qpath = ui->outputdirlineEdit->text() + "/";
                else qpath = ui->resourcedirlineEdit->text() + "/../Output/";

                QString resourcefolder = ui->resourcedirlineEdit->text() + "/";

                if (!QDir(qpath).exists())
                {
                    QDir().mkdir(qpath);
                }

                QFile file(qpath + "selected_fragments.dat");
                QFile xyzfile(qpath + "Fragments.dat");
                if (file.open(QIODevice::WriteOnly | QIODevice::Text))
                {
                    if (xyzfile.open(QIODevice::WriteOnly | QIODevice::Text))
                    {
                        QTextStream filestream(&file);
                        QTextStream xyzstream(&xyzfile);

                        QFile fragmentfile(resourcefolder + "Fragments.dat");
                        if (fragmentfile.open(QIODevice::ReadOnly | QIODevice::Text))
                        {
                            QTextStream textstream(&fragmentfile);
                            while (!textstream.atEnd())
                            {
                                QString line = textstream.readLine();

                                int ringNumber = line.at(4).digitValue();

                                int lastindex;
                                QString subString;
                                QToolButton *button;

                                //if (!QDir(qpath + "Fragments").exists())
                                //    QDir().mkdir(qpath+"Fragments");

                                switch (ringNumber)
                                {
                                case 0:
                                    filestream << line << endl;
                                    xyzstream << "./" + line + ".xyz" << endl;
                                    QFile::copy(resourcefolder + "/" + line + "/" + line + ".xyz", qpath + "/" + line + ".xyz");
                                    break;
                                case 1:
                                    lastindex = line.lastIndexOf("_");
                                    subString = line.mid(0, lastindex);
                                    button = ui->oneringScrollArea->findChild<QToolButton *>(subString);
                                    if ((button) && button->isChecked()){                                                                                 // Checks if the button was found, and if the button has been clicked
                                        filestream << line << endl;
                                        xyzstream << "./" + line + ".xyz" << endl;
                                        QFile::copy(resourcefolder + "/" + line + "/" + line + ".xyz", qpath + "/" + line + ".xyz");
                                    }
                                    break;
                                case 2:
                                    lastindex = line.lastIndexOf("_");
                                    subString = line.mid(0, lastindex);
                                    button = ui->tworingScrollArea->findChild<QToolButton *>(subString);
                                    if ((button) && button->isChecked()){                                                                                 // Checks if the button was found, and if the button has been clicked
                                        filestream << line << endl;
                                        xyzstream << "./" + line + ".xyz" << endl;
                                        QFile::copy(resourcefolder + "/" + line + "/" + line + ".xyz", qpath + "/" + line + ".xyz");
                                    }
                                    break;
                                default:    break;

                                }
                            }
                            fragmentfile.close();
                        }else{
                            ui->progressLabel->setText("Can't find fragment file: Fragments.dat");
                            return;
                        }
                        xyzfile.close();
                    }
                    file.close();
                }
                ui->progressLabel->setText("Fragment Files Created");
            }
        }else ui->progressLabel->setText("No database selected");
    }else ui->progressLabel->setText("No resource directory selected");
}


void MainWindow::on_runButton_clicked()
{
    if (ui->resourcedirlineEdit->text() != "")
    {
        if (ui->databaselineEdit->text() != "")
        {
            if (numNoRings || numOneRings || numTwoRings)
            {
                // Starts program to find firstfrag and transition probabilities
                generate();
            }
        }else ui->progressLabel->setText("No database selected");
    }else ui->progressLabel->setText("No resource directory selected");
}


void MainWindow::generate()
{
    flag = 0;
    ui->progressLabel->setText("Working...\n");

    QString program = "./SearchGUI.exe";
    QStringList arguments;
    QString database = ui->databaselineEdit->text();
    QString resourcefolder = ui->resourcedirlineEdit->text();
    QString outputfolder;
    if (ui->outputdirlineEdit->text() != "")
        outputfolder = ui->outputdirlineEdit->text() + "/";
    else outputfolder = ui->resourcedirlineEdit->text() + "/../Output/";

    string filename = outputfolder.toStdString() + "selected_fragments.dat";
    int numberFragments = sizeFile(filename.c_str());
    ui->progressBar->setMinimum(1);
    ui->progressBar->setValue(1);
    ui->progressBar->setMaximum(numberFragments*(numberFragments + 1));

    arguments  << "p" << database << resourcefolder << outputfolder << "q";

    myProcess = new QProcess(this);
    connect(myProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(processlog()));
    connect(myProcess, SIGNAL(readyReadStandardError()), this, SLOT(processerr()));
    myProcess->start(program, arguments);
}

void MainWindow::processlog()
{
    QString string = myProcess->readAllStandardOutput();
    qDebug() << string;
    if (string.contains("ping")){
        QByteArray flagdata = QByteArray::number(flag);
        flagdata.append("\n");
        myProcess->write(flagdata);
    }else if (string.contains("!"))
    {
        ui->progressLabel->setText(string);
    }else if (string.contains("#"))
    {
        string.remove(0,1);
        int number = string.toInt();
        ui->progressBar->setMaximum(number);
    }else{
        int number = string.toInt();
        ui->progressBar->setValue(number);
    }
    if (string.contains("Done"))
    {
        ui->progressLabel->setText("Done!");
        int number = ui->progressBar->maximum();
        ui->progressBar->setValue(number);
    }
}


void MainWindow::processerr()
{
    QString string = myProcess->readAllStandardError();
    ui->progressLabel->setText(string);
}

void MainWindow::findThreemers()
{
    flag = 0;
    ui->progressLabel->setText("Working...\n");

    QString program = "./SearchGUI.exe";
    QStringList arguments;
    QString database = ui->databaselineEdit->text();
    QString resourcefolder = ui->resourcedirlineEdit->text();
    QString outputfolder;
    if (ui->outputdirlineEdit->text() != "")
        outputfolder = ui->outputdirlineEdit->text() + "/";
    else outputfolder = ui->resourcedirlineEdit->text() + "/../Output/";

    arguments  << "t" << database << resourcefolder << outputfolder << "q";

    myProcess = new QProcess(this);
    connect(myProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(processlog()));
    connect(myProcess, SIGNAL(readyReadStandardError()), this, SLOT(processerr()));
    myProcess->start(program, arguments);
}

void MainWindow::delay()
{
    QTime dieTime= QTime::currentTime().addSecs(1);
    while( QTime::currentTime() < dieTime )
        QCoreApplication::processEvents(QEventLoop::AllEvents);
}

int MainWindow::sizeFile(char const fileName[])
{
    // Open the file and check it can be opened.
    ifstream listFile(fileName, ios::in);
    if(!listFile) {
        ui->progressLabel->setText("ERROR: The file is missing, or I can't read it.\n");
        exit (-1);
    }
    // If a read line is not empty, increase the counter
    int fragmentListSize=0;
    string line;
    while(getline(listFile, line)) {
        if(strcmp(line.c_str(),"")) { fragmentListSize++; }
    }
    listFile.close();
    return fragmentListSize;
}

void MainWindow::on_threemerpushButton_clicked()
{
    if (ui->resourcedirlineEdit->text() != "")
    {
        if (ui->databaselineEdit->text() != "")
        {
            if (numNoRings || numOneRings || numTwoRings)
            {
                // Starts program to find unfound threemers
                findThreemers();
            }
        }else ui->progressLabel->setText("No database selected");
    }else ui->progressLabel->setText("No resource directory selected");
}

void MainWindow::on_killButton_clicked()
{
    if(myProcess->pid())
    {
        flag = 1;
    }
}


void MainWindow::on_showIconButton_clicked()
{
    if (ui->showIconButton->isChecked())
    {
        QString iconPath = ui->resourcedirlineEdit->text() + "/" ;

        for(int i = 0; i < numOneRings; i++){
            // Find each button by its name (assigned when making the button)
            QToolButton *button = ui->oneringScrollArea->findChild<QToolButton *>(ring1[i].second);

            QFile iconfile(iconPath + ring1[i].second + "_1/" + ring1[i].second + ".svg");
            if(iconfile.exists())
            {
                button->setIcon(QIcon(iconPath + ring1[i].second + "_1/" + ring1[i].second + ".svg"));
            }else{
                QPixmap pm(80,80);
                pm.fill(Qt::transparent);
                QPainter *painter = new QPainter(&pm);
                painter->drawText(QPoint(20,40), ring1[i].second);
                delete painter;
                button->setIcon(QIcon(pm));
            }
        }

        for(int i = 0; i < numTwoRings; i++){
            QToolButton *button = ui->tworingScrollArea->findChild<QToolButton *>(ring2[i].second);

            QFile iconfile(iconPath + ring2[i].second + "_1/" + ring2[i].second + ".svg");
            if(iconfile.exists())
            {
                button->setIcon(QIcon(iconPath + ring2[i].second +"_1/" + ring2[i].second + ".svg"));
            }else{
                QPixmap pm(180,180);
                pm.fill(Qt::transparent);
                QPainter *painter = new QPainter(&pm);
                painter->drawText(QPoint(40,90), ring2[i].second);
                delete painter;
                button->setIcon(QIcon(pm));
            }

        }
    }else{
        for(int i = 0; i < numOneRings; i++){
            // Find each button by its name (assigned when making the button)
            QToolButton *button = ui->oneringScrollArea->findChild<QToolButton *>(ring1[i].second);
            QPixmap pm(80,80);
            pm.fill(Qt::transparent);
            QPainter *painter = new QPainter(&pm);
            painter->drawText(QPoint(20,40), ring1[i].second);
            delete painter;
            button->setIcon(QIcon(pm));
        }

        for(int i = 0; i < numTwoRings; i++){
            QToolButton *button = ui->tworingScrollArea->findChild<QToolButton *>(ring2[i].second);
            QPixmap pm(180,180);
            pm.fill(Qt::transparent);
            QPainter *painter = new QPainter(&pm);
            painter->drawText(QPoint(40,90), ring2[i].second);
            delete painter;
            button->setIcon(QIcon(pm));
        }
    }
}


void MainWindow::on_threemerlistpushButton_clicked()
{
    if (ui->resourcedirlineEdit->text() != "")
    {
        if (ui->databaselineEdit->text() != "")
        {
            if (numNoRings || numOneRings || numTwoRings)
            {
                QString qresourcepath = ui->resourcedirlineEdit->text() + "/";
                QString qoutputpath = ui->outputdirlineEdit->text();

                string resourcepath = qresourcepath.toStdString();
                string outputpath = qoutputpath.toStdString();

                if (outputpath == "")
                {
                    outputpath = resourcepath + "/../Output/";
                }else outputpath = outputpath + "/";

                string fragmentFileName = outputpath + "selected_fragments.dat";
                string probaTransitionFileName = outputpath + "Proba_Transition.dat";
                string threemersToSearch = outputpath + "3mers_to_search.dat";

                // Define some parameters.
                int numberFragments = sizeFile(fragmentFileName.c_str());
                vector<string> smartsAtomOnRight (numberFragments);
                vector<string> smartsAtomOnLeft (numberFragments);
                vector<string> nameFragment (numberFragments);
                vector<string> smartsMiddleFragment (numberFragments);

                // Open the fragment file, and check if it has been opened.
                ifstream fragmentFilestream(fragmentFileName.c_str());
                if(fragmentFilestream.is_open())
                {
                    for (int i=0 ; i < numberFragments ; i++)
                    {
                        fragmentFilestream >> nameFragment[i];
                    }
                }else{
                    ui->progressLabel->setText("Cant find fragment file!");
                    return;
                }

                // Read each fragment file
                for (int i = 0; i < numberFragments; i++)
                {
                    // Smarts
                    string fragFileName = resourcepath + "/" + nameFragment[i] + "/" + nameFragment[i] + ".smarts";
                    ifstream fragFileStream(fragFileName.c_str());
                    if (fragFileStream.is_open())
                    {
                        fragFileStream >> smartsAtomOnRight[i];
                        fragFileStream >> smartsAtomOnLeft[i];
                        fragFileStream >> smartsMiddleFragment[i];
                    }
                }

                // Read from transition probability file
                //float data[numberFragments][numberFragments];
                vector<vector<float> > data (numberFragments);
                for (int i = 0; i < numberFragments; i++)
                    data[i].resize(numberFragments);

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
                if (toSearch)
                {
                    // Add patterns that we know we don't want such as aketals.
                    toSearch << "[O;X2;!$(OC(=O));!$(OP(=O))]-[C;X4]-[O;X2;!$(OC(=O));!$(OP(=O))]"         << endl;
                    toSearch << "[O;X2;!$(OC(=O));!$(OP(=O))]-[C;X4]-[N;X3;!$(NS(=O)(=O));!$(NC(=O))]"     << endl;
                    toSearch << "[N;X3;!$(NS(=O)(=O));!$(NC(=O))]-[C;X4]-[N;X3;!$(NS(=O)(=O));!$(NC(=O))]" << endl;
                    // Then work with the other patterns.
                    for (int i=0 ; i < numberFragments ; i++)
                    {
                        string rightPattern(smartsAtomOnRight[i]);

                        // Calculate the number of connections.
                        for (int j=0 ; j < numberFragments ; j++)
                        {
                            string middlePattern(smartsMiddleFragment[j]);
                            int numberofSpots = std::count(middlePattern.begin(), middlePattern.end(), 'y');
                            for (int k=0; k < numberFragments; k++)
                            {
                                string leftPattern(smartsAtomOnLeft[k]);
                                for (int j2=0; j2 < numberofSpots; j2++)
                                {
                                    string newmiddlePattern = middlePattern;
                                    for (int j3 = 0; j3 < numberofSpots; j3++)
                                    {
                                        found = newmiddlePattern.find('y');
                                        if (j3 == j2)
                                        {
                                            newmiddlePattern.replace(found, 1, leftPattern);
                                        }else{
                                            newmiddlePattern.replace(found-1, 3, "");
                                        }
                                    }
                                    AB = data[i][j];
                                    string smartsPattern = rightPattern + "-" + newmiddlePattern;
                                    if (AB)
                                    {
                                        toSearch << smartsPattern << endl;
                                    }
                                }
                            }
                        }
                    }
                    if (toSearch)
                    {
                        ui->progressLabel->setText("3mer File Created!");
                    }
                }
            }
        }else ui->progressLabel->setText("No database selected");
    }else ui->progressLabel->setText("No resource directory selected");
}

void MainWindow::on_checkBox_clicked()
{
    if (numOneRings || numTwoRings)
    {
        bool selectall = ui->checkBox->isChecked();
        if (selectall)
            ui->checkBox_2->setChecked(0);
        for(int i = 0; i < numOneRings; i++){
            // Find each button by its name (assigned when making the button)
            QToolButton *button = ui->oneringScrollArea->findChild<QToolButton *>(ring1[i].second);
            button->setChecked(selectall);
        }
        for(int i = 0; i < numTwoRings; i++){
            QToolButton *button = ui->tworingScrollArea->findChild<QToolButton *>(ring2[i].second);
            button->setChecked(selectall);
        }
    }
}

void MainWindow::on_checkBox_2_clicked()
{
    if (numOneRings || numTwoRings)
    {
        bool selectall = ui->checkBox_2->isChecked();
        if (selectall)
            ui->checkBox->setChecked(0);
        for(int i = 0; i < numOneRings; i++){
            // Find each button by its name (assigned when making the button)
            QToolButton *button = ui->oneringScrollArea->findChild<QToolButton *>(ring1[i].second);
            button->setChecked(0);
            if (button->text().toInt() != 0)
                button->setChecked(selectall);
        }
        for(int i = 0; i < numTwoRings; i++){
            QToolButton *button = ui->tworingScrollArea->findChild<QToolButton *>(ring2[i].second);
            button->setChecked(0);
            if (button->text().toInt() != 0)
                button->setChecked(selectall);
        }
    }
}

