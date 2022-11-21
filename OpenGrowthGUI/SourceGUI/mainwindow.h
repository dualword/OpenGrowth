#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMultiMap>
#include <QProcess>
#include <ctime>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:

private slots:
    void processlog();
    void processerr();
    void loadprocesslog();
    void loadprocesserr();
    void killProcess();
    void on_selectDatabase_clicked();
    void on_updateDatabase_clicked();
    void on_runButton_clicked();
    void on_createFragmentFileButton_clicked();
    void on_outputdirpushButton_clicked();
    void on_resourcedirpushButton_clicked();
    void on_threemerpushButton_clicked();
    void on_killButton_clicked();
    void on_showIconButton_clicked();
    void on_checkBox_clicked();
    void on_threemerlistpushButton_clicked();
    void on_checkBox_2_clicked();

private:
    Ui::MainWindow *ui;

    int numNoRings;
    QMultiMap<int, QString> ring1map;                           // <number of times found in database, fragment ID>
    int numOneRings;
    QMultiMap<int, QString> ring2map;
    int numTwoRings;

    QList<QPair<QString, QString> > ring0;
    QList<QPair<QString, QString> > ring1;
    QList<QPair<QString, QString> > ring2;

    void searchForFragment(std::string);
    void findFragments();
    int  sizeFile(const char*);
    void showButtons();
    void deleteButtons();
    void generate();
    void delay();
    void loadFragments();
    void findThreemers();
    int  fragmentNo;
    int  ringType;
    int  flag;
    QProcess *myProcess;
};

#endif // MAINWINDOW_H

