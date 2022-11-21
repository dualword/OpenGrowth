#include "mainwindow.h"
#include <QApplication>
#include <QLibrary>
#include <QDebug>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    /*
    // Tells gui where to find Qt libraries
    QString libPath = QCoreApplication::applicationDirPath() + "/plugins";
    qApp->addLibraryPath(libPath);

    // Loads required libraries
    QLibrary icudt51("icudt51.dll");
    if (!icudt51.load())    { return 1; }

    QLibrary icuin51("icuin51.dll");
    if (!icuin51.load())    { return 1; }

    QLibrary icuuc51("icuuc51.dll");
    if (!icuuc51.load())    { return 1; }

    QLibrary libgcc("libgcc_s_dw2-1.dll");
    if (!libgcc.load())     { return 1; }

    QLibrary libstd("libstdc++-6.dll");
    if (!libstd.load())     { return 1; }

    QLibrary libwinp("libwinpthread-1.dll");
    if (!libwinp.load())    { return 1; }

    QLibrary Qt5Core("Qt5Core.dll");
    if (!Qt5Core.load())    { return 1; }

    QLibrary Qt5Gui("Qt5Gui.dll");
    if (!Qt5Gui.load())     { return 1; }

    QLibrary Qt5Widgets("Qt5Widgets.dll");
    if (!Qt5Widgets.load()) { return 1; }

    QLibrary Qt5Xml("Qt5Xml.dll");
    if (!Qt5Xml.load())     { return 1; }

    QLibrary Qt5Svg("Qt5Svg.dll");
    if (!Qt5Svg.load())     { return 1; }

    QLibrary svgicon("qsvgicon.dll");
    if (!svgicon.load())    { return 1; }

    QLibrary svg("qsvg.dll");
    if (!svg.load())        { return 1; }
    */

    MainWindow w;
    w.show();
    return a.exec();
}

