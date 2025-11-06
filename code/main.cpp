#include "mainwindow.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QtWidgets/QStyleFactory>
#include <QtCore/QLibraryInfo>

#ifdef _WIN32
#include <windows.h>
#endif

int main(int argc, char *argv[])
{
#ifdef _WIN32
    // Show console window on Windows for debug output
    AllocConsole();
    FILE* fp;
    freopen_s(&fp, "CONOUT$", "w", stdout);
    freopen_s(&fp, "CONOUT$", "w", stderr);
#endif

    QApplication app(argc, argv);

    QStringList styles = QStyleFactory::keys();
    if (styles.contains("windowsvista", Qt::CaseInsensitive)) {
        QApplication::setStyle("windowsvista");
    }
    else if (styles.contains("fusion", Qt::CaseInsensitive)) {
        QApplication::setStyle("fusion");
    }
    else {
        qDebug() << "Styles 'windowsvista' and 'fusion' not found!";
        qDebug() << "Available Qt Styles:";
        for (const QString& style : styles) {
            qDebug() << " " << style.toLower();
        }
        qDebug() << "The application's appearance may differ from the intended design.";
        qDebug() << "Plugins path:" << QLibraryInfo::path(QLibraryInfo::PluginsPath);
    }

    QSurfaceFormat f;
    f.setVersion(3, 3);
    f.setProfile(QSurfaceFormat::CoreProfile);
    QSurfaceFormat::setDefaultFormat(f);

    MainWindow mainWin;
    mainWin.showMaximized();

    return app.exec();
}
