#include "mainwindow.h"
#include <QApplication>


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    setlocale(LC_NUMERIC, "C"); // принудительно устанавливает разделитель между целой и дроброй частями нецелых чисел (точка)

    MainWindow w;
    try
    {
        //w.showFullScreen();
        //w.show();
        w.showMaximized();
        return a.exec();
    }
    catch(const char *err)
    {
        QMessageBox::about(0, "Main", err);
        return 1;
    }
}
