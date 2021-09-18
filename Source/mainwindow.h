#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <fstream>
#include <vector>

#include <QMainWindow>
#include <QMessageBox>
#include <QLineEdit>
#include <QPalette>
#include <QTableWidget>
#include <QtCharts/QtCharts>

#include "tests.h"

#include "solvingThread.h"
#include "visualization2d.h"

namespace Ui {
class MainWindow;
}



class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void save();
    void load();
    double QLineEdit_NumberToDouble(QLineEdit *le);
    void QLineEdit_NumberUpdate(QLineEdit *le, std::string &outPoint);
    bool QTableWidgetItem_NumberUpdate(QTableWidgetItem *item, std::string &outPoint);
    bool QTableWidgetItem_NameUpdate(QTableWidgetItem *item, std::string &outPoint);
        QString DoubleToText(double value);

    QChart *outChart;
    QChartView *outChartView;

    QLabel *out_graph;
    std::vector<QLabel *>progress_graph;

    void setMaximum(int maximum);
    void setValue(int value);
    void addMessage(const char *newMessage);
    void updateLog();
    void solvingFinished(int result);
    void updateStepsGraphs();
    void updateGraphs();
    bool started;
    bool mainwindowClosing;
    TimeIntervals::timeInterval processUpdateSteps_updateInterval;
    TimeIntervals::timeInterval processUpdate_updateLog;
    FILE *logfile = nullptr;
    long logfile_pos_write = 0;   // позиция в файле logfile для записи
    long logfile_pos_read = 0;    // позиция в файле logfile для чтения


    SlauSolving::SolverParameters slausolver_parameters;
    Visualization2d::Visualization2d_parameters vis2d_parameters;
    // логгер
    SolidThread::SolverLogger logger;
    // поток для решателя
    SolidThread::SolvingThread solvingThread;
    // входные данные для решателя
    Solid::Task task;
    Tests::Test_base *testBuilder = nullptr;
    // выходные данные от решателя
    Solid::OutData *out = nullptr;

private slots:
    void on_action_2_triggered();
    void on_action_3_triggered();
    void on_action_4_triggered();
    void on_action_triggered();

    virtual void closeEvent(QCloseEvent *event) override;
    void on_RunSolverButton_clicked();
    void workOnReceivedMessage(Threads::Message m);

    void on_slau_type_currentIndexChanged(int index);
    void on_slau_preconditioning_currentIndexChanged(int index);
    void on_slau_accuracy_editingFinished();
    void on_slau_maxiter_editingFinished();

    void on_valueIndex_currentIndexChanged(int index);
    void on_colorModeIndex_currentIndexChanged(int index);
    void on_globalStepIndex_valueChanged(int arg1);
    void on_z_index_valueChanged(int arg1);
    void on_savePictureButton_clicked();
    void on_pictureSize_valueChanged(int arg1);
    void on_textFontSize_valueChanged(int arg1);
    void on_resGraph_updateAll_PushButton_clicked();
    void on_testIndex_comboBox_currentIndexChanged(int index);
    //on_resGraph_Index_comboBox_currentIndexChanged
    //on_testIndex_comboBox_currentIndexChanged
    void on_resGraph_Index_comboBox_currentIndexChanged(int index);
    void on_processUpdate_pushButton_clicked();


    void on_slau_blocks_currentIndexChanged(int index);

public:
    Ui::MainWindow *ui;
    //ProgressBarDialog *ui_ProgressBarDialog;
};

#endif // MAINWINDOW_H
