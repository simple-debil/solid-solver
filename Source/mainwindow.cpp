#include <iostream>

#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "console.h"
#include "interpolation.h"

const int GRAPHLABELSX = 3;
const int GRAPHLABELSY = 3;
Visualization2d::myWidget_out *out_picture;

void test_lagrange1Unregular()
{
    using namespace Elementary;
    // инициализация интерполянта
    Interpolation::Interpolant1D_base *it;
    it = new Interpolation::Interpolant1D_Lagrange1Unregular;
    // построение интерполянта по значениям ф-и в узлах
    std::vector<POINT1> coordinates;    //набора узлов интерполянта
    coordinates.clear();
    double h = 10;  // длина отрезка на котором задана ф-я
    double n = 100; // точек
    for(size_t i = 0; i <= n; i++)
    {
        POINT1 x = 1. + i*i/n/n*h;
        double nodeValue = (3*x + 1);
        it->addPoint(x, nodeValue);
        coordinates.push_back(x);
    }
    it->buildInterpolant();

    size_t t = 10; // количество точек на КЭ
    for(size_t i = 0; i < coordinates.size() - 1; i++)
    {
        for(size_t j = 0; j <= t; j++)
        {
            POINT1 x = coordinates[i] + j*(coordinates[i + 1] - coordinates[i])/t;
            std::cerr << x << "\t" << it->fun(x) - 1*(3*x + 1) << std::endl;
        }
    }

}
void test_lagrange1()
{
    using namespace Elementary;
    // инициализация интерполянта
    Interpolation::Interpolant1D_base *it;
    it = new Interpolation::Interpolant1D_Lagrange1;
    Grid::GridRectangleRegular1D it_grid;
    it_grid.init(1,11,100);
    double CoefAlpha = 0.0000000001;
    it->init(it_grid, CoefAlpha);
    // получение набора узлов интерполянта
    std::vector<POINT1> coordinates;
    it->getNodesCoordinates(coordinates);
    // построение интерполянта по значениям ф-и в узлах
    Vector nodeValue;
    nodeValue.resize(coordinates.size());
    // для лагранжего сплайна
    for(size_t i = 0; i < coordinates.size(); i++)
    {
        POINT1 x = coordinates[i];

        if(i == 0)
            x = coordinates[i] + 0.0001;
        else
            x = coordinates[i] - 0.0001;
        nodeValue[i] = (3*x + 1);
        //if(i != 3)
        {
            it->addPoint(x, nodeValue[i]);
        }
    }
    //it->buildInterpolantByAllNodes(nodeValue);
    it->buildInterpolant();

    int t = 10;
    for(size_t i = 0; i <= (coordinates.size() - 1)*t; i++)
    {
        POINT1 x = it_grid.rect[0] + i*it->h/t;
        std::cerr << x << "\t" << it->fun(x) - (3*x + 1) << std::endl;
    }
}
void test_lagrange3()
{
    // тест 1D интерполянта
    using namespace Elementary;
    // инициализация интерполянта
    Interpolation::Interpolant1D_base *it;
    it = new Interpolation::Interpolant1D_Lagrange3;
    Grid::GridRectangleRegular1D it_grid;
    it_grid.init(1,11,4);
    double CoefAlpha = 0.0000000001;
    it->init(it_grid, CoefAlpha);
    // получение набора узлов интерполянта
    std::vector<POINT1> coordinates;
    it->getNodesCoordinates(coordinates);
    // построение интерполянта по значениям ф-и в узлах
    Vector nodeValue;
    nodeValue.resize(coordinates.size());
    // для лагранжего сплайна
    for(size_t i = 0; i < coordinates.size(); i++)
    {
        POINT1 x;

        if(i == 0)
            x = coordinates[i] + 0.0001;
        else
            x = coordinates[i] - 0.0001;
        nodeValue[i] = (x*x*x + x*x + x + 1);
        if(i != 3)
        {
            it->addPoint(x, nodeValue[i]);
        }
    }
    //it->buildInterpolantByAllNodes(nodeValue);
    it->buildInterpolant();

    int t = 10;
    for(size_t i = 0; i <= coordinates.size()/3*t; i++)
    {
        POINT1 x = it_grid.rect[0] + i*it->h/t;
        std::cerr << x << "\t" << it->fun(x) - (x*x*x + x*x + x + 1) << std::endl;
    }
}
void test_hermite()
{
    // тест 1D интерполянта
    using namespace Elementary;
    // инициализация интерполянта
    Interpolation::Interpolant1D_base *it;
    it = new Interpolation::Interpolant1D_Hermite3;
    Grid::GridRectangleRegular1D it_grid;
    it_grid.init(1,11,20);
    double CoefAlpha = 0.000001;
    it->init(it_grid, CoefAlpha);
    // получение набора узлов интерполянта
    std::vector<POINT1> coordinates;
    it->getNodesCoordinates(coordinates);
    // построение интерполянта по значениям ф-и в узлах
    Vector nodeValue;
    nodeValue.resize(coordinates.size());
    for(size_t i = 0; i < coordinates.size(); i++)
    {
        POINT1 x;
        // buildInterpolant

        if(i != 0)
        {
            x = coordinates[i] - 0.0001;
            nodeValue[i] = (x*x*x + x*x + x + 1);
            it->addPoint(x, nodeValue[i]);
        }
        if(i != coordinates.size() - 1)
        {
            x = coordinates[i] + 0.0001;
            nodeValue[i] = (x*x*x + x*x + x + 1);
            it->addPoint(x, nodeValue[i]);
        }

        // buildInterpolantByAllNodes
        //x = coordinates[i];
        //nodeValue[i] = (x*x*x + x*x + x + 1);
    }
    //it->buildInterpolantByAllNodes(nodeValue);
    it->buildInterpolant();

    int t = 10;
    for(size_t i = 0; i <= (coordinates.size() - 1)*t; i++)
    {
        POINT1 x = it_grid.rect[0] + i*it->h/t;
        std::cerr << x << "\t" << it->fun(x) - (x*x*x + x*x + x + 1) << std::endl;
    }
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    logger.init();
    logger.start();

    QObject::connect(&logger, SIGNAL(workOnReceivedMessage(Threads::Message)), this, SLOT(workOnReceivedMessage(Threads::Message)), Qt::BlockingQueuedConnection);


    /*std::stringstream out;
    for(int dig = 0; dig < 16; dig++)
    {
        double x = -2.e-100;
        double y = Elementary::roundToDigits(x, dig);
        fprintf(stderr, "dig = %d: %le -> %le\n", dig, x, y);
    }*/
    //ui->progress_textEdit->insertPlainText(out.);

    //test_lagrange3();
    //test_hermite();
    //test_lagrange1();
    //test_lagrange1Unregular();









    // загрузка
    load();

    // график
    outChart = new QChart();
    outChartView = new QChartView(outChart);
    outChartView->setRenderHint(QPainter::Antialiasing);
    outChartView->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
    ////ui->outGraphic_layout->addWidget(outChartView);
    // виджет для рисования картинки
    out_picture = new Visualization2d::myWidget_out;
    QScrollArea* scrollArea_out = new QScrollArea;
    scrollArea_out->setWidget(out_picture);
    //scrollArea_out->setBackgroundRole(QPalette::HighlightedText);
        //scrollArea_out->setAutoFillBackground(true);
        //scrollArea_out->show();
    ui->picture_layout->addWidget(scrollArea_out);

    // виджет для рисования графиков
    out_graph = new QLabel;
    out_graph->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    QScrollArea* scrollArea_out_graph = new QScrollArea;
    scrollArea_out_graph->setWidget(out_graph);
    ui->graph_grid_Layout->addWidget(scrollArea_out_graph);


    progress_graph.resize(GRAPHLABELSX*GRAPHLABELSY);
    for(size_t ind = 0; ind < progress_graph.size(); ind++)
    {
        progress_graph[ind] = new QLabel;
        progress_graph[ind]->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    }

    // виджеты для графиков с данными о шагах и итерациях
    QGridLayout *l = new QGridLayout;
    int ind = 0;
    for(int ly = 0; ly < GRAPHLABELSY; ly++)
    {
        for(int lx = 0; lx < GRAPHLABELSX; lx++)
        {
            l->addWidget(progress_graph[ind], ly, lx);
            ind++;
        }
    }


    for(size_t ind = 0; ind < GRAPHLABELSX; ind++)
        l->addWidget(progress_graph[ind], 0, (int)(ind));
    for(size_t ind = 3; ind < progress_graph.size(); ind++)
        l->addWidget(progress_graph[ind], 1, (int)(ind) - 3);

    //for(size_t ind = 0; ind < progress_graph.size(); ind++)
    //    l->addWidget(progress_graph[ind], (int)(ind/2), (int)(ind%2));

    QWidget *w = new QWidget;
    w->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    w->setLayout(l);

    QScrollArea* scrollArea_progress_graph = new QScrollArea;
    scrollArea_progress_graph->setWidget(w);
    scrollArea_progress_graph->setWidgetResizable(true);
    ui->processGraph_gridLayout->addWidget(scrollArea_progress_graph);




    //progress_graph2 = new QLabel;
    //QScrollArea* scrollArea_progress_graph = new QScrollArea;
    //scrollArea_progress_graph->setWidget(progress_graph1);
    //ui->processGraph_gridLayout->addWidget(scrollArea_progress_graph);

    // тесты
    ui->testIndex_comboBox->clear();
    //ui->testIndex_comboBox->addItem("полый шар", QVariant((int)Tests::Test_base::Type::Sphere));
    ui->testIndex_comboBox->addItem("полый шар: теплопроводность", QVariant((int)Tests::Test_base::Type::Sphere_T));
    ui->testIndex_comboBox->addItem("полый шар: упруго-пластичность", QVariant((int)Tests::Test_base::Type::Sphere_ep));
    ui->testIndex_comboBox->addItem("полый шар: ползучесть", QVariant((int)Tests::Test_base::Type::Sphere_creep));
    ui->testIndex_comboBox->addItem("полый шар: термо-упругость", QVariant((int)Tests::Test_base::Type::Sphere_Te));
    ui->testIndex_comboBox->addItem("формовка", QVariant((int)Tests::Test_base::Type::Formovka));
    ui->testIndex_comboBox->addItem("контакт цилиндра и полупространства", QVariant((int)Tests::Test_base::Type::ContactCylinder));
    ui->testIndex_comboBox->addItem("контакт шара и полупространства", QVariant((int)Tests::Test_base::Type::ContactSphere));
    ui->testIndex_comboBox->addItem("ползучесть", QVariant((int)Tests::Test_base::Type::Creep));
    ui->testIndex_comboBox->addItem("задача Кирша", QVariant((int)Tests::Test_base::Type::Kirsch));
    ui->testIndex_comboBox->addItem("неоднородная пластинка с трещиной(все бф интерполянты)", QVariant((int)Tests::Test_base::Type::Crack_plate_gen));
    // теплопроводность
    //ui->testIndex_comboBox->setCurrentIndex(0);
    // пластичность
    //ui->testIndex_comboBox->setCurrentIndex(1); // шар
    // термо-упругость
    //ui->testIndex_comboBox->setCurrentIndex(3); // шар
    // ползучесть
    //ui->testIndex_comboBox->setCurrentIndex(7); // балка
    //ui->testIndex_comboBox->setCurrentIndex(2); // шар
    // формовка
    //ui->testIndex_comboBox->setCurrentIndex(4);
    // контакт цилиндра и полупространства
    //ui->testIndex_comboBox->setCurrentIndex(5);
    // контакт шара и полупространства
    ui->testIndex_comboBox->setCurrentIndex((int)Tests::Test_base::Type::Crack_plate_gen);


    // поле по умолчанию - эквивалентные деформации
    ui->valueIndex->setCurrentIndex(1);

    started = false;
    mainwindowClosing = false;
    ui->progressBar->setEnabled(true);
    ui->progressBar->setHidden(true);
}
MainWindow::~MainWindow()
{
    logger.signalToRelease();
    logger.release();
    delete ui;
}

void MainWindow::load()
{
    using namespace std;
    // параметры решателя
    fstream bin("_autosave.sav", ios::binary | ios::in);
    bin.read((char *)&slausolver_parameters, sizeof(SlauSolving::SolverParameters));
    bin.close();
    // Заполнение форм
    ui->slau_type->setCurrentIndex((int)slausolver_parameters.type);
    ui->slau_preconditioning->setCurrentIndex((int)slausolver_parameters.preconditioning);
    ui->slau_blocks->setCurrentIndex((int)slausolver_parameters.blocks);
    ui->slau_accuracy->setText(DoubleToText(slausolver_parameters.eps));
    ui->slau_maxiter->setValue(slausolver_parameters.maxIter);
    //ui->integration_type->setCurrentIndex((int)task.parameters0.integration_type_cube);
}
void MainWindow::save()
{
    using namespace std;
    // параметры решателя
    fstream bout("_autosave.sav", ios::binary | ios::out);
    bout.write((char *)&slausolver_parameters, sizeof(SlauSolving::SolverParameters));
    bout.close();
}




double MainWindow::QLineEdit_NumberToDouble(QLineEdit *le)
{
    QPalette p;
    QString s = le->text();
    bool ok;
    double value = s.toDouble(&ok);
    if(!ok)
    {
        p = le->palette();
        p.setColor(QPalette::Text, Qt::red);
        le->setPalette(p);
        le->setFocus();
        ui->RunSolverButton->setEnabled(false);
    }
    else
    {
        p = QApplication::palette();
        le->setPalette(p);
        //if(ui_ProgressBarDialog->isHidden())
        //    ui->RunSolverButton->setEnabled(true);
    }
    return value;
}
void MainWindow::QLineEdit_NumberUpdate(QLineEdit *le, std::string &out)
{
    QPalette p;
    std::string s = le->text().toStdString();
    if(!FunParser::strToDoubleTest(s))
    {
        p = le->palette();
        p.setColor(QPalette::Text, Qt::red);
        le->setPalette(p);
        le->setFocus();
        ui->RunSolverButton->setEnabled(false);
        // out не меняется, т.к. формат не соответствует
    }
    else
    {
        p = QApplication::palette();
        le->setPalette(p);
        //if(ui_ProgressBarDialog->isHidden())
        //    ui->RunSolverButton->setEnabled(true);
        out = s;
    }
}
bool MainWindow::QTableWidgetItem_NumberUpdate(QTableWidgetItem *item, std::string &out)
{
    std::string s = item->text().toStdString();
    if(FunParser::strToDoubleTest(s))
    {
        item->setForeground(QApplication::palette().foreground());
        //if(ui_ProgressBarDialog->isHidden())
        //    ui->RunSolverButton->setEnabled(true);
        out = s;
        return true;
    }
    else
    {
        item->setForeground(Qt::red);
        item->setSelected(true);
        ui->RunSolverButton->setEnabled(false);
        return false;
    }
}
bool MainWindow::QTableWidgetItem_NameUpdate(QTableWidgetItem *item, std::string &out)
{
    //char c = item->text().toStdString()[0];
    if(true)//(c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_')
    {
        item->setForeground(QApplication::palette().foreground());
        //if(ui_ProgressBarDialog->isHidden())
        //    ui->RunSolverButton->setEnabled(true);
        out = item->text().toStdString();
        return true;
    }
    else
    {
        item->setForeground(Qt::red);
        item->setSelected(true);
        ui->RunSolverButton->setEnabled(false);
        return false;
    }
}
QString MainWindow::DoubleToText(double value)
{
    QString s;
    for(int i = 0; i <= 16; i++)
    {
        s = QString::number(value, 'g', i);
        if(s.toDouble() == value) break;
    }
    return s;
}

// Нажатия на кнопки

void MainWindow::on_action_2_triggered()
{
    load();
}
void MainWindow::on_action_3_triggered()
{
    save();
}
void MainWindow::on_action_4_triggered()
{
    close();
}
void MainWindow::on_action_triggered()
{
    QMessageBox::about(0, "О программе", "версия 0.0.1");
}


void MainWindow::setMaximum(int maximum)
{
    ui->progressBar->setMaximum(maximum);
}
void MainWindow::setValue(int value)
{
    ui->progressBar->setValue(value);
}
void MainWindow::addMessage(const char *newMessage)
{
    if(logfile != NULL)
    {
        // установка позиции на конец файла
        fseek(logfile, logfile_pos_write, SEEK_SET);
        // определение размера порции для записи
        size_t writeSize = strlen(newMessage);
        // запись порции в файл
        fwrite(newMessage, 1, writeSize, logfile);
        // обновление позиции конца файла
        logfile_pos_write += writeSize;
        //fprintf(logfile_write, "%s", newMessage);
    }
}
void MainWindow::updateLog()
{
    if(logfile != NULL)
    {
        // цикл порционного чтения из файла
        const long buf_size = 1025;
        char buf[buf_size];
        //QTextCursor cursor0 = ui->progress_textEdit->textCursor();
        //QCursor cursor0 = ui->progress_textEdit->cursor();
        for(;;)
        {
            // определение размера следующей порции для чтения
            size_t readSize = MIN(buf_size - 1, logfile_pos_write - logfile_pos_read);
            if(readSize <= 0)
                break;
            // установка позиции на начало непрочитанных данных
            fseek(logfile, logfile_pos_read, SEEK_SET);
            // чтение порции из файла
            fread(buf, 1, readSize, logfile);
            // добавление нулевого символа в конец порции
            buf[readSize] = 0;
            // обновление позиции начала непрочитанных данных
            logfile_pos_read += readSize;
            // вывод порции на виджет
            QTextCursor cursor = ui->progress_textEdit->textCursor();
            cursor.movePosition(QTextCursor::End);
            cursor.insertText(buf);
            /*
            QTextCursor cursor = ui->progress_textEdit->textCursor();
            cursor.movePosition(QTextCursor::End);//cursor.setPosition(QTextCursor::End);
            ui->progress_textEdit->setTextCursor(cursor);
            ui->progress_textEdit->insertPlainText(buf);
            */
        }
        //ui->progress_textEdit->setTextCursor(cursor0);
        //ui->progress_textEdit->setCursor(cursor0);
    }
}


void MainWindow::closeEvent(QCloseEvent *event)
{
    if(started)
    {
        event->ignore();
        mainwindowClosing = true;
        solvingThread.signalToRelease();
    }
}

void MainWindow::on_RunSolverButton_clicked()
{
    if(started)
    {
        // клик на "Остановить"
        ui->RunSolverButton->setEnabled(false);
        solvingThread.signalToRelease();
        //solvingThread.kill();
    }else
    {
        // клик на "Запустить"
        ui->testIndex_comboBox->setEnabled(false);
        ui->RunSolverButton->setText("Остановить");
        ui->progress_textEdit->clear();
        setValue(0);
        ui->progressBar->setEnabled(true);
        ui->progressBar->setHidden(false);
        ui->processStep1_spinBox->setEnabled(false);
        ui->processStep2_spinBox->setEnabled(false);
        ui->processUpdate_pushButton->setEnabled(false);
        // запуск решателя
        processUpdateSteps_updateInterval.begin();
        processUpdate_updateLog.begin();
        if(logfile != NULL)
            fclose(logfile);
        logfile = fopen(testBuilder->get_logFileName()->c_str(), "w+b");
        logfile_pos_write = 0;
        logfile_pos_read = 0;
        out = nullptr;
        solvingThread.setTestParameters(&task, testBuilder, &logger);
        solvingThread.start();
        started = true;
    }

}
void MainWindow::solvingFinished(int result)
{
    if(result == 0)
    {
        solvingThread.release();
        solvingThread.getResults(out);
        int stepsNumberMaximum = 0;
        if(out->mechOut != nullptr)
        {
            stepsNumberMaximum = (int)(*out->mechOut).size() - 2;
        }
        else
        if(out->thermOut != nullptr)
        {
            stepsNumberMaximum = (int)(*out->thermOut).size() - 2;
            //stepsNumberMaximum = 2;
        }
        ui->globalStepIndex->setMaximum(stepsNumberMaximum);
        ui->globalStepIndex->setValue(stepsNumberMaximum);

        out_picture->setResults(&task, testBuilder, out, &vis2d_parameters);
        out_picture->resize(1280, 1024);
        out_picture->repaint();
        // графики строятся после каждого глобального шага
        //testBuilder->makeResultGraphs();

        // loggen.sync()?##
        //if(logfile != NULL)
        //    fclose(logfile);
        //logfile = fopen(testBuilder->get_logFileName()->c_str(), "a+b");

        ui->resGraph_Index_comboBox->currentIndexChanged(ui->resGraph_Index_comboBox->currentIndex());
        ui->progressBar->setEnabled(true);
        ui->progressBar->setHidden(true);
    }
    else
    {
        solvingThread.release();
        ui->progressBar->setEnabled(false);
        ui->progressBar->setHidden(false);
        //QMessageBox::about(0, "Сообщение", "поток преждевременно завершен");
    }
    started = false;
    ui->RunSolverButton->setText("Запустить");
    ui->RunSolverButton->setEnabled(true);
    ui->testIndex_comboBox->setEnabled(true);
    ui->processStep2_spinBox->setMaximum(testBuilder->get_stepsInfLastSavedStepNomber());
    ui->processStep2_spinBox->setValue(ui->processStep2_spinBox->maximum());
    ui->processStep1_spinBox->setEnabled(true);
    ui->processStep2_spinBox->setEnabled(true);
    ui->processUpdate_pushButton->setEnabled(true);
    int stepIndex2 = ui->processStep2_spinBox->maximum();
    int stepIndex1 = stepIndex2 - ui->processUpdateSteps_spinBox->value() + 1;
    if(stepIndex1 < 0)
        stepIndex1 = 0;
    logger.updateGraphs(stepIndex1, stepIndex2);
    //testBuilder->makeStepsGraphs(stepIndex1, stepIndex2);
    //updateStepsGraphs();

    if(mainwindowClosing)
    {
        close();
    }
}
void MainWindow::updateStepsGraphs()
{
    // обновление графиков
    const Tests::GraphInfs *stepsGraph_all_iterations;
    const Tests::GraphInfs *stepsGraph_results_of_iterations;
    testBuilder->get_oldStepsGraphs(stepsGraph_all_iterations, stepsGraph_results_of_iterations);
    size_t countInd = 0;
    // первая строчка
    //countInd = 0*GRAPHLABELSX + 0;
    for(size_t lx = 0; lx < (*stepsGraph_all_iterations).size(); lx++)
    {
        QPixmap picture((*stepsGraph_all_iterations)[lx].fn_pic.c_str());
        if(!picture.isNull())
        {
            progress_graph[countInd]->setPixmap(picture);
            progress_graph[countInd]->resize(progress_graph[countInd]->sizeHint());
            progress_graph[countInd]->show();
            countInd++;
        }
    }
    for(; countInd < 1*GRAPHLABELSX + 0; countInd++)
    {
        progress_graph[countInd]->clear();
    }
    // вторая строчка
    //countInd = 1*GRAPHLABELSX + 0;
    for(size_t lx = 0; lx < (*stepsGraph_results_of_iterations).size(); lx++)
    {
        QPixmap picture((*stepsGraph_results_of_iterations)[lx].fn_pic.c_str());
        if(!picture.isNull())
        {
            progress_graph[countInd]->setPixmap(picture);
            progress_graph[countInd]->resize(progress_graph[countInd]->sizeHint());
            progress_graph[countInd]->show();
            countInd++;
        }
    }
    for(; countInd < 2*GRAPHLABELSX + 0; countInd++)
    {
        progress_graph[countInd]->clear();
    }
    // обновление лога
    updateLog();
}
void MainWindow::workOnReceivedMessage(Threads::Message m)
{
    if(m.type == Threads::Message::Type::string || m.type == Threads::Message::Type::stringstream)
    {
        addMessage(m.data);
        processUpdate_updateLog.end();
        if(processUpdate_updateLog.getDuration() >= 5.0)
        {
            updateLog();
            processUpdate_updateLog.begin();
        }
    }
    if(m.type == Threads::Message::Type::maximumChanged)
    {
        setMaximum(*((int *)m.data));
    }
    if(m.type == Threads::Message::Type::valueChanged)
    {
        setValue(*((int *)m.data));
    }
    if(m.type == Threads::Message::Type::solvingFinished)
    {
        int result = *((int *)m.data);
        solvingFinished(result);
    }
    if(m.type == Threads::Message::Type::stepInf)
    {
        Solid::OutDataForLogger stepsData;
        memcpy(&stepsData, m.data, m.size);
        testBuilder->saveLastStepInf(stepsData);
        // обновляем не слишком часто
        processUpdateSteps_updateInterval.end();
        if(processUpdateSteps_updateInterval.getDuration() >= ui->processUpdateTime_spinBox->value())
        {
            // отправка логгеру запроса на обновление графиков (возврат возможен т.к. логгер сейчас висит)
            int stepIndex2 = stepsData.stepNumber;
            int stepIndex1 = stepIndex2 - ui->processUpdateSteps_spinBox->value() + 1;
            if(stepIndex1 < 0)
                stepIndex1 = 0;
            logger.updateGraphs(stepIndex1, stepIndex2);
            processUpdateSteps_updateInterval.begin();
            processUpdateSteps_updateInterval.beginTime += std::chrono::hours(1000);
        }
    }
    if(m.type == Threads::Message::Type::updateGraphs)
    {
        // обновление графиков
        int stepIndex1 = ((int *)m.data)[0];
        int stepIndex2 = ((int *)m.data)[1];
        testBuilder->makeStepsGraphs(stepIndex1, stepIndex2);
        updateStepsGraphs();
        processUpdateSteps_updateInterval.begin();
    }
    if(m.data != nullptr)
        delete []m.data;
    //std::cout << "[logger]message(size = " << newMessage.size << "): " << std::string(newMessage.data);
}


// Данные формы

// Параметры решателя

// параметры решателя СЛАУ
void MainWindow::on_slau_type_currentIndexChanged(int index)
{
    slausolver_parameters.type = (SlauSolving::SolverType)index;
    if(testBuilder != nullptr)
        testBuilder->slausolver_parameters = slausolver_parameters;
}
void MainWindow::on_slau_preconditioning_currentIndexChanged(int index)
{
    slausolver_parameters.preconditioning = (SlauSolving::Preconditioning)(index);
    if(testBuilder != nullptr)
        testBuilder->slausolver_parameters = slausolver_parameters;
}
void MainWindow::on_slau_blocks_currentIndexChanged(int index)
{
    if(index == 0)
        slausolver_parameters.blocks = false;
    if(index == 1)
        slausolver_parameters.blocks = true;
    if(testBuilder != nullptr)
        testBuilder->slausolver_parameters = slausolver_parameters;
}
void MainWindow::on_slau_maxiter_editingFinished()
{
    slausolver_parameters.maxIter = ui->slau_maxiter->value();
    if(testBuilder != nullptr)
        testBuilder->slausolver_parameters = slausolver_parameters;
}
void MainWindow::on_slau_accuracy_editingFinished()
{
    slausolver_parameters.eps = QLineEdit_NumberToDouble(ui->slau_accuracy);
    if(testBuilder != nullptr)
        testBuilder->slausolver_parameters = slausolver_parameters;
}

// Параметры входных данных
void MainWindow::on_valueIndex_currentIndexChanged(int index)
{
    vis2d_parameters.valueIndex = index;
    out_picture->repaint();
}
void MainWindow::on_colorModeIndex_currentIndexChanged(int index)
{
    vis2d_parameters.colorModeIndex = index;
    out_picture->repaint();
}
void MainWindow::on_globalStepIndex_valueChanged(int arg1)
{
    vis2d_parameters.globalStepIndex = arg1;
    out_picture->repaint();
}
void MainWindow::on_z_index_valueChanged(int arg1)
{
    vis2d_parameters.z_index = arg1;
    out_picture->repaint();
}
void MainWindow::on_savePictureButton_clicked()
{
    // сохранение картинки в файл
    QString fileName = QFileDialog::getSaveFileName(this, tr("Сохранить график"), "", tr("Изображения (*.png)"));
    //out_picture->setEnabled(true);
    //out_picture->
    //out_picture->hide();
    //out_picture->repaint();
    //out_picture->show();
    //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    //QPixmap screen = QPixmap::grabWindow(out_picture->winId());
    QPixmap screen = out_picture->grab();
    screen.save(fileName, "png");
}
void MainWindow::on_pictureSize_valueChanged(int arg1)
{
    vis2d_parameters.pictureSize = arg1;
    out_picture->repaint();
}
void MainWindow::on_textFontSize_valueChanged(int arg1)
{
    vis2d_parameters.textFontSize = arg1;
    out_picture->repaint();
}



/*
    // Отрисовка интерполянта
    p.setBrush(QBrush(QColor(240, 240, 240, 255), Qt::SolidPattern));
    p.setPen(QColor(0, 0, 0, 0));
int sizeX = 100;
int sizeY = 100;
for(int y = 0; y < sizeY; y++)
{
    for(int x = 0; x < sizeX; x++)
    {
        double alpha = 2.*PI*x/sizeX;
        double zz = 5.*y/sizeY;
        double valuex = itx.fun({alpha, zz, 0});    // -R..R
        valuex = (valuex + R) / (2*R);
        double valuey = ity.fun({alpha, zz, 0});    // -R..R
        valuey = (valuey + R) / (2*R);
        double valuez = itx.fun({alpha, zz, 0});    // 0..5
        valuez = (valuez) / (5.);
        QPointF point;
        point.setX(x);
        point.setX(y);
        int c = 255 - (int)(valuex*255);
        p.setPen(QColor(c, c, c, 255));
        p.drawPoint(point);
    }
}
*/

void MainWindow::on_resGraph_updateAll_PushButton_clicked()
{
    testBuilder->makeResultGraphs();
    ui->resGraph_Index_comboBox->currentIndexChanged(ui->resGraph_Index_comboBox->currentIndex());
}

void MainWindow::on_testIndex_comboBox_currentIndexChanged(int) //index
{
    // выбор теста, задание параметров теста
    //ui->testIndex_comboBox->resize(ui->testIndex_comboBox->sizeHint());
    Tests::Test_base::Type testType = (Tests::Test_base::Type)ui->testIndex_comboBox->currentData().toInt();
    testBuilder = Tests::Test_base::gen(testType);
    testBuilder->slausolver_parameters = slausolver_parameters;
    // графики
    const Tests::GraphInfs *resultGraph = testBuilder->get_oldResultGraphs();
    ui->resGraph_Index_comboBox->clear();
    for(size_t ind = 0; ind < (*resultGraph).size(); ind++)
        ui->resGraph_Index_comboBox->addItem((*resultGraph)[ind].fn_pic_short.c_str(), QVariant((int)ind));
    ui->processStep2_spinBox->setMaximum(testBuilder->get_stepsInfLastSavedStepNomber());
    ui->processStep2_spinBox->setValue(ui->processStep2_spinBox->maximum());
    // лог
    ui->progress_textEdit->clear();
    if(logfile != NULL)
        fclose(logfile);
    logfile = fopen(testBuilder->get_logFileName()->c_str(), "a+b");
    if(logfile != NULL)
    {
        fseek(logfile, 0L, SEEK_END);// определение позиции конца файла (последнего байта файла + 1)
        long lastpos = ftell(logfile);
        logfile_pos_write = lastpos;
        logfile_pos_read = 0;
    }
    updateStepsGraphs();
}

void MainWindow::on_resGraph_Index_comboBox_currentIndexChanged(int index)
{
    //ui->graphIndex_comboBox->resize(ui->graphIndex_comboBox->sizeHint());
    const Tests::GraphInfs *resultGraph = testBuilder->get_oldResultGraphs();

    /*for(int i = 0; i < resultGraph->size(); i++)
    {
        std::cout << (*resultGraph)[i].fn_pic << std::endl;
        std::cout << (*resultGraph)[i].fn_graph << std::endl;
        std::cout << (*resultGraph)[i].fn_pic_short << std::endl;
    }
    std::cout << std::endl<< std::endl<< std::endl<< std::endl;*/

    if((index < 0 || index >= (int)(*resultGraph).size()))
    {
        out_graph->clear();
        return;
    }
    QPixmap picture((*resultGraph)[index].fn_pic.c_str());
    // вывод графика
    //out_graph->setScaledContents(true);
    out_graph->setPixmap(picture);
    out_graph->resize(out_graph->sizeHint());
    out_graph->show();
}

void MainWindow::on_processUpdate_pushButton_clicked()
{
    //const Tests::GraphInfs *stepsGraph = testBuilder->getOldStepsGraphs();
    int stepIndex1 = ui->processStep1_spinBox->value();
    int stepIndex2 = ui->processStep2_spinBox->value();
    logger.updateGraphs(stepIndex1, stepIndex2);
    //testBuilder->makeStepsGraphs(stepIndex1, stepIndex2);
    //updateStepsGraphs();
}


