#ifndef SOLVINGTHREAD_H
#define SOLVINGTHREAD_H

#include "solving.h"
#include "tests.h"
#include "logger.h"
//#include "progressbardialog.h"
#include <QObject>


namespace SolidThread
{
// поток для решателя
class SolvingThread: public Threads::Thread
{
public:
    //SolvingThread(QObject *parent = nullptr) : QThread(parent){}
    // задание параметров теста
    void setTestParameters(Solid::Task *set_task, Tests::Test_base *set_testBuilder, Threads::Logger_base *set_logger);
    // процесс
    virtual void run() override;
    // получение результата теста
    void getResults(Solid::OutData *&_out);
protected:
    Threads::Logger_base *logger;
    Solid::Task *task;
    Tests::Test_base *testBuilder;
    Solid::Solver solver;
    Solid::OutData *out;
};

// логгер для решателя (реализация обработки сообщений)
class SolverLogger: public QObject, public Threads::Logger
{
    Q_OBJECT
public:
    void init();
signals:
    void workOnReceivedMessage(Threads::Message m) override;    // virtual
};
} // namespace SolidThread

#endif // SOLVINGTHREAD_H
