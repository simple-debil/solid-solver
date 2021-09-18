#include "solvingThread.h"

namespace SolidThread
{
void SolvingThread::setTestParameters(Solid::Task *set_task, Tests::Test_base *set_testBuilder, Threads::Logger_base *set_logger)
{
    task = set_task;
    testBuilder = set_testBuilder;
    logger = set_logger;
}
void SolvingThread::run()
{
    logger->sendString("Started\n");
    //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    testBuilder->initTask(*task);
    task->stepResultsWriter = testBuilder;// для записи рузкльтатов после каждого шага по времени
    logger->sendString("Task inited\n");
    //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    out = solver.solve(*task, logger, signal);
    // проверка сигнала о необходимости завершения процесса
    if(signal->get_needToReleaseState() == false)
    {
        // процесс завершился нормально
        //testBuilder->writeResults(*task, *out);
        //task->stepResultsWriter->writeResults(*task, *out);
        logger->sendString("\nFinished\n");
        logger->solvingFinished(0);
    }
    else
    {
        // процесс прерван
        logger->sendString("\nTerminated\n");
        logger->solvingFinished(1);
    }
    //std::this_thread::sleep_for(std::chrono::milliseconds(5000));
}
void SolvingThread::getResults(Solid::OutData *&_out)
{
    _out = out;
}

void SolverLogger::init()
{
}
}   // namespace SolidThread
