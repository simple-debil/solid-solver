#include <iostream>
#include <string.h>
#include <functional>

#include "logger.h"

// в случае отсутствия новых сообщений, отдых в течение 1000 мс перед следующей проверкой
#define SleepMilliseconds 1000

namespace Threads
{
bool SignalToThread::get_needToReleaseState()
{
    mutex_needToRelease.lock();
    if(needToReleaseState)
    {
        mutex_needToRelease.unlock();
        return true;
    }
    else
    {
        mutex_needToRelease.unlock();
        return false;
    }
}
void SignalToThread::set_needToReleaseState(const bool new_needToReleaseState)
{
    mutex_needToRelease.lock();
    needToReleaseState = new_needToReleaseState;
    mutex_needToRelease.unlock();
}

void Thread::start()
{
    thread = new std::thread(std::bind(&Thread::run, this));
    signal = new SignalToThread;
}
void Thread::signalToRelease()
{
    signal->set_needToReleaseState(true);
}
void Thread::release()
{
    thread->join();
    delete thread;
    delete signal;
}
void Thread::kill()
{
    thread->detach();
    thread->~thread();  // осталось поймать исключение внутри потока!
    //if(thread->joinable())
    //    std::terminate();
    //delete thread;
}
/*
void Thread_base::run_catch()
{
    try
    {
        run();
    }
    catch (std::string err)
    {
        std::cout << err;
        return;
    }
}
*/

void Logger::send(const Message newMessage)
{
    mutex_queue.lock();
    messagesQueue.push(newMessage);
    cv_queue.notify_one();
    mutex_queue.unlock();
}
void Logger::sendString(const std::string &str)
{
    Message m;
    m.type = Message::Type::string;
    int size = str.size();
    m.size = size + 1;  // + символ "\0"
    m.data = new char[m.size];
    memcpy(m.data, str.c_str(), size);
    m.data[m.size - 1] = 0;
    send(m);
    //str.clear();
    /*outStream << next_iterNumber
              << " ci = " << next_contact_iterNumber
              << " pi = " << next_plastic_iterNumber
              << " genGMb..";
    logger->sendBuffer(outStream);
    */
}
void Logger::sendBuffer(std::stringstream &buffer)
{
    Message m;
    m.type = Message::Type::stringstream;
    buffer.seekg(0, std::ios::end);
    int size = buffer.tellg();
    buffer.seekg(0, std::ios::beg);
    m.size = size + 1;  // + символ "\0"
    m.data = new char[m.size];
    buffer.read(m.data, size);
    m.data[m.size - 1] = 0;
    send(m);
    buffer.str("");
}
void Logger::maximumChanged(const int newMaximum)
{
    Message m;
    m.type = Message::Type::maximumChanged;
    m.size = sizeof(newMaximum);
    m.data = new char[m.size];
    memcpy(m.data, &newMaximum, m.size);
    send(m);
}
void Logger::valueChanged(const int newValue)
{
    Message m;
    m.type = Message::Type::valueChanged;
    m.size = sizeof(newValue);
    m.data = new char[m.size];
    memcpy(m.data, &newValue, m.size);
    send(m);
}
void Logger::solvingFinished(const int result)
{
    Message m;
    m.type = Message::Type::solvingFinished;
    m.size = sizeof(result);
    m.data = new char[m.size];
    memcpy(m.data, &result, m.size);
    send(m);
}

void Logger::updateGraphs(const int stepIndex1, const int stepIndex2)
{
    Message m;
    m.type = Message::Type::updateGraphs;
    m.size = sizeof(stepIndex1) + sizeof(stepIndex2);
    m.data = new char[m.size];
    memcpy(m.data, &stepIndex1, sizeof(stepIndex1));
    memcpy(m.data + sizeof(stepIndex1), &stepIndex2, sizeof(stepIndex2));
    send(m);
}
void Logger::run()
{
    for(;;)
    {
        // ждём накопления сообщений
        //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        // ожидание сообщения
        std::unique_lock<std::mutex> locker_queue(mutex_queue);
        for(;;)
        {
            if(!messagesQueue.empty())
                break;
            if(cv_queue.wait_for(locker_queue, std::chrono::milliseconds(SleepMilliseconds)) == std::cv_status::timeout)
            {
                //std::cout << "[logger]listen..." << std::endl;
                if(signal->get_needToReleaseState() == true)
                {
                    goto exit;
                }
            }
        }
        for(;;)
        {
            if(messagesQueue.empty())
                break;
            Message &m = messagesQueue.front();
            locker_queue.unlock();
            // обработка
            workOnReceivedMessage(m);
            // обработка
            locker_queue.lock();
            messagesQueue.pop();
        }
        locker_queue.unlock();
    }
exit:;
}
}   //namespace Threads
