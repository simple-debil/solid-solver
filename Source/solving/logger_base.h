#ifndef LOGGER_BASE_H
#define LOGGER_BASE_H

#include <cstddef>
#include <sstream>

#define ARGS(...) __VA_ARGS__
#define PRINT(format, args){ char strtemp[1000]; sprintf(strtemp, format, args); logger->sendString(strtemp);}
#define PRINT1(format) { char strtemp[1000]; sprintf(strtemp, format); logger->sendString(strtemp);}


namespace Threads
{
// сигнал потоку
class SignalToThread_base
{
public:
    // получение значения сигнала (true/false)
    virtual bool get_needToReleaseState() = 0;
    // изменение значения сигнала
    virtual void set_needToReleaseState(const bool new_needToReleaseState) = 0;
    virtual ~SignalToThread_base(){};
};
// сообщение логгеру
struct Message
{
    enum class Type
    {
        solvingFinished = 0,
        string = 1,
        stringstream = 2,
        stepInf = 3,
        valueChanged = 4,
        maximumChanged = 5,
        updateGraphs = 6,
    };
    Type type;
    size_t size;
    char *data;
};
// логгер, который принимает и обрабатывает сообщения
class Logger_base
{
public:
    // отправка сообщения логгеру
    // сообщение будет находится в очереди, пока не обработается, поэтому лучше выделить память для копии данных в newMessage
    virtual void send(const Message newMessage) = 0;
    // заворачивание в объект Message строки
    // отправка посредством send, и очистка buffer
    virtual void sendString(const std::string &str) = 0;
    // заворачивание в объект Message текстового сообщения, сформированного в buffer,
    // отправка посредством send, и очистка buffer
    virtual void sendBuffer(std::stringstream &buffer) = 0;
    // заворачивание в объект Message сообщения об изменении максимального значения шкалы прогресса
    // и отправка посредством send
    virtual void maximumChanged(const int newMaximum) = 0;
    // заворачивание в объект Message сообщения об изменении значения шкалы прогресса
    // и отправка посредством send
    virtual void valueChanged(const int newValue) = 0;
    // заворачивание в объект Message сообщения о завершении прогресса
    // и отправка посредством send
    virtual void solvingFinished(const int result) = 0;
    // функция-обработчик вызывается логгером для каждого полученного сообщения
    // посылается копия данных
    // освобождение памяти на стороне обработчика сообщений
    virtual void updateGraphs(const int stepIndex1, const int stepIndex2) = 0;
    virtual void workOnReceivedMessage(Message newMessage) = 0;
    virtual ~Logger_base(){};
};
// содержит компоненты для логгирования
class Logging
{
public:
    void initLogging(Threads::Logger_base *set_logger, Threads::SignalToThread_base *set_signal)
    {
        logger = set_logger;
        signal = set_signal;
        outStream.str("");
    }
    Threads::Logger_base *logger;           // логгер
    Threads::SignalToThread_base *signal;   // для передачи сигналов решателю
    std::stringstream outStream;            // поток вывода текста, который отправляется логгеру
};
}   //namespace Threads

#endif // LOGGER_BASE_H
