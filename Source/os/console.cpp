#include "console.h"

namespace OS
{
bool Console::isEnabledCheck()
{
    if(isEnabled())
    {
        std::cerr << "error: command processor is disabled" << std::endl;
        return false;
    }
    else
    {
        return true;
    }
}
Console::Console():
    enabled((std::system(nullptr) == 0))
{
}
bool Console::isEnabled()
{
    return enabled;
}
void Console::exec(const char *command)
{
    if(!isEnabledCheck())
        return;
    if(std::system(command) != 0)
    {
        std::cerr << "error of execution command \"" << command <<  "\"" << std::endl;
    }
}
void Console::exec(const char *command, std::vector<char> &out)
{
    if(!isEnabledCheck())
        return;
    FILE *pipe;
    // запуск команды и открытие трубы для чтения
#ifdef WIN32
    pipe = _popen(command, "r");
#else
    pipe  = popen(command, "r");
#endif
    if(!pipe)
    {
        std::cerr << "error: cant open pipe to execute a command \"" << command << "\"" << std::endl;
        return;
    }
    // чтение вывода из трубы
    for(;;)
    {
        char c;
        if(fscanf(pipe, "%c", &c) == EOF)
            break;
        else
            out.push_back(c);
    }
    // закрытие трубы
#ifdef WIN32
    _pclose(pipe);
#else
    pclose(pipe);
#endif
}

void Gnuplot::exec(const char *fileName)
{
    std::string command;
    command = "gnuplot ";
    command += fileName;
    console.exec(command.c_str());
}
void Gnuplot::exec(const char *fileName, std::vector<char> &out)
{
    std::string command;
#ifdef WIN32
    command = "pgnuplot -persist ";
    command += fileName;
#else
    command = "gnuplot -persist ";
    command += fileName;
#endif
    //
    // запуск gnuplot через pipe и чтение вывода
    console.exec(command.c_str(), out);
}
}   // namespace OS
