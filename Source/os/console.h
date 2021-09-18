#ifndef CONSOLE_H
#define CONSOLE_H

#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

namespace OS
{
// командный процессор
class Console
{
protected:
    bool enabled;
    // проверяет доступность командного процессора с выводом ошибки, если не доступен
    bool isEnabledCheck();
public:
    Console();
    // возвращает true если командный процессор доступен
    bool isEnabled();
    // выполнение команды
    void exec(const char *command);
    // выполнение команды через pipe и возврат выходного потока
    void exec(const char *command, std::vector<char> &out);
};

// Gnuplot (вызывается через консоль)
class Gnuplot
{
protected:
    Console console;
public:
    // натравить gnuplot на файл fileName
    // результат выполнения не проверяется
    void exec(const char *fileName);
    // натравить gnuplot на файл fileName (через pipe) и вернуть выходной поток stdout в вектор out
    // скрипт должен выдавать файл с картинкой через выходной поток,
    // то есть вывод должен задаваться без параметра: "set output"
    void exec(const char *fileName, std::vector<char> &out);
};
}   // namespace OS


#endif // CONSOLE_H
