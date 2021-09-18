#ifndef TIMEINTERVAL_H
#define TIMEINTERVAL_H

#include <stdio.h>
#include <chrono>

namespace TimeIntervals
{
// временной интервал
struct timeInterval
{
    std::chrono::steady_clock::time_point beginTime; // начало временного интервала
    double duration;                                 // продолжительностей последнего из измеренных промежутков
    double sumDurations;                             // сумма продолжительностей всех измеренных промежутков
    int intervalsNumber;                             // общее количество интервалов
    // инициализация суммарной продолжительности и количества интервалов нулём
    timeInterval();
    // начало отсчёта времени
    void begin();
    // конец отсчёта времени
    // вычисление продолжительности временного промежутка, добавление её к суммарной продолжительности, увеличение счётчика измеренных промежутков
    void end();
    // конец отсчёта времени
    // возвращает прошедшее время после вызова begin() до вызова end(), без подсчёта статистики
    double getCurrentDuration();
    // возвращает продолжительность последнего временного промежутка, в секундах
    double getDuration()const;
    // возвращает среднюю продолжительность временного промежутка
    double getAverageDuration()const;
    // вывод в файл
    void print(FILE *f);
};
}  // namespace TimeIntervals

#endif // TIMEINTERVAL_H
