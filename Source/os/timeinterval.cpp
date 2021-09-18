#include "timeinterval.h"
namespace TimeIntervals
{
timeInterval::timeInterval()
{
    sumDurations = 0;        //sumTime = std::chrono::steady_clock::duration{0};
    intervalsNumber = 0;
}
void timeInterval::begin()
{
    beginTime = std::chrono::steady_clock::now();
}
void timeInterval::end()
{
    duration = getCurrentDuration();
    sumDurations += duration;
    intervalsNumber++;
}

double timeInterval::getCurrentDuration()
{
    std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();   // конец временного интервала
    return (double)std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - beginTime).count() * 1.e-9;;
}

double timeInterval::getDuration()const
{
    return duration;
}
double timeInterval::getAverageDuration()const
{
    return sumDurations / intervalsNumber;
}
void timeInterval::print(FILE *f)
{
    fprintf(f, "duration = %le, sumDurations = %le, averageDuration = %le\n", duration, sumDurations, getAverageDuration());
}
}  // namespace TimeIntervals
