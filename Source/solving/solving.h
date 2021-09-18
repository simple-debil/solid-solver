/* --------------------------------------------------------- */
// РЕШАТЕЛЬ
/* --------------------------------------------------------- */

#ifndef SOLVING_H
#define SOLVING_H

#include "heat.h"
#include "solid.h"
#include "logger_base.h"

namespace Solid
{
// МДТТ + теплопроводность

struct Task;
struct OutData;
// записывальщик результатов шага
struct StepResultsWriter_base
{
    // вывод результата теста и построение графиков
    virtual void writeResults(const Solid::Task &task, const Solid::OutData &out) = 0;
};

// входные данные
struct Task
{
    Grid::Grid3D *grid;             // сетка
    std::vector<Grid::GlobalStep> *step;  // информация о глобальных шагах по времени
    Heat::ThermTask thermTask;      // входные данные для решателя задачи теплопроводности
    Solid::MechTask mechTask;       // входные данные для решателя задачи МДТТ
    StepResultsWriter_base *stepResultsWriter = nullptr;
    // сохранение состояния в файл
    void save(const std::string &subdir) const;
};

// выходные данные
struct OutData
{
    std::vector<Heat::ThermOutGlobalStepData> *thermOut;   // выходные данные из решателя задачи теплопроводности
    std::vector<Solid::MechOutGlobalStepData> *mechOut;    // выходные данные из решателя задачи МДТТ
};

// решатель
class Solver: public Threads::Logging
{
public:
    OutData *solve(Task &task, Threads::Logger_base *set_logger, Threads::SignalToThread_base *set_signal);
    OutData *out;                // выходные данные
protected:
    Heat::ThermSolver ts;
    Solid::MechSolver ms;
};

}   // namespace Solving

#endif  // SOLVING_H
