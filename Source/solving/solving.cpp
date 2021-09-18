#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "stdio.h"
#include <string.h>

#include "solving.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
using namespace Elementary::Operations;

namespace Solid
{

void Task::save(const std::string &subdir) const
{
    grid->save(subdir);
    mechTask.save(subdir);
}

// решатель
OutData *Solver::solve(Task &task, Threads::Logger_base *set_logger, Threads::SignalToThread_base *set_signal)
{
    initLogging(set_logger, set_signal);
    ms.stepInf.clear();
    TimeIntervals::timeInterval time1;
    time1.begin();
    int progress = 0;
    logger->maximumChanged(100);
    logger->valueChanged(progress);
    // временные слои
    Grid::TimeLayers tl;
    int globalStepNumber = 0;      // номер глобального шага

    // подготовка к первому глобальному шагу
    {
        Grid::GlobalStep &step0 = (*task.step)[0];
        tl.dt = step0.dt0;              // текущий шаг по времени
        tl.t0 = step0.t_start;               // время конца шага
        tl.t_prev1 = step0.t_start - step0.dt0;   // время начала шага
        tl.t_prev2 = 0;                 // время 2 шага назад(не важно, т.к. после первого 3-точечного шага затрется)
        tl.t_prev3 = 0;                 // время 3 шага назад(не важно, т.к. сразу затрется)
    }


    ts.init(task.thermTask, set_logger, set_signal);
    ms.init(task.mechTask, set_logger, set_signal, tl);
    // Шаги по времени
    for (int stepNumber = 0;;stepNumber++)
    {
        Grid::GlobalStep &step_el = (*task.step)[globalStepNumber];
        ms.stepInf.push_back({});
        ms.stepInf.back().globalStepNumber = globalStepNumber;
        ms.stepInf.back().stepNumber = stepNumber;

        ms.stepInf.back().prepareForStep.begin();
        ms.prepareForStep((*ms.mechStep)[globalStepNumber]);
        ms.stepInf.back().prepareForStep.end();
        ts.prepareForStep();

        tl.t_prev3 = tl.t_prev2;
        tl.t_prev2 = tl.t_prev1;
        tl.t_prev1 = tl.t0;
retryAdaptive_t:;
        if(ms.stepInf.size() <= 1)
        {
            ms.stepInf.back().iterInf.back() = MechIterInf({});
        }
        else
        {
            ms.stepInf.back().iterInf.back() = ms.stepInf[ms.stepInf.size() - 2].iterInf.back();
        }
        //ms.stepInf.back().iterInf.back() = MechIterInf({});

        // вычисление времени конца шага (начало шага tl.t1)
        // если почти достигли конца промежутка времени то шагаем ровно к концу, иначе просто делаем шаг
        // (шаг увеличится не более чем на 10%)
        bool saveDetailedInf = false;
        if(fabs(tl.t_prev1 + tl.dt - step_el.t_finish)/fabs(tl.dt) < 0.1)
        {
            tl.dt = step_el.t_finish - tl.t_prev1;
            tl.t0 = step_el.t_finish;
            saveDetailedInf = true;
        }
        else
        {
            tl.t0 = tl.t_prev1 + tl.dt;
        }
        PRINT("\nGlobalStep = %d, stepNumber = %d, t = %lf, dt = %lf\n", ARGS(globalStepNumber, stepNumber, tl.t0, tl.dt));
        // шаг по времени
        ts.tryStep(globalStepNumber, tl);

        // пересылка решения из температурного решателя в механический
        if(task.mechTask.enabled && task.thermTask.enabled)
        {
            for(size_t feInd = 0; feInd < task.thermTask.fe->size(); feInd++)
            {
                Solid::MechFeData_base *feDataEl = (*task.mechTask.fe)[feInd];
                Heat::ThermFeData &feDataElT = (*task.thermTask.fe)[feInd];
                Heat::ThermFePointData &fePointDataElT = feDataElT.pd[0];
                feDataEl->setT(fePointDataElT.T,
                               fePointDataElT.newT - fePointDataElT.T);
            }
        }

        ms.stepInf.back().prepareForIter.begin();
        ms.prepareForIter((*ms.mechStep)[globalStepNumber], tl);
        ms.stepInf.back().prepareForIter.end();

        // шаг по времени
        MechIterationsGeneralResult mechIterationsGeneralResult;
        ms.stepInf.back().tryStep.begin();
        ms.tryStep((*ms.mechStep)[globalStepNumber], tl,
                                     mechIterationsGeneralResult);
        ms.stepInf.back().tryStep.end();

        // проверка результата итераций
        switch (mechIterationsGeneralResult)
        {
        case MechIterationsGeneralResult::Interrupted_Ok:
        {
            // ok
        }break;
        case MechIterationsGeneralResult::Interrupted_OutOfIterationsLimit:
        {
            // желаемая точность не достигнута!
            logger->sendString("WARNING: tochnost ne dostignuta!");
            int controlMode = (*ms.mechStep)[globalStepNumber].controlMode;
            if(controlMode == 0)
            {
                logger->sendString(" Ignoriruu\n");
            }
            if(controlMode == 1)
            {
                // адаптация шага по времени
                logger->sendString(" dt umenshau v 2 raza\n");
                tl.dt /= 2;
                goto retryAdaptive_t;
            }
            if(controlMode == 2)
            {
                logger->sendString(" Prejdevremenno zavershau\n");
                goto endIters;
            }
        }break;
        case MechIterationsGeneralResult::Interrupted_Terminated:
        case MechIterationsGeneralResult::Continue_NeedMoreIterations:
        {
            logger->sendString("Solver terminating...\n");
            goto endIters;
        }break;
        }

        // добавление приращений и сдвиг сетки
        ts.finalizeStep();
        ms.stepInf.back().finalizeStep.begin();
        ms.finalizeStep((*ms.mechStep)[globalStepNumber], tl);
        ms.stepInf.back().finalizeStep.end();

        // сохранение результатов шага по времени
        ts.saveResultsOfStep(globalStepNumber, tl);
        ms.stepInf.back().saveResultsOfStep.begin();
        ms.saveResultsOfStep(saveDetailedInf, globalStepNumber, stepNumber, tl);
        ms.stepInf.back().saveResultsOfStep.end();

        ms.stepInf.back().writeResults.begin();
        if(saveDetailedInf)
        {
            // запуск обработки результатов
            {
                OutData ttt;
                ttt.thermOut = ts.out;
                ttt.mechOut  = ms.out;
                task.stepResultsWriter->writeResults(task, ttt);
            }

            // переход к следующему глобальному шагу
            globalStepNumber++;
            if(globalStepNumber < (int)task.step->size())
            {
                Grid::GlobalStep &new_step_el = (*task.step)[globalStepNumber];
                tl.dt = new_step_el.dt0;
            }
            else
            {
                progress = 100;
                break;
            }
        }
        else
        {
            /*
            MechGlobalStep &mech_step_el = (*task.mechTask.mechStep)[globalStepNumber];
            if(mech_step_el.controlMode == 1)
            {
                // адаптация шага по времени
                if(nonlinearIter <= 8)
                {
                    if(abs(dproportion*1.5) <= dproportionMax)
                        dproportion *= 1.5;
                    else
                        dproportion = direction * dproportionMax;
                }
                if(nonlinearIter >= 32)
                    dproportion /= 1.5;
            }*/
        }
        ms.stepInf.back().writeResults.end();
        ms.stepInf.back().print();


        // примерный процент проделанной работы
        double t_start; // t_start = время конца прошлого шага, т.к. начало необходимо задавать только для 1-го шага
        if(globalStepNumber == 0)
            t_start = step_el.t_start;
        else
            t_start = (*task.step)[globalStepNumber - 1].t_finish;
        double partOfGlobalStep;
        if(step_el.t_finish - t_start == 0)
            partOfGlobalStep = 0;
        else
            partOfGlobalStep = fabs((tl.t0 - t_start) / (step_el.t_finish - t_start));
        progress = (int)round(((double)globalStepNumber + partOfGlobalStep)/task.step->size()*100);
        logger->valueChanged(progress);
    }
endIters:
    PRINT("Time = %.1lf\n", time1.getCurrentDuration());
    logger->valueChanged(progress);
    //fclose(ff);
    out = new OutData;
    out->thermOut = ts.out;
    out->mechOut  = ms.out;
    return out;
}

}   // namespace Solving
