/* --------------------------------------------------------- */
// МАТЕРИАЛЫ, КРАЕВЫЕ УСЛОВИЯ ДЛЯ МДТТ
/* --------------------------------------------------------- */

#ifndef SOLID_H
#define SOLID_H

#include "logger_base.h"
#include "timeinterval.h"
#include "solidPlastic.h"
#include "solidContact.h"
#include "slau3x3.h"

namespace Solid
{
using namespace Elementary;

// затраты времени на различные этапы вычислений на некоторой итерации
struct MechTimeIterInf
{
    TimeIntervals::timeInterval buildGMb;
     TimeIntervals::timeInterval buildGMb_buildNonContactLocalGMb;
     TimeIntervals::timeInterval buildGMb_addNonContactLocalGMb_to_global;
     TimeIntervals::timeInterval buildGMb_addContact;
    TimeIntervals::timeInterval SSCMaddBoundaryCondition1;
    TimeIntervals::timeInterval solveSLAU;
     TimeIntervals::timeInterval solveSLAU_Gpreconditioner_init;
     TimeIntervals::timeInterval solveSLAU_Gsolver_solve;
    TimeIntervals::timeInterval calcIterResults;
     TimeIntervals::timeInterval calcIterResults_fe;
     TimeIntervals::timeInterval calcIterResults_contact;
     void print(const int numIter, FILE *f);
};

// информация о решении СЛАУ
struct MechSlauIterInf
{
    // что изменилось
    bool matrixG_elementsChanged;       // true если требуется обновить элементы матрицы G
    bool matrixG_portraitChanged;       // true если требуется обновить портрет матрицы G
    bool matrixM_elementsChanged;       // true если требуется обновить элементы матрицы M
    bool matrixM_portraitChanged;       // true если требуется обновить портрет матрицы M
    bool precChanged;                   // true если требуется обновить предобусловливатель
    int GFirstStrChanged;               // номер первой изменённой строки
    bool GsolverChanged;                // true если требуется инициализация решателя СЛАУ
    // требуется ли вывести отладочную информацию о матрицах(параметры и портреты)
    bool debugMatrixesInformationNeedToWrite;
    // информация о решении СЛАУ на текущем шаге
    double slauResidual;
    double slauRelativeResidual;
    int slauIterations;
    double slauTime;
    int slauSolvedWithoutPreconditioner;        // количество раз решали СЛАУ без предобусловливания
    double slauTimeSumWithoutPreconditioner;    // общее время, затраченное на решение СЛАУ (+ предобусловливание), начиная с последнего построения предобусловливателя
    // невязки
    double slauResidualWithLastq;
    // невязки предыдущего шага
    double lastSlauResidualWithLastq;

    MechSlauIterInf();
    // завершение шага по времени (добавка приращений, обновление значений величин)
    void finalizeStep();
};
// информация о итерации
struct MechIterInf
{
    int iterNumber;                     // номер итерации (по всем нелинейностям)
    bool accuracyAchieved;              // true если итерации можно завершить на текущей итерации
     bool improvingAccuracy;                // true если итерации продолжают улучшать невязки (для дожимания)
    MechSlauIterInf slau;
    MechPlasticIterInf plastic;
    MechContactIterInf contact;
    MechTimeIterInf timeIterInf;
    bool BC1Changed;                    // true если требуется обновить первые краевые
    bool BC2Changed;                    // true если требуется обновить вторые краевые
    MechIterInf();
    // подготовительные процедуры перед итерациями
    void prepareForIter(const MechGlobalStepParameters &step_el);
    // обнуление счётчиков (вызывается в начале calcIterResults)
    void clearCounts(const int matrixSize);
    // вычисление результатов итерации
    void calcIterResults(const MechGlobalStepParameters &step_el, const bool needToReleaseState, const size_t matrixSize, MechIterationsGeneralResult &mechIterationsGeneralResult);
    // завершение шага по времени (добавка приращений, обновление значений величин)
    void finalizeStep();
    void save(FILE *f);
    void load(FILE *f);
};
// информация о шаге по времени
struct MechStepInf
{
    std::vector<MechIterInf> iterInf;
    int globalStepNumber;     // номер глобального шага
    int stepNumber;        // номер шага
    // затраты времени на различные этапы вычислений
    TimeIntervals::timeInterval prepareForStep;
    TimeIntervals::timeInterval prepareForIter;
     TimeIntervals::timeInterval prepareForIter_prepareBC1;
     TimeIntervals::timeInterval prepareForIter_prepareBC2;
     TimeIntervals::timeInterval prepareForIter_prepareContact;
    TimeIntervals::timeInterval tryStep;
    TimeIntervals::timeInterval finalizeStep;
     TimeIntervals::timeInterval finalizeStep_fe;
     TimeIntervals::timeInterval finalizeStep_vertex;
     TimeIntervals::timeInterval finalizeStep_contact;
     TimeIntervals::timeInterval finalizeStep_movingGrid;
    TimeIntervals::timeInterval saveResultsOfStep;//##
    TimeIntervals::timeInterval writeResults;//##
    // вывод структуры в файлы
    void print();
};
// характеристики шага по времени
struct MechOutStepData
{
    int globalStepNumber;  // номер глобального шага
    int stepNumber;        // номер шага
    double t0;           // текущее время(после шага)
    MechStepInf nlInf;
};
// данные о шагах для логгера
struct OutDataForLogger
{
    std::vector<MechStepInf> *stepInf;
    int stepNumber;
    int iterNumber;
};

// входные данные
struct MechTask
{
    bool enabled;            // инициализирован?
    std::vector<Grid::GlobalStep> *step;               // общие параметры глобальных шагов
    std::vector<MechGlobalStepParameters> *mechStep;   // последовательность нагружений

    Grid::Grid3D *grid;                             // сетка

    std::vector<MechVertexData> *vertex;            // начальные перемещения
        std::vector<MechVertexData> *vertexForCurvature;            // для криволинейных границ ###
    std::vector<MechFeData_base *> *fe;             // начальные состояния в конечных элементах
    Vector *V0;                                     // начальные скорости
    Vector *dV0;                                    // начальные ускорения



    // поверхности, контакты между ними (КЭ-поверхности, заданные наборами граней сетки, заданы в grid)
    std::vector<Grid::Surface_base *> *rigidSurface;            // жёсткие поверхности для механического контакта
    std::vector<ContactCondition_FE_Rigid> *CC_FE_Rigid;        // контакты КЭ-поверхность - жёсткая поверхность
    //std::vector<MechContactConditionSource_FE_FE> *CCSource_FE_FE;          // контакт КЭ-поверхность - КЭ-поверхность

    // информация в вершинах для учёта реакции опоры
    std::vector<ContactSurfaceVertexData_FE_Rigid> contactVertex_FE_Rigid;

    // инициализация нулями данных: fe (метод интегрирования set_integrationType), vertex, vertexForCurvature, V0, dV0
    // для инициализации fe используются КЭ MechFeData_LinearHexagon_homogeny_E_NLE_EP
    void initNull(const Integration::IntegrationType set_integrationType);
    // сохранение состояния в файл ##начато
    void save(const std::string &subdir) const;
};

// выходные данные
struct MechOutGlobalStepData
{
    std::vector<MechOutStepData> step;
    std::vector<MechOutVertexData> vertex;
    std::vector<MechOutFeData> fe;
};

// решатель задачи МДТТ
struct MechSolver: public MechTask, public Threads::Logging
{
    // размер матрицы
    size_t matrixSize;
    // инфа о шагах по времени
    std::vector<MechStepInf> stepInf;
    //std::vector<MechNonlinearInf> step_nlInf;
    // матрицы
    SlauSolving::SSCM Gmatrix;                          // матрица G
    SlauSolving::SSCMPreconditioner_base *Gpreconditioner;    // предобусловливатель для матрицы G
    SlauSolving::SSCMSolver_base *Gsolver;               // решатель СЛАУ
     SlauSolving::SSCM3x3 Gmatrix3x3;                          // матрица G
     Vector3 b3x3;               // правая часть
     Vector3 dq3x3;              // решение СЛАУ на текущем шаге
     Vector3 dqLastIter3x3;      // перемещения в узлах на предыдущей итерации
     Vector3 r3x3;               // для подсчета невязки
     SlauSolving::SSCM3x3Preconditioner_base *Gpreconditioner3x3;    // предобусловливатель для матрицы G
     SlauSolving::SSCM3x3Solver_base *Gsolver3x3;               // решатель СЛАУ
    SlauSolving::SSCM Mmatrix;                          // матрица M









    // векторы
    Vector b;               // правая часть
     Vector dq;              // решение СЛАУ на текущем шаге - используется для 1х1 и для 3х3
    Vector q;               // перемещения в узлах на текущем шаге
    Vector q1;              // перемещения в узлах 1 шаг назад
    Vector q2;              // перемещения в узлах 2 шага назад
    Vector q3;              // перемещения в узлах 3 шага назад
    Vector dqLastIter;      // перемещения в узлах на предыдущей итерации
    //Vector R;           // R <- Aq - b_other + R
    //Vector R;           // R <- integral(sumSigma) - integral(sumP) - невязка с прошлого шага
    // вспомогательные векторы
    Vector r;               // для подсчета невязки
     Vector Mddq;
     Vector qForM;
    // вспомогательные данные для учёта первых краевых условий
     std::vector<double> bc1_u0;
     std::vector<bool> bc1_state;
    // выходные данные
    std::vector<MechOutGlobalStepData> *out;   // результаты решения задачи

    // вспомогательные процедуры
    // подготовительная процедура для учёта первых краевых условий
    void prepareBC1(const MechGlobalStepParameters &step_el);
    // подготовительная процедура для учёта вторых краевых условий
    void prepareBC2(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl);
    //
    void prepareContact(const Grid::TimeLayers &tl);
    // сборка матрицы и вектора СЛАУ
    void buildGMb(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl);
    // расчёт локальных матриц и векторов (кроме контакта)
    void buildNonContactLocalGMb(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl);
    // расчёт локальных векторов невязки
    void buildLocalR(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl);
    // занесение локальных матриц и векторов в глобальные матрицу и вектор (кроме контакта)
    void addNonContactLocalGMb_to_global(const MechGlobalStepParameters &step_el);
    // сборка локальных константных слагаемых и добавление к глобальным
    void addContact(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl);
    // решение СЛАУ
    void solveSLAU(MechGlobalStepParameters &step_el);
    // проверка условий завершения итераций по пластичности
    //void checkTerminationConditions(const MechGlobalStep &step_el);
    // проверка условий завершения итераций контактной задачи
    //bool checkContactTerminationConditions();

    // основные процедуры решателя
    // инициализация матриц, векторов
    void init(MechTask &mechTask, Threads::Logger_base *set_logger, Threads::SignalToThread_base *set_signal, const Grid::TimeLayers &tl);
    // подготовительные процедуры перед шагом по времени (для динамики)
    void prepareForStep(const MechGlobalStepParameters &step_el);
    // подготовительные процедуры перед итерациями (вычисление функций для вторых краевых условий и начального приближения контактной задачи)
    void prepareForIter(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl);
    // шаг по времени (собственно, итерации)
    void tryStep(MechGlobalStepParameters &step_el, Grid::TimeLayers &tl,
                MechIterationsGeneralResult &mechIterationsGeneralResult);
    // завершение шага по времени (добавка приращений, обновление значений величин)
    void finalizeStep(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl);
    // сохранение найденных значений величин после шага по времени
    void saveResultsOfStep(const bool saveDetailedInf, const int globalStepNumber, const int stepNumber, const Grid::TimeLayers &tl);


    void calcIterResults(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl);
};

}   // namespace Solid
#endif  // SOLID_H
