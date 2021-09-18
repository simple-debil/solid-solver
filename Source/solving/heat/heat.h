/* --------------------------------------------------------- */
// МАТЕРИАЛЫ, КРАЕВЫЕ УСЛОВИЯ ДЛЯ ЗАДАЧИ ТЕПЛОПРОВОДНОСТИ
/* --------------------------------------------------------- */

#ifndef HEAT_H
#define HEAT_H

#include <vector>

#include "elementary.h"
#include "slausolving.h"
#include "integration.h"
#include "grid.h"
#include "fem.h"
 #include "funParser.h"
#include "logger_base.h"


namespace Heat
{
using namespace Elementary;

// тип базисных функций
#define BasisType_1L                                    1

#define Solid_FeData_HomogenyMode_homogeny              0
#define Solid_FeData_HomogenyMode_not_homogeny          1

// типы вторых краевых условий
enum class ThermBoundaryCondition2Type
{
    none,               // отсутствует
    scalar,             // скаляр
};

// материал
struct ThermMaterialSource
{
    double ro;  // плотность
    MATR3x3 L;  // тензор теплопроводности
    double c;   // удельная теплоёмкость
     double f;   // мощность внутренних объёмных источников (стоков) тепла
};

// 1-е краевое условие
struct ThermBoundaryCondition1Source
{
    int mode;           // -1 - отсутствует
    double T0;          // температура
};

// 2-е краевое условие
struct ThermBoundaryCondition2Source
{
    ThermBoundaryCondition2Type mode;
    double q;           // текущая плотность набегающего теплового потока
    double hi;          // коэффициент конвективного теплообмена
    double Ta;          // температура окружающей среды
};

// параметры глобального шага
struct ThermGlobalStep
{
    // решение СЛАУ
    SlauSolving::SolverParameters slausolverParameters;        // параметры решателя СЛАУ
    // метод дискретизации по времени
    int timeMode;  // 0 - квазистатический, 1 - двухточечный неявный
    std::vector<ThermMaterialSource> *material;             // параметры материалов
    std::vector<ThermBoundaryCondition1Source> *bc1SourceT; // 1-е краевые
    std::vector<ThermBoundaryCondition2Source> *bc2SourceT; // 2-е краевые
};

// данные для вершин
struct ThermVertexData
{
    double T;           // текущая температура
    double newT;        // температура после шага по времени
};

// данные для точки внутри конечного элемента
struct ThermFePointData
{
    double T;           // текущая температура
    double newT;        // температура после шага по времени
};

// данные в точках Гаусса внутри конечного элемента
struct ThermFeData
{
    int basisType;        // тип базисных функций
                          // 1 - линейные базисные функции
    int homogenyMode;     // однороден ли элемент
                          // 0 - температура расчитывается только в центрах
                          // 1 - температура расчитывается в точках Гаусса
    Integration::IntegrationType integrationType;   // способ интегрирования
    std::vector<ThermFePointData> pd;               // данные для точек внутри конечного элемента
    void init(const int setBasisType, const int setHomogenyMode, const Integration::IntegrationType setIntegrationType, double T0)
    {
        basisType = setBasisType;
        homogenyMode = setHomogenyMode;
        integrationType = setIntegrationType;
        size_t size = 0;
        if(homogenyMode == Solid_FeData_HomogenyMode_homogeny)
        {
            size = 1;
        }
        if(homogenyMode == Solid_FeData_HomogenyMode_not_homogeny)
        {
            switch (integrationType)
            {
            case Integration::IntegrationType::Gauss2:
                size = Integration::Integration3D_Size_Gauss2;
                break;
            case Integration::IntegrationType::Gauss3:
                size = Integration::Integration3D_Size_Gauss3;
                break;
            case Integration::IntegrationType::Gauss5:
                size = Integration::Integration3D_Size_Gauss5;
                break;
            }
        }
        pd.resize(size);
        for (size_t pInd = 0; pInd < pd.size(); pInd++) // индекс точки Гаусса внутри конечного элемента
        {
            ThermFePointData &fePD = pd[pInd];
            fePD.T = T0;
            fePD.newT = T0;
        }
    }
};


// информация об итерации нелинейной задачи
struct ThermNonlinearInf
{
    // информация о решении СЛАУ
    double slauResidual;
    double slauRelativeResidual;
    int slauIterations;
    double slauTime;
};

// общая информация о шаге по времени
struct ThermOutStepData
{
    double t;           // текущее время (после шага)
    double time;        // время расчёта
};

// результаты в узлах сетки
struct ThermOutVertexData
{
    double T;      // температура
};

// результаты в точке шестигранника
struct ThermOutFePointData
{
    POINT3 p;      // координаты точки
    double T;      // температура
};

// результаты в точках Гаусса шестигранника
struct ThermOutFeData
{
    std::vector<ThermOutFePointData> pd;
    void init(const int setSize)
    {
        pd.resize(setSize);
    }
};


// входные данные
struct ThermTask
{
    bool enabled;            // инициализирован?
    Grid::Grid3D *grid;                                     // сетка
    std::vector<Grid::GlobalStep> *step;                          // общие параметры глобальных шагов
    std::vector<ThermGlobalStep> *thermStep;                // параметры глобальных шагов
    std::vector<ThermVertexData> *vertex;                   // начальная температура
    std::vector<ThermFeData> *fe;                           // начальная температура
};

// выходные данные
struct ThermOutGlobalStepData
{
    std::vector<ThermOutStepData> step;
    std::vector<ThermOutVertexData> vertex;
    std::vector<ThermOutFeData> fe;
};

// решатель задачи теплопроводности
struct ThermSolver: public ThermTask, public Threads::Logging
{
    bool firstStep;                             // true если матрицы ниразу не строились
    // размер матрицы
    int matrixSizeT;
    // инфа для итераций
    ThermNonlinearInf nlInf;
    // матрицы
    SlauSolving::SSCM GmatrixT;                          // матрица G
    SlauSolving::SSCMBulder GbulderT;                    // сборщик
    //SlauSolving::SSCMPreconditioner GpreconditionerT;    // предобусловливатель
    //SlauSolving::SSCMSolver GsolverT;                    // решатель СЛАУ
    SlauSolving::SSCMPreconditioner_base *GpreconditionerT;    // предобусловливатель
    SlauSolving::SSCMSolver_base *GsolverT;                    // решатель СЛАУ

    SlauSolving::SSCM MmatrixT;                          // матрица M
    SlauSolving::SSCMBulder MbulderT;
    // векторы
    Vector bT;
    Vector qT;
    Vector q1T;
    Vector MqT;
    // вспомогательные данные для учёта первых краевых условий
    std::vector<double> bc1_u0T;
    std::vector<bool> bc1_stateT;
    // выходные данные
    std::vector<ThermOutGlobalStepData> *out;   // результаты решения задачи

    void init(ThermTask &thermTask, Threads::Logger_base *set_logger, Threads::SignalToThread_base *set_signal);
    void prepareForStep();
    void tryStep(const int globalStepNumber, Grid::TimeLayers &tl);
    void finalizeStep();
    void saveResultsOfStep(const int globalStepNumber, const Grid::TimeLayers &tl);
};

// вспомогательные процедуры
void buldGMb_T(const Grid::Grid3D *grid, const std::vector<ThermFeData> *fe, std::vector<ThermMaterialSource> *material,
               SlauSolving::SSCMBulder &GbulderT, SlauSolving::SSCMBulder &MbulderT, Vector &bT);
void addBc2_T(const Grid::Grid3D *grid, const std::vector<ThermBoundaryCondition2Source> *bc2Source,
              Vector &b, SlauSolving::SSCMBulder &GbulderT);
}   // namespace Solid
#endif // HEAT_H
