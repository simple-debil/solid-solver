#ifndef SOLID_PLASTIC_H
#define SOLID_PLASTIC_H

#include "elementary.h"
#include "slausolving.h"
#include "slau3x3.h"
#include "grid.h"

#include "solid_base.h"

namespace Solid
{
using namespace Elementary;
struct MechFeData_base;

// информация о нелинейном процессе для учёта пластичности/ползучести
struct MechPlasticIterInf
{
    // общее состояние итерационного процесса
    IterationMode iterationMode;
    int iterNumber;                     // номер итерации по пластичности
    bool accuracyAchieved;              // true если итерации можно завершить на текущей итерации
    bool improvingAccuracy;                // true если итерации продолжают улучшать невязки (для дожимания)
    // информация о различных изменениях
    int indexesChangedNumber;
    int tablesChangedNumber;
     int preliminarily_GLocalChangedNumber;    // количество КЭ, для которых матрица точно меняется (временная переменная)
    int GLocalChangedNumber;    // количество КЭ, для которых матрица меняется
     int GGlobalChangedState;    // = 1 если GLocalChangedNumber>=1 - то есть глобальную матрицу таки пришлось изменить из-за разгрузок
    int bLocalChangedNumber;
    int nonlinearStateFENumber;     // количество КЭ, для которых идёт итерационный процесс
    int unloadStateFENumber;     // количество КЭ, которые разгружаются
    int GFirstStrChanged;   // номер первой изменённой строки
    // невязки на последней итерации
    // ## сходимость по изменению tan2?
    MaxSigmaEpsResidual maxPlasticResidual;
    // невязки на предпоследней итерации
    MaxSigmaEpsResidual prevMaxPlasticResidual;
    // обнуление счётчиков
    void clearCounts(const int matrixSize);
    // добавление информации из конечного элемента в счётчики (вызывается из Fe->calcIterResults)
    void addFeIterResultsToCounts(const MechFeData_base *feEl);
    // расчёт общих результов итерации (accuracyAchieved, improvingState), данные для итераций itData и itData_new не нужны
    void calcIterResults(const MechGlobalStepParameters &step_el);
};

struct LoadingBehavior
{
    bool isActiveLoading;  // активное нагружение
    bool isElasticLoading; // упругое нагружение
    bool isNeutralLoading; // нейтральное нагружение (легендарный случай)
    bool isUnloading;      // разгрузка(если dW < 0, но предел есть превышение предела упругости, то считается, что типа происходит активное нагружение)
    bool readyToUnload; // true если нагружение происходит далеко не перпендикулярно поверхности текучести (cos угла < cos_min)
    void clear()
    {
        isActiveLoading = false;
        isElasticLoading = false;
        isNeutralLoading = false;
        isUnloading = false;
        readyToUnload = false;
    }
};

struct MechFePointIterationData
{
    // параметры итерации
    TENSOR4 C_pl;       // тензор ОС, для алгоритма Соловейчика
    TanAngle_Yeld A; // параметр для построения матрицы определяющих соотношений
                     // A.is_elastic => D_pl = D
    TENSOR2 dsigma0;          // начальные напряжения
    // результаты итерации: искомые приращения, режим нагружения и прочее
    double dq;             // приращение параметра Одквиста
    double dW;             // работа приращения деформации
    TENSOR2 depsElastic;   // приращение упругих деформации
    TENSOR2 depsPlastic;   // приращение пластичной деформации
    TENSOR2 dsigma;        // приращение напряжения
    TENSOR2 dsigma_trial;  // "пробное" приращение напряжения (если считать полное приращение деформации упругой)
    TENSOR2 h;             // направляющий тензор для начальных напряжений, единичный в смысле эвклидовой нормы
    TENSOR2 z;             // направляющий тензор для матрицы C_pl, единичный в смысле эквивалентной деформации
    MaxSigmaEpsResidual residual; // невязки
    LoadingBehavior lb;     // режим нагружения
    double cosTetta;        // угол между девиатором напряжений и "пробным" напряжений
    bool preliminarily_needUpdateMatrix; // true если точно требуется пересобрать локальную матрицу жёсткости для следующей итерации
    bool needIterations; // true если требуются итерации (не зависит от невязок, т.е. указывает на наличие итерационного процесса)
    bool needUpdateMatrix; // true если требуется пересобрать локальную матрицу жёсткости для следующей итерации
    // needUpdateMatrix = preliminarily_needUpdateMatrix или (хотя бы для 1 КЭ придётся менять матрицу и lb.readyToUnload)
    void clear()
    {
        A.is_elastic = true;
        dsigma0.clear();
    }
};

// данные для точек внутри конечного элемента
struct MechFePointData
{
    // изменения, связанные с изменением температуры (вычисляются или пересылаются извне перед 0-й итерацией)
    double T;       // температура (пересылается на вход из температурного решателя)
    double deltaT;  // изменение температуры на текущем шаге по времени (пересылается на вход из температурного решателя)
    TENSOR2 depsTerm;   // приращение температурной деформации
    ElasticParameters elasticParameters_T; // упругие параметры материала при температуре T + deltaT
    // текущее состояние (вычисляется в конце шага, не меняются в процессе итераций)
    double q;              // параметр Одквиста - сумма эквивалентных приращений пластических деформаций
    TENSOR2 sumEpsElastic; // упругая деформация
    TENSOR2 sumEpsPlastic; // пластическая деформация
    TENSOR2 sumEpsTerm;    // температурная деформация
    TENSOR2 sumSigma;      // напряжение
    bool isYelding;  // пластичное состояние, т.е. предел упругости был превышен и после последнего активного нагружения не было разгрузки
     // (вспомогательные данные текущего состояния, вычисляются перед 0-й итерацией - чтобы повторно не вычислять)
     TENSOR2 h0; // направляющий тензор для начальных напряжений, единичный в смысле эвклидовой нормы
     TENSOR2 z0; // направляющий тензор для матрицы C_pl, единичный в смысле эквивалентной деформации
     double sumSigma_eqv; // эквивалентное напряжение в конце прошлого шага
    // параметры текущей итерации (заменяется на it_new перед сборкой матрицы, если итерации активны)
    MechFePointIterationData it;
    // параметры для следующей итерации (вычисляются в конце итерации)
    MechFePointIterationData it_new;
public:
    // основные функции
    // инициализация нулями
    void init()
    {
        q = 0;
        sumEpsElastic.clear();
        sumEpsPlastic.clear();
        sumEpsTerm.clear();
        sumSigma.clear();
        isYelding = false;
    }
    // инициализировать первую итерацию, либо установить данные для итераций (которые были получены в результате предыдущей итерации)
    void initNewIteration(MechMaterialSource &m0, const Grid::TimeLayers &tl, const bool firstIter, const IterationMode iterationMode);
    // расчёт результатов итерации
    void calcIterResults(MechMaterialSource &m0, const VECTOR3 *du, const Grid::TimeLayers &tl);
    // расчёт нового приближения для следующей итерации
    void prepareForNextIter(MechMaterialSource &m0, const int iterNumber, const Grid::TimeLayers &tl, int preliminarily_GLocalChangedNumber);
    // завершение шага по времени - прибавление приращений
    void finalizeStep(const MechMaterialSource &m0, const int movingGridMode);
    // сохранение результатов шага по времени
    void saveResultsOfStep(MechOutFePointData &outFePD);

    //вспомогательные функции
    // проекция на поверхность текучести с учётом упрочнения
    void project(MechMaterialSource &m0, const Grid::TimeLayers &tl,
                 TENSOR2 &dsigma0_corr, double &d, bool &search_crashed);
    void debug(MechMaterialSource &m0) const;
};

// данные конечного элемента
// базовый класс
struct MechFeData_base
{
    // константные таблицы
    // Гаусс-2
    static double Gauss2_integrationwSource[8];
    static double Gauss2_basCubeSource[8*8];
    static VECTOR3 Gauss2_dLinearBasCubeSource[8*8];
    static VECTOR3 Gauss2_dQuadraticBasCubeSource[8*27];
    // Гаусс-3
    static double Gauss3_integrationwSource[27];
    static double Gauss3_basCubeSource[27*8];
    static VECTOR3 Gauss3_dLinearBasCubeSource[27*8];
    static VECTOR3 Gauss3_dQuadraticBasCubeSource[27*27];
    // инициализация константных таблиц для КЭ расчётов
    static void initConstantTables();
    Integration::IntegrationType integrationType;   // способ интегрирования
    bool indexesChanged;    // true если изменилась нумерация узлов сетки и нужно перестраивать таблицы индексов
    bool tablesChanged;     // true если сетка сместилась и нужно перестраивать таблицы множителей для расчёта элементов матриц и векторов
    bool GLocalChanged;     // true если нужно перестраивать локальную матрицу
    bool bLocalChanged;     // true если нужно перестраивать локальный вектор
    //int GFirstStrChanged;   // номер первой изменённой строки
    // пластичность
    MaxSigmaEpsResidual maxPlasticResidual;
    bool plasticNeedIterations;
    bool isUnloading;
    bool needUpdateMatrix;
    // инициализация
    virtual void init(Integration::IntegrationType set_integrationType, const Grid::FE_base *feEl) = 0;
    // задание температуры и прироста температуры, однородно (из температурного решателя)
    virtual void setT(const double T, const double deltaT) = 0;
    // обновление значений локальной матрицы и вектора
    // на 0-й итерации также задаётся eps1Eqv и sigma1Eqv и первое приближение по касательной (даже если iterationMode = ReadOnly)
    // на последующих итерациях данные для итераций itData обновляются только если iterationMode = Free
    virtual void updateLocalGMb(const Grid::Grid3D *grid, const Grid::GlobalDOFs *DOFs, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0, const MechPlasticIterInf &nlInfPlastic, const bool firstIter, const IncForsesMode incForsesMode) = 0;
    // расчёт невязки
     virtual void calcLocalR(const Grid::Grid3D *grid, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0, const MechPlasticIterInf &nlInfPlastic, const IncForsesMode incForsesMode) = 0;
    // добавление значений от второго краевого условия в локальный вектор
    virtual void add_BC2_to_localb(const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechBoundaryCondition2Source_base *bc2SourceEl, const int faceIndex, const Integration::Integrator &integrationFoursquare) = 0;
    // добавление к невязке краевых условий
     virtual void add_BC2_to_localR(const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechBoundaryCondition2Source_base *bc2SourceEl, const int faceIndex, const Integration::Integrator &integrationFoursquare, const IncForsesMode incForsesMode) = 0;



    // добавление глобальных индексов конечного элемента в портрет матрицы
    virtual void addGlobalIndexesToGlobalPortrait1x1(SlauSolving::SSCMPortraitBulder &pb) = 0;
    // обновление индексов (прямых ссылок) для добавления локальной матрицы к глобальной
    virtual void updateSSCM_indexes1x1(const SlauSolving::SSCMPortrait &portrait) = 0;
    // добавление локальный матрицы к глобальной
    virtual void addLocalGM_to_global1x1(SlauSolving::SSCMElements &elements) = 0;
    // добавление локального вектора к глобальному
    virtual void addLocalb_to_global1x1(Vector &b) = 0;


    // добавление глобальных индексов конечного элемента в портрет матрицы
    virtual void addGlobalIndexesToGlobalPortrait3x3(SlauSolving::SSCMPortraitBulder &pb) = 0;
    // обновление индексов (прямых ссылок) для добавления локальной матрицы к глобальной
    virtual void updateSSCM_indexes3x3(const SlauSolving::SSCMPortrait &portrait) = 0;
    // добавление локальный матрицы к глобальной
    virtual void addLocalGM_to_global3x3(SlauSolving::SSCM3x3Elements &elements) = 0;
    // добавление локального вектора к глобальному
    virtual void addLocalb_to_global3x3(Vector3 &b) = 0;

    // проверка равен ли локальный вектор 0
    virtual bool localb_forces_is_null() = 0;
    // получение минимального глобального индекса базисных ф-й
    virtual size_t getMinGlobalIndex()const = 0;



    // расчёт результатов итерации
    // расчитывается новое приближение (новые данные для итераций itData_new), погрешности и общие результаты итерации
    virtual void calcIterResults1(const Vector &q, const Grid::Grid3D *grid, const Grid::GlobalDOFs *DOFs, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0,
                                  std::vector<MechVertexData> *vertex, std::vector<MechVertexData> *vertexForCurvature, MechPlasticIterInf &nlInfPlastic) = 0;
    virtual void calcIterResults2(const Grid::TimeLayers &tl, MechMaterialSource &m0,
                                   MechPlasticIterInf &nlInfPlastic) = 0;
    // завершение шага по времени (добавка приращений, обновление значений величин)
    virtual void finalizeStep(const MechMaterialSource &m0, const int movingGridMode) = 0;
    // сохранение найденных значений величин после шага по времени
    virtual void saveResultsOfStep(const Vector &q, const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechOutFeData &outFeDataEl) = 0;

    virtual void debug(MechMaterialSource &)
    {};
};
// линейный базисные функции, линейное или квадратичное отображение
struct MechFeData_LinearHexagon:
        public MechFeData_base
{
    // переменные таблицы для данного конечного элемента, меняются только при перестройке сетки
    int local_mn[24];       // отсортированные локальные номера базисных функций по возрастанию глобальных номеров
    int local_ik[24];
    int globalIndex[24];    // m*3+i-й элемент соответствует номеру глобальной строки/глобального столбца для m-й базисной функции(0..7) и i-й координате(0..2)
    size_t SSCM_index[300]; // индексы для быстрого занесения локальных матриц в глобальную
    // локальные матрица и вектор
    double GLocalL[300];
    double bLocal[24];
    double RLocal[24];
    inline double &G_local_el(const int &ind)
    {
        return GLocalL[ind];
    }
    inline double &G_local_el(const int &i, const int &j)
    {
        return GLocalL[((i+1)*i/2) + j];
    }
    inline double &b_local_el(const int &i)
    {
        return bLocal[i];
    }
    inline double &R_local_el(const int &i)
    {
        return RLocal[i];
    }
    void G_local_clear();
    void b_local_clear();
    void R_local_clear();
    void updateIndexes(const Grid::FE_base *feEl, const Grid::GlobalDOFs *DOFs);
    virtual void add_BC2_to_localb(const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechBoundaryCondition2Source_base *bc2SourceEl, const int faceIndex, const Integration::Integrator &integrationFoursquare) override;
    virtual void add_BC2_to_localR(const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechBoundaryCondition2Source_base *bc2SourceEl, const int faceIndex, const Integration::Integrator &integrationFoursquare, const IncForsesMode incForsesMode) override;

    virtual void addGlobalIndexesToGlobalPortrait1x1(SlauSolving::SSCMPortraitBulder &pb) override;
    virtual void updateSSCM_indexes1x1(const SlauSolving::SSCMPortrait &portrait) override;
    virtual void addLocalGM_to_global1x1(SlauSolving::SSCMElements &elements) override;
    virtual void addLocalb_to_global1x1(Vector &b) override;

    virtual void addGlobalIndexesToGlobalPortrait3x3(SlauSolving::SSCMPortraitBulder &pb) override;
    virtual void updateSSCM_indexes3x3(const SlauSolving::SSCMPortrait &portrait) override;
    virtual void addLocalGM_to_global3x3(SlauSolving::SSCM3x3Elements &elements) override;
    virtual void addLocalb_to_global3x3(Vector3 &b) override;

    virtual bool localb_forces_is_null() override;
    virtual size_t getMinGlobalIndex()const override;

};
// линейный базисные функции, линейное или квадратичное отображение, однородный материал на КЭ
struct MechFeData_LinearHexagon_homogeny:
        public MechFeData_LinearHexagon
{
    MechFePointData pd;                 // данные только в центре конечного элемента, т.к. он однородный
    double intForG_mnjl[8][8][3][3];    // множители для расчёта элементов матрицы G
    double intForb1_m[8];               // множители для расчёта 1-го слагаемого вектора b
    double intForb23_mj[8][3];          // множители для расчёта 2-го и 3-го слагаемых вектора b
    //double intForM_mn[8][8];
    void updateTables(const Grid::Grid3D *grid, const Grid::FE_base *feEl, const int numPoints, const double *integrationw, const double *basCube, const VECTOR3 *dLinearBasCube, const VECTOR3 *dQuadraticBasCube);
    virtual void init(Integration::IntegrationType set_integrationType, const Grid::FE_base *feEl) override;
    virtual void setT(const double T, const double deltaT) override;
    virtual void finalizeStep(const MechMaterialSource &m0, const int movingGridMode) override;
    virtual void saveResultsOfStep(const Vector &q, const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechOutFeData &outFeDataEl) override;
};
// линейный базисные функции, линейное или квадратичное отображение, однородный материал на КЭ, упругий/нелинейно-упругий или упруго-пластичный материал
struct MechFeData_LinearHexagon_homogeny_E_NLE_EP:
        public MechFeData_LinearHexagon_homogeny
{
    virtual void updateLocalGMb(const Grid::Grid3D *grid, const Grid::GlobalDOFs *DOFs, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0, const MechPlasticIterInf &nlInfPlastic, const bool firstIter, const IncForsesMode incForsesMode) override;
    virtual void calcLocalR(const Grid::Grid3D *grid, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0, const MechPlasticIterInf &nlInfPlastic, const IncForsesMode incForsesMode) override;
    virtual void calcIterResults1(const Vector &q, const Grid::Grid3D *grid, const Grid::GlobalDOFs *DOFs, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0,
                                  std::vector<MechVertexData> *vertex, std::vector<MechVertexData> *vertexForCurvature, MechPlasticIterInf &nlInfPlastic) override;
    virtual void calcIterResults2(const Grid::TimeLayers &tl, MechMaterialSource &m0,
                                   MechPlasticIterInf &nlInfPlastic) override;
    virtual void debug(MechMaterialSource &m0) override;
};


// линейный базисные функции, линейное отображение (все базисные ф-и задаются интерполянтами)
// ##buildLocalR стала пустой, R строится в updateLocalGMb (поэтому работает только режим IncForsesMode::MinusIntegral)
struct MechFeData_LinearXFEM_gen:
        public MechFeData_base
{
    // переменные таблицы для данного конечного элемента, меняются только при перестройке сетки
    std::vector<Grid::LocalDOF_ID> local_DOF_ID;  // локальная нумерация степеней свободы
    std::vector<int> local_mn;       // отсортированные локальные номера базисных функций по возрастанию глобальных номеров
    std::vector<int> local_ik;
    std::vector<int> globalIndex;    // m*3+i-й элемент соответствует номеру глобальной строки/глобального столбца для m-й базисной функции(0..7) и i-й координате(0..2)
    std::vector<size_t> SSCM_index;  // индексы для быстрого занесения локальных матриц в глобальную
    // множители для интегрирования по 6-гранным подобластям
    //std::vector<Grid::HexagonIntTable> hexagonIntTable;
    // локальные матрица и вектор
    std::vector<double> GLocalL;
    std::vector<double> bLocal;
    std::vector<double> RLocal;
    // данные в центрах 6-гранных подобластей КЭ
    std::vector<MechFePointData> pd;
    inline double &G_local_el(const int &ind)
    {
        return GLocalL[ind];
    }
    inline double &G_local_el(const int &i, const int &j)
    {
        return GLocalL[((i+1)*i/2) + j];
    }
    inline double &b_local_el(const int &i)
    {
        return bLocal[i];
    }
    inline double &R_local_el(const int &i)
    {
        return RLocal[i];
    }
    void G_local_clear();
    void b_local_clear();
    void R_local_clear();
    void updateIndexes(const Grid::FE_base *feEl, const Grid::GlobalDOFs *DOFs);
    virtual void add_BC2_to_localb(const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechBoundaryCondition2Source_base *bc2SourceEl, const int subSurfaceIndex, const Integration::Integrator &integrationFoursquare) override;
    virtual void add_BC2_to_localR(const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechBoundaryCondition2Source_base *bc2SourceEl, const int faceIndex, const Integration::Integrator &integrationFoursquare, const IncForsesMode incForsesMode) override;

    virtual void addGlobalIndexesToGlobalPortrait1x1(SlauSolving::SSCMPortraitBulder &pb) override;
    virtual void updateSSCM_indexes1x1(const SlauSolving::SSCMPortrait &portrait) override;
    virtual void addLocalGM_to_global1x1(SlauSolving::SSCMElements &elements) override;
    virtual void addLocalb_to_global1x1(Vector &b) override;

    virtual void addGlobalIndexesToGlobalPortrait3x3(SlauSolving::SSCMPortraitBulder &pb) override;
    virtual void updateSSCM_indexes3x3(const SlauSolving::SSCMPortrait &portrait) override;
    virtual void addLocalGM_to_global3x3(SlauSolving::SSCM3x3Elements &elements) override;
    virtual void addLocalb_to_global3x3(Vector3 &b) override;

    virtual bool localb_forces_is_null() override;
    virtual size_t getMinGlobalIndex()const override;

    void updateTables(const Grid::FE_base *feEl);
    virtual void init(Integration::IntegrationType set_integrationType, const Grid::FE_base *feEl) override;
    virtual void setT(const double T, const double deltaT) override;
    virtual void finalizeStep(const MechMaterialSource &m0, const int movingGridMode) override;
    virtual void saveResultsOfStep(const Vector &q, const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechOutFeData &outFeDataEl) override;

    virtual void updateLocalGMb(const Grid::Grid3D *grid, const Grid::GlobalDOFs *DOFs, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0, const MechPlasticIterInf &nlInfPlastic, const bool firstIter, const IncForsesMode incForsesMode) override;
    virtual void calcLocalR(const Grid::Grid3D *grid, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0, const MechPlasticIterInf &nlInfPlastic, const IncForsesMode incForsesMode) override;
    virtual void calcIterResults1(const Vector &q, const Grid::Grid3D *grid, const Grid::GlobalDOFs *DOFs, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0,
                                  std::vector<MechVertexData> *vertex, std::vector<MechVertexData> *vertexForCurvature, MechPlasticIterInf &nlInfPlastic) override;
    virtual void calcIterResults2(const Grid::TimeLayers &tl, MechMaterialSource &m0,
                                   MechPlasticIterInf &nlInfPlastic) override;
    void moveSubVertexes(const Grid::Grid3D *grid, const Grid::GlobalDOFs *DOFs, const Vector &dq, Grid::FE_base *feEl);
};
}   // namespace Solid

#endif // SOLID_GENERAL_H
