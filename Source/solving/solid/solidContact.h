#ifndef SOLID_CONTACT_H
#define SOLID_CONTACT_H

#include "elementary.h"
#include "slausolving.h"
#include "slau3x3.h"
#include "grid.h"

#include "solid_base.h"

namespace Solid
{
using namespace Elementary;
struct ContactSurfaceVertexData_FE_Rigid;

// информация о нелинейном процессе для учёта контакта
struct MechContactIterInf
{
    // общее состояние итерационного процесса
    IterationMode iterationMode;
    int iterNumber;                     // номер итерации по контакту
    bool accuracyAchieved;              // true если итерации можно завершить на текущей итерации
     bool improvingAccuracy;                // true если итерации продолжают улучшать невязки (для дожимания)
    // информация о различных изменениях
    int contactChangedNumber;           // количество узлов, которые изменили своё состояние в результате последней итерации (нет контакта -> контакт или контакт -> нет контакта)
    int contactNumber;                  // количество узлов, которые находятся в состоянии контакта (F>0)
    int separationNumber;               // количество узлов, которые выходят из состояния контакта (F>0 -> F=0)
    int GFirstStrChanged;   // номер первой изменённой строки
    // невязки на последней итерации
    double max_endPoint_residual;
    double max_deltaF_residual;
    // обнуление счётчиков
    void clearCounts(const int matrixSize);
    // добавление информации из контактноой вершины в счётчики (вызывается из ContactSurfaceVertexData_FE_Rigid->calcIterResults)
    // учитываются новые данные для итераций itData_new
    void addContactVertexIterResultsToCounts(const ContactSurfaceVertexData_FE_Rigid &cvd, const int iterNumber, const bool constantNormal);
    // расчёт общих результов итерации (accuracyAchieved, improvingState)
    void calcIterResults(const MechGlobalStepParameters &step_el);
};

// способ учёта контакта
enum class ContactType
{
    AugmentedLagrange = 0, // расширенный метод Лагранжа (многократное применение метода штрафа)
    Penalty = 1,           // метод штрафа (1 раз)
};


// контакт КЭ-поверхность - жёсткая поверхность
struct ContactCondition_FE_Rigid
{
    ContactType method;       // метод учёта контакта
    bool constantNormal;      // true - нормаль вычисляется к исходной точке, false - нормаль вычисляется к точке контакта, жёсткость не обновляется
    double w_stiffness;     // множитель, на который домножается жёсткость, взятая из диагонали локальной матрицы жёсткости узла
    int FESurfaceInd;       // индекс КЭ поверхности
    int RigidSurfaceInd;    // индекс жёсткой поверхности
    bool noContactRadiusOptimization;   // использовать оптимизацию
};
// контакт КЭ-поверхность - КЭ-поверхность
struct ContactCondition_FE_FE
{
    double stiffness0;       // начальная контактная жёсткость
    int FEsurfaceInd1;      // индекс 1-й КЭ поверхности
    int FEsurfaceInd2;      // индекс 2-й КЭ поверхности
};





/*
Дано:
набор RigidSurface_base поверхностей (жёстких),
набор FESurface поверхностей,
набор контактов (пара поверхность-поверхность, коэффициент жёсткости)
*/

// по заданному массиву граней решатель составляет массив данных для вершин контактирующих граней

// данные, которые влияют на сборку матриц и векторов
struct ContactVertexIterationData_FE_Rigid
{
    // новое приближение - расчитывается в результате итерации
    bool GLocalChanged;         // true если нужно перестраивать локальную матрицу
    bool contactNew;            // новое состояние
    bool contactNewChanged;     // состояние (contactNew) изменилось (дополнительная информация)
    bool separation;            // true если contact = true и contactNew = false (дополнительная информация)
    double stiffness;           // жёсткость
    POINT3 endPoint;        // точка, в которую предположительно сместится узел при контакте
    VECTOR3 a;              // вектор из начального положения в endPoint
    VECTOR3 normal;         // соответствующая нормаль к поверхности опоры в этой точке - направление реакции опоры
    VECTOR3 F;              // новое приближение - суммарная реакция опоры в узле
    void init()
    {
        GLocalChanged = false;
        contactNew = false;
        contactNewChanged = false;
        separation = false;
        stiffness = 0;
        endPoint.clear();
        a.clear();
        normal.clear();
        F.clear();
    }
    void updateEndPoint(const POINT3 &set_endPoint, const POINT3 &x1)
    {
        endPoint = set_endPoint;
        a = endPoint - x1;
    }
    void updateNormal(const VECTOR3 &set_normal, const std::vector<bool> &bc1_state, const VECTOR3_int DOFindex, const VECTOR3 &localNodeMatrixDiag, const double w_stiffness)
    {
        normal = set_normal;
        // учёт 1-х краевых
        {
            // зануление компонент нормалей в соответствии с 1-ми краевыми
            for(int k = 0; k < 3; k++)
            {
                int global_mk = DOFindex[k];
                if(bc1_state[global_mk]) // есть первое краевое для компоненты k вектора перемещений
                {
                    normal[k] = 0;
                }
            }
            // нормализация нормали
            double new_normal_abs = normal.abs();
            if(new_normal_abs != 0.)
            {
                normal /= new_normal_abs;
            }
        }
        // контактная жёсткость
        {
            double sum = 0;
            for(int k = 0; k < 3; k++)
                sum += fabs(localNodeMatrixDiag[k]*normal[k]);
            stiffness = sum*w_stiffness;
        }
    }
};

// данные в вершине, которая контактирует с поверхностями (расширенный метод Лагранжа / метод штрафа)
struct ContactSurfaceVertexData_FE_Rigid
{
    size_t vertexIndex;            // индекс вершины
    VECTOR3_int DOFindex;    // индексы базисных функций, привязанных к узлу
    size_t si;                     // индекс контактного условия с поверхностью (не может быть -1)
    // данные для итерационного процесса
    ContactVertexIterationData_FE_Rigid it;        // исходные данные текущей итерации
    ContactVertexIterationData_FE_Rigid it_new;    // обновлённые данные для следующей итерации
    VECTOR3 localNodeMatrixDiag;  // элементы локальной узловой матрицы жёсткости после занесения локальных матриц(кроме контактных) в глобальную
     Grid::SurfacePositionData ssd1;     // данные для поиска ближайшей точки (начальное приближение)
     Grid::SurfacePositionData ssd2;
    double h;                        // (дополнительная информация) расстояние до поверхности
    double endPoint_residual;        // (дополнительная информация) = |endPointNew - endPoint|, вычисляебся только для контактирующих узлов
    double deltaF_residual;          // = (deltaF_new - deltaF).abs() / (F_sum.abs() + deltaF_new.abs() + deltaF.abs()), вычисляебся только для контактирующих узлов
    // данные, не меняющиеся при итерациях
    bool contact;           // предыдущее состояние: true - есть контакт, false - нет контакта
    VECTOR3 F_sum;          // суммарная реакция опоры в узле, действовавшая на предыдущем шаге
    double noContactRadius; // внутри сферы с таким радиусом и центром в узле отсутствует контактная поверхность, = -1 если в данный момент есть контакт или отлипание
    // локальная матрица
    double G_local[3][3];
    // локальный вектор
    double b_local[3];
    //int SSCM_index[3][3];
    // инициализация
    void init(const size_t set_vertexIndex, const size_t set_si);
    // обновление индексов базисных функций, привязанных к контактному узлу
    void updateIndexes(const Grid::GlobalDOFs *DOFs);
    // подготовительная процедура (пустая)
    void prepareForIter();
    // обновление значений локальной матрицы и вектора
    // данные для итераций it обновляются только если iterationMode = Free
    void updateLocalGMb(const Grid::Grid3D *grid, const std::vector<Grid::Surface_base *> *rigidSurface, const std::vector<ContactCondition_FE_Rigid> *CCSource_FE_Rigid, const std::vector<bool> &bc1_state, const Grid::GlobalDOFs *DOFs, const MechContactIterInf &nlInf, const IncForsesMode incForsesMode);



    // добавление локальный матрицы к глобальной
    void addLocalGM_to_global1x1(SlauSolving::SSCMElements &elements, SlauSolving::SSCMPortrait &portrait);
    // добавление локального вектора к глобальному
    void addLocalb_to_global1x1(Vector &b);
    // сохранение копий диагональных элементов локальной узловой матрицы жёсткости
    // вызывается после занесения локальных матриц(кроме контактных) в глобальную
    void copyDiagStiffnessesFromG1x1(const SlauSolving::SSCMElements &Gelements);

    // добавление локальный матрицы к глобальной
    void addLocalGM_to_global3x3(SlauSolving::SSCM3x3Elements &elements, SlauSolving::SSCMPortrait &portrait);
    // добавление локального вектора к глобальному
    void addLocalb_to_global3x3(Vector3 &b);
    // сохранение копий диагональных элементов локальной узловой матрицы жёсткости
    // вызывается после занесения локальных матриц(кроме контактных) в глобальную
    void copyDiagStiffnessesFromG3x3(const SlauSolving::SSCM3x3Elements &Gelements);

    // проверка равен ли локальный вектор 0
    bool localb_forces_is_null();



    // расчёт результатов итерации
    // расчитываются ближайшие точки поверхности
    // расчитываются нормали с учётом 1-х краевых условий
    // расчитывается новое приближение (новые данные для итераций itData_new), погрешности и общие результаты итерации
    void calcIterResults1(const Grid::Grid3D *grid, const std::vector<MechVertexData> *vertex, const std::vector<Grid::Surface_base *> *rigidSurface, const std::vector<ContactCondition_FE_Rigid> *CC_FE_Rigid, const int iterNumber, const std::vector<bool> &bc1_state);
    // сбор общей информации
    void calcIterResults2(MechContactIterInf &nlInf, const std::vector<Grid::Surface_base *> *rigidSurface, const std::vector<ContactCondition_FE_Rigid> *CC_FE_Rigid);
    // завершение шага по времени - вычислениеление суммарной силы реакции опоры и статуса наличия контакта
    void finalizeStep(const std::vector<MechVertexData> *vertex);
    // сохранение найденных значений величин после шага по времени
    void saveResultsOfStep(std::vector<MechOutVertexData> &outVertex);
};

}   // namespace Solid

#endif // SOLID_CONTACT_H
