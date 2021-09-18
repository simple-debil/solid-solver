/* --------------------------------------------------------- */
// СЕТКИ
/* --------------------------------------------------------- */

#ifndef GRID_H
#define GRID_H

#include <set>

#include "elementary.h"
#include "integration.h"

#include "grid_base.h"
#include "regionIndexation.h"


// структуры и классы для использования сеток
namespace Grid
{
using namespace Elementary;

// параметры глобального шага
struct GlobalStep
{
    double t_start;      // начальное время (имеет значение только для первого глобального шага)
    double t_finish;     // конечное время
    double dt0;          // начальный  шаг по времени
};

// временные слои
struct TimeLayers
{
    double dt;      // текущий шаг по времени
    double t0;      // время в конце текущего шага
    double t_prev1;      // время в начале текущего шага (1 шаг назад)
     double t_prev2;      // время 2 шага назад
     double t_prev3;      // время 3 шага назад
    // шаг по времени: t_prev1 -> t0
};

// тип конечного элемента
enum class FEType
{
    LinearHexagon = 0,
    QuadraticHexagon = 1,
    LinearHexagon_XFEM = 2,
};

// шестигранный конечный элемент, заданный индексами вершин и индексом материала
struct FE_base
{
    int mi;         // индекс материала

    virtual FEType get_FEType()const = 0;
    // получение глобальных индексов вершин КЭ
    // если порядок отображения выше порядка базисных функций, то индексы дополнительных вершин ТОЖЕ ВОЗВРАЩАЮТСЯ,
    // но индексы этих вершин не соответствуют индексам глобальных базисных функций
    virtual void getGeomVertexIndexes(int *vertexIndexes)const = 0;
    // обновление глобальных индексов вершин КЭ
    // если порядок отображения выше порядка базисных функций, то индексы дополнительных вершин ТОЖЕ ОБНОВЛЯЮТСЯ,
    // но индексы этих вершин не соответствуют индексам глобальных базисных функций
    virtual void setGeomVertexIndexes(const int *vertexIndexes) = 0;
    // получение вершин КЭ
    // если порядок отображения выше порядка базисных функций, то дополнительные вершины ТОЖЕ ВОЗВРАЩАЮТСЯ
    virtual void getGeomVertexes(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature,
                             POINT3 *v)const = 0;

    // получение глобальных индексов вершин КЭ
    // если порядок отображения выше порядка базисных функций, то дополнительные вершины игнорируются
    virtual void getVertexIndexes(int *vertexIndexes)const = 0;
    // обновление глобальных индексов вершин КЭ
    // если порядок отображения выше порядка базисных функций, то дополнительные вершины игнорируются
    virtual void setVertexIndexes(const int *new_vertexIndexes) = 0;


    // получение индексов вершин некоторой грани, заданной индексом
    // если порядок отображения выше порядка базисных функций, то дополнительные вершины игнорируются
    virtual int getFaceVertexIndexes(const int faceInd,
                                      int *vi)const = 0;
    // получение локальных и глобальных индексов вершин некоторой грани, заданной индексом
    // если порядок отображения выше порядка базисных функций, то дополнительные вершины игнорируются
    virtual int getFaceVertexIndexes(const int faceInd,
                                      int *vi, int *vi_local)const = 0;
    // получение вершин некоторой грани, заданной индексом
    // если порядок отображения выше порядка базисных функций, то дополнительные вершины игнорируются
    virtual void getFaceVertexes(const std::vector<POINT3> &vertex, const int faceInd,
                                 POINT3 *v)const = 0;
    // получение сферы, внутри которой находится грань с индексом faceInd
    virtual void getFaceOuterSphere(const std::vector<POINT3> &vertex, const int faceInd,
                                    Sphere &s)const = 0;
    // получение нормали в некоторой точке заданной грани
    virtual void getFaceNormal(const POINT3 *v, const int faceIndex, const POINT2_CUBE p,
                               VECTOR3 &normal)const = 0;
    // p - заданная точка
    // nearestPoint - ближайшая к точке p точка грани
    // normal - нормаль наружу в точке nearestPoint
    // h - высота заданной точки над гранью(|p-nearestPoint|), положительная если точка находится с внешней стороны поверхности
    // возвращает true, если nearestPoint внутри грани (т.е. normal || p-nearestPoint)
    virtual bool project(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const int faceInd, const POINT3 p,
                         POINT3 &nearestPoint, VECTOR3 &normal, double &h)const = 0;

    // заполнение таблиц окончательных множителей для численного интегрирования
    // и производных базисных функций в точках интегрирования
    // входные данные:
    // numPoints - количество точек численного интегрирования,
    // integrationw - множителями для численного интегрирования,
    // dLinearBasCube - производные базисных функций на шаблонном КЭ,
    // если порядок отображения выше порядка базисных функций, то dQuadraticBasCube - производные базисных функций отображения на шаблонном КЭ
    // выходные данные:
    // w - множители для численного интегрирования, домноженные на значение модулей детерминанта отображения в точках интегрирования
    // dbas - производные базисных функций в точках интегрирования
    virtual void calcBasisFuncValues(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const int numPoints, const double *integrationw, const VECTOR3 *dLinearBasCube, const VECTOR3 *dQuadraticBasCube,
                                      double *w, double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3])const = 0;
    // заполнение таблиц окончательных множителей для 2-го краевого (реализация для шестигранников)
    virtual void calcBc2BasisFuncValues(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const Integration::Integrator &integrationFoursquare, const int faceIndex,
                                         double *basCube, double *hexagonCoef, VECTOR3 *hexagonNormal)const = 0;
    // получение координат point(x,y,z) точки из заданного шестигранника по координатам cubePoint(X,Y,Z) шаблонного куба
    // берутся производные dif
    // (массив v должен быть получен при помощи getGeomVertexes)
    virtual void cubeToHexagon(const POINT3 *v, const POINT3_CUBE &cubePoint, const DIF_STATE3 &dif,
                               POINT3 &point)const = 0;
};
struct FE_LinearHexagon:
    virtual public FE_base
{
    int vi[8];      // массив глобальных индексов вершин
    virtual FEType get_FEType()const override;
    virtual void getGeomVertexIndexes(int *vertexIndexes)const override;
    virtual void setGeomVertexIndexes(const int *new_vertexIndexes) override;
    virtual void getGeomVertexes(const std::vector<POINT3> &vertex, const std::vector<POINT3> &,
                             POINT3 *v)const override;

    virtual void getVertexIndexes(int *vertexIndexes)const override;
    virtual void setVertexIndexes(const int *new_vertexIndexes) override;

    virtual void calcBasisFuncValues(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const int numPoints, const double *integrationw, const VECTOR3 *dLinearBasCube, const VECTOR3 *,
                                      double *w, double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3])const override;
    virtual void calcBc2BasisFuncValues(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const Integration::Integrator &integrationFoursquare, const int faceIndex,
                                         double *basCube, double *hexagonCoef, VECTOR3 *hexagonNormal)const override;
    virtual void cubeToHexagon(const POINT3 *v, const POINT3_CUBE &cubePoint, const DIF_STATE3 &dif, POINT3 &point)const override;

    virtual int getFaceVertexIndexes(const int faceInd,
                                      int *vertexIndexes)const override;
    virtual int getFaceVertexIndexes(const int faceInd,
                                      int *vertexIndexes_global, int *vertexIndexes_local)const override;
    virtual void getFaceVertexes(const std::vector<POINT3> &vertex, const int faceInd,
                                 POINT3 *v)const override;
    virtual void getFaceOuterSphere(const std::vector<POINT3> &vertex, const int faceInd,
                                    Sphere &s)const override;
    virtual void getFaceNormal(const POINT3 *v, const int faceIndex, const POINT2_CUBE p,
                               VECTOR3 &normal)const override;
    virtual bool project(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const int faceInd, const POINT3 p,
                         POINT3 &nearestPoint, VECTOR3 &normal, double &h)const override;

};
struct FE_QuadraticHexagon:
    virtual public FE_LinearHexagon
{
    int vi_curvature[19]; // массив дополнительных индексов вершин для квадратичного отображения
    virtual FEType get_FEType()const override;
    virtual void getGeomVertexIndexes(int *vertexIndexes)const override;
    virtual void setGeomVertexIndexes(const int *new_vertexIndexes)override;
    virtual void getGeomVertexes(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature,
                             POINT3 *v)const override;
    virtual void calcBasisFuncValues(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const int numPoints, const double *integrationw, const VECTOR3 *dLinearBasCube, const VECTOR3 *dQuadraticBasCube,
                              double *w, double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3])const override;
    virtual void cubeToHexagon(const POINT3 *v, const POINT3_CUBE &cubePoint, const DIF_STATE3 &dif,
                               POINT3 &point)const override;
};


// грань конечного элемента (индексы)
struct FEFace
{
    int feIndex;        // индекс кэ, которому принадлежит грань
    int faceIndex;      // индекс грани (см. HexagonFace[][])
};

// 1-е краевое условие (индексы)
struct BoundaryCondition1
{
    int vertexIndex;        // индекс вершины
    int bc1SourceIndex;     // индекс ресурса 10го краевого условия в этой вершине
};

// 2-е краевое условие (индексы):
struct MechBoundaryCondition2
{
    int FEsurfaceInd;       // индекс КЭ поверхности (массив: номер КЭ + номер грани)
    int bc2SourceIndex;     // индекс ресурса 2-го краевого (поверхностные силы)
};

// 2-е краевое условие для поверхности трещины (индексы):
struct MechBoundaryCondition2Crack
{
    size_t FEIndex;            // индекс КЭ (FE_LinearHexagon_XFEM)
    size_t subSurfaceIndex;    // индекс поверхности внутри КЭ (0 или 1)
    size_t bc2SourceIndex;     // индекс ресурса 2-го краевого (поверхностные силы)
};

// КЭ-поверхность, заданная набором граней сетки
struct FiniteElementSurface:
        public Surface_base
{
    RegionWithSpheres *r0;           // корневая область
    std::vector<FEFace> face;
    const Grid3D *grid;   // сетка
    bool buildRegionsOnce;// строить декомпозицию пространства только 1 раз (обновлять не нужно), то есть поверхность неподвижна
    void init(const bool set_buildRegionsOnce, const Grid3D *set_grid);
    // обновление декомпозиции пространства
    void buildRegions();
    // нахождение индексов тел, которые могут пересекаться с заданной сферой
    void findIntersectingBodies(const Sphere &s, std::vector<size_t> &bodyInd);
    virtual void update(const double, const RegionDecompositionParameters *set_param)override;
    virtual void findNearestPoint(const POINT3 &p, const VECTOR3 &u, const SurfacePositionData &prevSolutionData,
                                  SurfacePositionData &solutionData, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const override;
    // ## не доделано
    virtual void findIntersection(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &prevSolutionData1, const Grid::SurfacePositionData &prevSolutionData2,
                                  Grid::SurfacePositionData &solutionData1, Grid::SurfacePositionData &solutionData2, POINT3 &intersectionPoint, POINT3 &normal, int &side, bool &found)const override;
};




// тип базисной ф-и
enum class FuncType: size_t
{
    linear = 0,
    H = 1,
    tip1 = 2,
    tip2 = 3,
    tip3 = 4,
    tip4 = 5,
    tip_ep1 = 6,
    tip_ep2 = 7,
    tip_ep3 = 8,
    tip_ep4 = 9,
    tip_ep5 = 10,
    tip_ep6 = 11,
    _SIZE = 12
};

// идентификатор базисной ф-и (тип + привязка к конкретной поверхности трещины)
struct FuncID
{
    FuncType type; // тип базисной ф-и
    size_t crackIndex;   // индекс трешины
    inline bool operator==(const FuncID &f)const
    {
        return type == f.type && crackIndex == f.crackIndex;
    }
};
// данные в вершине подобласти
struct SubVertexData
{
    double lagrFuncsValue[8];        // значения ф-й формы Лагранжа в вершинах
    std::vector<double> funcValue;   // значения базисных ф-й, заданных в вершине, funcValue[funcsIDIndex] - значение базисной ф-и с индексом funcsID[funcsIDIndex]
};
// 6-гранная подобласть, заданная индексами вершин
struct SubHex
{
    size_t vi[8];   // индексы вершин (из subVertex)
};
// типы поверхностей
enum class CrackSurfaceType
{
    none = -1,
    Analitical_vertical_x0 = 0,
};
// поверхность трещины
// базовый класс
struct CrackSurface_base
{
    CrackSurfaceType CracksurfaceType; // тип поверхности
    // тангенсальная ф-я уровня
    virtual double ls_tan(const POINT3 &p)const = 0;
    // нормальная ф-я уровня
    virtual double ls_normal(const POINT3 &p)const = 0;
    // расчёт обогащающих ф-й с индексами funcsSource[i].funcIndex, результаты помещаются в subVertexData_el.enrFuncValue[funcsSource[i].funcIndex]
    virtual void solveFuncValues(const POINT3 &p, const std::vector<FuncID> &funcsID, const size_t crackIndex,
                                 SubVertexData &subVertexData_el)const = 0;
};
// вертикальная трещина x = x0, y = y1...y2
struct CrackSurface_Analitical_vertical_x0:
    virtual public CrackSurface_base
{
    double x0;
    double y1;
    double y2;
    double n;
    void init(const double &set_x0, const double &set_y1, const double &set_y2, const double set_n);
    virtual double ls_tan(const POINT3 &p)const override final;
    virtual double ls_normal(const POINT3 &p)const override final;
    virtual void solveFuncValues(const POINT3 &p, const std::vector<FuncID> &funcsID, const size_t crackIndex,
                                 SubVertexData &subVertexData_el)const override final;
};
// Локальная нумерация степеней свободы
struct LocalDOF_ID
{
    size_t m;     // локальный индекс узла, в котором задана базисная ф-я (0...7)
    size_t funcID_index;    // индекс интерполянта, funcsID[funcIDIndex] - ID базисной ф-и
    size_t DOFs_index;      // DOFs[vi[m]][DOFs_index].DOF_index - глобальные степени свободы (по x, y, z), (DOFs[vi[m]][DOFs_index].funcID - ID базисной ф-и)
};

// Нумерация базисных ф-й
struct VertexDOF
{
    Grid::FuncID funcID;     // идентификатор базисной ф-и
    VECTOR3_int DOF_index;   // индексы степеней свободы для базисной ф-и funcID (по x, y, z)
};
// Для каждого узла задаются степени свободы
// первыми должны добавляться ф-и формы Лагранжа
struct GlobalDOFs
{
    size_t DOF_count;    // количество степеней свободы в сумме
    std::vector<std::vector<VertexDOF>> DOFs; // DOFs[vertexIndex][i].DOF_index - индексы степеней свободы для базисной ф-и DOFs[vertexIndex][i].type
    void init(const size_t set_size)
    {
        DOFs.resize(set_size);
        for(size_t i = 0; i < set_size; i++)
        {
            DOFs[i].clear();
        }
        DOF_count = 0;
    }
    // добавление степеней свободы (по x, y, z) для базисной ф-и type в узле vertexIndex
    void addDOFs(const size_t vertexIndex, const Grid::FuncID &funcID)
    {
        size_t sz = DOFs[vertexIndex].size();
        // дважды одну и ту же базисную ф-ю не добавляем
        for(size_t funcSourceIndex = 0; funcSourceIndex < sz; funcSourceIndex++)
        {
            if(DOFs[vertexIndex][funcSourceIndex].funcID == funcID)
                return;
        }
        // добавление степеней свободы для узла
        DOFs[vertexIndex].resize(sz + 1);
        DOFs[vertexIndex][sz].funcID = funcID;
        for(int j = 0; j < 3; j++)
        {
            DOFs[vertexIndex][sz].DOF_index[j] = DOF_count;
            DOF_count++;
        }
    }
    const VECTOR3_int findDOFIndex(const size_t vertexIndex, const Grid::FuncID &funcID) const
    {
        for(size_t i = 0; i < DOFs[vertexIndex].size(); i++)
        {
            if(DOFs[vertexIndex][i].funcID == funcID)
            {
                return DOFs[vertexIndex][i].DOF_index;
            }
        }
        return VECTOR3_int(-1, -1, -1); //# ошибка
    }
    const VECTOR3_int &findDOFIndex_Lagr1(const size_t &vertexIndex) const
    {
        return DOFs[vertexIndex][0].DOF_index;
    }
};

// грань конечного элемента (индексы)
struct SubHexFace
{
    int hexInd;      // индекс 6-гранной подобласти, которой принадлежит грань
    int faceIndex;   // индекс грани (см. HexagonFace[][])
};


struct HexagonIntTable;

struct FE_LinearHexagon_XFEM:
    virtual public FE_LinearHexagon
{
    // при построении сетки задаются crackIndex и массивы funcsID, subVertex, subHex
    static const size_t DOF_max = 256;//64;
    static const size_t Vertexes_max = 1024;
    // идентификаторы базисных ф-й, которые заданы на КЭ (будут строится их интерполянты в subVertexData)
    std::vector<FuncID> funcsID;       // funcsID[funcsIDIndex] - ID базиснойой ф-и
    // глобальные координаты вершин подобластей
    std::vector<POINT3> subVertex;     // subVertex[subHex[hexInd].vi[i]] - глобальные координаты i-й вершины 6-гранника hexInd
    // перемещения глобальных координат вершин подобластей
    std::vector<POINT3> dsubVertex;
    // индексы вершин 6-гранных подобластей (из subVertex)
    std::vector<SubHex> subHex;
    // данные в вершинах подобластей
    std::vector<SubVertexData> subVertexData;   // subVertexData[subHex[hexInd].vi[i]] - данные в i-й вершине 6-гранника hexInd (т.е. в вершине с глобальными координатами subVertex[subHex[hexInd].vi[i]])
    // индексы вершин подобластей, которые соответствуют вершинам КЭ (т.е. в вершиных, координаты которых выдаёт ф-я getGeomVertexes)
    size_t FESubVertexIndex[8];
    // поверхности трещины
    std::vector<SubHexFace> subHexFace[2]; // subHexFace[0/1][i] - грань 6-гранной подобласти


    // 3) универсальный случай (все базисные ф-и задаются интерполянтами)
    // инициализация координат вершин шестигранных подобластей будет задаваться при генерации сетки препроцессором
    // subVertex_XYZ[subHex[hexInd].vi[i]] - локальные координаты i-й вершины 6-гранника hexInd
    void gn_initSubVertexes(const std::vector<CrackSurface_base *> &cs, const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const std::vector<POINT3_CUBE> &subVertex_XYZ);
    // получение глобальных координат вершин шестигранной подобласти с индексом hexInd
    void gn_getSubVertexes(const size_t hexInd,
                           POINT3 (&subv)[8]) const;
    // получение вершин грани faceInd шестигранной подобласти hexInd
    void gn_getSubFaceVertexes(const size_t hexInd, const int faceInd,
                               POINT3 (&v)[4]) const;

    void gn_calcBc2BasisFuncValues(const std::vector<LocalDOF_ID> &DOF_ID, const Integration::Integrator &integrationFoursquare, const size_t hexInd, const int faceIndex,
                                   double *basCube, double *hexagonCoef, VECTOR3 *hexagonNormal)const;
    // рассчёт таблицы для численного интегрирования по подобласти
    void gn_calcIntegralsTable(const double (&Gauss3_integrationwSource)[27], const VECTOR3 (&Gauss3_dLinearBasCubeSource)[27*8], const std::vector<LocalDOF_ID> &DOF_ID, const size_t hexInd, Grid::HexagonIntTable &hexagonIntTable) const;
    // сдвиг сетки (q[m][i] - значение коэффициента перед базисной ф-й m, координата i)
    void gn_moveSubVertexes(const std::vector<LocalDOF_ID> &DOF_ID, const VECTOR3 (&q_local)[FE_LinearHexagon_XFEM::DOF_max]);
    // (вспомогательная) рассчёт значений базисных функций в точке cubePoint подобласти hexInd
    void gn_calcSubBasFuncs(const std::vector<LocalDOF_ID> &DOF_ID, const size_t hexInd, const POINT3_CUBE &cubePoint,
                            double (&bas)[FE_LinearHexagon_XFEM::DOF_max]) const;
    // рассчёт значений производных базисных функций в точке cubePoint подобласти hexInd
    void gn_calcSubDifBasFuncs(const std::vector<LocalDOF_ID> &DOF_ID, const size_t hexInd, const POINT3_CUBE &cubePoint,
                               double (&dbas)[FE_LinearHexagon_XFEM::DOF_max][3]) const;

    virtual FEType get_FEType()const override;
};

// множители для интегрирования (произвольное количество базисных ф-й)
struct HexagonIntTable
{
    double intForG_mnjl[FE_LinearHexagon_XFEM::DOF_max][FE_LinearHexagon_XFEM::DOF_max][3][3]; // множители для расчёта элементов матрицы G
    double intForb23_mj[FE_LinearHexagon_XFEM::DOF_max][3];        // множители для расчёта 2-го и 3-го слагаемых вектора b
};






// область-шестигранник
struct HexagonArea
{
    std::array<int, 8> vi;      // индексы вершин
    std::array<int, 3> N;       // количество промежутков на ребре (т.е. узлов на 1 больше)
    double condensation_q;      // коэффициент сгущения
    int condensation_coord;     // локальная координата сгущения (0|1|2) или -1 если сгущение не требуется
    //std::array<double, 3> q;
    std::array<int, 6> surfaceIndex;// индексы поверхностей для каждой грани (-1 - поверхность не задана)
    std::array<int, 6> bc1Index;// индексы 1-х краевых для каждой грани (-1 - не заданы)
    int mi;         // индекс материала
    // узел внутри 6-гранника по его координатам (v - массив из 8 элементов)
    void calcNodeCoord(const POINT3 *v, const VECTOR3_uint &i, POINT3 &x)const;
    // вычисление состояний нахождения узла на каждой из 6-ти поверхностей
    // (faceState - массив из 6 элементов)
    void calcFaceStates(const VECTOR3_uint &i, bool *faceState)const;
    // вычисление состояний совпадения каждой из 6-ти поверхностей 6-гранного КЭ на соответствубщей грани 6-гранной области
    // faceState - массив из 6 элементов, i - координаты 0-й вершины 6-гранного КЭ
    void calcFEFaceStates(const VECTOR3_uint &i0, bool *faceState)const;
};



// параметры сетки - ползучесть
struct CreepParameters
{
    double E;               // модуль упругости
    double Nu;              // коэффициент Пуассона
    double sigma0;
    double alpha;
    double timeCreep;
    double timeRelaxation1;
    double timeRelaxation2;
    double A;
    double n;
    double m;

    CUBE cube;              // параллелепипед
    VECTOR3_uint N;         // разбиение сетки
};

// параметры сетки - вдавливание жёсткого цилиндра в упругое полупространство
struct ContactCylinderParameters
{
    SurfaceType surfaceType;                // способ представления поверхности
    RegionDecompositionType decompositionType;    // способ декомпозиции поверхности

    double d;   // перемещение цилиндра по вертикали
    double y0;  // начальное положение по y
    double R0;  // радиус(константа)

    Interpolation::Interpolant1D_base *R_fun;
    Interpolation::Interpolant1D_base *x_fun;
    Interpolation::Interpolant1D_base *y_fun;

    double w_stiffness; // множитель для вычисления контактных жёсткостей
    //double stiffness_min;   // минимальная контактная жёсткость
    //double stiffness_max;   // максимальная контактная жёсткость
    //double h_max;           // ограничение сверху расстояния от узла до поверхности в случае контакта
    double E;              // модуль упругости
    double Nu;             // коэффициент Пуассона

    int STEPS_PUSH;
    int CORR;

    std::vector<POINT3> av;
    std::vector<HexagonArea> ah;

    double L;   // длина цилиндра(для рассчёта аналитического решения)
    double Ly;  // высота полупространства(для рассчёта аналитического решения)
    int Nz; // количество разбиений по z(для визуализации)
    bool halhGrid;  // true - строить половинку сетки
};

// параметры сетки - вдавливание жёсткого шара в упругое полупространство
struct ContactSphereParameters
{
    SurfaceType surfaceType;                // способ представления поверхности
    RegionDecompositionType decompositionType;    // способ декомпозиции поверхности

    bool plasticOut;
    // параметры теста для определения глубина - сила
    double bd_max;
    double bYeldCoeff;
    double n = 0;
    double P_y;
    double d_y;
    double A_y;
     double bP_max_numb; // численное значение безразмерного давления в конце нагружения
     double bA_max_numb_shamanstvo; // численное значение безразмерного радиуса зоны контакта в конце нагружения
     double unload_bd_res_numb = -1; // численное значение глубины при которой происходит отлипание
     // максимальные значения за всю историю нагружения
     double max_sigma_eqv; // максимальное эквивалентное напряжение
     double max_epsElastoPlastic_eqv; // максимальная эквивалентная деформация
     double max_q; // максимальное значение параметра Одквиста



    int NN;
    double c1;
    double c2;
    double c3;
    double grid_mnojitel1;   // =1.5 - магический коэффициент, на который домнажается sqrt() для расчёта размера образца
    double grid_mnojitel2;   // =1.5 - магический коэффициент, на который домнажается sqrt() для расчёта размера образца
    double q;
    double d0;   // перемещение сферы по вертикали
    double z0;  // начальное положение по z
    double R0;  // радиус(константа)


    int STEPS_PUSH;
    int CORR;
    bool sortNodes; // переупорядочивать узлы сетки

    Interpolation::Interpolant1D_base *R_fun;
    Interpolation::Interpolant1D_base *x_fun;
    Interpolation::Interpolant1D_base *y_fun;
    Interpolation::Interpolant1D_base *z_fun;

    double w_stiffness; // множитель для вычисления контактных жёсткостей
    //double stiffness_min;   // минимальная контактная жёсткость
    //double stiffness_max;   // максимальная контактная жёсткость
    //double h_max;       // ограничение сверху расстояния от узла до поверхности в случае контакта

    int plasticityMethodType;   // метод решения задачи пластичности
    int incForsesMode;          // 1 - безумство, 2 - без безумства
    int constantNormal;         // 1 - нормаль вычисляется к исходной точке, 0 - нормаль вычисляется к точке контакта, жёсткость не обновляется

    double w_project;   // коэфициент регуляризации
    double w_midPoint;
    double cosTettaMin;   // минимальный косинус угла для разгрузки заранее

    double E;           // модуль упругости
    double Nu;          // коэффициент Пуассона
    double elasticSigmaLimit;   // предел упругости
    double sigmaResidualLimit;
    double epsResidualLimit;
    double contactDeltaFResidualLimit;
    double contactEndPointResidualLimit;


    std::vector<POINT3> av;
    std::vector<HexagonArea> ah;
};

// параметры сетки формовка
struct FormovkaParameters
{
    POINT3 p1, p2;  // 2 крайние вершины параллелепипеда
    VECTOR3_uint N; // разбиения на КЭ
    int N_p;        // размер области, на которую давят (количество КЭ)
    int N_react;    // размер области, которая закреплена
    SurfaceType surfaceType;                // способ представления поверхности
    RegionDecompositionType decompositionType;    // способ декомпозиции поверхности

    double y0;  // положение по y(константа)
    double R0;  // радиус(константа)
    Interpolation::Interpolant1D_base *R_fun;
    Interpolation::Interpolant1D_base *x_fun;
    Interpolation::Interpolant1D_base *y_fun;

    // кривая ползучести
    double B;
    double n;
};

// параметры сетки толстостенной сферы 3D
struct SphereParameters
{
    //POINT3 p0;      // центр сферы
    double r1;      // внутренний радиус
    double r2;      // внешний радиус
    int Nparts;     // количество частей, на который делится сфера осоэдром
    int N;          // четверть от количества узлов на границе одной части сферы
    double q;       // коэффициент разрядки
    int curvilinear;// 0 - линейное отображение, 1 - квадратичное
    int buildingMethod; // способ построения (0 - делим дуги, 1 - отображение на куб)

    // кривая ползучести
    double B;
    double n;
};

// параметры сетки для задачи Кирша
struct KirschParameters
{
    bool XFEM;      // true для теста XFEM
    double NU;      // коэффициент Пуассона
    double E;       // модуль Юнга
    double Lx;   // размер пластины x
    double Ly;   // размер пластины y
    double Lz;   // размер пластины z
    double r;   // радиус ответстия в центре
    double Px;      // поверхностная сила по x
    double Py;      // поверхностная сила по y
    int N_r;    // разбиение по угловой координате
    int N_fi;   // разбиение по радиальной координате
    int N_z;    // разбиение по толщине пластинки
    double q;   // коэффициент сгущения сетки к отверстию
};

// параметры сетки, неоднородная пластинка с трещиной
struct CrackPlateGenParameters
{
    int grid_mode;      // 1 - XFEM, 0 - FEM
    POINT3 p1, p2;  // 2 крайние вершины параллелепипеда
    POINT2 matetial_surface_y; // координаты границ между материалами
    double NU;      // коэффициент Пуассона
    double E[3];
    bool fix_z;     // true = фиксировать по z везде
    double Px;      // поверхностная сила по x
    double P_in_left;    // давление на левую часть трещины
    double P_in_right;   // давление на левую часть трещины
    VECTOR3_uint N; // разбиения на КЭ
    int fem_multipler;  // во сколько раз в мкэ плотнее сетка
    POINT2 crack_y; // координаты вершин трещины
     double crack_x_local;   // локальная координата трещины в КЭ (-1...+1)
     double crack_x_global;  // глобальная координата трещины в КЭ
     VECTOR3_uint XFEM_N;    // разбиение КЭ для интегрирования (XFEM)
     int XFEM_crack_Nx_left;      // количество элементов подсетки слева от трещины
     int XFEM_crack_Nx_right;     // количество элементов подсетки справа от трещины
    std::vector<int> centr_ind;
    std::vector<std::set<size_t>> subVertex_ind;
    double minSigmaEqv;
    double maxSigmaEqv;
    bool vertex_tip_add;
    bool tip_elastic;
    bool tip_elastoplastic;
    double pl_hardening_n;  // коэффициент упрочнения
};

// сетка 3D
struct Grid3D
{
    // построение регулярной 2D сетки
    //void genRegular2D(const Regular2DParameters &gp);


    // вывод сетки в формате openSCAD
    // построение модели тела в формате OpenSCAD
    void buldOpenSCADModel(const char *file_scad);
    // строит сетку для тестирования ползучести
    void genCreep(const CreepParameters &cp);
    // строит сетку для статической задачи герца: вдавливание жёсткого цилиндра в полупространство
    void genContactCylinder(const ContactCylinderParameters &gp);
    // строит сетку для статической задачи герца: вдавливание жёсткого шара в полупространство
    void genContactSphere(const ContactSphereParameters &gp);
    // строит сетку для формовки
    void genFormovka(const FormovkaParameters &gp);
    // строит толстостенную сферу
    // на x=0, y=0, z=0 первые краевые с индексами 0,1,2 соответственно
    // на внутренней, внешней поверхности вторые краевые с индексами 0,1 соответственно
    void genSphere(const SphereParameters &gp);
    // построение сетки для задачи Кирша
    void genKirch(const KirschParameters &gp);
    // неоднородная пластинка с трещиной (обобщение: базисные ф-и задаются интерполянтами)
    void genCrack_plate_gen_XFEM(CrackPlateGenParameters &gp);
    void genCrack_plate_gen_FEM(CrackPlateGenParameters &gp);


    // построение сетки на основе набора шестигранников с гранями находящимися в одной плоскости
    void build(const std::vector<POINT3> &areaVertex, const std::vector<HexagonArea> &areaHexagon);
    // сортировка вершин, контактные вершины в конце
    void cuthillMcKee(const std::vector<int> &contactFESurfaceIndex);
    // переупорядочивание ущлов методом вложенных сечений
    void nestedDissection(const std::vector<int> &contactFESurfaceIndex);
    // сохранение состояния в файл
    void save(const std::string &subdir) const;

    //void genRectangle(const GridRectangle3D &gp);
    std::vector<FE_base *> fe;                      // индекс материала + функции для обращения к КЭ
    std::vector<POINT3> vertex;                     // вершины (соответствуют базисным функциям)
        std::vector<POINT3> vertexForCurvature ;        // дополнительные вершины для КЭ, в которых порядок отображения выше порядка базисных функций
    std::vector<BoundaryCondition1> bc1;            // индекс вершины + индекс первого краевого условия
    std::vector<FiniteElementSurface> FESurface;    // массив поверхностей: индекс КЭ + индекс грани - будет использоваться для 2-х краевых и для контакта
    std::vector<Interpolation::InterpolantSurface_base *> analiticalSurface;  // массив поверхностей, заданных аналитически
    std::vector<Interpolation::InterpolantSurface_base *> ISurface;    // массив поверхностей, заданных интерполянтами(сплайны)
    RegionDecompositionParameters regionIndexationParameters;      // параметры индексации пространства
    // XFEM
    std::vector<Grid::CrackSurface_base *> cs; // поверхности трещин

    GlobalDOFs *DOFs;                               // нумерация узловых базисных функций

    std::vector<Grid::MechBoundaryCondition2Crack> bc2Crack;       // 2-е краевые для поверхностей трещин в КЭ
};

}   // namespace Grid

#endif  // GRID_H
