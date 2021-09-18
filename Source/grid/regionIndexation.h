#ifndef REGIONINDEXATION_H
#define REGIONINDEXATION_H

#include  <unordered_map>

#include "elementary.h"
#include "grid_base.h"

namespace Grid
{
using namespace Elementary;

// индекс фигуры + сфера, которая содержит данное тело
struct SpheredBody
{
    size_t bodyInd;                 // индекс тела
    Sphere s;                       // сфера, которая содержит данное тело
};

typedef std::vector<size_t> SpheredBodyArray;
// индексация пространства
// область, которой соответствует набор шаров
struct RegionWithSpheres
{
    CUBE q0;                        // куб, внутри которого содержатся центры тел
    RegionWithSpheres *down;                   // ссылка на массив из 8ми подобластей или null (если узел листовой)
    int bodysNomber;                // количество тел в области
    Sphere s0;                      // сфера, внутри которой содержатся тела (не имеет значения если bodysNomber = 0)
    SpheredBody *sb;                // null или массивиз MAXSIZE тел, центры которых находятся внутри куба q0 (если узел листовой)
    // инициализация области
    void init(const CUBE &set_q0);
    // освобождение памяти
    void release();
    // добавление тела (индекс + шар)
    void insert(const SpheredBody &sbEl);
    // дополняет массив bodyInd[] индексами тел, которые могут (не обязательно) пересекаться с шаром s
    void findIntersectingBodies(const Sphere &s, SpheredBodyArray &body)const;
    // дополняет массив bodyInd[] индексами тел, которые могут (не обязательно) пересекаться с отрезком p1-p2
    void findIntersectingBodies(const POINT3 &p1, const POINT3 &p2, SpheredBodyArray &body)const;
};

// индексация куба
// есть множество кубов, для которых известны defined, side, onBorder,   (q0,h,...)
// есть множество вершин кубов, для которых известны R, side, onBorder
// идентификатор вершины ID = (X, Y, Z)
// (X, Y, Z = 0..2^H0 - индексы по координатам при наиболее мелком разбиении глубины h_max)

// куб делится, если
// 1) h < h_min
// 2) h < h_max и состояния вершин куба не одинаковые
// 3) h < h_max и состояния вершин куба одинаковые, но куб не находится внутри какой-то из 8 сфер

// поиск
// если defined=true, то вернуть значение, иначе смотреть ветки
// если в листике defined=false, то вернуть первое приближение

// состояние точки
// считаем состояния различными, не обращая внимания на onBorder
struct RegionIndexation_State
{
    int side;
    bool onBorder;
    friend inline bool operator==(const RegionIndexation_State &l, const RegionIndexation_State &r)
    {
        //return l.side == r.side && l.onBorder == r.onBorder;
        return l.side == r.side;
    }
    friend inline bool operator!=(const RegionIndexation_State &l, const RegionIndexation_State &r)
    {
        //return l.side != r.side || l.onBorder != r.onBorder;
        return l.side != r.side;
    }
};

// информация о вершине куба
struct RegionIndexation_vertexData
{
    double R;                       // расстояние до поверхности
    RegionIndexation_State state;   // положение относительно поверхности
    Grid::SurfacePositionData ssd;  // информация о том, где ближайшая точка находится на интерполянте
};

// индекс вершины куба
// 0..2^H0 - индексы по координатам при наиболее мелком разбиении глубины h_max
//typedef VECTOR3_uint RegionIndexation_vertexIndex;
struct RegionIndexation_vertexIndex:
        public VECTOR3_uint
{
};

// функтор равенства индексов вершин куба
struct RegionIndexationVertexIndex_equal
{
    bool operator()(const RegionIndexation_vertexIndex &vIndex1, const RegionIndexation_vertexIndex &vIndex2) const;
};

// хеш функция от индекса вершины куба
struct RegionIndexation_vertexIndex_hash
{
    size_t operator()(const RegionIndexation_vertexIndex &vIndex) const;
};

// контейнер для хранения информации о вершинах куба с доступом по индексу вершины куба
typedef std::unordered_map<RegionIndexation_vertexIndex, RegionIndexation_vertexData, RegionIndexation_vertexIndex_hash, RegionIndexationVertexIndex_equal> RegionIndexationVertexDataMap;

// индексация области
struct RegionIndexation
{
    RegionIndexation_vertexIndex v0Index;   // индекс вершины куба с наименьшими координатами
    unsigned int h;                         // высота узла дерева (h_min <= h <= h_max), по которому определяетмя размер куба и индексы его 8-ми вершин
    RegionIndexation *down;                 // ссылка на массив из 8-ми подобластей или null (если узел листовой)
    bool defined;
    RegionIndexation_State state;
    // defined = true если внутри куба state одинаковые
    // значение state определяют состояния точек внутри куба, если defined = true
    // если defined = false, то state не имеют значения
    // если defined = true, то state определяют состояния любой точки внутри куба

    // инициализация области
    void init(const RegionDecompositionParameters &param, const RegionIndexation_vertexIndex &set_v0Index, const unsigned int set_h, const Grid::Surface_base *s, const Grid::SurfacePositionData &set_prevSolutionData,
              RegionIndexationVertexDataMap *m);
    // освобождение памяти
    void release();
    // получение информации о заданной точке
    // если p_defined = true, то состояние однозначно определено и точке p соответствуют значения p_side и p_onBorder
    // если p_defined = false, то состояние однозначно не определено и точке p соответствует значение approximateSolutionData
    // approximateSolutionData - первое приближение для подробного поиска ближайшей точки поверхности
    void getInformation(const RegionDecompositionParameters &param, const RegionIndexationVertexDataMap *m, const POINT3 &p,
                        bool &p_defined, int &p_side, bool &p_onBorder, Grid::SurfacePositionData &approximateSolutionData)const;

    // вспомогательные процедуры
    // расчёт координат вершины по индексу
    void getVertex(const RegionDecompositionParameters &param, const RegionIndexation_vertexIndex &vIndex,
                   POINT3 &v)const;
    // расчёт координат куба
    void getCube(const RegionDecompositionParameters &param,
                 CUBE &cube)const;
    // расчёт длины диагонали куба
    double getDiag(const RegionDecompositionParameters &param, const unsigned int h)const;
    // нахождение информации информации о вершине куба (предполагается что в контейнере эта информация отсутствует)
    // поиск готовых данных в вершине vIndex или их расчёт и сохранение
    void initOrGetvertexData(const RegionDecompositionParameters &param, const Grid::Surface_base *s, const RegionIndexation_vertexIndex &vIndex, const SurfacePositionData &prevSolutionData,
                             Grid::RegionIndexationVertexDataMap *m, RegionIndexation_vertexData &vd)const;

    //RegionIndexationVertexDataMap::const_iterator initOrGetCubeVertexInformation(const RegionIndexationParameters &param, const Grid::Surface_base *s, const RegionIndexation_vertexIndex &vIndex, const SurfaceSolutionData &prevSolutionData,
    //                                                                             Grid::RegionIndexationVertexDataMap *m)const;
    //void findIntersectingBodies(const Sphere &s, std::vector<int> &body)const;
    //RegionIndexation *adjacentCell[6];      // массив из 6ти указателей на соседние ячейки (null указывает на край индексируемой области)
};

}   // namespace Grid
#endif // REGIONINDEXATION_H
