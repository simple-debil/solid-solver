/* --------------------------------------------------------- */
// ИНТЕРПОЛИРОВАНИЕ
/* --------------------------------------------------------- */

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "elementary.h"

#include "grid_base.h"
#include "regionIndexation.h"

#define UNKNOWN_FE_INDEX        -1

namespace Interpolation
{
using namespace Elementary;

// рисунок 32 бит x*y
typedef struct
{
    unsigned int *c;    // буффер
    int x;              // ширина
    int y;              // высота
} BMP_INF;

// вспомогательная структура 1
typedef struct
{
    unsigned short  tip;
    unsigned int sizef;
    unsigned short rez1;
    unsigned short rez2;
    unsigned int gotobits;
} BMP_INF1;

// вспомогательная структура 2
typedef struct
{
    unsigned int lenstruct;
    unsigned short x;
    unsigned short y;
    unsigned short planes;
    unsigned short bits;
} BMP_INF2;

void calc_rgb(const int colorMode, const double value, const double valueMin, const double valueNull, const double valueMax, int &r, int &g, int &b);
//возвращает индекс цвета (a, red,green,blue), каждая компонента <=2^8
inline unsigned int RGB32(int a, int r, int g, int b)
{
    return ((b)+((g)<<8)+((r)<<16)+((a)<<24));
}
void bmp_save(const char *file_name, const BMP_INF &inf);
void bmp_load(const char *file_name, BMP_INF &inf);



// линейный лагранжев интерполянт скалярной функции на регулярной равномерной сетке
struct Interpolant1D_Lagrange1:
        public Interpolant1D_base
{
    static double G[2][2];		// локальная матрица жесткости
    virtual void buildInterpolant() override;
    virtual void getNodesCoordinates(std::vector<POINT1> &coordinates)const override;
    virtual void buildInterpolantByAllNodes(const Vector &nodeValue) override;
    virtual double fun(const POINT1 &p)const override;
    virtual double difFun(const POINT1 &p, const DIF_STATE1 &dif)const override;
    // вспомогательные функции
    // вычисление глобального индекса базисной ф-и с локальным индексом localBusFunIndex, задан индекс КЭ feIndex
    int getBasFunGlobaIndex(const int feIndex, const int localBusFunIndex)const;
    void buildLagrange1DLocalMatrix();
};
// линейный лагранжев интерполянт скалярной функции на неравномерной сетке
// gr и alpha не используются
struct Interpolant1D_Lagrange1Unregular:
        public Interpolant1D_base
{
    std::vector<POINT1_VALUE_UNREGULAR> mu;// массив добавленных точек, индекс = номер КЭ (кроме последней точки)
    virtual void addPoint(const POINT1 &p, const double value) override;
    virtual void findFe(const POINT1 &p,
                        int &i)const override;
    virtual void buildInterpolant() override;
    // пустышка
    virtual void getNodesCoordinates(std::vector<POINT1> &coordinates)const override;
    // пустышка
    virtual void buildInterpolantByAllNodes(const Vector &nodeValue) override;
    virtual double fun(const POINT1 &p)const override;
    virtual double difFun(const POINT1 &p, const DIF_STATE1 &dif)const override;
};

// кубический эрмитов интерполянт скалярной функции на регулярной равномерной сетке
struct Interpolant1D_Hermite3:
        public Interpolant1D_base
{
    static double G[4][4];		// локальная матрица жесткости
    virtual void buildInterpolant() override;
    virtual void getNodesCoordinates(std::vector<POINT1> &coordinates)const override;
    virtual void buildInterpolantByAllNodes(const Vector &nodeValue) override;
    virtual double fun(const POINT1 &p)const override;
    virtual double difFun(const POINT1 &p, const DIF_STATE1 &dif)const override;
    // вспомогательные функции
    // вычисление глобального индекса базисной ф-и с локальным индексом localBusFunIndex, задан индекс КЭ feIndex
    int getBasFunGlobaIndex(const int feIndex, const int localBusFunIndex)const;
    void buildHermite1DLocalMatrix();
};
// кубический лагранжев интерполянт скалярной функции на регулярной равномерной сетке
struct Interpolant1D_Lagrange3:
        public Interpolant1D_base
{
    static double G[4][4];		// локальная матрица жесткости
    virtual void buildInterpolant() override;
    virtual void getNodesCoordinates(std::vector<POINT1> &coordinates)const override;
    virtual void buildInterpolantByAllNodes(const Vector &nodeValue) override;
    virtual double fun(const POINT1 &p)const override;
    virtual double difFun(const POINT1 &p, const DIF_STATE1 &dif)const override;
    // вспомогательные функции
    // вычисление глобального индекса базисной ф-и с локальным индексом localBusFunIndex, задан индекс КЭ feIndex
    int getBasFunGlobaIndex(const int feIndex, const int localBusFunIndex)const;
    void buildLagrange1DLocalMatrix();
};

// бикубический эрмитов интерполянт скалярной функции на регулярной равномерной сетке
struct Interpolant2D_Hermite3:
        public Interpolant2D_base
{
    static double G[16][16];		// локальная матрица жесткости
    virtual void buildInterpolant() override;
    virtual void getNodesCoordinates(std::vector<POINT2> &coordinates)const override;
    virtual void buildInterpolantByAllNodes(const Vector &nodeValue) override;
    virtual double fun(const POINT2 &p)const override;
    virtual double difFun(const POINT2 &p, const DIF_STATE2 &dif)const override;
    // вспомогательные функции
    // вычисление глобального индекса базисной ф-и с локальным индексом localBusFunIndex, задан индекс КЭ feIndex
    int getBasFunGlobaIndex(const int feIndex, const int localBusFunIndex)const;
    void buildHermite2DLocalMatrix();
};
// бикубический лагранжев интерполянт скалярной функции на регулярной равномерной сетке
struct Interpolant2D_Lagrange3:
        public Interpolant2D_base
{
    static double G[16][16];		// локальная матрица жесткости
    virtual void buildInterpolant() override;
    virtual void getNodesCoordinates(std::vector<POINT2> &coordinates)const override;
    virtual void buildInterpolantByAllNodes(const Vector &nodeValue) override;
    virtual double fun(const POINT2 &p)const override;
    virtual double difFun(const POINT2 &p, const DIF_STATE2 &dif)const override;
    // вспомогательные функции
    // вычисление глобального индекса базисной ф-и с локальным индексом localBusFunIndex, задан индекс КЭ feIndex
    int getBasFunGlobaIndex(const int feIndex, const int localBusFunIndex)const;
    void buildLagrange2DLocalMatrix();
};

// отсутствие декомпозиции пространства
struct InterpolantSurfaceDecomposition_none:
        public InterpolantSurfaceDecomposition_base
{
    virtual void buildRegions(const Grid::RegionDecompositionParameters *, const InterpolantSurface_base *) override;
    virtual void excludeIntersectionCase(const POINT3 &, const POINT3 &, const Grid::SurfacePositionData &, const Grid::SurfacePositionData &,
                                         Grid::SurfacePositionData &, Grid::SurfacePositionData &, bool &mayBeFound)const override;
};
// декомпозиция пространства на вложенные кубы, некоторые из которых имеют определённое положение относительно поверхности
struct InterpolantSurfaceDecomposition_cubeVertexIndexation:
        public InterpolantSurfaceDecomposition_base
{
    Grid::RegionIndexationVertexDataMap map;
    Grid::RegionIndexation ri;             // корневая область
    virtual void buildRegions(const Grid::RegionDecompositionParameters *set_param, const InterpolantSurface_base *surface) override;
    virtual void excludeIntersectionCase(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &prevSolutionData1, const Grid::SurfacePositionData &prevSolutionData2,
                                         Grid::SurfacePositionData &solutionData1, Grid::SurfacePositionData &solutionData2, bool &mayBeFound)const override;
};
// представление поверхности множеством сфер
struct InterpolantSurfaceDecomposition_surfaceAsSpheres:
        public InterpolantSurfaceDecomposition_base
{
    Grid::RegionWithSpheres ri;             // корневая область
    virtual void buildRegions(const Grid::RegionDecompositionParameters *set_param, const InterpolantSurface_base *surface) override;
    virtual void excludeIntersectionCase(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &, const Grid::SurfacePositionData &,
                                         Grid::SurfacePositionData &, Grid::SurfacePositionData &, bool &mayBeFound)const override;
};

// цилиндрическая поверхность
struct AnaliticalSurface_Cylinder:
        public InterpolantSurface_base
{
    Interpolant1D_base *R_fun;
    Interpolant1D_base *x_fun;
    Interpolant1D_base *y_fun;
    POINT3 C;
    double R;
    AnaliticalSurface_Cylinder(const Grid::RegionDecompositionType set_decompositionType, const Grid::GridRectangleRegular2D &set_gr, Interpolant1D_base *set_R, Interpolant1D_base *set_x, Interpolant1D_base *set_y, const int set_maxIter, const double set_eps_uv_2d, const double set_eps_uv_1d, const double set_eps_grad, const double set_eps_inerception);
    virtual void findNearestPoint(const POINT3 &p, const VECTOR3 &, const Grid::SurfacePositionData &,
                                  Grid::SurfacePositionData &, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const override;
    virtual void findIntersection(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &prevSolutionData1, const Grid::SurfacePositionData &prevSolutionData2,
                                  Grid::SurfacePositionData &solutionData1, Grid::SurfacePositionData &solutionData2, POINT3 &intersectionPoint, POINT3 &normal, int &side, bool &found)const override;
    virtual void update(const double time, const Grid::RegionDecompositionParameters *set_param) override;
    virtual POINT3 s(const POINT2 &) const override;
    virtual POINT3 difs(const POINT2 &, const DIF_STATE2 &) const override;
    void findIntersection_analitical(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &prevSolutionData1, const Grid::SurfacePositionData &prevSolutionData2,
                                     Grid::SurfacePositionData &solutionData1, Grid::SurfacePositionData &solutionData2, POINT3 &intersectionPoint, POINT3 &normal, int &side, bool &found)const;
};
// сферическая поверхность
struct AnaliticalSurface_Sphere:
        public InterpolantSurface_base
{
    Interpolant1D_base *R_fun;
    Interpolant1D_base *x_fun;
    Interpolant1D_base *y_fun;
    Interpolant1D_base *z_fun;
    POINT3 C;
    double R;
    AnaliticalSurface_Sphere(const Grid::RegionDecompositionType set_decompositionType, const Grid::GridRectangleRegular2D &set_gr, Interpolant1D_base *set_R, Interpolant1D_base *set_x, Interpolant1D_base *set_y, Interpolant1D_base *set_z, const int set_maxIter, const double set_eps_uv_2d, const double set_eps_uv_1d, const double set_eps_grad, const double set_eps_inerception);
    virtual void findNearestPoint(const POINT3 &p, const VECTOR3 &, const Grid::SurfacePositionData &,
                                  Grid::SurfacePositionData &, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const override;
    virtual void findIntersection(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &prevSolutionData1, const Grid::SurfacePositionData &prevSolutionData2,
                                  Grid::SurfacePositionData &solutionData1, Grid::SurfacePositionData &solutionData2, POINT3 &intersectionPoint, POINT3 &normal, int &side, bool &found)const override;
    virtual void update(const double time, const Grid::RegionDecompositionParameters *set_param) override;
    virtual POINT3 s(const POINT2 &) const override;
    virtual POINT3 difs(const POINT2 &, const DIF_STATE2 &) const override;
};
// интерполянт поверхности эрмита
struct InterpolantSurface_Hermite3:
        public InterpolantSurface_base
{
    InterpolantSurface_Hermite3(const Grid::RegionDecompositionType set_decompositionType, const Interpolant2D_Hermite3 *set_it_x, const Interpolant2D_Hermite3 *set_it_y, const Interpolant2D_Hermite3 *set_it_z, const int set_maxIter, const double set_eps_uv_2d, const double set_eps_uv_1d, const double set_eps_grad, const double set_eps_inerception);
    virtual void findNearestPoint(const POINT3 &p, const VECTOR3 &u, const Grid::SurfacePositionData &prevSolutionData,
                                  Grid::SurfacePositionData &solutionData, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const override;
    virtual POINT3 s(const POINT2 &uv) const override;
    virtual POINT3 difs(const POINT2 &uv, const DIF_STATE2 &dif) const override;
};
// интерполянт поверхности лагранжа
struct InterpolantSurface_Lagrange3:
        public InterpolantSurface_base
{
    InterpolantSurface_Lagrange3(const Grid::RegionDecompositionType set_decompositionType, const Interpolant2D_Lagrange3 *set_it_x, const Interpolant2D_Lagrange3 *set_it_y, const Interpolant2D_Lagrange3 *set_it_z, const int set_maxIter, const double set_eps_uv_2d, const double set_eps_uv_1d, const double set_eps_grad, const double set_eps_inerception);
    virtual void findNearestPoint(const POINT3 &p, const VECTOR3 &u, const Grid::SurfacePositionData &prevSolutionData,
                                  Grid::SurfacePositionData &solutionData, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const override;
    virtual POINT3 s(const POINT2 &uv) const override;
    virtual POINT3 difs(const POINT2 &uv, const DIF_STATE2 &dif) const override;
};
}   // namespace Interpolation

#endif  // INTERPOLATION_H
