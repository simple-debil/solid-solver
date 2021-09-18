#ifndef GRID_BASE_H
#define GRID_BASE_H

#include "elementary.h"

namespace Grid
{
using namespace Elementary;
struct Grid3D;
struct Surface_base;

// регулярная сетка 1D
struct GridRectangleRegular1D
{
    INTERVAL rect;  // промежуток
    int N;          // количество разбиений
    void init(const double x1, const double x2, const int N)
    {
        rect[0] = x1;
        rect[1] = x2;
        this->N = N;
    }
};

// регулярная прямоугольная сетка 2D
struct GridRectangleRegular2D
{
    INTERVAL rect[2];    // прямоугольная область (интервалы по x, y)
    int N[2];       // количество разбиений по координатам x, y
    void init(const double x1, const double x2, const double y1, const double y2, const int Nx, const int Ny)
    {
        rect[0][0] = x1;
        rect[0][1] = x2;
        rect[1][0] = y1;
        rect[1][1] = y2;
        N[0] = Nx;
        N[1] = Ny;
    }
};

// типы поверхностей
enum class SurfaceType
{
    none = -1,
    AnaliticalSurface_Cylinder = 0,
    InterpolantSurface_Hermite3 = 1,
    InterpolantSurface_Lagrange3 = 2,
    FiniteElementSurface = 3,
    AnaliticalSurface_Sphere = 4,
};

// данные о положении точки на поверхности
// если поверхность задана интерполянтом, то index = -1 если данные отсутствуют, index = 0 если данные присутствуют
// если поверхность задана КЭ сеткой, то index = -1 если данные отсутствуют, index = индекс КЭ >= 0 если данные присутствуют
struct SurfacePositionData
{
    POINT2 uv0;
    int index;
};

// типы декомпозиции
enum class RegionDecompositionType
{
    none = -1,
    cubeVertexIndexation = 0,
    surfaceAsSpheres = 1,
};

// общие параметры индексации
struct RegionDecompositionParameters
{
    CUBE q0;                     // корневой куб
    unsigned int h_min;          // минимальная высота дерева
    unsigned int h_max;          // максимальная высота дерева
    unsigned int div;            // сколько раз дробить каждый КЭ при составлении массива сфер
};

// поверхность
// базовый класс
struct Surface_base
{
    SurfaceType surfaceType;                // тип поверхности
    RegionDecompositionType decompositionType;    // тип декомпозиции пространства
    // процедура, вызываемая перед шагом по времени (изменение поверхности в зависимости от времени и/или обновление декомпозиции пространства)
    virtual void update(const double time, const RegionDecompositionParameters *set_param) = 0;
    // p - заданная точка
    // nearestPoint - ближайшая к точке p точка поверхности
    // normal - нормаль наружу в точке nearestPoint (|| p-nearestPoint)
    // ###h - высота заданной точки над поверхностью(|p-nearestPoint|), положительная если точка находится с внешней стороны поверхности
    virtual void findNearestPoint(const POINT3 &p, const VECTOR3 &u, const SurfacePositionData &prevSolutionData,
                                  SurfacePositionData &solutionData, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const = 0;
    // пересечение отрезка p1-p2 с поверхностью цилиндра
    // nearestPoint - точка пересечения
    // normal - нормаль к поверхности в точке nearestPoint
    // direction = +1, если p2 находится с внешней стороны поверхности
    // direction = -1, если p2 находится с внутренней стороны поверхности
    // возвращает true, если пересечение есть
    virtual void findIntersection(const POINT3 &p1, const POINT3 &p2, const SurfacePositionData &prevSolutionData1, const SurfacePositionData &prevSolutionData2,
                                  SurfacePositionData &solutionData1, SurfacePositionData &solutionData2, POINT3 &intersectionPoint, POINT3 &normal, int &side, bool &found)const = 0;
    /*
    1) findNearestPoint
    входные данные:
    точка p и вектор перемещения u
    выходные данные:
    nearestPoint - ближайшая к точке p точка поверхности
    normal - внешняя нормаль к поверхности в точке nearestPoint*
    side = +1, если вектор p - nearestPoint нулевой или направлен в сторону нормали (т.е. p находится с внешней стороны поверхности)
    side = -1, если вектор p - nearestPoint направлен против нормали (т.е. p находится с внутренней стороны поверхности)
    onBorder = true, если nearestPoint находится на границе поверхности

    *Если nearestPoint находится на изломе, то среди смежных участков поверхности выбирается нормаль, которая направлена против перемещения u

    2) findIntersection
    входные данные:
    точки p1 и p2
    выходные данные:
    intersectionPoint - первая (ближайшая к p1) точка пересечения отрезка p1_p2 с поверхностью
    normal - внешняя нормаль к поверхности в точке intersectionPoint
    side = +1, если вектор p2 - intersectionPoint нулевой или направлен в сторону нормали (т.е. p1 находится с внутренней стороны поверхности)
    side = -1, если вектор p2 - intersectionPoint направлен против нормали (т.е. p1 находится с внешней стороны поверхности)
    found = true если пересечение найдено
    */
};
}   // namespace Grid


namespace Interpolation
{
using namespace Elementary;
// 1D точка, значение ф-и в точке и индекс КЭ
struct POINT1_VALUE
{
    POINT1 p;       // точка
    double value;   // значение скаляра в этой точке
    int feIndex;      // идентификатор КЭ, содержащего точку p
};

// 1D точка, значение ф-и в точке
struct POINT1_VALUE_UNREGULAR
{
    POINT1 p;       // точка
    double value;   // значение скаляра в этой точке
};

// 2D точка, значение ф-и в точке и индекс КЭ
struct POINT2_VALUE
{
    POINT2 p;       // точка
    double value;   // значение скаляра в этой точке
    int feIndex;      // идентификатор КЭ, содержащего точку p
};

// одномерный интерполянт скалярной функции на регулярной равномерной сетке
// базовый класс
struct Interpolant1D_base
{
    Grid::GridRectangleRegular1D gr;   // сетка
    double alpha;            // параметр, влияющий на точность интерполяции. Чем меньше, тем более точный, но тем больше точек необходимо для построения
    double l;                // длина отрезка
    double h;                // размер ячеек
    std::vector<POINT1_VALUE> m;// массив добавленных точек
    Vector res;                 // решение СЛАУ (после вызова buildInterpolant)
    double valueMin;
    double valueMax;
    double valueNull;
    // инициализация регулярной равномерной сеткой и параметром регуляризации
    virtual void init(const Grid::GridRectangleRegular1D &set_grid, const double set_alpha);
    // освобождение ресурсов
    virtual void release();
    // добавление точки и значения функции в этой точке
    virtual void addPoint(const POINT1 &p, const double value);
    // добавленное множество точек и значений превращается в интерполянт
    virtual void buildInterpolant() = 0;
    // getNodesCoordinates возвращает массив координат узлов coordinates
    // если поместить в массив nodeValue соответствующие значения интерполируемой ф-и в этих узлах,
    // то buildInterpolantByAllNodes построит интерполянт
    virtual void getNodesCoordinates(std::vector<POINT1> &coordinates)const = 0;
    virtual void buildInterpolantByAllNodes(const Vector &nodeValue) = 0;
    // вычисление значения ф-и в точке (должен быть построен интерполянт ф-ми buildInterpolant или buildInterpolantByAllNodes)
    virtual double fun(const POINT1 &p)const = 0;
    // вычисление производной ф-и в точке по параметрам-локальным координатам (должен быть построен интерполянт ф-ми buildInterpolant или buildInterpolantByAllNodes)
    virtual double difFun(const POINT1 &p, const DIF_STATE1 &dif)const = 0;
    // сохранение изображения интерполянта в bmp-файл
    //virtual void save(const char *file_name, const int bmp_x, const int bmp_y, const int colorMode, const int pictureMode);
protected:
    // вспомогательные функции
    // нахождение координат КЭ, которому принадлежит точка
    // i - номер КЭ
    virtual void findFe(const POINT1 &p,
                        int &i)const;
};

// двумерный интерполянт скалярной функции на регулярной равномерной сетке
// базовый класс
struct Interpolant2D_base
{
    Grid::GridRectangleRegular2D gr;   // сетка
    double alpha;               // параметр, влияющий на точность интерполяции. Чем меньше, тем более точный, но тем больше точек необходимо для построения
    double l[2];                // размерности прямоугольной регулярной сетки
    double h[2];                // размерности ячеек прямоугольной регулярной сетки
    std::vector<POINT2_VALUE> m;// массив добавленных точек
    Vector res;                 // решение СЛАУ (после вызова buildInterpolant)
    double valueMin;
    double valueMax;
    double valueNull;
    // инициализация регулярной равномерной сеткой и параметром регуляризации
    virtual void init(const Grid::GridRectangleRegular2D &set_grid, const double set_alpha);
    // освобождение ресурсов
    virtual void release();
    // добавление точки и значения функции в этой точке
    virtual void addPoint(const POINT2 &p, const double value);
    // добавленное множество точек и значений превращается в интерполянт
    virtual void buildInterpolant() = 0;
    // getNodesCoordinates возвращает массив координат узлов coordinates
    // если поместить в массив nodeValue соответствующие значения интерполируемой ф-и в этих узлах,
    // то buildInterpolantByAllNodes построит интерполянт
    virtual void getNodesCoordinates(std::vector<POINT2> &coordinates)const = 0;
    virtual void buildInterpolantByAllNodes(const Vector &nodeValue) = 0;
    // вычисление значения ф-и в точке (должен быть построен интерполянт ф-ми buildInterpolant или buildInterpolantByAllNodes)
    virtual double fun(const POINT2 &p)const = 0;
    // вычисление производной ф-и в точке по параметрам-локальным координатам (должен быть построен интерполянт ф-ми buildInterpolant или buildInterpolantByAllNodes)
    virtual double difFun(const POINT2 &p, const DIF_STATE2 &dif)const = 0;
    // сохранение изображения интерполянта в bmp-файл
    virtual void save(const char *file_name, const int bmp_x, const int bmp_y, const int colorMode, const int pictureMode);
    // вспомогательные функции
    // нахождение координат КЭ, которому принадлежит точка
    // i - номер КЭ по y, j - номер КЭ по горизонтали
    virtual void findFe_ij(const POINT2 &p,
                          int &i, int &j)const;
    // индекс КЭ по координатам
    virtual int getFeIndex(const int i, const int j)const;
};

struct InterpolantSurface_base;

// декомпозиция пространства для быстрого определения отсутствия пересечения отрезка с поверхностью
// базовый класс
struct InterpolantSurfaceDecomposition_base
{
    static InterpolantSurfaceDecomposition_base *newEl(const Grid::RegionDecompositionType elType);
    Grid::RegionDecompositionParameters param;     // параметры индексации
    // построение декомпозиции пространства
    virtual void buildRegions(const Grid::RegionDecompositionParameters *set_param, const InterpolantSurface_base *surface) = 0;
    // быстрая проверка на отсутствие пересечения
    // mayBeFound = false если пересечения точно нет
    // mayBeFound = true если пересечение не исключено
    virtual void excludeIntersectionCase(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &prevSolutionData1, const Grid::SurfacePositionData &prevSolutionData2,
                                         Grid::SurfacePositionData &solutionData1, Grid::SurfacePositionData &solutionData2, bool &mayBeFound)const = 0;
};

// интерполянт поверхности
struct InterpolantSurface_base:
        public Grid::Surface_base
{
    InterpolantSurfaceDecomposition_base *decomposition;
    bool regionIndexationInited;
    Grid::GridRectangleRegular2D gr;    // равномерная регулярная сетка, общая для 3-х сплайнов
    int maxIter;
    double eps_uv_2d;   // абсолютная точность по изменению сплайновых координат точки поверхности
    double eps_uv_1d;   // относительная точность по изменению параметра при одномерном поиске
    double eps_grad;    // абсолютная точность по модулю градиента
    double eps_inerception; // абсолютная точность по изменению положения точки при поиске пересечения
    double l[2];                // размерности прямоугольной регулярной сетки
    double h[2];                // размерности ячеек прямоугольной регулярной сетки
    Vector interpolantRes;      // коэффициенты интерполянта
    InterpolantSurface_base(const Grid::RegionDecompositionType set_decompositionType, const Grid::GridRectangleRegular2D &set_gr, const int set_maxIter, const double set_eps_uv_2d, const double set_eps_uv_1d, const double set_eps_grad, const double set_eps_inerception);
    virtual void findIntersection(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &prevSolutionData1, const Grid::SurfacePositionData &prevSolutionData2,
                                  Grid::SurfacePositionData &solutionData1, Grid::SurfacePositionData &solutionData2, POINT3 &intersectionPoint, POINT3 &normal, int &side, bool &found)const override;
    virtual void update(const double, const Grid::RegionDecompositionParameters *) override;
    // точка поверхности в зависимости от параметров u и v
    virtual POINT3 s(const POINT2 &uv) const = 0;
    // производная точки поверхности по параметрам u и v
    virtual POINT3 difs(const POINT2 &uv, const DIF_STATE2 &dif) const = 0;
    // поиск КЭ (i,j), которому принадлежит точка uv и узла этого КЭ с наименьшими координатами
    inline void findFe(const POINT2 &uv,
                       int &i, int &j, POINT2 &uv0)const;
    // одномерный поиск ближайшей точки
    double findNearestPoint_find1d_zolot_sechenie (const POINT3 &p, const POINT2 &uv0, const int coordinateIndex, double a, double b, double eps)const;
    double findNearestPoint_find1d_zolot_sechenie_vector (const POINT3 &p, const POINT2 &uv0, const VECTOR2 &duv, double a, double b, double eps)const;
    // проверка находится ли точка на границе поверхности
    bool pointOnBorder(const POINT2 &uv, const double eps_uv_2d)const;
    // расчёт интервара a <= t <= b, чтобы точки uv0 + duv_1*t не выходили за пределы сетки gr
    void calcInterval(const POINT2 &uv0, const POINT2 duv_1,
                       double &a, double &b)const;
    // восстановление решения
    void loadSolutionData(const Grid::SurfacePositionData &solutionData,
                          POINT2 &uv0)const;
    // сохранение данных о решении
    void saveSolutionData(const POINT2 &uv0,
                          Grid::SurfacePositionData &solutionData)const;
};

}   // namespace Interpolation

#endif // GRID_BASE_H
