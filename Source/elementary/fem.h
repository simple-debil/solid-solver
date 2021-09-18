/* --------------------------------------------------------- */
// МЕТОД КОНЕЧНЫХ ЭЛЕМЕНТОВ
/* --------------------------------------------------------- */

#ifndef FEM_H
#define FEM_H

#include "elementary.h"
#include "integration.h"

// типы данных и функции для мкэ
namespace Fem
{
using namespace Elementary;

//  1-мерные стандартные базисные функции и их производные на шаблонном отрезке
//  cubePoint - координата точки на шаблонном отрезке [-1, 1]
//  i - номер базисной функции
//  dif - количество  производных
// стандартная скалярная лагранжева функция 1-го порядка, i=0|1
double cubeLagrange1_1D(const POINT1_CUBE cubePoint, const int i);
double cubeLagrange1_1D(const POINT1_CUBE cubePoint, const int i, const int dif);
// стандартная скалярная лагранжева функция 2-го порядка, i=0|1|2
double cubeLagrange2_1D(const POINT1_CUBE cubePoint, const int i);
double cubeLagrange2_1D(const POINT1_CUBE cubePoint, const int i, const int dif);

//  3-мерные стандартные базисные функции и их производные на шаблонном отрезке
//  cubePoint - координаты точки на шаблонном кубе [-1, 1]^3
//  i - номер базисной функции
//  dif - количество  производных по каждой координате
// стандартная скалярная лагранжева функция 1-го порядка, 0<=i<8
double cubeLagrange1_3D(const POINT3_CUBE &cubePoint, const int i, const DIF_STATE3 &dif);
// стандартная скалярная лагранжева функция 1-го порядка, 0<=i<27
double cubeLagrange2_3D(const POINT3_CUBE &cubePoint, const int i, const DIF_STATE3 &dif);

//  расчёт координат точки шестигранника по координатам шаблонного куба
//  вычисляется как линейная комбинация стандартных базисных функций, поэтому можно брать производные
//  v - массив вершин шестигранника
//  cubePoint = {X,Y,Z} - координаты шаблонного куба [-1, 1]^3
//  dif - количество производных по каждой координате
//  point = {x,y,z} - координаты искомой точки шестигранника
// линейное отображение шаблонного куба в шестигранник, линейные базисные функции
// 8 вершин шестигранника (в v)
void cubeToLagrange1Hexagon(const POINT3 *v, const POINT3_CUBE &cubePoint, const DIF_STATE3 &dif,
                            POINT3 &point);
// квадратичное отображение шаблонного куба в шестигранник, линейные базисные функции
// 27 вершин шестигранника (в v)
void cubeToLagrange2Hexagon(const POINT3 *v, const POINT3_CUBE &cubePoint, const DIF_STATE3 &dif,
                            POINT3 &point);

//  расчёт матрицы отображения шаблонного куба в шестигранник и детерминанта
//  v - массив вершин шестигранника
//  cubePoint - координаты шаблонного куба
//  J - матрица преобразования координат, переводящего шаблонный куб в шестигранник
//  detJ - детерминант J
// отображение линейное, 8 вершин шестигранника (в v)
void lagrange1Hexagon_mappingJ(const POINT3 *v, const POINT3_CUBE &cubePoint,
                               MATR3x3 &J, double &detJ);
// отображение квадратичное, 27 вершин шестигранника (в v)
void lagrange2Hexagon_mappingJ(const POINT3 *v, const POINT3_CUBE &cubePoint,
                               MATR3x3 &J, double &detJ);
//  расчёт детерминанта отображения шаблонного куба в шестигранник
//  v - массив вершин шестигранника
//  cubePoint - координаты шаблонного куба
//  возврашает детерминант J
// отображение линейное, 8 вершин шестигранника (в v)
double lagrange1Hexagon_mappingDetJ(const POINT3 *v, const POINT3_CUBE &cubePoint);
// отображение квадратичное, 27 вершин шестигранника (в v)
double lagrange2Hexagon_mappingDetJ(const POINT3 *v, const POINT3_CUBE &cubePoint);

//  производные базисной функции 1-го порядка в точке (координаты шаблонного куба)
//  v - массив вершин шестигранника
//  cubePoint - координаты шаблонного куба
//  i - номер базисной функции
//  difBasisFunction = {dFi/dx, dFi/dy, dFi/dz} - производные по обычным координатам
//  возвращает число det(J), J - матрица отображения
// отображение линейное, 8 вершин шестигранника (в v)
double lagrange1Hexagon_difLagrange1(const POINT3 *v, const POINT3_CUBE &cubePoint, const int i,
                                     VECTOR3 &difBasisFunction);
// отображение квадратичное, 27 вершин шестигранника (в v)
double lagrange2Hexagon_difLagrange1(const POINT3 *v, const POINT3_CUBE &cubePoint, const int i,
                                     VECTOR3 &difBasisFunction);

//  1-мерные скалярные базисные функции и их производные на отрезке
//  p0 - координата начала отрезка
//  h - длина отрезка
//  p - координата точки на отрезке
//  i - номер базисной функции
//  dif - количество  производных
// скалярная лагранжева базисная функция 1-го порядка, i=0|1
double lagrange1_1D(const POINT1 &p0, const double h, const POINT1 p, const int i, const DIF_STATE1 &dif);
// скалярная лагранжева базисная функция 3-го порядка, i=0|1|2|3
double lagrange3_1D(const POINT1 &p0, const double h, const POINT1 p, const int i, const DIF_STATE1 &dif);
// скалярная эрмитова базисная функция 3-го порядка, i=0|1|2|3
double hermite_1D(const POINT1 &p0, const double h, const POINT1 p, const int i, const DIF_STATE1 &dif);

//  2-мерные скалярные базисные функции и их производные на прямоугольнике
//  p0 - координаты начала прямоугольника
//  hx, hy - стороны прямоугольника (p0.x, p0.y) - (p0.x + hx, p0.y + hy)
//  p - координаты точки на прямоугольнике
//  i - номер базисной функции
//  dif - количество  производных по каждой координате
// скалярная лагранжева базисная функция 3-го порядка (i = 0..15)
double lagrange3_2D(const POINT2 &p0, const double hx, const double hy, const POINT2 &p, const int i, const DIF_STATE2 &dif);
// скалярная эрмитова базисная функция 3-го порядка (i = 0..15)
double hermite_2D(const POINT2 &p0, const double hx, const double hy, const POINT2 &p, const int i, const DIF_STATE2 &dif);









// заполнение таблиц точек Гаусса, предварительных коэффициентов для численного интегрирования,
// значений базисных функций и их производных на шаблонном кубе
void calcCubeLinearBasisFuncValues(const Integration::IntegrationType integrationType,
                                    double *integrationw, double *linearBasCube, VECTOR3 *dLinearBasCube);
void calcCubeLinearBasisFuncValues_ttt(const Integration::IntegrationType integrationType,
                          std::vector<double> &integrationw, std::vector<double> &linearBasCube, std::vector<VECTOR3> &dLinearBasCube);


void calcCubeQuadraticBasisFuncValues(const Integration::IntegrationType integrationType,
                          VECTOR3 *dQuadraticBasCube);
void calcCubeQuadraticBasisFuncValues_ttt(const Integration::IntegrationType integrationType,
                          std::vector<VECTOR3> &dQuadraticBasCube);


void calcLinearHexagonLinearBasisFuncValues(const POINT3 *v, const int numPoints, const double *integrationw, const VECTOR3 *dLinearBasCube,
                          double *w, double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3]);

void calcQuadraticHexagonLinearBasisFuncValues(const POINT3 *v, const int numPoints, const double *integrationw, const VECTOR3 *dLinearBasCube, const VECTOR3 *dQuadraticBasCube,
                          double *w, double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3]);




void calcIntForG(const int numPoints, const double *w, const double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3],
                          double (&intForG_mnjl)[8][8][3][3]);
void calcIntForM_ttt(const int numPoints, const double *w, const std::vector<double> &basCube,
                          double (&intForM_mn)[8][8]);
void calcIntForb_f(const int numPoints, const double *w, const double *basCube,
                          double (&intForb1_m)[8]);
void calcIntForb_df(const int numPoints, const double *w, const double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3],
                          double (&intForb2_mj)[8][3]);
//double linear1D(const double a0, const double a1, const double a, const int i);
//double linear1D(const double a0, const double a1, const double a, const int i, const int dif);
//double rectangleScalarLinear3D(const POINT3 &p0, const POINT3 &p, const int i, const double hx, const double hy, const double hz);
//void solveTriangleMatrix(const double x1, const double x2, const double x3, const double y1, const double y2, const double y3, MATR3x3 &inverse_D, double &det);
//double triangleScalarLinear(const MATR3x3 &inverse_D, const double x, const double y, const int i, const int difx, const int dify);
//double prism3ScalarLinear(const MATR3x3 &inverse_D, const double z1, const double z2, const POINT3 &p, const int i, const DIF_STATE3 &dif);

}   // namespace Fem

#endif  // FEM_H
