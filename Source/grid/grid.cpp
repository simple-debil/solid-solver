#define _CRT_SECURE_NO_WARNINGS

#include <map>
#include <unordered_set>
#include <vector>
#include <deque>
#include "stdio.h"

#include "timeinterval.h"
#include "slausolving.h"
#include "fem.h"
#include "grid.h"
#include "regionIndexation.h"

// округляет до d знаячащих цифр
double roundToDigits(const double x, const int d)
{
    if(x == 0.0)
        return 0;
    double scale = pow(10.0, -ceil(log10(fabs(x))) + (double)d);
    return round(x * scale) / scale;
}

// индексы вершин граней (HexagonArea привязан к такой же нумерации)
const int HexagonFace[6][4] =
{   // если переставить 3-й и 4-й индексы то направление (по часовой или против часовой) смотрится буравчиком изнутри 6-гранника
    {0, 1, 2, 3},   // (0) Z = -1   сначала меняется X потом меняется Y - против часовой
    {4, 5, 6, 7},   // (1) Z = +1   сначала меняется X потом меняется Y - по часовой
    {0, 2, 4, 6},   // (2) X = -1   сначала меняется Y потом меняется Z - против часовой
    {1, 3, 5, 7},   // (3) X = +1   сначала меняется Y потом меняется Z - по часовой
    {0, 4, 1, 5},   // (4) Y = -1   сначала меняется Z потом меняется X - против часовой
    {2, 6, 3, 7},   // (5) Y = +1   сначала меняется Z потом меняется X - по часовой
};
const int HexagonFaceCoordinates[6][4] =
{
    {0, 1, 2, -1},  // 0: переменные X, Y, константа Z = -1
    {0, 1, 2, +1},  // 1: переменные X, Y, константа Z = +1
    {1, 2, 0, -1},  // 2: переменные Y, Z, константа X = -1
    {1, 2, 0, +1},  // 3: переменные Y, Z, константа X = +1
    {2, 0, 1, -1},  // 4: переменные Z, X, константа Y = -1
    {2, 0, 1, +1},  // 5: переменные Z, X, константа Y = +1
};

int edge4[6][4][2] =
{
    {{ 0,1 },{ 0,2 },{ 1,3 },{ 2,3 }}, //-3 Z = -1
    {{ 4,5 },{ 4,6 },{ 5,7 },{ 6,7 }}, //+3 Z = +1
    {{ 0,2 },{ 0,4 },{ 2,6 },{ 4,6 }}, //-1 X = -1
    {{ 1,3 },{ 1,5 },{ 3,7 },{ 5,7 }}, //+1 X = +1
    {{ 0,1 },{ 0,4 },{ 1,5 },{ 4,5 }}, //-2 Y = -1
    {{ 2,3 },{ 2,6 },{ 3,7 },{ 6,7 }}, //+2 Y = +1
};

int edgeHexogen[12][2] =
{
    { 0,1 },{ 2,3 },{ 4,5 },{ 6,7 },// ||x
    { 0,2 },{ 1,3 },{ 4,6 },{ 5,7 },// ||y
    { 0,4 },{ 1,5 },{ 2,6 },{ 3,7 },// ||z
    //{ 0,1 },{ 0,2 },{ 0,4 },
    //{ 3,1 },{ 3,2 },{ 3,7 },
    //{ 5,1 },{ 5,4 },{ 5,7 },
    //{ 6,2 },{ 6,4 },{ 6,7 },      // пары локальных номеров вершин 6-гренника - ребра
};

int edgeHexogen_extended[16][2] =
{
    { 0,1 },{ 2,3 },{ 4,5 },{ 6,7 },// ||x
    { 0,2 },{ 1,3 },{ 4,6 },{ 5,7 },// ||y
    { 0,4 },{ 1,5 },{ 2,6 },{ 3,7 },// ||z
    { 0,7 },{ 1,6 },{ 2,5 },{ 3,4 },
    //{ 0,1 },{ 0,2 },{ 0,4 },
    //{ 3,1 },{ 3,2 },{ 3,7 },
    //{ 5,1 },{ 5,4 },{ 5,7 },
    //{ 6,2 },{ 6,4 },{ 6,7 },      // пары локальных номеров вершин 6-гренника - ребра
};

namespace Grid
{

FEType FE_LinearHexagon::get_FEType()const
{
    return FEType::LinearHexagon;
}
void FE_LinearHexagon::getGeomVertexIndexes(int *vertexIndexes) const
{
    for (int t = 0; t < 8; t++) // t - локальный номер вершины
        vertexIndexes[t] = vi[t];
}
void FE_LinearHexagon::setGeomVertexIndexes(const int *new_vertexIndexes)
{
    for (int t = 0; t < 8; t++) // t - локальный номер вершины
        vi[t] = new_vertexIndexes[t];

}
void FE_LinearHexagon::getGeomVertexes(const std::vector<POINT3> &vertex, const std::vector<POINT3> &,
                                              POINT3 *v)const
{
    for (int t = 0; t < 8; t++) // t - локальный номер вершины
    {
        v[t] = vertex[vi[t]];
    }
}
void FE_LinearHexagon::getVertexIndexes(int *vertexIndexes)const
{
    for (int t = 0; t < 8; t++) // t - локальный номер вершины
        vertexIndexes[t] = vi[t];
}
void FE_LinearHexagon::setVertexIndexes(const int *new_vertexIndexes)
{
    for (int t = 0; t < 8; t++) // t - локальный номер вершины
        vi[t] = new_vertexIndexes[t];
}
void FE_LinearHexagon::calcBasisFuncValues(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const int numPoints, const double *integrationw, const VECTOR3 *dLinearBasCube, const VECTOR3 *,
                                  double *w, double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3])const
{
    POINT3 v[8];
    getGeomVertexes(vertex, vertexForCurvature, v);
    Fem::calcLinearHexagonLinearBasisFuncValues(v, numPoints, integrationw, dLinearBasCube,
                                                w, dbas);
}
void FE_LinearHexagon::calcBc2BasisFuncValues(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const Integration::Integrator &integrationFoursquare, const int faceIndex, double *basCube, double *hexagonCoef, VECTOR3 *hexagonNormal) const
{
    POINT3 v[27];
    // копируем вершины шестигранника
    getGeomVertexes(vertex, vertexForCurvature, v);
    POINT3_CUBE cubePoint;
    POINT3 dxyz[2];
    int c1 = HexagonFaceCoordinates[faceIndex][0]; // первый индекс свободной координаты
    int c2 = HexagonFaceCoordinates[faceIndex][1]; // второй индекс свободной координаты
    int cc = HexagonFaceCoordinates[faceIndex][2]; // индекс зафиксированной координаты
    int ccValue = HexagonFaceCoordinates[faceIndex][3];  // значение зафиксированной координаты

    for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
    {
        cubePoint[c1] = integrationFoursquare.p2[valIndex][0];
        cubePoint[c2] = integrationFoursquare.p2[valIndex][1];
        cubePoint[cc] = ccValue;
        // линейное либо квадратичное отображение
        cubeToHexagon(v, cubePoint, Fem::dif_XYZ3[c1], dxyz[0]);
        cubeToHexagon(v, cubePoint, Fem::dif_XYZ3[c2], dxyz[1]);
        // нормаль
        hexagonNormal[valIndex] = vector3Mul(dxyz[0], dxyz[1]);
        hexagonNormal[valIndex] = ccValue*hexagonNormal[valIndex]/hexagonNormal[valIndex].abs();
        hexagonCoef[valIndex] =
                integrationFoursquare.w[valIndex] * integrationFoursquare.detJ *
               sqrt(
               (SQR(dxyz[0][0]) +
                SQR(dxyz[0][1]) +
                SQR(dxyz[0][2]))
                *
               (SQR(dxyz[1][0]) +
                SQR(dxyz[1][1]) +
                SQR(dxyz[1][2]))
                -
               SQR(dxyz[0][0]*dxyz[1][0] +
                   dxyz[0][1]*dxyz[1][1] +
                   dxyz[0][2]*dxyz[1][2]));
        for (int m = 0; m < 8; m++)
            basCube[m*integrationFoursquare.size + valIndex] = Fem::cubeLagrange1_3D(cubePoint, m, Fem::dif_NULL3);
    }
}
void FE_LinearHexagon::cubeToHexagon(const POINT3 *v, const POINT3_CUBE &cubePoint, const DIF_STATE3 &dif,
                                     POINT3 &point)const
{
    Fem::cubeToLagrange1Hexagon(v, cubePoint, dif, point);
}
int FE_LinearHexagon::getFaceVertexIndexes(const int faceInd,
                                           int *vertexIndexes)const
{
    for(int i = 0; i < 4; i++)
    {
        vertexIndexes[i] = vi[HexagonFace[faceInd][i]];
    }
    return 4;
}
int FE_LinearHexagon::getFaceVertexIndexes(const int faceInd,
                                           int *vertexIndexes_global, int *vertexIndexes_local)const
{
    for(int i = 0; i < 4; i++)
    {
        vertexIndexes_global[i] = vi[HexagonFace[faceInd][i]];
        vertexIndexes_local[i] = HexagonFace[faceInd][i];
    }
    return 4;
}
void FE_LinearHexagon::getFaceVertexes(const std::vector<POINT3> &vertex, const int faceInd,
                                       POINT3 *v) const
{
    int vertexIndexes[4];
    getFaceVertexIndexes(faceInd, vertexIndexes);
    v[0] = vertex[vertexIndexes[0]];
    v[1] = vertex[vertexIndexes[1]];
    v[2] = vertex[vertexIndexes[2]];
    v[3] = vertex[vertexIndexes[3]];
}
void FE_LinearHexagon::getFaceOuterSphere(const std::vector<POINT3> &vertex, const int faceInd,
                                 Sphere &s)const
{
    POINT3 v[4];
    getFaceVertexes(vertex, faceInd,
                    v);
    s.init0(v[0]);
    s.expand(v[1]);
    s.expand(v[2]);
    s.expand(v[3]);
}
void FE_LinearHexagon::getFaceNormal(const POINT3 *v, const int faceIndex, const POINT2_CUBE p,
                                     VECTOR3 &normal) const
{
    POINT3_CUBE cubePoint;
    POINT3 dxyz[2];
    int c1 = HexagonFaceCoordinates[faceIndex][0]; // первый индекс свободной координаты
    int c2 = HexagonFaceCoordinates[faceIndex][1]; // второй индекс свободной координаты
    int cc = HexagonFaceCoordinates[faceIndex][2]; // индекс зафиксированной координаты
    int ccValue = HexagonFaceCoordinates[faceIndex][3];  // значение зафиксированной координаты
    cubePoint[c1] = p.x[0];
    cubePoint[c2] = p.x[1];
    cubePoint[cc] = ccValue;
    // линейное либо квадратичное отображение
    cubeToHexagon(v, cubePoint, Fem::dif_XYZ3[c1], dxyz[0]);
    cubeToHexagon(v, cubePoint, Fem::dif_XYZ3[c2], dxyz[1]);
    normal = vector3Mul(dxyz[0], dxyz[1]);
    normal = ccValue*normal/normal.abs();
}
bool FE_LinearHexagon::project(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const int faceInd, const POINT3 p, POINT3 &nearestPoint, VECTOR3 &normal, double &h) const
{
    // 1) находим нормаль в центре грани
    {
        POINT3 v[27];
        getGeomVertexes(vertex, vertexForCurvature,
                    v);
        getFaceNormal(v, faceInd, {0, 0},
                      normal);
    }

    // 2) разбиваем грань 6-гранника на 4 треугольника
    TRIANGLE tr[4];  // треугольники
    {
        POINT3 v[4];
        getFaceVertexes(vertex, faceInd,
                        v);
        POINT3 C = (v[0] + v[1] + v[2] + v[3]) / 4;    // центр грани
        if(HexagonFaceCoordinates[faceInd][3] == -1)
        {
            // вершины v[0], v[1], v[3], v[2] расположены против часовой стрелки
            tr[0].v[0] = v[1];
            tr[0].v[1] = v[0];
            tr[0].v[2] = C;

            tr[1].v[0] = v[3];
            tr[1].v[1] = v[1];
            tr[1].v[2] = C;

            tr[2].v[0] = v[2];
            tr[2].v[1] = v[3];
            tr[2].v[2] = C;

            tr[3].v[0] = v[0];
            tr[3].v[1] = v[2];
            tr[3].v[2] = C;
        }
        else
        {
            // вершины v[0], v[1], v[3], v[2] расположены по часовой стрелки
            tr[0].v[0] = v[0];
            tr[0].v[1] = v[1];
            tr[0].v[2] = C;

            tr[1].v[0] = v[1];
            tr[1].v[1] = v[3];
            tr[1].v[2] = C;

            tr[2].v[0] = v[3];
            tr[2].v[1] = v[2];
            tr[2].v[2] = C;

            tr[3].v[0] = v[2];
            tr[3].v[1] = v[0];
            tr[3].v[2] = C;
        }
        // вершины треугольников упорядочены так, чтобы буравчик вкручивался в треугольник в направлении внешней нормали
    }

    // 3) проецируем точку на каждый треугольник
    for(int trInd = 0; trInd < 4; trInd++)
    {
        POINT3 O;
        double t0;
        bool res = tr[trInd].intersect(p, normal,
                                       O, t0);
        nearestPoint = O;
        h = -t0;
        if(res)
        {
            return true;    // пересечение прямой n*t с плоскостью треугольника лежит внутри треугольника
            // то есть "проекция" найдена
        }
    }
    return false;   // пересечение прямой n*t с плоскостью треугольника лежит вне треугольника или прямая n*t параллельна плоскости треугольника
}


FEType FE_QuadraticHexagon::get_FEType()const
{
    return FEType::QuadraticHexagon;
}
void FE_QuadraticHexagon::getGeomVertexIndexes(int *vertexIndexes) const
{
    // узлы 6-гранника
    vertexIndexes[0] = vi[0];
    vertexIndexes[2] = vi[1];
    vertexIndexes[6] = vi[2];
    vertexIndexes[8] = vi[3];
    vertexIndexes[18] = vi[4];
    vertexIndexes[20] = vi[5];
    vertexIndexes[24] = vi[6];
    vertexIndexes[26] = vi[7];
    // дополнительные вершины для квадратичного отображения
    vertexIndexes[1] = vi_curvature[0];
    vertexIndexes[3] = vi_curvature[1];
    vertexIndexes[4] = vi_curvature[2];
    vertexIndexes[5] = vi_curvature[3];
    vertexIndexes[7] = vi_curvature[4];
    vertexIndexes[9] = vi_curvature[5];
    vertexIndexes[10] = vi_curvature[6];
    vertexIndexes[11] = vi_curvature[7];
    vertexIndexes[12] = vi_curvature[8];
    vertexIndexes[13] = vi_curvature[9];
    vertexIndexes[14] = vi_curvature[10];
    vertexIndexes[15] = vi_curvature[11];
    vertexIndexes[16] = vi_curvature[12];
    vertexIndexes[17] = vi_curvature[13];
    vertexIndexes[19] = vi_curvature[14];
    vertexIndexes[21] = vi_curvature[15];
    vertexIndexes[22] = vi_curvature[16];
    vertexIndexes[23] = vi_curvature[17];
    vertexIndexes[25] = vi_curvature[18];
}
void FE_QuadraticHexagon::setGeomVertexIndexes(const int *new_vertexIndexes)
{
    // узлы 6-гранника
    vi[0] = new_vertexIndexes[0];
    vi[1] = new_vertexIndexes[2];
    vi[2] = new_vertexIndexes[6];
    vi[3] = new_vertexIndexes[8];
    vi[4] = new_vertexIndexes[18];
    vi[5] = new_vertexIndexes[20];
    vi[6] = new_vertexIndexes[24];
    vi[7] = new_vertexIndexes[26];
    // дополнительные вершины для квадратичного отображения
    vi_curvature[0] = new_vertexIndexes[1];
    vi_curvature[1] = new_vertexIndexes[3];
    vi_curvature[2] = new_vertexIndexes[4];
    vi_curvature[3] = new_vertexIndexes[5];
    vi_curvature[4] = new_vertexIndexes[7];
    vi_curvature[5] = new_vertexIndexes[9];
    vi_curvature[6] = new_vertexIndexes[10];
    vi_curvature[7] = new_vertexIndexes[11];
    vi_curvature[8] = new_vertexIndexes[12];
    vi_curvature[9] = new_vertexIndexes[13];
    vi_curvature[10] = new_vertexIndexes[14];
    vi_curvature[11] = new_vertexIndexes[15];
    vi_curvature[12] = new_vertexIndexes[16];
    vi_curvature[13] = new_vertexIndexes[17];
    vi_curvature[14] = new_vertexIndexes[19];
    vi_curvature[15] = new_vertexIndexes[21];
    vi_curvature[16] = new_vertexIndexes[22];
    vi_curvature[17] = new_vertexIndexes[23];
    vi_curvature[18] = new_vertexIndexes[25];
}
void FE_QuadraticHexagon::getGeomVertexes(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature,
                                                 POINT3 *v)const
{
    // узлы 6-гранника
    v[0] = vertex[vi[0]];
    v[2] = vertex[vi[1]];
    v[6] = vertex[vi[2]];
    v[8] = vertex[vi[3]];
    v[18] = vertex[vi[4]];
    v[20] = vertex[vi[5]];
    v[24] = vertex[vi[6]];
    v[26] = vertex[vi[7]];
    // дополнительные вершины для квадратичного отображения
    v[1] = vertexForCurvature[vi_curvature[0]];
    v[3] = vertexForCurvature[vi_curvature[1]];
    v[4] = vertexForCurvature[vi_curvature[2]];
    v[5] = vertexForCurvature[vi_curvature[3]];
    v[7] = vertexForCurvature[vi_curvature[4]];
    v[9] = vertexForCurvature[vi_curvature[5]];
    v[10] = vertexForCurvature[vi_curvature[6]];
    v[11] = vertexForCurvature[vi_curvature[7]];
    v[12] = vertexForCurvature[vi_curvature[8]];
    v[13] = vertexForCurvature[vi_curvature[9]];
    v[14] = vertexForCurvature[vi_curvature[10]];
    v[15] = vertexForCurvature[vi_curvature[11]];
    v[16] = vertexForCurvature[vi_curvature[12]];
    v[17] = vertexForCurvature[vi_curvature[13]];
    v[19] = vertexForCurvature[vi_curvature[14]];
    v[21] = vertexForCurvature[vi_curvature[15]];
    v[22] = vertexForCurvature[vi_curvature[16]];
    v[23] = vertexForCurvature[vi_curvature[17]];
    v[25] = vertexForCurvature[vi_curvature[18]];
}
void FE_QuadraticHexagon::calcBasisFuncValues(const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const int numPoints, const double *integrationw, const VECTOR3 *dLinearBasCube, const VECTOR3 *dQuadraticBasCube,
                          double *w, double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3])const
{
    POINT3 v[27];
    getGeomVertexes(vertex, vertexForCurvature, v);
    Fem::calcQuadraticHexagonLinearBasisFuncValues(v, numPoints, integrationw, dLinearBasCube, dQuadraticBasCube,
                                                 w, dbas);
}

void FE_QuadraticHexagon::cubeToHexagon(const POINT3 *v, const POINT3_CUBE &cubePoint, const DIF_STATE3 &dif,
                                                   POINT3 &point)const
{
    Fem::cubeToLagrange2Hexagon(v, cubePoint, dif, point);
}






void FiniteElementSurface::init(const bool set_buildRegionsOnce, const Grid3D *set_grid)
{
    buildRegionsOnce = set_buildRegionsOnce;
    grid = set_grid;
    r0 = nullptr;
    grid = nullptr;
    surfaceType = SurfaceType::FiniteElementSurface;
}
void FiniteElementSurface::buildRegions()   // построение декомпозиции пространства
{
    std::vector<SpheredBody> sb;
    sb.clear();
    // 0) составления массива сферированных тел (поверхностей)
    for(size_t faceInd = 0; faceInd < face.size(); faceInd++)
    {
        SpheredBody sbEl;
        sbEl.bodyInd = faceInd;    // индекс тела = индекс грани
        grid->fe[face[faceInd].feIndex]->getFaceOuterSphere(grid->vertex, face[faceInd].faceIndex,
                                                            sbEl.s);
        sb.push_back(sbEl);
    }
    // 1) поиск параллелепипеда, который содержит все грани и инициализация корневой области
    CUBE q0;
    q0.initBySphere(sb[0].s);
    for(size_t i = 1; i < sb.size(); i++)
    {
        // каждый новый шар может расширить параллелепипед
        q0.expand(sb[i].s);
    }
    // нашли параллелепипед, который содержит все тела
    // преобразуем его в куб
    q0.toCube();
    // 2) инициализация корневой области
    r0 = new RegionWithSpheres;
    r0->init(q0);
    // 3) добавление сферированных тел (индекс + шар) в корневую область r0
    for(size_t i = 0; i < sb.size(); i++)
    {
        r0->insert(sb[i]);
    }
}
void FiniteElementSurface::findIntersectingBodies(const Sphere &s, std::vector<size_t> &bodyInd)
{
    bodyInd.clear();
    r0->findIntersectingBodies(s,
                               bodyInd);
}
void FiniteElementSurface::update(const double, const RegionDecompositionParameters *)
{
    if(buildRegionsOnce && r0 != nullptr)
    {
        // декомпозиция пространства уже построена и обновлять её не нужно
        return;
    }
    else
    {
        if(r0 != nullptr)
        {
            // удаление декомпозиции
            r0->release();
            delete r0;
        }
        // построение декомпозиции занова
        buildRegions();
    }
}
void FiniteElementSurface::findNearestPoint(const POINT3 &p, const VECTOR3 &u, const SurfacePositionData &prevSolutionData,
                                            SurfacePositionData &solutionData, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const
{
    // ищем, на какие грани можно "спроецировать"(по направлению нормали к грани) вершину p
    Sphere s;
    s.init0(p);     // сфера 0-го радиуса, соответствующая вершине подвижной поверхности
    SpheredBodyArray body;  // индексы граней, с которыми возможно пересечение шара s (с центром в p нулевого радиуса)
    body.clear();
    r0->findIntersectingBodies(s, body);
    for(size_t bodyInd = 0; bodyInd < body.size(); bodyInd++)
    {
        size_t fesInd = body[bodyInd];    // индекс грани, на которую может быть получится "спроецировать" узел p
        int feInd = face[fesInd].feIndex;          // индекс конечного элемента
        int faceLocalInd = face[fesInd].faceIndex; // локальный индекс грани
        POINT3 nearestPoint0;
        VECTOR3 normal0;
        double h0;
        bool res = grid->fe[feInd]->project(grid->vertex, grid->vertexForCurvature, faceLocalInd, p,
                                            nearestPoint0, normal0, h0);
        if(res)
        {
            nearestPoint = nearestPoint0;
            normal = normal0;
            //h = h0;
            if(h0 >= 0)
                side = 1;
            else
                side = -1;
            onBorder = false;
            return;
        }
    }
    return;   // неопредеоённое поведение ###
    // в некоторых случаях "проекция" может не найтись, но она должна существовать всегда!
}
void FiniteElementSurface::findIntersection(const POINT3 &p1, const POINT3 &p2, const SurfacePositionData &prevSolutionData1, const SurfacePositionData &prevSolutionData2, SurfacePositionData &solutionData1, SurfacePositionData &solutionData2,
                                            POINT3 &intersectionPoint, POINT3 &normal, int &side, bool &found)const
{
    // ищем, с какими гранями может быть пересечение отрезка p1-p2
    Sphere s;
    s.init0(p1);
    s.expand(p2);
    // шар s сожержит точки p1 и p2
    SpheredBodyArray body;  // индексы граней, с которыми возможно пересекается шар s
    body.clear();
    r0->findIntersectingBodies(s, body);
    for(size_t bodyInd = 0; bodyInd < body.size(); bodyInd++)
    {/*
        int fesInd = body[bodyInd];       // индекс грани, с которой возможено пересечение отрезка p1-p2
        int feInd = face[fesInd].feIndex; // индекс конечного элемента
        int faceLocalInd = face[fesInd].faceIndex; // индекс грани
        POINT3 nearestPoint0;
        VECTOR3 normal0;
        double direction0;*/
        /*
        bool res = grid->fe[feInd]->intersect(grid->vertex, grid->vertexForCurvature, faceLocalInd, p1, p2,
                                            nearestPoint0, normal0, direction0);
        if(res)
        {
            nearestPoint = nearestPoint0;
            normal = normal0;
            direction = direction0;
            return true;
        }
        */
    }
    found = false;   // пересечения нет
}


#define MALENKOE_CHISLO (1.e-20L)
void Grid3D::buldOpenSCADModel(const char *file_scad)
{
    printf("%s: ", file_scad);
    FILE *f = fopen(file_scad, "w");
    //if(mode.projection == 1)	// проекция на плоскость z=0
    //	fprintf(f, "projection(cut = false){");

    // ## поворот для сферы
    //fprintf(f, "rotate([45,-35,0]){\n");
    fprintf(f, "rotate([0,0,0]){\n");

    // отрисовка ребер
    for(size_t k = 0; k < fe.size(); k++)	// k - индекс конечного элемента
    {
        int vi[8];                          // глобальные индексы вершин шестигранника
        // копии глобальных индексов и координат вершин 6-гранника
        fe[k]->getVertexIndexes(vi);
        for(int i = 0; i < 12; i++)	// i - индекс ребра rebroHexogen[]
        {
            POINT3 p1, p2;
            VECTOR3 dp;
            double len, b, c;
            p1 = vertex[vi[edgeHexogen[i][0]]];
            p2 = vertex[vi[edgeHexogen[i][1]]];	// задано ребро p1-p2
            dp[0] = p2[0] - p1[0];
            dp[1] = p2[1] - p1[1];
            dp[2] = p2[2] - p1[2];	// вектор dp из точки p1 совпадет с ребром p1-p2
            len = sqrt(SQR(dp[0]) + SQR(dp[1]) + SQR(dp[2]));	// длина вектора
            b = acos(dp[2]/len)*180.L/PI;
            c = (fabs(dp[0]) < MALENKOE_CHISLO) ? (SIGN(dp[1])*90) : ((dp[0]>0) ? atan(dp[1]/dp[0])*180.L/PI : atan(dp[1]/dp[0])*180.L/PI+180);
            fprintf(f, "color([0.0,0.0,0.0,1])translate([%.16le,%.16le,%.16le])rotate([0,%.16le,%.16le])translate([0,0,%.16le])cube(size=[%.16le,%.16le,%.16le],center=true);\n",
                p1[0], p1[1], p1[2],
                b, c,
                len/2,
                0.03, 0.03, len);
/*            if(k == T_Nparts*(T_N*T_N - T_N)/2)
                fprintf(f, "color([0.0,0.0,1.0,1])translate([%.16le,%.16le,%.16le])rotate([0,%.16le,%.16le])translate([0,0,%.16le])cube(size=[%.16le,%.16le,%.16le],center=true);\n",
                    p1[0], p1[1], p1[2],
                    b, c,
                    len/2,
                    0.1L, 0.1L, len);
            else
            fprintf(f, "color([0.0,0.0,0.0,1])translate([%.16le,%.16le,%.16le])rotate([0,%.16le,%.16le])translate([0,0,%.16le])cube(size=[%.16le,%.16le,%.16le],center=true);\n",
                p1[0], p1[1], p1[2],
                b, c,
                len/2,
                0.03L, 0.03L, len);*/
            //fprintf(f, "color([0.0,0.0,0.0,1])translate([%lf,%lf,%lf])rotate([0, %lf, %lf])cylinder(h=%lf, r=%lf);\n", p1[0], p1[1], p1[2], b, c, len, mode.getCount()_cylinder);
        }
    }
    /*
    // отрисовка вершин
    for (size_t i = 0; i < vertex.size(); i++)
    {
        POINT3 xyz = vertex[i];	// координаты вершины
        // вершина
        fprintf(f, "color([0.0,0.0,0.0,1])translate([%.16le,%.16le,%.16le])cube(size=[%.16le,%.16le,%.16le],center=true);\n",
            xyz[0], xyz[1], xyz[2],
            0.03L, 0.03L, 0.03L);
    }
    // отрисовка дополнительных вершин для криволинейных шестигранников
    for (size_t i = 0; i < vertexForCurvature.size(); i++)
    {
        POINT3 xyz = vertexForCurvature[i];	// координаты вершины
        // вершина
        fprintf(f, "color([0.0,0.0,0.0,1])translate([%.16le,%.16le,%.16le])cube(size=[%.16le,%.16le,%.16le],center=true);\n",
            xyz[0], xyz[1], xyz[2],
            0.03L, 0.03L, 0.03L);
    }
    */
    // отрисовка первых краевых условий
    for (size_t i = 0; i < bc1.size(); i++)
    {
        if(bc1[i].bc1SourceIndex <= 1)
        //if(bc1[i].si <= 2)
        {
            POINT3 xyz = vertex[bc1[i].vertexIndex];	// координаты вершины
            // вершина
            fprintf(f, "color([1.0,0.0,0.0,1])translate([%.16le,%.16le,%.16le])cube(size=[%.16le,%.16le,%.16le],center=true);\n",
                xyz[0], xyz[1], xyz[2],
                0.1, 0.1, 0.1);
        }
    }
    // отрисовка вторых краевых условий
    for (size_t FEsurfaceInd = 0; FEsurfaceInd < FESurface.size(); FEsurfaceInd++)
    {
        FiniteElementSurface &s = FESurface[FEsurfaceInd];
        for (size_t faceInd = 0; faceInd < s.face.size(); faceInd++)
        {
            FEFace &face = s.face[faceInd];
            //if(bc2_el.si == 0)    //##
            //if(bc2_el.si == -10)
            {
                int vi[8];                  // глобальные индексы вершин шестигранника
                // копии глобальных индексов и координат вершин 6-гранника
                fe[face.feIndex]->getVertexIndexes(vi);
                for(int k = 0; k < 4; k++)  // k - индекс ребра rebro4[]
                {
                    POINT3 p1, p2;
                    VECTOR3 dp;
                    double len, b, c;
                    p1 = vertex[vi[edge4[face.faceIndex][k][0]]];
                    p2 = vertex[vi[edge4[face.faceIndex][k][1]]];        // задано ребро p1-p2
                    dp[0] = p2[0] - p1[0];
                    dp[1] = p2[1] - p1[1];
                    dp[2] = p2[2] - p1[2];	// вектор dp из точки p1 совпадет с ребром p1-p2
                    len = sqrt(SQR(dp[0]) + SQR(dp[1]) + SQR(dp[2]));	// длина вектора
                    b = acos(dp[2]/len)*180.L/PI;
                    c = (fabs(dp[0]) < MALENKOE_CHISLO) ? (SIGN(dp[1])*90) : ((dp[0]>0) ? atan(dp[1]/dp[0])*180.L/PI : atan(dp[1]/dp[0])*180.L/PI+180);
                    fprintf(f, "color([0.0,0.0,1.0,1])translate([%.16le,%.16le,%.16le])rotate([0,%.16le,%.16le])translate([0,0,%.16le])cube(size=[%.16le,%.16le,%.16le],center=true);\n",
                            p1[0], p1[1], p1[2],
                            b, c,
                            len/2,
                            0.1, 0.1, len);
                    //fprintf(f, "color([0.0,0.0,0.0,1])translate([%lf,%lf,%lf])rotate([0, %lf, %lf])cylinder(h=%lf, r=%lf);\n", p1[0], p1[1], p1[2], b, c, len, mode.getCount()_cylinder);
                }
            }

        }
    }
    fprintf(f, "}\n");

    fclose(f);
    printf("\n");
}

void HexagonArea::calcNodeCoord(const POINT3 *v, const VECTOR3_uint &i, POINT3 &x) const
{
    using namespace Operations;
    const int vi[3][4][2] =
    {
        {{0, 2}, {4, 6},  // пара рёбер при x = -1, вершины в порядке возрастания y(1), рёбра в порядке возрастания z(2)
         {1, 3}, {5, 7}}, // пара рёбер при x = +1, вершины в порядке возрастания y(1), рёбра в порядке возрастания z(2)
        {{0, 4}, {1, 5},  // пара рёбер при y = -1, вершины в порядке возрастания z(2), рёбра в порядке возрастания x(0)
         {2, 6}, {3, 7}}, // пара рёбер при y = +1, вершины в порядке возрастания z(2), рёбра в порядке возрастания x(0)
        {{0, 1}, {2, 3},  // пара рёбер при z = -1, вершины в порядке возрастания x(0), рёбра в порядке возрастания y(1)
         {4, 5}, {6, 7}}, // пара рёбер при z = +1, вершины в порядке возрастания x(0), рёбра в порядке возрастания y(1)
    };
    const int ci[3][3] =
    {
        {1, 2, 0}, // y, z, x
        {2, 0, 1}, // z, x, y
        {0, 1, 2}, // x, y, z
    };
    POINT3 p1;
    POINT3 p2;
    POINT3 p1_2;
    POINT3 p3;
    POINT3 p4;
    POINT3 p3_4;
    int c = condensation_coord;
    double q = condensation_q;
    if(c == -1)
    {
        c = 0;
        q = 1;
    }
    findPointOnTheLine_3d(v[vi[c][0][0]], v[vi[c][0][1]], N[ci[c][0]], 1, i[ci[c][0]], p1);
    findPointOnTheLine_3d(v[vi[c][1][0]], v[vi[c][1][1]], N[ci[c][0]], 1, i[ci[c][0]], p2);
    findPointOnTheLine_3d(v[vi[c][2][0]], v[vi[c][2][1]], N[ci[c][0]], 1, i[ci[c][0]], p3);
    findPointOnTheLine_3d(v[vi[c][3][0]], v[vi[c][3][1]], N[ci[c][0]], 1, i[ci[c][0]], p4);
    findPointOnTheLine_3d(p1, p2, N[ci[c][1]], 1, i[ci[c][1]], p1_2);
    findPointOnTheLine_3d(p3, p4, N[ci[c][1]], 1, i[ci[c][1]], p3_4);
    findPointOnTheLine_3d(p1_2, p3_4, N[ci[c][2]], q, i[ci[c][2]], x);
}
void HexagonArea::calcFaceStates(const VECTOR3_uint &i, bool *faceState) const
{
    faceState[0] = (i[2] == 0);
    faceState[1] = (i[2] == N[2]);
    faceState[2] = (i[0] == 0);
    faceState[3] = (i[0] == N[0]);
    faceState[4] = (i[1] == 0);
    faceState[5] = (i[1] == N[1]);
}
void HexagonArea::calcFEFaceStates(const VECTOR3_uint &i0, bool *faceState) const
{
    faceState[0] = (i0[2] == 0);
    faceState[1] = (i0[2] + 1 == N[2]);
    faceState[2] = (i0[0] == 0);
    faceState[3] = (i0[0] + 1 == N[0]);
    faceState[4] = (i0[1] == 0);
    faceState[5] = (i0[1] + 1 == N[1]);
}

// узел сетки(координаты)
struct NODE: public POINT3
{
    static const int digits = 8;
};
// функтор равенства узлов
struct NODE_equal
{
    bool operator()(const NODE &n1, const NODE &n2) const
    {
        return
                roundToDigits(n1[0], NODE::digits) == roundToDigits(n2[0], NODE::digits) &&
                roundToDigits(n1[1], NODE::digits) == roundToDigits(n2[1], NODE::digits) &&
                roundToDigits(n1[2], NODE::digits) == roundToDigits(n2[2], NODE::digits);
    }
};
// хеш функция от индекса вершины куба
struct NODE_hash
{
    size_t operator()(const NODE &n) const
    {
        size_t h0 = std::hash<double>()(roundToDigits(n[0], NODE::digits));
        size_t h1 = std::hash<double>()(roundToDigits(n[1], NODE::digits));
        size_t h2 = std::hash<double>()(roundToDigits(n[2], NODE::digits));
        return h0 ^ ((h1 ^ (h2 << 1)) << 1);
    }
};
// хэш таблица узлов сетки
typedef std::unordered_map<NODE, int, NODE_hash, NODE_equal> NodeMap;

int findNode(const NODE &x, NodeMap &m)
{
    NodeMap::const_iterator it = m.find(x);
    if(it == m.end())
    {
        // новый узел
        return -1;
    }
    else
    {
        return (*it).second;
    }
}
void Grid3D::build(const std::vector<POINT3> &areaVertex, const std::vector<HexagonArea> &hexagon)
{
    using namespace Operations;
    fe.clear();
    vertex.clear();
    vertexForCurvature.clear();
    bc1.clear();
    FESurface.clear();
    analiticalSurface.clear();
    ISurface.clear();

    // контейнер для хранения узлов
    NodeMap m;
    // количество поверхностей (максимальный индекс + 1)
    int maxSurfaceIndex = 0;

    // заполнение хэш таблицы и массива узлов сетки, задание 1-х краевых
    for (size_t hexagonInd = 0; hexagonInd < hexagon.size(); hexagonInd++)//индекс 6-гранной области
    {
        const HexagonArea &h = hexagon[hexagonInd];
        for (int faceInd = 0; faceInd < 6; faceInd++) // t - локальный номер грани
        {
            if(h.surfaceIndex[faceInd] > maxSurfaceIndex)
                maxSurfaceIndex = h.surfaceIndex[faceInd];
        }
        POINT3 v[8];
        // массив из 8-ми вершин области-шестигранника
        for (int t = 0; t < 8; t++) // t - локальный номер вершины
            v[t] = areaVertex[h.vi[t]];
        VECTOR3_uint i;
        for(i[2] = 0; i[2] <= h.N[2]; i[2]++)
        for(i[1] = 0; i[1] <= h.N[1]; i[1]++)
        for(i[0] = 0; i[0] <= h.N[0]; i[0]++)
        {
            NODE x;
            h.calcNodeCoord(v, i,
                             x);
            int vIndex = findNode(x, m);
            if(vIndex == -1)
            {
                // новая вершина
                // добавление вершины в массив вершин
                vertex.push_back(x);
                vIndex = (int)vertex.size() - 1;
                // добавление вершины в хеш таблицу вершин
                std::pair<const NODE, int> m_value = {x, vIndex};
                m.insert(m_value);
            }
            // первые краевые(возможны повторения но зато правильно)
            bool faceState[6];
            h.calcFaceStates(i, faceState);
            for (int faceInd = 0; faceInd < 6; faceInd++) // t - локальный номер грани
            {
                if(faceState[faceInd] && h.bc1Index[faceInd] != -1)
                {
                    // узел vIndex принадлежит грани faceInd и на этой грани задано 1-е краевое
                    BoundaryCondition1 bc1_el;
                    bc1_el.bc1SourceIndex = h.bc1Index[faceInd];
                    bc1_el.vertexIndex = vIndex;
                    bc1.push_back(bc1_el);
                }
            }

        }
    }

    // построение конечных элементов и поверхностей
    FESurface.resize(maxSurfaceIndex + 1);
    for (size_t FEsurfaceInd = 0; FEsurfaceInd < FESurface.size(); FEsurfaceInd++)
        FESurface[FEsurfaceInd].face.clear();
    for (size_t hexagonInd = 0; hexagonInd < hexagon.size(); hexagonInd++)//индекс 6-гранной области
    {
        const HexagonArea &h = hexagon[hexagonInd];
        POINT3 v[8];
        // массив из 8-ми вершин области-шестигранника
        for (int t = 0; t < 8; t++) // t - локальный номер вершины
            v[t] = areaVertex[h.vi[t]];
        VECTOR3_uint i0;
        for(i0[2] = 0; i0[2] < h.N[2]; i0[2]++)
        for(i0[1] = 0; i0[1] < h.N[1]; i0[1]++)
        for(i0[0] = 0; i0[0] < h.N[0]; i0[0]++)
        {
            // индексы вершин 6-гранника
            int vi[8];
            VECTOR3_uint di;
            int vi_index = 0;
            for(di[2] = 0; di[2] <= 1; di[2]++)
            for(di[1] = 0; di[1] <= 1; di[1]++)
            for(di[0] = 0; di[0] <= 1; di[0]++)
            {
                VECTOR3_uint i = i0 + di;
                NODE x;
                h.calcNodeCoord(v, i,
                                 x);
                int vIndex = findNode(x, m);
                // информация о вершине с индексом fe_i должна существовать в хеш таблице m
                if(vIndex == -1)
                {
                    for(;;);
                }
                //POINT3 v = vertex[vIndex];
                vi[vi_index] = vIndex;
                vi_index++;
            }
            // создание шестигранного КЭ
            FE_LinearHexagon *fe_el = new FE_LinearHexagon;   // шестигранник, линейное отображение
            fe_el->mi = h.mi;                                 // индекс материала
            for (int t = 0; t < 8; t++)
                fe_el->vi[t] = vi[t];
            fe.push_back(fe_el);
            // поверхности
            bool faceState[6];
            h.calcFEFaceStates(i0, faceState);
            for (int faceInd = 0; faceInd < 6; faceInd++) // t - локальный номер грани
            {
                if(faceState[faceInd] && h.surfaceIndex[faceInd] != -1)
                {
                    Grid::FEFace face;
                    face.feIndex = (int)fe.size() - 1;
                    face.faceIndex = faceInd;
                    FESurface[h.surfaceIndex[faceInd]].face.push_back(face);
                }
            }
        }
    }
    // нумерация базисных функций
    DOFs = new Grid::GlobalDOFs;
    DOFs->init(vertex.size());
    // ф-и первого порядка
    for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
    {
        Grid::FuncID funcID;
        funcID.type = Grid::FuncType::linear;
        funcID.crackIndex = -1;
        DOFs->addDOFs(vertexIndex, funcID);
    }
}
void Grid3D::save(const std::string &subdir) const
{
    std::string fn_vertex = subdir + "vertex.dat";
    std::string fn_fe = subdir + "fe.dat";
    FILE *f_vertex = fopen(fn_vertex.c_str(), "w");
    FILE *f_fe = fopen(fn_fe.c_str(), "w");
    // узлы
    fprintf(f_vertex, "%d\n", (int)vertex.size());
    for(int vertexInd = 0; vertexInd < (int)vertex.size(); vertexInd++)
    {
        for(int i = 0; i < 3; i++)
        {
            fprintf(f_vertex, "%le ", vertex[vertexInd][i]);
        }
        fprintf(f_vertex, "\n");
    }
    // КЭ
    fprintf(f_fe, "%d\n", (int)fe.size());
    for(int feInd = 0; feInd < (int)fe.size(); feInd++)
    {
        const FE_base *fe_el = fe[feInd];
        int vi_8[8];
        fe_el->getVertexIndexes(vi_8);
        for(int i = 0; i < 8; i++)
        {
            fprintf(f_fe, "%d ", vi_8[i]);
        }
        fprintf(f_fe, "\n");
    }
    fclose(f_vertex);
    fclose(f_fe);
}









FILE *debug_tree;
// узел
struct NDNode
{
    std::vector<size_t> edge;   // рёбра: массив индексов соседних узлов (в исходной нумерации узлов)
    int level;
    // level >= 0 - номер уровня если
    // level = -1 - узел не размечен
    // level = -2 - узел удалён
    NDNode()
    {
        edge.reserve(32);
        level = -1;
    }
    void addEdge(const size_t newEdgeVertexIndex)
    {
        for(size_t edgeInd = 0; edgeInd < edge.size(); edgeInd++)
        {
            if(edge[edgeInd] == newEdgeVertexIndex)
                return;// такое ребро уже есть
        }
        edge.push_back(newEdgeVertexIndex);
    }
    // возвращает true если узел связан с узлами уровня level + 1
    bool associatedWithHigherLevel(const std::vector<NDNode> &node)const
    {
        for(size_t edgeInd = 0; edgeInd < edge.size(); edgeInd++)
        {
            if(node[edge[edgeInd]].level == level + 1)
                return true;
        }
        return false;
    }
    size_t degree(const std::vector<NDNode> &node)
    {
        size_t E = 0;
        for(size_t i = 0; i < edge.size(); i++)
        {
            if(node[edge[i]].level != -2)
                E++;
        }
        return E;
        //return edge.size();
    }
};


/*
// очистка разметки всех узлов всех уровней (исходная разметка не важна)
void clearAllAreaLevels(std::vector<NDNode> &node, const size_t levelsNumber, const std::vector<size_t> &ind, const std::vector<size_t> &vi)
{
    for(size_t i = 0; i < ind[levelsNumber]; i++)
    {
        size_t vertexIndex = vi[i];     // индекс вершины подобласти
        node[vertexIndex].level = -1;   // очистка разметки
    }
}*/

// разметка нулями связанной подобласти с корнем в root
// подобласть состоит из узлов с разметкой level >= 0
void clearArea(const size_t root, std::vector<NDNode> &node, std::deque<size_t> &vi_temp, const int levelValue)
{
    vi_temp.clear();
    node[root].level = -1; // очистка разметки корневого узла
    vi_temp.push_back(root); // добавление узла в очередь
    // обход в ширину
    for(;;)
    {
        if(vi_temp.empty())
            break;
        // взятие из начала очереди очередного узла
        size_t vertexIndex = vi_temp.front();
        vi_temp.pop_front();
        NDNode &nodeEl = node[vertexIndex]; // узел с очищенной разметкой
        // добавление в конец очереди соседних узлов, которые ещё не очищены
        for(size_t edgeInd = 0; edgeInd < nodeEl.edge.size(); edgeInd++)
        {
            size_t associatedVertexIndex = nodeEl.edge[edgeInd];
            if(node[associatedVertexIndex].level >= 0)
            {
                node[associatedVertexIndex].level = levelValue; // очистка разметки связанного узла
                vi_temp.push_back(associatedVertexIndex); // добавление связанного узла в очередь
            }
        }
    }
}

// построение уровней из корня root
// возвращает levelsNumber - количество уровней
// vi[ind[l]]...vi[ind[l + 1] - 1] - индексы узлов уровня l
void buildLevels(const size_t root, std::vector<NDNode> &node,
                 size_t &levelsNumber, std::vector<size_t> &ind, std::vector<size_t> &vi)
{
    size_t l = 0;
    node[root].level = 0;
    ind[0] = 0;
    ind[1] = 1;
    vi[0] = root;
    for(;;)
    {
        l++;
        ind[l + 1] = ind[l];
        for(size_t i = ind[l - 1]; i < ind[l]; i++)
        {
            size_t vertexIndex = vi[i];     // индекс вершины уровня (l - 1)
            NDNode nodeEl = node[vertexIndex];
            for(size_t edgeInd = 0; edgeInd < nodeEl.edge.size(); edgeInd++)
            {
                size_t associatedVertexIndex = nodeEl.edge[edgeInd];
                if(node[associatedVertexIndex].level == -1)
                {
                    node[associatedVertexIndex].level = l;
                    vi[ind[l + 1]] = associatedVertexIndex;
                    ind[l + 1]++;
                }
            }
        }
        if(ind[l] == ind[l + 1])
            break;
    }
    levelsNumber = l;
    //fprintf(stderr, "l = %d\n", (int)l);
}

// поиск псевдопереферийного узла подобласти root + построение уровней с корнем в этом узле
// возвращает bestRoot - псевдоперефирийный узел и строит уровни с корнем в bestRoot
// levelsNumber - количество уровней
// vi[ind[l]]...vi[ind[l + 1] - 1] - индексы узлов уровня l
void findBestRoot(const size_t root, std::vector<NDNode> &node,
                  size_t &bestRoot, size_t &levelsNumber, std::vector<size_t> &ind, std::vector<size_t> &vi)
{
    fprintf(debug_tree, "findBestRoot..");
    buildLevels(root, node, levelsNumber, ind, vi);
    for(;;)
    {
        size_t lastLevel_start = ind[levelsNumber - 1];
        size_t lastLevel_end = ind[levelsNumber];
        // поиск узла с минимальной степенью из последнего уровня (levelsNumber - 1)
        size_t i = lastLevel_start;
        size_t minDegree_vertexIndex = vi[i];
        size_t minDegree = node[minDegree_vertexIndex].degree(node);
        for(;;)
        {
            i++;
            if(!(i < lastLevel_end))
                break;
            size_t vertexIndex = vi[i];     // индекс вершины подобласти
            if(node[vertexIndex].degree(node) < minDegree)
            {
                minDegree_vertexIndex = vertexIndex;
                minDegree = node[vertexIndex].degree(node);
            }
        }
        // теперь корень minDegree_vi
        size_t prev_levelsNumber = levelsNumber;
        for(size_t i = 0; i < ind[levelsNumber]; i++)
        {
            size_t vertexIndex = vi[i];     // индекс вершины подобласти
            node[vertexIndex].level = -1;   // очистка разметки
        }
        //clearAllAreaLevels(node, levelsNumber, ind, vi);
        fprintf(debug_tree, "%d ", (int)levelsNumber);
        buildLevels(minDegree_vertexIndex, node, levelsNumber, ind, vi);
        if(levelsNumber <= prev_levelsNumber)
        {
            bestRoot = minDegree_vertexIndex;
            fprintf(debug_tree, "%d ", (int)levelsNumber);
            return;
        }
    }
}

struct SeparationAreas
{
    std::vector<size_t> root; // индексы корней подобластей, образовавшихся после удаления разделителя
};


// разделение подобласти root разделителем по центру из псевдоперефирийного узла
// если разделитель найден, то
//  возвращается true и узлы полученных после разделения подобластей: root_1 и root_2,
//  индексы узлов разделителя добавляются в newOrder, разделитель удаляется из графа
// если разделитель не найден, то
//  возвращается false
//  индексы узлов всей подобласти добавляются в newOrder, узлы удаляются из графа
void selectSeparator(const size_t root, std::vector<NDNode> &node, std::vector<size_t> &newOrder_seporator, std::vector<size_t> &newOrder_nonseporator,
                           std::vector<size_t> &ind, std::vector<size_t> &vi, std::deque<size_t> &vi_temp, SeparationAreas &sa)
{
    fprintf(debug_tree, "sep...");
    // 1) поиск псевдопереферийного узла
    size_t bestRoot;
    size_t levelsNumber;
    findBestRoot(root, node,
                 bestRoot, levelsNumber, ind, vi);
    fprintf(debug_tree, " size = %d\n", (int)ind[levelsNumber]);
    int debug_count = 0;
    // 2) поиск разделителя по центру
    /*
    size_t nodesNomber = ind[levelsNumber];
    size_t sepLevel = 1;
    for(; sepLevel < levelsNumber; sepLevel++)
    {
        //fprintf(debug_tree, "size(%d) = %d\n", (int)sepLevel, (int)ind[sepLevel]);
        if(ind[sepLevel] > nodesNomber/2)
        {
            sepLevel--;
            break;
        }
    }
    if(sepLevel == levelsNumber && levelsNumber >= 3)
    {
        sepLevel -= 2;
    }
    fprintf(debug_tree, "sepLevel = %d\n", (int)sepLevel);
    // 3)
    if(sepLevel == 0 || sepLevel >= levelsNumber - 1)
    {
        for(size_t i = 0; i < ind[levelsNumber]; i++)
        {
            size_t vertexIndex = vi[i];     // индекс вершины подобласти
            node[vertexIndex].level = -2;   // удаление
            newOrder_nonseporator.push_back(vertexIndex);// сбрасывание
            debug_count++;
        }
        fprintf(debug_tree, "|R| = %d\n", debug_count);
        return;
    }*/
    // 2) поиск разделителя по центру
    size_t sepLevel = (levelsNumber - 1) / 2;
    // 3) если уровней <= 2 то разделитель найти нельзя - сбрасываеется в newOrder вся область
    if(levelsNumber <= 2)
    {
        debug_count = 0;
        for(size_t i = 0; i < ind[levelsNumber]; i++)
        {
            size_t vertexIndex = vi[i];     // индекс вершины подобласти
            node[vertexIndex].level = -2;   // удаление
            newOrder_nonseporator.push_back(vertexIndex);// сбрасывание
            debug_count++;
        }
        fprintf(debug_tree, "|R| = %d\n", debug_count);
        return;
    }
    // 4) сбрасывание разделителя в newOrder
    // можно убирать из разделителя узлы, которые не связаны с узлами на 1 большего уровня
    // убранные из разделителя узлы помечаются level = 0
    debug_count = 0;
    for(size_t i = ind[sepLevel]; i < ind[sepLevel + 1]; i++)
    {
        size_t vertexIndex = vi[i];     // индекс вершины уровня sepLevel
        if(!node[vertexIndex].associatedWithHigherLevel(node))
        {
            node[vertexIndex].level = 0;   // узел перестал принадлежать разделителю
        }
        else
        {
            node[vertexIndex].level = -2;   // удаление
            newOrder_seporator.push_back(vertexIndex);// сбрасывание
            debug_count++;
        }
    }
    fprintf(debug_tree, "|S| = %d\n", debug_count);
    // 5) очистка разметки уровней 0...sepLevel
    for(size_t i = 0; i < ind[sepLevel + 1]; i++)
    {
        size_t vertexIndex = vi[i];     // индекс вершины уровня <= sepLevel
        if(node[vertexIndex].level >= 0)
        {
            node[vertexIndex].level = -1;   // очистка разметки
        }
    }
    // 6) поиск связных компонент исходной подобласти, которые получились после удаления разделителя
    //    и добавление корней каждой из связанной подобласти в очередь
    sa.root.push_back(vi[ind[0]]); // уровни 0...(sepLevel-1)+часть уровня sepLevel связаны (по построению уровней) и эта область задаётся корнем vi[ind[0]]=bestRoot
    debug_count = 1;
    for(size_t i = ind[sepLevel + 1]; i < ind[sepLevel + 2]; i++)
    {
        size_t vertexIndex = vi[i];     // индекс вершины уровня (sepLevel + 1)
        if(node[vertexIndex].level >= 0)
        {
            // поиск связной подобласти с корнем в vertexIndex
            // разметка узлов подобласти с корнем vertexIndex очищается (level = -1)
            clearArea(vertexIndex, node, vi_temp, -1);
            sa.root.push_back(vertexIndex);    // добавление корня связанной подобласти в очередь
            debug_count++;
        }
    }
    fprintf(debug_tree, " |P| = %d\n", debug_count);

}

void Grid3D::nestedDissection(const std::vector<int> &contactFESurfaceIndex)
{
debug_tree = fopen("__tree.txt", "w");
TimeIntervals::timeInterval t;
t.begin();
    const size_t size = vertex.size();
    // инициализация
    std::vector<NDNode> node; // узлы со связями
    node.resize(size);
    std::vector<size_t> newOrder_seporator; // изменения индексов узлов: newOrder[i] - старый индекс вершины i
    std::vector<size_t> newOrder_nonseporator; // изменения индексов узлов: newOrder[i] - старый индекс вершины i
    //newOrder.reserve(size); // (изначально пуст)
    std::vector<SeparationAreas> area;  // множество несвязных подграфов, задаются одним узлом
    area.reserve(size*100);//##
    //area.reserve(size);
    std::vector<size_t> ind;
    std::vector<size_t> vi;
    ind.resize(size + 1);
    vi.resize(size);
    std::deque<size_t> vi_temp;
    //SeparationNode *sn = new SeparationNode;

    // обход всех КЭ и составление массивов рёберных связей для каждого узла
    for (size_t feInd = 0; feInd < fe.size(); feInd++)   // индекс конечного элемента
    {
        int vi_8[8];
        fe[feInd]->getVertexIndexes(vi_8);
        // каждая пара вершин 6-гранника связана между собой
        for(int i = 0; i < 8; i++)	// i - локальный индекс вершины
        {
            for(int j = 0; j < 8; j++)	// i - локальный индекс вершины
            {
                size_t vertexIndex1 = vi_8[i];
                size_t vertexIndex2 = vi_8[j]; // vertexIndex1 - vertexIndex2 - индексы вершин связи
                node[vertexIndex1].addEdge(vertexIndex2);
                node[vertexIndex2].addEdge(vertexIndex1);
            }
        }
    }

    // нумерация контактных вершин
    // все контактные узлы, и которые связаны с ними
    for(size_t csIndex = 0; csIndex < contactFESurfaceIndex.size(); csIndex++)
    {
        size_t FESurfaceInd = contactFESurfaceIndex[csIndex]; // индекс КЭ поверхности, которая будет учавствовать в контакте
        std::vector<Grid::FEFace> &face = FESurface[FESurfaceInd].face;    // набор граней поверхности
        // сами контактные вершины
        for(size_t faceInd = 0; faceInd < face.size(); faceInd++)
        {
            int vi_4[4];   // индексы вершин, принадлежащие данной грани
            size_t surfaceVertexesNumber = fe[face[faceInd].feIndex]->getFaceVertexIndexes(face[faceInd].faceIndex, vi_4);
            for(size_t i = 0; i < surfaceVertexesNumber; i++)
            {
                size_t vertexIndex = vi_4[i]; // индекс вершины грани
                if(node[vertexIndex].level != -2)
                {
                    node[vertexIndex].level = -2;   // удаление
                    newOrder_seporator.push_back(vertexIndex);// сбрасывание
                }
            }
        }
        // связанные вершины
        for(size_t faceInd = 0; faceInd < face.size(); faceInd++)
        {
            int vi_4[4];   // индексы вершин, принадлежащие данной грани
            size_t surfaceVertexesNumber = fe[face[faceInd].feIndex]->getFaceVertexIndexes(face[faceInd].faceIndex, vi_4);
            // сами контактные вершины
            for(size_t i = 0; i < surfaceVertexesNumber; i++)
            {
                size_t vertexIndex = vi_4[i]; // индекс вершины грани
                NDNode nodeEl = node[vertexIndex];
                for(size_t edgeInd = 0; edgeInd < nodeEl.edge.size(); edgeInd++)
                {
                    size_t associatedVertexIndex = nodeEl.edge[edgeInd];
                    if(node[associatedVertexIndex].level != -2)
                    {
                        node[associatedVertexIndex].level = -2; // удаление
                        newOrder_seporator.push_back(associatedVertexIndex);// сбрасывание
                    }
                }
            }
        }
    }

    // изначально задана 1 область с любым корнем
    for(size_t i = 0; i < size; i++)
    {
        if(node[i].level != -2)
        {
            area.push_back({});
            area.back().root.push_back(i);
            break;
        }
    }
    for(;;)
    {
         // I) охуенный способ
        if(area.empty())
            break;
        // извлечение корней подобластей последнего уровня
        SeparationAreas &sa = area.back();
        if(sa.root.empty())
        {
            area.pop_back();
        }
        else
        {
            // стандартный вариант
            size_t root = sa.root.back();
            sa.root.pop_back();
            area.push_back({});
            // разделение подобласти root разделителем и добавление корней оставшихся после разделения подобластей
            selectSeparator(root, node, newOrder_seporator, newOrder_seporator,
                            ind, vi, vi_temp, area.back());

            /*
            // самый норм вариант(сдвинуты сечения)
            for(size_t rootIndex = 0; rootIndex < sa.root.size(); rootIndex++)
            {
                size_t root = sa.root[rootIndex];
                area.push_back({});
                // разделение подобласти root разделителем и добавление корней оставшихся после разделения подобластей
                selectSeparator(root, node, newOrder_seporator, newOrder_seporator,
                                ind, vi, vi_temp, area.back());
            }
            sa.root.clear();
            */

            /*
            // ещё вариант: сечения сдвинуты, но наоборот - какашка: медленное символическое разложение, медленное численное разложение
            if(sa.root.size() != 0)
            for(size_t rootIndex = sa.root.size() - 1;;)
            {
                size_t root = sa.root[rootIndex];
                area.push_back({});
                // разделение подобласти root разделителем и добавление корней оставшихся после разделения подобластей
                selectSeparator(root, node, newOrder_seporator, newOrder_seporator,
                                ind, vi, vi_temp, area.back());
                if(rootIndex == 0)
                    break;
                rootIndex--;
            }
            sa.root.clear();*/

        }
        // II) стандартный способ
        /*
        if(area.empty())
            break;
        // извлечение корней подобластей последнего уровня
        SeparationAreas &sa = area.back();
        if(sa.root.empty())
        {
            area.pop_back();
        }
        else
        {
            size_t &root = sa.root.back();
            sa.root.pop_back();
            area.push_back({});
            // разделение подобласти root разделителем и добавление корней оставшихся после разделения подобластей
            selectSeparator(root, node, newOrder_seporator, newOrder_seporator,
                            ind, vi, vi_temp, area.back());
        }*/
        /*
        if(area.empty())
            break;
        // 1) извлечение корня очередной подобласти
        size_t root = area.back();
        area.pop_back();
        // 2) разделение подобласти root разделителем и добавление корней оставшихся после разделения подобластей
        //selectSeparator(root, node, newOrder_seporator, newOrder_nonseporator,
        //                   ind, vi, vi_temp, area);
        selectSeparator(root, node, newOrder_seporator, newOrder_seporator,
                           ind, vi, vi_temp, area);
        */
    }

    /*for(int ll = 0;;ll++)
    {
        if(area.empty())
            break;
        size_t r = area.front();
        area.pop_front();
        fprintf(stderr, "r = %d\n", (int)r);
        fprintf(stderr, "r level = %d\n", (int)node[r].level);
        size_t ln;
        buildLevels(r, node, ln, ind, vi);
        fprintf(stderr, "size(r) = %d\n", (int)ind[ln]);
        for(size_t i = 0; i < ind[ln]; i++)
        {
            size_t vertexIndex = vi[i];     // индекс вершины подобласти
            node[vertexIndex].level = ll;   // удаление
            newOrder.push_back(vertexIndex);// сбрасывание
        }
    }
    for(size_t i = 0; i < size; i++)
    {
        NDNode &nodeEl = node[i];
        int l1 = nodeEl.level;
        if(l1 != -2)
        for(size_t edgeInd = 0; edgeInd < nodeEl.edge.size(); edgeInd++)
        {
            size_t associatedVertexIndex = nodeEl.edge[edgeInd];
            int l2 = node[associatedVertexIndex].level;
            if(l2 != -2 && l2 != l1)
            {
                fprintf(stderr, "vi: %d -> %d level: %d -> %d\n", (int)i, (int)associatedVertexIndex, l1, l2);
            }
        }
    }*/


    fprintf(debug_tree, "size = %d\n", (int)size);
    fprintf(stderr, "size = %d\n", (int)size);
    fprintf(stderr, "newOrder_seporator.size() = %d\n", (int)newOrder_seporator.size());
    fprintf(stderr, "newOrder_nonseporator.size() = %d\n", (int)newOrder_nonseporator.size());
    fprintf(stderr, "newOrder.size = %d\n", (int)newOrder_seporator.size() + (int)newOrder_nonseporator.size());

    // заполнение искомого массива перенумировки
    std::vector<size_t> reversedNewOrder; // изменения индексов узлов в обратном виде и обратном порядке: newOrder[vi] - новый индекс вершины vi
    reversedNewOrder.resize(size);
    for(size_t i = 0; i < newOrder_nonseporator.size(); i++)
    {
        newOrder_seporator.push_back(newOrder_nonseporator[i]);
    }
    for(size_t i = 0; i < size; i++)
        reversedNewOrder[newOrder_seporator[i]] = size - 1 - i;
        //reversedNewOrder[i] = newOrder[i];


    // перенумерация индексов вершин: vertexIndex -> node[vertexIndex].newVertexIndex
    // fe
    for (size_t feInd = 0; feInd < fe.size(); feInd++)   // индекс конечного элемента
    {
        int vi_8[8];
        fe[feInd]->getVertexIndexes(vi_8);
        int new_vi_8[8];
        for(int i = 0; i < 8; i++)
        {
            new_vi_8[i] = reversedNewOrder[vi_8[i]];
        }
        fe[feInd]->setVertexIndexes(new_vi_8);
    }
    // bc1
    for (size_t bc1Ind = 0; bc1Ind < bc1.size(); bc1Ind++)
    {
        bc1[bc1Ind].vertexIndex = reversedNewOrder[bc1[bc1Ind].vertexIndex];
    }
    // перестановка самих вершин
    std::vector<POINT3> vertex_copy;
    vertex_copy.resize(vertex.size());
    for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
    {
        vertex_copy[vertexIndex] = vertex[vertexIndex];
    }
    for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
    {
        vertex[reversedNewOrder[vertexIndex]] = vertex_copy[vertexIndex];
    }
fprintf(stderr, "nestedDissection time = %le\n", t.getCurrentDuration());
    fclose(debug_tree);
}


struct GraphNode
{
    //size_t vertexIndex;       // индекс узла в исходной нумерации
    size_t newVertexIndex;      // индекс узла после сортировки
    bool isSorted;              // = true если узлу присвоен новый индекс
    std::vector<size_t> edge;   // рёбра: массив индексов соседних узлов (в исходной нумерации)
    GraphNode()
    {
        edge.reserve(32);
        isSorted = false;
    }
    void addEdge(const size_t newEdgeVertexIndex)
    {
        for(size_t edgeInd = 0; edgeInd < edge.size(); edgeInd++)
        {
            if(edge[edgeInd] == newEdgeVertexIndex)
                return;// такое ребро уже есть
        }
        edge.push_back(newEdgeVertexIndex);
    }
    /*
    void setNumberAndDel(size_t &currentMaxVertexIndex)
    {
        newVertexIndex = currentMaxVertexIndex;
        currentMaxVertexIndex--;
        isSorted = true;
    }*/
    void setNewVertexIndex(const size_t push_vertexIndex, std::deque<size_t> &activeVertexIndex, size_t &currentMaxVertexIndex)
    {
        if(!isSorted)
        {
            newVertexIndex = currentMaxVertexIndex;
            currentMaxVertexIndex--;
            isSorted = true;
            activeVertexIndex.push_back(push_vertexIndex);
        }
    }
    // возвращает степень узла
    size_t degree(const std::vector<GraphNode> &node)
    {
        size_t E = 0;
        for(size_t i = 0; i < edge.size(); i++)
        {
            if(!node[edge[i]].isSorted)
                E++;
        }
        return E;
    }
};

// индекс + степень узла (для сортировки)
struct NodeWithDegree
{
    size_t vertexIndex;
    size_t degree;
};

void Grid3D::cuthillMcKee(const std::vector<int> &contactFESurfaceIndex)
{
TimeIntervals::timeInterval t;
t.begin();
    std::vector<GraphNode> node;
    node.resize(vertex.size());

    // обход всех КЭ и составление массивов рёберных связей для каждого узла
    for (size_t feInd = 0; feInd < fe.size(); feInd++)   // индекс конечного элемента
    {
        int vi_8[8];
        fe[feInd]->getVertexIndexes(vi_8);
        // каждая пара вершин 6-гранника связана между собой
        for(int i = 0; i < 8; i++)	// i - локальный индекс вершины
        {
            for(int j = 0; j < 8; j++)	// i - локальный индекс вершины
            {
                size_t vertexIndex1 = vi_8[i];
                size_t vertexIndex2 = vi_8[j]; // vertexIndex1 - vertexIndex2 - индексы вершин связи
                node[vertexIndex1].addEdge(vertexIndex2);
                node[vertexIndex2].addEdge(vertexIndex1);
            }
        }
    } // feInd

    // нумерация контактных вершин самыми большими индексами (начиная с vertex.size() - 1)
    // и инициализация массива активных узлов
    size_t currentMaxVertexIndex = vertex.size() - 1;
    std::deque<size_t> activeVertexIndex;       // массив активных узлов (индексы узлов в старой нумерации)
    for(size_t csIndex = 0; csIndex < contactFESurfaceIndex.size(); csIndex++)
    {
        size_t FESurfaceInd = contactFESurfaceIndex[csIndex]; // индекс КЭ поверхности, которая будет учавствовать в контакте
        std::vector<Grid::FEFace> &face = FESurface[FESurfaceInd].face;    // набор граней поверхности
        for(size_t faceInd = 0; faceInd < face.size(); faceInd++)
        {
            // контактные узлы и узлы соседних КЭ
            /*
            int vi_8[8];
            fe[face[faceInd].feIndex]->getVertexIndexes(vi_8);
            for(int i = 0; i < 8; i++)
            {
                node[vi_8[i]].setNewVertexIndex(vi_8[i], activeVertexIndex, currentMaxVertexIndex);
            }*/
            // только сами контактные узлы
            int vi_4[4];   // индексы вершин, принадлежащие данной грани
            int surfaceVertexesNumber = fe[face[faceInd].feIndex]->getFaceVertexIndexes(face[faceInd].faceIndex, vi_4);
            for(int i = 0; i < surfaceVertexesNumber; i++)
            {
                //vi_4[i] - индекс вершины грани
                node[vi_4[i]].setNewVertexIndex(vi_4[i], activeVertexIndex, currentMaxVertexIndex);
            }
        }
    }
    // если контактных узлов нет, то инициализируем узлом с индексом vertex.size() - 1
    if(activeVertexIndex.empty())
    {
        size_t firstVertexIndex = vertex.size() - 1;
        node[firstVertexIndex].setNewVertexIndex(firstVertexIndex, activeVertexIndex, currentMaxVertexIndex);
    }

    std::vector<NodeWithDegree> nwd;
    nwd.reserve(64);
    // обход в ширину
    for(;;)
    {
        if(activeVertexIndex.empty())
            break;
        // взятие из начала очереди очередного активного узла
        size_t activeVertexIndex_el = activeVertexIndex.front();
        activeVertexIndex.pop_front();
        // составление массива неудалённых соседей узла activeVertexIndex_el
        GraphNode &nodeEl = node[activeVertexIndex_el]; // активный узел
        for(size_t edgeInd = 0; edgeInd < nodeEl.edge.size(); edgeInd++)
        {
            size_t associatedVertexIndex = nodeEl.edge[edgeInd];
            GraphNode &agn = node[associatedVertexIndex]; // связанный узел
            if(!agn.isSorted)
            {
                // ещё не пронумерован
                NodeWithDegree nwd_el;
                nwd_el.vertexIndex = associatedVertexIndex;
                nwd_el.degree = agn.degree(node);
                nwd.push_back(nwd_el);
            }
        }
        // сортировка в порядке возрастания степени
        for(size_t i = 0; i < nwd.size(); i++)
            for(size_t j = i + 1; j < nwd.size(); j++)
                if(nwd[i].degree > nwd[j].degree)
                    std::swap(nwd[i], nwd[j]);
        // удаление из графа и добавление в конец очереди соседних узлов (в порядке возрастания степени)
        for(size_t i = 0; i < nwd.size(); i++)
        {
            size_t vi = nwd[i].vertexIndex;
            node[vi].setNewVertexIndex(vi, activeVertexIndex, currentMaxVertexIndex);
        }
        nwd.clear();
    }
    fprintf(stderr, "currentMaxVertexIndex = %d\n", (int)currentMaxVertexIndex);

    // перенумерация индексов вершин: vertexIndex -> node[vertexIndex].newVertexIndex
    // fe
    for (size_t feInd = 0; feInd < fe.size(); feInd++)   // индекс конечного элемента
    {
        int vi_8[8];
        fe[feInd]->getVertexIndexes(vi_8);
        int new_vi_8[8];
        for(int i = 0; i < 8; i++)
        {
            new_vi_8[i] = node[vi_8[i]].newVertexIndex;
        }
        fe[feInd]->setVertexIndexes(new_vi_8);
    } // feInd
    // bc1
    for (size_t bc1Ind = 0; bc1Ind < bc1.size(); bc1Ind++)
    {
        bc1[bc1Ind].vertexIndex = node[bc1[bc1Ind].vertexIndex].newVertexIndex;;
    }

    // перестановка самих вершин
    std::vector<POINT3> vertex_copy;
    vertex_copy.resize(vertex.size());
    for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
    {
        vertex_copy[vertexIndex] = vertex[vertexIndex];
    }
    for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
    {
        //vertex[vertexIndex] = vertex_copy[node[vertexIndex].newVertexIndex];
        vertex[node[vertexIndex].newVertexIndex] = vertex_copy[vertexIndex];
    }
fprintf(stderr, "cuthillMcKee time = %le\n", t.getCurrentDuration());
}





// a[4] - 4 точки 4-угольной области
// i - 3 индекса
void getAreaPoint(const POINT2 *a, const double z0, const double z1, const VECTOR3_uint N, const VECTOR3_uint i,
               POINT3 &v)
{
    using namespace Operations;
    double z;
    POINT2 v_bot;
    POINT2 v_top;
    // 1) z
    findPointOnTheLine_1d(z0, z1, N[2], 1, i[2], z);
    // 2) точка на нижней стороне области (v_bot)
    findPointOnTheLine_1d(a[0][0], a[1][0], N[0], 1, i[0], v_bot[0]);
    findPointOnTheLine_1d(a[0][1], a[1][1], N[0], 1, i[0], v_bot[1]);
    // 3) точка на верхней стороне области (v_top)
    findPointOnTheLine_1d(a[2][0], a[3][0], N[0], 1, i[0], v_top[0]);
    findPointOnTheLine_1d(a[2][1], a[3][1], N[0], 1, i[0], v_top[1]);
    // 4) искомая точка на отрезке v_bot - v_top
    findPointOnTheLine_1d(v_bot[0], v_top[0], N[1], 1, i[1], v[0]);
    findPointOnTheLine_1d(v_bot[1], v_top[1], N[1], 1, i[1], v[1]);
    v[2] = z;
}
// идентификатор вершины области
/*struct CC_vertexIndex
{
    unsigned int areaIndex; // индекс подобласти
    VECTOR3_uint i;         // индексы координат узла подобласти
    void getVertex(const ContactCylinderParameters &ccp,
                  POINT3 &v) const
    {
        getAreaPoint(ccp.a[areaIndex], ccp.z0, ccp.z1, ccp.N[areaIndex], i,
                     v);
    }

};
// функтор равенства индексов вершин куба
struct CC_vertexIndex_equal
{
    bool operator()(const CC_vertexIndex &e1, const CC_vertexIndex &e2) const
    {
        return e1.areaIndex == e2.areaIndex && e1.i == e2.i;
    }
};
// хеш функция от индекса вершины куба
struct CC_vertexIndex_hash
{
    size_t operator()(const CC_vertexIndex &e) const
    {
        size_t h0 = std::hash<unsigned int>()(e.areaIndex);
        size_t h1 = std::hash<unsigned int>()(e.i[0]);
        size_t h2 = std::hash<unsigned int>()(e.i[1]);
        size_t h3 = std::hash<unsigned int>()(e.i[2]);
        return h0 ^ ((h1 ^ ((h2 ^ (h3 << 1)) << 1)) << 1);
    }
};
*/
// контейнер для хранения узлов сетки с доступом по индексу узла
// данные узла - индекс в массиве узлов
//typedef std::unordered_map<CC_vertexIndex, int, CC_vertexIndex_hash, CC_vertexIndex_equal> CC_vertexMap;
struct SpherePointKey
{
    #define LOWER(a,b){if(a < b) return true;if(a > b) return false;}
    int r;
    int x[3];
    friend inline bool operator<(const SpherePointKey &l, const SpherePointKey &r)
    {
        LOWER(l.r, r.r);
        LOWER(l.x[0], r.x[0]);
        LOWER(l.x[1], r.x[1]);
        LOWER(l.x[2], r.x[2]);
        return false;
    }
};
typedef std::map<SpherePointKey, int> PointsMap;
void Grid3D::genSphere(const SphereParameters &gp)
{
    //FILE *out = fopen("R.txt", "w");
    // индексы вершин (r,z,y,x)->ind
    PointsMap vind;
    vind.clear();
    // построение куба, вписанного во внутреннюю сферу
    POINT3 cube[8]; // вершины куба, вписанного в сферу
    POINT3_int cube0[8]; // индексация вершин на сфере через вершины на кубе
    double r0 = 1;
    double A = r0/sqrt(3.); // половина стороны куба, вписанного во внутреннюю сферу
    int count = 0;
    for(int i2 = -1; i2 <= 1; i2 += 2)              // z
        for(int i1 = -1; i1 <= 1; i1 += 2)          // y
            for(int i0 = -1; i0 <= 1; i0 += 2)      // x
            {
                cube[count] = POINT3(i0*A, i1*A, i2*A);
                int x0 = (i0+1)/2 * gp.N;
                int x1 = (i1+1)/2 * gp.N;
                int x2 = (i2+1)/2 * gp.N;
                cube0[count] = POINT3_int(x0, x1, x2);
                //printf("%lf %lf %lf\n",cube[count][0],cube[count][1],cube[count][2]);
                //printf("%lf\n",cube[count][0]*cube[count][0]);
                count++;
            }
    int indexes[6][4] =
    {
        {0,2,1,3},
        {1,3,5,7},
        {1,5,0,4},
        {0,4,2,6},
        {2,6,3,7},
        {4,5,6,7},        //{6,7,4,5},
    };

    fe.clear();
    vertex.clear();
    vertexForCurvature.clear();
    bc1.clear();
    FESurface.clear();

    // построение вершин и первых краевых условий

    for(int t = 0; t <= gp.Nparts; t++) // тиражирование
    for(int k = 0; k < 6; k++)  // k - номер грани
    {
        double angle;
        if(gp.buildingMethod == 0)
        {
            // подготовка для первого способа
            angle = cube[indexes[k][0]]^cube[indexes[k][1]]; // угол, под которым видно ребро вписанного куба
        }
        POINT3_int dpc0 = (cube0[indexes[k][2]]-cube0[indexes[k][0]])/gp.N;   // прирост pc0 - внутренний цикл
        POINT3_int dpc1 = (cube0[indexes[k][1]]-cube0[indexes[k][0]])/gp.N;   // прирост pc1 - внешний цикл
        POINT3_int pc1 = cube0[indexes[k][0]];    // точка на кубе для внешнего цикла
        for(int i = 0; i <= gp.N; i++)
        {
            // первый способ
            POINT3 arcPoint[2];  // arcPoint[0]-arcPoint[1], arcPoint[2]-arcPoint[3] - дуги, между которыми строятся четырехугольники
            double ArcAngle;
            if(gp.buildingMethod == 0)
            {
                Operations::sphereArcPoint(r0, cube[indexes[k][0]], cube[indexes[k][1]], angle*i/gp.N, arcPoint[0]);
                Operations::sphereArcPoint(r0, cube[indexes[k][2]], cube[indexes[k][3]], angle*i/gp.N, arcPoint[1]);
                ArcAngle = arcPoint[0]^arcPoint[1]; // уголы, под которыми видна дуга arcPoint[0]-arcPoint[1]
            }

            // второй способ
            // arcPoint[0]-arcPoint[1], arcPoint[2]-arcPoint[3] - рёбра вписанного куба
            if(gp.buildingMethod == 1)
            {
                arcPoint[0] = cube[indexes[k][0]] + (cube[indexes[k][1]] - cube[indexes[k][0]])*i/gp.N;
                arcPoint[1] = cube[indexes[k][2]] + (cube[indexes[k][3]] - cube[indexes[k][2]])*i/gp.N;
            }
            POINT3_int pc0 = pc1;
            for(int j = 0; j <= gp.N; j++)
            {
                if(pc0[0] <= gp.N/2 &&
                   pc0[1] <= gp.N/2 &&
                   pc0[2] <= gp.N/2)
                {
                    POINT3 p0;    // вершина на единичной сфере

                    // первый способ
                    if(gp.buildingMethod == 0)
                    {
                        Operations::sphereArcPoint(r0, arcPoint[0], arcPoint[1], ArcAngle*j/gp.N, p0);
                    }

                    // второй способ
                    if(gp.buildingMethod == 1)
                    {
                        p0 = arcPoint[0] + (arcPoint[1] - arcPoint[0])*j/gp.N;
                        p0 = p0 / p0.abs();
                    }
                    // заведомо нулевые координаты
                    if(pc0[0] == gp.N/2) p0[0] = 0;
                    if(pc0[1] == gp.N/2) p0[1] = 0;
                    if(pc0[2] == gp.N/2) p0[2] = 0;
                    {
                        double R;
                        Operations::findPointOnTheLine_1d(gp.r1, gp.r2, gp.Nparts, gp.q, t, R);
                        //double R = gp.r1 + (gp.r2 - gp.r1)*t/(gp.Nparts);
                        POINT3 p = p0*R;    // вершина на сфере радиуса R
                        SpherePointKey pkey = {t,{(int)pc0[0],(int)pc0[1],(int)pc0[2]}};
                        auto it = vind.find(pkey);
                        if (it == vind.end())   // новая точка
                        {
                            vertex.push_back(p);
                            size_t p_ind = vertex.size() - 1;  // индекс вершины p
                            BoundaryCondition1 bc1_el;
                            vind.insert(PointsMap::value_type(pkey, (int)p_ind));
                            // первые краевые условия
                            if(pc0[0] == gp.N/2)
                            {
                                // x = 0
                                vertex[p_ind][0] = 0;
                                bc1_el.bc1SourceIndex = 0;
                                bc1_el.vertexIndex = (int)p_ind;
                                bc1.push_back(bc1_el);
                            }
                            if(pc0[1] == gp.N/2)
                            {
                                // y = 0
                                vertex[p_ind][1] = 0;
                                bc1_el.bc1SourceIndex = 1;
                                bc1_el.vertexIndex = (int)p_ind;
                                bc1.push_back(bc1_el);
                            }
                            if(pc0[2] == gp.N/2)
                            {
                                // z = 0
                                vertex[p_ind][2] = 0;
                                bc1_el.bc1SourceIndex = 2;
                                bc1_el.vertexIndex = (int)p_ind;
                                bc1.push_back(bc1_el);
                            }
                            if(t == 0)
                            {
                                // вершина на внутренней границе
                                bc1_el.bc1SourceIndex = 3;
                                bc1_el.vertexIndex = (int)p_ind;
                                bc1.push_back(bc1_el);
                            }
                            if(t == gp.Nparts)
                            {
                                // вершина на внешней границе
                                bc1_el.bc1SourceIndex = 4;
                                bc1_el.vertexIndex = (int)p_ind;
                                bc1.push_back(bc1_el);
                            }
                            //printf("(%d; %.0lf, %.0lf, %.0lf)\n", t, pc0[0],pc0[1],pc0[2]);
                        }
                    }
                }
                pc0 += dpc0;
            }
            pc1 += dpc1;
        }
    }


    // построение шестигранников и вторых краевых условий
    // будет 2 поверхности
    FESurface.resize(2);
    for (size_t FEsurfaceInd = 0; FEsurfaceInd < FESurface.size(); FEsurfaceInd++)
        FESurface[FEsurfaceInd].face.clear();

    for(int t = 0; t < gp.Nparts; t++) // тиражирование
    for(int k = 0; k < 6; k++)  // k - номер грани
    {
        POINT3_int dpc0 = (cube0[indexes[k][2]]-cube0[indexes[k][0]])/gp.N;   // прирост pc0 - внутренний цикл
        POINT3_int dpc1 = (cube0[indexes[k][1]]-cube0[indexes[k][0]])/gp.N;   // прирост pc1 - внешний цикл
        POINT3_int pc1 = cube0[indexes[k][0]];    // точка на кубе для внешнего цикла
        for(int i = 0; i < gp.N; i++)
        {
            POINT3_int pc0 = pc1;
            for(int j = 0; j < gp.N; j++)
            {
                POINT3_int c[4];
                c[0] = pc0;
                c[1] = pc0 + dpc1;
                c[2] = pc0 + dpc0;
                c[3] = pc0 + dpc0 + dpc1;
                // максимальные "координаты" вершин 4-угольника
                POINT3_int cMax = c[0];
                for(int k = 0; k < 4; k++)
                    for(int kk = 0; kk < 3; kk++)
                        if(cMax[kk] < c[k][kk]) cMax[kk] = c[k][kk];
                if(cMax[0] <= gp.N/2 &&
                   cMax[1] <= gp.N/2 &&
                   cMax[2] <= gp.N/2)
                {
                    {
                        // добавление шестигранника
                        //if(t >= 0)
                        // квадратичное отображение
                        if(gp.curvilinear == 1)
                        {
                            FE_QuadraticHexagon *fe_el = new FE_QuadraticHexagon;
                            fe_el->mi = 0;
                            FE_LinearHexagon fe_el_t; // временный шестигранник из 8 вершин будет хранить индексы 8-ми вершин
                            for(int n = 0; n < 8; n++)
                            {
                                int nn, tt;
                                if(n < 4)
                                {
                                    nn = n;
                                    tt = t;
                                }
                                else
                                {
                                    nn = n - 4;
                                    tt = t + 1;
                                }
                                SpherePointKey pkey = {tt,{(int)c[nn][0],(int)c[nn][1],(int)c[nn][2]}};
                                auto it = vind.find(pkey);
                                fe_el_t.vi[n] = (*it).second;
                            }
                            // искомые 27 точек на искривлённом шестиграннике
                            int vi_27[27];
                            // 8 вершин
                            vi_27[0] = fe_el_t.vi[0];
                            vi_27[2] = fe_el_t.vi[1];
                            vi_27[6] = fe_el_t.vi[2];
                            vi_27[8] = fe_el_t.vi[3];
                            vi_27[18] = fe_el_t.vi[4];
                            vi_27[20] = fe_el_t.vi[5];
                            vi_27[24] = fe_el_t.vi[6];
                            vi_27[26] = fe_el_t.vi[7];
                            // остальные точки записываются в отдельный массив
                            POINT3 p;
                            double R1 = vertex[fe_el_t.vi[0]].abs();
                            double R2 = vertex[fe_el_t.vi[4]].abs();
                            double Rc = (R1 + R2) / 2;
                            // центр шестигранника
                            p = VECTOR3_NULL;
                            for(int n = 0; n < 8; n++)
                                p += vertex[fe_el_t.vi[n]];
                            p /= 8;
                            p = p*(Rc/p.abs());
                            vi_27[13] = (int)vertexForCurvature.size();
                            vertexForCurvature.push_back(p);
                            // 6 центров граней
                            {
                                int gvi[6] =     {12,        14,        10,        16,        4,         22,};
                                int g4vi[6][4] = {{0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, {0,1,2,3}, {4,5,6,7},};
                                for(int ind = 0; ind < 6; ind++)
                                {
                                    p = VECTOR3_NULL;
                                    for(int n = 0; n < 4; n++)
                                        p += vertex[fe_el_t.vi[g4vi[ind][n]]];
                                    p /= 4;
                                    if(gvi[ind] == 4)
                                        p = p*(R1/p.abs());
                                    if(gvi[ind] == 22)
                                        p = p*(R2/p.abs());
                                    if(gvi[ind] == 10 ||
                                       gvi[ind] == 12 ||
                                       gvi[ind] == 14 ||
                                       gvi[ind] == 16)
                                        p = p*(Rc/p.abs());
                                    vi_27[gvi[ind]] = (int)vertexForCurvature.size();
                                    vertexForCurvature.push_back(p);
                                }
                            }
                            // 12 середин рёбер
                            {
                                int gvi[12] =     {1,     19,    9,     11,    7,     25,    15,    17,    3,     21,    5,     23,};
                                int g4vi[12][2] = {{0,1}, {4,5}, {0,4}, {1,5}, {2,3}, {6,7}, {2,6}, {3,7}, {0,2}, {4,6}, {1,3}, {5,7},};
                                for(int ind = 0; ind < 12; ind++)
                                {
                                    p = VECTOR3_NULL;
                                    for(int n = 0; n < 2; n++)
                                        p += vertex[fe_el_t.vi[g4vi[ind][n]]];
                                    p /= 2;
                                    if(gvi[ind] == 1 ||
                                       gvi[ind] == 3 ||
                                       gvi[ind] == 5 ||
                                       gvi[ind] == 7)
                                        p = p*(R1/p.abs());
                                    if(gvi[ind] == 19 ||
                                       gvi[ind] == 21 ||
                                       gvi[ind] == 23 ||
                                       gvi[ind] == 25)
                                        p = p*(R2/p.abs());
                                    if(gvi[ind] == 9 ||
                                       gvi[ind] == 11 ||
                                       gvi[ind] == 15 ||
                                       gvi[ind] == 17)
                                        p = p*(Rc/p.abs());
                                    vi_27[gvi[ind]] = (int)vertexForCurvature.size();
                                    vertexForCurvature.push_back(p);
                                }
                            }
                            fe_el->setGeomVertexIndexes(vi_27);
                            fe.push_back(fe_el);
                        }
                        // линейное отображение
                        if(gp.curvilinear == 0)
                        {
                            FE_LinearHexagon *fe_el = new FE_LinearHexagon;
                            fe_el->mi = 0;
                            for(int n = 0; n < 8; n++)
                            {
                                int nn, tt;
                                if(n < 4)
                                {
                                    nn = n;
                                    tt = t;
                                }
                                else
                                {
                                    nn = n - 4;
                                    tt = t + 1;
                                }
                                SpherePointKey pkey = {tt,{(int)c[nn][0],(int)c[nn][1],(int)c[nn][2]}};
                                auto it = vind.find(pkey);
                                fe_el->vi[n] = (*it).second;
                            }
                            fe.push_back(fe_el);
                        }
                        // вторые краевые условия
                        FEFace face;
                        if(t == 0)
                        {
                            // r=a
                            face.feIndex = (int)fe.size() - 1;
                            face.faceIndex = 0;   // Z = -1
                            FESurface[0].face.push_back(face);
                        }
                        if(t == gp.Nparts - 1)
                        {
                            // r=b
                            face.feIndex = (int)fe.size() - 1;
                            face.faceIndex = 1;   // Z = +1
                            FESurface[1].face.push_back(face);
                        }
                    }
                }
                pc0 += dpc0;
            }
            pc1 += dpc1;
        }
    }
    //fclose(out);

    {
        //std::vector<int> contactFESurfaceIndex;
        //contactFESurfaceIndex.push_back(0); // КЭ поверхность с индексом 0 - контактсная
        //cuthillMcKee(contactFESurfaceIndex);
        //nestedDissection(contactFESurfaceIndex);
    }
    // нумерация базисных функций
    DOFs = new Grid::GlobalDOFs;
    DOFs->init(vertex.size());
    // ф-и первого порядка
    for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
    {
        Grid::FuncID funcID;
        funcID.type = Grid::FuncType::linear;
        funcID.crackIndex = -1;
        DOFs->addDOFs(vertexIndex, funcID);
    }
}

void CrackSurface_Analitical_vertical_x0::init(const double &set_x0, const double &set_y1, const double &set_y2, const double set_n)
{
    x0 = set_x0;
    y1 = set_y1;
    y2 = set_y2;
    n = set_n;
}
double CrackSurface_Analitical_vertical_x0::ls_tan(const POINT3 &p) const
{
    if(p[1] > (y1 + y2) / 2)
        return p[1] - y2;
    else
        return -p[1] + y1;
}
double CrackSurface_Analitical_vertical_x0::ls_normal(const POINT3 &p) const
{
    return p[0] - x0;
}
/*
void CrackSurface_Analitical_vertical_x0::solveFuncValues_t(const POINT3 &p, const std::vector<FuncID> &funcsID, SubVertexData &subVertexData_el) const
{
    for(size_t funcsIDIndex = 0; funcsIDIndex < funcsID.size(); funcsIDIndex++)
    {
        FuncType type = funcsID[funcsIDIndex].type;
        double func_value;
        double t = ls_tan(p);
        double n = ls_normal(p);
        double r;
        double a;
        // вычисление угла
        {
            //return atan(n/t);
            if(t == 0)
                a = -PI/2 * SIGN(n);
            else
            {
                a = atan(n/t);
                if(t < 0)
                {
                    if(a < 0)
                        a += PI;
                    else
                        a -= PI;
                }
            }
        }
        // вычисление расстояния до фронта трещины
        {
            r = sqrt(t*t + n*n);
        }
        // вычисление функций
        switch (type)
        {
        // Лагранжа
        case FuncType::linear:
        {
            func_value = 1;
            break;
        }
            // ф-и для сильного разрыва (Хевисайда)
        case FuncType::H:
        {
            func_value = SIGN(n);
            break;
        }
            // ф-и для вершины трещины (асимптотическое решение Вильямса)
        case FuncType::tip1:
        {
            func_value = sqrt(r)*sin(a/2);
            break;
        }
        case FuncType::tip2:
        {
            func_value = sqrt(r)*sin(a/2)*sin(a);
            break;
        }
        case FuncType::tip3:
        {
            func_value = sqrt(r)*cos(a/2);
            break;
        }
        case FuncType::tip4:
        {
            func_value = sqrt(r)*cos(a/2)*sin(a);
            break;
        }
        }
        subVertexData_el.funcValue[funcsIDIndex] = func_value;
    }
}
*/
void CrackSurface_Analitical_vertical_x0::solveFuncValues(const POINT3 &p, const std::vector<FuncID> &funcsID, const size_t crackIndex, SubVertexData &subVertexData_el) const
{
    for(size_t funcsIDIndex = 0; funcsIDIndex < funcsID.size(); funcsIDIndex++)
    {
        if(funcsID[funcsIDIndex].crackIndex == crackIndex)
        {
            FuncType type = funcsID[funcsIDIndex].type;
            double func_value;
            double t = ls_tan(p);
            double n = ls_normal(p);
            double r;
            double a;
            // вычисление угла
            {
                //return atan(n/t);
                if(t == 0)
                    a = -PI/2 * SIGN(n);
                else
                {
                    a = atan(n/t);
                    if(t < 0)
                    {
                        if(a < 0)
                            a += PI;
                        else
                            a -= PI;
                    }
                }
            }
            // вычисление расстояния до фронта трещины
            {
                r = sqrt(t*t + n*n);
            }
            // вычисление функций
            switch (type)
            {
            // Лагранжа
            case FuncType::linear:
            {
                func_value = 1;
                break;
            }
                // ф-и для сильного разрыва (Хевисайда)
            case FuncType::H:
            {
                func_value = SIGN(n);
                break;
            }
                // ф-и для вершины трещины (асимптотическое решение Вильямса)
            case FuncType::tip1:
            {
                func_value = sqrt(r)*sin(a/2);
                break;
            }
            case FuncType::tip2:
            {
                func_value = sqrt(r)*sin(a/2)*sin(a);
                break;
            }
            case FuncType::tip3:
            {
                func_value = sqrt(r)*cos(a/2);
                break;
            }
            case FuncType::tip4:
            {
                func_value = sqrt(r)*cos(a/2)*sin(a);
                break;
            }
                // упруго-пластичность (3.20b)
            case FuncType::tip_ep1:
            {
                func_value = pow(r, 1./(n + 1))*sin(a/2);
                break;
            }
            case FuncType::tip_ep2:
            {
                func_value = pow(r, 1./(n + 1))*cos(a/2);
                break;
            }
            case FuncType::tip_ep3:
            {
                func_value = pow(r, 1./(n + 1))*sin(a/2)*sin(a);
                break;
            }
            case FuncType::tip_ep4:
            {
                func_value = pow(r, 1./(n + 1))*cos(a/2)*sin(a);
                break;
            }
            case FuncType::tip_ep5:
            {
                func_value = pow(r, 1./(n + 1))*sin(a/2)*sin(3*a);
                break;
            }
            case FuncType::tip_ep6:
            {
                func_value = pow(r, 1./(n + 1))*cos(a/2)*sin(3*a);
                break;
            }
                // ## ошибка
            case FuncType::_SIZE:
            {
                func_value = 0;
                break;
            }
            }
            subVertexData_el.funcValue[funcsIDIndex] = func_value;
        }
    }
}

double Gauss3_integrationwSource[27];
//double Gauss3_basCubeSource[27*8];
VECTOR3 Gauss3_dLinearBasCubeSource[27*8];


void FE_LinearHexagon_XFEM::gn_initSubVertexes(const std::vector<CrackSurface_base *> &cs, const std::vector<POINT3> &vertex, const std::vector<POINT3> &vertexForCurvature, const std::vector<POINT3_CUBE> &subVertex_XYZ)
{
    // глобальные координаты вершин 6-гранного КЭ
    POINT3 v[8];
    getGeomVertexes(vertex, vertexForCurvature, v);
    // глобальные координаты вершин шестигранных подобластей
    for(size_t sub_vi = 0; sub_vi < subVertex.size(); sub_vi++) // индекс вершины подобласти
    {
        // локальные координаты
        POINT3_CUBE XYZ = subVertex_XYZ[sub_vi];
        // глобальные координаты
        POINT3 xyz;
        Fem::cubeToLagrange1Hexagon(v, XYZ, dif_NULL3, xyz);
        subVertex[sub_vi] = xyz;
        dsubVertex[sub_vi] = POINT3_CUBE(0, 0, 0);
    }
    // значения ф-й формы Лагранжа в узлах подобластей
    for(size_t sub_vi = 0; sub_vi < subVertex.size(); sub_vi++) // индекс вершины подобласти
    {
        // локальные координаты
        POINT3_CUBE XYZ = subVertex_XYZ[sub_vi];
        for(size_t m = 0; m < 8; m++) // индекс узла
        {
            double lagr_func_value = Fem::cubeLagrange1_3D(XYZ, m, dif_NULL3);
            subVertexData[sub_vi].lagrFuncsValue[m] = lagr_func_value;
        }
    }
    // значения базисных ф-й в узлах подобластей
    /*
    for(size_t sub_vi = 0; sub_vi < subVertex.size(); sub_vi++) // индекс вершины подобласти
    {
        POINT3 xyz = subVertex[sub_vi];
        cs[crackIndex]->solveFuncValues(xyz, funcsID,
                                        subVertexData[sub_vi]);
    }
    */
    for(size_t sub_vi = 0; sub_vi < subVertex.size(); sub_vi++) // индекс вершины подобласти
    {
        POINT3 xyz = subVertex[sub_vi];
        // ф-я Лагранжа
        subVertexData[sub_vi].funcValue[0] = 1;
        // остальные ф-и
        bool crackIndex_used[1024];
        size_t crackIndex_max = 0;
        for(size_t funcsIDIndex = 1; funcsIDIndex < funcsID.size(); funcsIDIndex++)
        {
            if(funcsID[funcsIDIndex].crackIndex + 1 > crackIndex_max)
            {
                crackIndex_max = funcsID[funcsIDIndex].crackIndex + 1;
            }
        }
        for(size_t crackIndex = 0; crackIndex < crackIndex_max; crackIndex++)
        {
            crackIndex_used[crackIndex] = false;
        }
        for(size_t funcsIDIndex = 1; funcsIDIndex < funcsID.size(); funcsIDIndex++)
        {
            crackIndex_used[funcsID[funcsIDIndex].crackIndex] = true;
        }
        for(size_t crackIndex = 0; crackIndex < crackIndex_max; crackIndex++)
        {
            if(crackIndex_used[crackIndex])
            {
                cs[crackIndex]->solveFuncValues(xyz, funcsID, crackIndex,
                                                subVertexData[sub_vi]);
            }
        }
    }
}
void FE_LinearHexagon_XFEM::gn_getSubVertexes(const size_t hexInd, POINT3 (&subv)[8]) const
{
    for(size_t i = 0; i < 8; i++) // индекс вершины 6-гранной подобласти
    {
        subv[i] = subVertex[subHex[hexInd].vi[i]];
    }
}
void FE_LinearHexagon_XFEM::gn_getSubFaceVertexes(const size_t hexInd, const int faceInd, POINT3 (&v)[4]) const
{
    for(size_t i = 0; i < 4; i++)
    {
        v[i] = subVertex[subHex[hexInd].vi[HexagonFace[faceInd][i]]];
    }
}

void FE_LinearHexagon_XFEM::gn_calcBc2BasisFuncValues(const std::vector<LocalDOF_ID> &DOF_ID, const Integration::Integrator &integrationFoursquare, const size_t hexInd, const int faceIndex, double *basCube, double *hexagonCoef, VECTOR3 *hexagonNormal) const
{
    // копируем вершины шестигранника
    POINT3 subv[8];
    gn_getSubVertexes(hexInd, subv);
    POINT3_CUBE cubePoint;
    POINT3 dxyz[2];
    int c1 = HexagonFaceCoordinates[faceIndex][0]; // первый индекс свободной координаты
    int c2 = HexagonFaceCoordinates[faceIndex][1]; // второй индекс свободной координаты
    int cc = HexagonFaceCoordinates[faceIndex][2]; // индекс зафиксированной координаты
    int ccValue = HexagonFaceCoordinates[faceIndex][3];  // значение зафиксированной координаты

    for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
    {
        cubePoint[c1] = integrationFoursquare.p2[valIndex][0];
        cubePoint[c2] = integrationFoursquare.p2[valIndex][1];
        cubePoint[cc] = ccValue;
        // линейное либо квадратичное отображение
        cubeToHexagon(subv, cubePoint, Fem::dif_XYZ3[c1], dxyz[0]);
        cubeToHexagon(subv, cubePoint, Fem::dif_XYZ3[c2], dxyz[1]);
        // нормаль
        hexagonNormal[valIndex] = vector3Mul(dxyz[0], dxyz[1]);
        hexagonNormal[valIndex] = ccValue*hexagonNormal[valIndex]/hexagonNormal[valIndex].abs();
        hexagonCoef[valIndex] =
                integrationFoursquare.w[valIndex] * integrationFoursquare.detJ *
               sqrt(
               (SQR(dxyz[0][0]) +
                SQR(dxyz[0][1]) +
                SQR(dxyz[0][2]))
                *
               (SQR(dxyz[1][0]) +
                SQR(dxyz[1][1]) +
                SQR(dxyz[1][2]))
                -
               SQR(dxyz[0][0]*dxyz[1][0] +
                   dxyz[0][1]*dxyz[1][1] +
                   dxyz[0][2]*dxyz[1][2]));
        double bas_t[FE_LinearHexagon_XFEM::DOF_max];
        gn_calcSubBasFuncs(DOF_ID, hexInd, cubePoint,
                           bas_t);
        for(size_t DOF_ID_index = 0; DOF_ID_index < DOF_ID.size(); DOF_ID_index++)
        {
            basCube[DOF_ID_index*integrationFoursquare.size + valIndex] = bas_t[DOF_ID_index];
        }
        //for (int m = 0; m < 8; m++)
        //    basCube[m*integrationFoursquare.size + valIndex] = Fem::cubeLagrange1_3D(cubePoint, m, Fem::dif_NULL3);
    }
}
void FE_LinearHexagon_XFEM::gn_calcIntegralsTable(const double (&Gauss3_integrationwSource)[27], const VECTOR3 (&Gauss3_dLinearBasCubeSource)[27*8], const std::vector<LocalDOF_ID> &DOF_ID, const size_t hexInd, HexagonIntTable &hexagonIntTable) const
{
    // расчёт искомых множителей для интегрирования
    double w[27];
    double dbas[FE_LinearHexagon_XFEM::DOF_max][27][3];   // значения производных базисных функций в каждой точке интегрирования 6-гранной подобласти hexInd
    // значения производных базисных функций, детерминантов и коэффициентов численного интегрирования
    for (size_t valIndex = 0; valIndex < 27; valIndex++) // индекс точки интегрирования
    {
        // якобиан отображения
        VECTOR3 J[3] = {};
        for (size_t i = 0; i < 8; i++)     // индекс вершины 6-гранной подобласти
        {
            POINT3 vertex_i = subVertex[subHex[hexInd].vi[i]];
            J[0] += Gauss3_dLinearBasCubeSource[valIndex*8 + i][0] * vertex_i; // производная по X
            J[1] += Gauss3_dLinearBasCubeSource[valIndex*8 + i][1] * vertex_i; // производная по Y
            J[2] += Gauss3_dLinearBasCubeSource[valIndex*8 + i][2] * vertex_i; // производная по Z
        }
        // детерминант
        double det = calcDet3x3(
                    J[0][0], J[0][1], J[0][2],
                    J[1][0], J[1][1], J[1][2],
                    J[2][0], J[2][1], J[2][2]);
        // множитель для интегрирования
        w[valIndex] = Gauss3_integrationwSource[valIndex]*fabs(det);
        // производные базисных ф-й по x, y, z
        for(size_t DOF_ID_index = 0; DOF_ID_index < DOF_ID.size(); DOF_ID_index++)
        {
            const LocalDOF_ID &DOF_ID_el = DOF_ID[DOF_ID_index];
            const FuncID &funcID_el = funcsID[DOF_ID_el.funcID_index];
            VECTOR3 dBasCube = VECTOR3_NULL;

            for (int i = 0; i < 8; i++) // индекс вершины 6-гранной подобласти
            {
                const SubVertexData &vd = subVertexData[subHex[hexInd].vi[i]];
                const SubVertexData &FEvd = subVertexData[FESubVertexIndex[DOF_ID_el.m]];
                double funcValue_i;
                if(funcID_el.type == FuncType::linear)
                {
                    funcValue_i = vd.lagrFuncsValue[DOF_ID_el.m];
                    //fprintf(stderr, "funcValue_i[%d][%d][%d] = %le\n", (int)valIndex, (int)DOF_ID_index, (int)i, funcValue_i);
                }
                else
                {
                    funcValue_i = (vd.funcValue[DOF_ID_el.funcID_index] - FEvd.funcValue[DOF_ID_el.funcID_index])
                            *vd.lagrFuncsValue[DOF_ID_el.m];   // значение ф-и funcsID[DOF_ID_el.funcID_index] в узле подобласти i (которая привязана к узлу КЭ DOF_ID_el.m]), домноженной на ф-ю Лагранжа DOF_ID_el.m
                }
                dBasCube[0] += Gauss3_dLinearBasCubeSource[valIndex*8 + i][0] * funcValue_i; // производная по X
                dBasCube[1] += Gauss3_dLinearBasCubeSource[valIndex*8 + i][1] * funcValue_i; // производная по Y
                dBasCube[2] += Gauss3_dLinearBasCubeSource[valIndex*8 + i][2] * funcValue_i; // производная по Z
            }
            // dbasm/dx
            dbas[DOF_ID_index][valIndex][0] = calcDet3x3(
                    dBasCube[0], J[0][1], J[0][2],
                    dBasCube[1], J[1][1], J[1][2],
                    dBasCube[2], J[2][1], J[2][2]) / det;
            // dbasm/dy
            dbas[DOF_ID_index][valIndex][1] = calcDet3x3(
                    J[0][0], dBasCube[0], J[0][2],
                    J[1][0], dBasCube[1], J[1][2],
                    J[2][0], dBasCube[2], J[2][2]) / det;
            // dbasm/dz
            dbas[DOF_ID_index][valIndex][2] = calcDet3x3(
                    J[0][0], J[0][1], dBasCube[0],
                    J[1][0], J[1][1], dBasCube[1],
                    J[2][0], J[2][1], dBasCube[2]) / det;
        }
    }
    for (size_t m = 0; m < DOF_ID.size(); m++)
        for (size_t n = 0; n < DOF_ID.size(); n++) // m, n - локальные номера базисных функций
            for (size_t j = 0; j < 3; j++)
                for (size_t l = 0; l < 3; l++)
                {
                    double E = 0;
                    for (size_t valIndex = 0; valIndex < 27; valIndex++)
                    {
                        E += w[valIndex] *
                             dbas[m][valIndex][j] * dbas[n][valIndex][l];
                    }
                    hexagonIntTable.intForG_mnjl[m][n][j][l] = E;
                }
    for (size_t m = 0; m < DOF_ID.size(); m++)
        for (size_t j = 0; j < 3; j++)
        {
            double E = 0;
            for (size_t valIndex = 0; valIndex < 27; valIndex++)
                E += w[valIndex] *
                     dbas[m][valIndex][j];
            hexagonIntTable.intForb23_mj[m][j] = E;
        }
}
void FE_LinearHexagon_XFEM::gn_moveSubVertexes(const std::vector<LocalDOF_ID> &DOF_ID, const VECTOR3 (&q_local)[FE_LinearHexagon_XFEM::DOF_max])
{
    // #значения базисных ф-й в узлах подобластей известны, их можно не вычислять
    double bas[FE_LinearHexagon_XFEM::DOF_max];
    bool subVertex_moved[FE_LinearHexagon_XFEM::Vertexes_max];
    for(size_t subVertexInd = 0; subVertexInd < subVertex.size(); subVertexInd++)
    {
        subVertex_moved[subVertexInd] = false;
    }
    /*
    for(size_t m = 0; m < DOF_ID.size(); m++)    // m - индекс степени свободы
    {
        fprintf(stderr, "q_local[%d]=(%le, %le, %le) %d\n", (int)m, q_local[m][0], q_local[m][1], q_local[m][2], (int)DOF_ID[m].funcID_index);
    }
    fprintf(stderr, "\n");
    */
    for(size_t hexInd = 0; hexInd < subHex.size(); hexInd++) // индекс 6-гранной подобласти КЭ
    {
        //fprintf(stderr, "hexInd = %d:\n", (int)hexInd);
        for(size_t i = 0; i < 8; i++)  // индекс вершины 6-гранной подобласти hexInd
        {
            if(subVertex_moved[subHex[hexInd].vi[i]] == false)
            {
                POINT3_CUBE XYZ = CUBE_VERTEX[i];    // координаты вершины шаблонного куба (-1 или +1 каждая координата)
                // в вершине XYZ вычисляем сумму значений базисных функций, домноженных на коэффициенты
                gn_calcSubBasFuncs(DOF_ID, hexInd, XYZ,
                                   bas);
                POINT3 U(0, 0, 0);
                for(size_t m = 0; m < DOF_ID.size(); m++)    // m - индекс степени свободы
                {
                    U += bas[m] * q_local[m];
                    //fprintf(stderr, "bas[%d]=%le\n", (int)m, bas[m]);
                }
                //fprintf(stderr, "\n");
                subVertex[subHex[hexInd].vi[i]] += U;
                dsubVertex[subHex[hexInd].vi[i]] += U;
                subVertex_moved[subHex[hexInd].vi[i]] = true;
                //fprintf(stderr, "U[%d]=(%le, %le, %le)\n", (int)subHex[hexInd].vi[i], U[0], U[1], U[2]);
            }
        }
    }

}

void FE_LinearHexagon_XFEM::gn_calcSubBasFuncs(const std::vector<LocalDOF_ID> &DOF_ID, const size_t hexInd, const POINT3_CUBE &cubePoint,
                                               double (&bas)[FE_LinearHexagon_XFEM::DOF_max]) const
{
    for(size_t DOF_ID_index = 0; DOF_ID_index < DOF_ID.size(); DOF_ID_index++)
    {
        const LocalDOF_ID &DOF_ID_el = DOF_ID[DOF_ID_index];
        const FuncID &funcID_el = funcsID[DOF_ID_el.funcID_index];
        bas[DOF_ID_index] = 0;
        for (int i = 0; i < 8; i++) // индекс вершины 6-гранной подобласти
        {
            const SubVertexData &vd = subVertexData[subHex[hexInd].vi[i]];
            const SubVertexData &FEvd = subVertexData[FESubVertexIndex[DOF_ID_el.m]];
            double funcValue_i;
            if(funcID_el.type == FuncType::linear)
            {
                funcValue_i = vd.lagrFuncsValue[DOF_ID_el.m];
            }
            else
            {
                funcValue_i = (vd.funcValue[DOF_ID_el.funcID_index] - FEvd.funcValue[DOF_ID_el.funcID_index])
                        *vd.lagrFuncsValue[DOF_ID_el.m];   // значение ф-и funcsID[DOF_ID_el.funcID_index] в узле подобласти i (которая привязана к узлу КЭ DOF_ID_el.m]), домноженной на ф-ю Лагранжа DOF_ID_el.m
            }
            bas[DOF_ID_index] += Fem::cubeLagrange1_3D(cubePoint, i, Fem::dif_NULL3) * funcValue_i;
        }
    }
}
void FE_LinearHexagon_XFEM::gn_calcSubDifBasFuncs(const std::vector<LocalDOF_ID> &DOF_ID, const size_t hexInd, const POINT3_CUBE &cubePoint,
                                                  double (&dbas)[FE_LinearHexagon_XFEM::DOF_max][3]) const
{
    // производные базисных ф-й 1-го порядка по локальным координатам
    VECTOR3 dLinearBasCube[8];
    for (int i = 0; i < 8; i++)
    {
        dLinearBasCube[i][0] = Fem::cubeLagrange1_3D(cubePoint, i, Fem::dif_XYZ3[0]);
        dLinearBasCube[i][1] = Fem::cubeLagrange1_3D(cubePoint, i, Fem::dif_XYZ3[1]);
        dLinearBasCube[i][2] = Fem::cubeLagrange1_3D(cubePoint, i, Fem::dif_XYZ3[2]);
    }
    // якобиан отображения
    VECTOR3 J[3] = {};
    for (size_t i = 0; i < 8; i++)     // индекс вершины 6-гранной подобласти
    {
        POINT3 vertex_i = subVertex[subHex[hexInd].vi[i]];
        J[0] += dLinearBasCube[i][0] * vertex_i; // производная по X
        J[1] += dLinearBasCube[i][1] * vertex_i; // производная по Y
        J[2] += dLinearBasCube[i][2] * vertex_i; // производная по Z
    }
    // детерминант
    double det = calcDet3x3(
                J[0][0], J[0][1], J[0][2],
                J[1][0], J[1][1], J[1][2],
                J[2][0], J[2][1], J[2][2]);
    // производные базисных ф-й по x, y, z
    for(size_t DOF_ID_index = 0; DOF_ID_index < DOF_ID.size(); DOF_ID_index++)
    {
        const LocalDOF_ID &DOF_ID_el = DOF_ID[DOF_ID_index];
        const FuncID &funcID_el = funcsID[DOF_ID_el.funcID_index];
        VECTOR3 dBasCube = VECTOR3_NULL;
        for (int i = 0; i < 8; i++) // индекс вершины 6-гранной подобласти
        {
            const SubVertexData &vd = subVertexData[subHex[hexInd].vi[i]];
            const SubVertexData &FEvd = subVertexData[FESubVertexIndex[DOF_ID_el.m]];
            double funcValue_i;
            if(funcID_el.type == FuncType::linear)
            {
                funcValue_i = vd.lagrFuncsValue[DOF_ID_el.m];
            }
            else
            {
                funcValue_i = (vd.funcValue[DOF_ID_el.funcID_index] - FEvd.funcValue[DOF_ID_el.funcID_index])
                        *vd.lagrFuncsValue[DOF_ID_el.m];   // значение ф-и funcsID[DOF_ID_el.funcID_index] в узле подобласти i (которая привязана к узлу КЭ DOF_ID_el.m]), домноженной на ф-ю Лагранжа DOF_ID_el.m
            }
            dBasCube[0] += dLinearBasCube[i][0] * funcValue_i; // производная по X
            dBasCube[1] += dLinearBasCube[i][1] * funcValue_i; // производная по Y
            dBasCube[2] += dLinearBasCube[i][2] * funcValue_i; // производная по Z
        }
        // dbasm/dx
        dbas[DOF_ID_index][0] = calcDet3x3(
                dBasCube[0], J[0][1], J[0][2],
                dBasCube[1], J[1][1], J[1][2],
                dBasCube[2], J[2][1], J[2][2]) / det;
        // dbasm/dy
        dbas[DOF_ID_index][1] = calcDet3x3(
                J[0][0], dBasCube[0], J[0][2],
                J[1][0], dBasCube[1], J[1][2],
                J[2][0], dBasCube[2], J[2][2]) / det;
        // dbasm/dz
        dbas[DOF_ID_index][2] = calcDet3x3(
                J[0][0], J[0][1], dBasCube[0],
                J[1][0], J[1][1], dBasCube[1],
                J[2][0], J[2][1], dBasCube[2]) / det;
    }
}
FEType FE_LinearHexagon_XFEM::get_FEType() const
{
    return FEType::LinearHexagon_XFEM;
}

}   // namespace Grid


