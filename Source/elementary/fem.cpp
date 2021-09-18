#define _CRT_SECURE_NO_WARNINGS

#include "fem.h"

namespace Fem
{

double cubeLagrange1_1D(const POINT1_CUBE cubePoint, const int i)
{
    if (i == 0) return (1 - cubePoint) / 2;
    if (i == 1) return (1 + cubePoint) / 2;
    return 0;
}
double cubeLagrange1_1D(const POINT1_CUBE cubePoint, const int i, const int dif)
{
    if (dif == 0)
    {
        if (i == 0) return (1 - cubePoint) / 2;
        if (i == 1) return (1 + cubePoint) / 2;
        return 0;
    }
    else
        if (dif == 1)
        {
            if (i == 0) return -0.5L;
            if (i == 1) return 0.5L;
            return 0;
        }
        else
            if (dif >= 2)
                return 0;
    return 0;
}
double cubeLagrange2_1D(const POINT1_CUBE cubePoint, const int i)
{
    if (i == 0) return cubePoint*(cubePoint - 1) / 2;
    if (i == 1) return 1 - cubePoint*cubePoint;
    if (i == 2) return cubePoint*(cubePoint + 1) / 2;
    return 0;
}
double cubeLagrange2_1D(const POINT1_CUBE cubePoint, const int i, const int dif)
{
    switch (dif)
    {
    case 0:
        if (i == 0) return cubePoint*(cubePoint - 1) / 2;
        if (i == 1) return 1 - cubePoint*cubePoint;
        if (i == 2) return cubePoint*(cubePoint + 1) / 2;
        return 0;
        break;
    case 1:
        if (i == 0) return cubePoint - 1. / 2.;
        if (i == 1) return -2*cubePoint;
        if (i == 2) return cubePoint + 1. / 2.;
        return 0;
        break;
    case 2:
        if (i == 0) return 1;
        if (i == 1) return -2;
        if (i == 2) return 1;
        return 0;
        break;
    default:
        return 0;
        break;
    }
    return 0;
}
double cubeLagrange1_3D(const POINT3_CUBE &cubePoint, const int i, const DIF_STATE3 &dif)
{
    int k1 = i % 2;
    int k2 = (i / 2) % 2;
    int k3 = i / 4;
    return cubeLagrange1_1D(cubePoint[0], k1, dif[0])*
            cubeLagrange1_1D(cubePoint[1], k2, dif[1])*
            cubeLagrange1_1D(cubePoint[2], k3, dif[2]);
}
double cubeLagrange2_3D(const POINT3_CUBE &cubePoint, const int i, const DIF_STATE3 &dif)
{
    int k1 = i % 3;
    int k2 = (i / 3) % 3;
    int k3 = i / 9;
    return cubeLagrange2_1D(cubePoint[0], k1, dif[0])*
            cubeLagrange2_1D(cubePoint[1], k2, dif[1])*
            cubeLagrange2_1D(cubePoint[2], k3, dif[2]);
}

void cubeToLagrange1Hexagon(const POINT3 *v, const POINT3_CUBE &cubePoint, const DIF_STATE3 &dif,
                         POINT3 &point)
{
    point[0] = 0;
    point[1] = 0;
    point[2] = 0;
    for (int i = 0; i < 8; i++)
        point += cubeLagrange1_3D(cubePoint, i, dif)*v[i];
}
void cubeToLagrange2Hexagon(const POINT3 *v, const POINT3_CUBE &cubePoint, const DIF_STATE3 &dif,
                            POINT3 &point)
{
    point[0] = 0;
    point[1] = 0;
    point[2] = 0;
    for (int i = 0; i < 27; i++)
        point += cubeLagrange2_3D(cubePoint, i, dif)*v[i];
}

void lagrange1Hexagon_mappingJ(const POINT3 *v, const POINT3_CUBE &cubePoint,
                          MATR3x3 &J, double &detJ)
{
    VECTOR3 str;
    // матрица J
    cubeToLagrange1Hexagon(v, cubePoint, dif_XYZ3[0], str);   // производная по X
    J.m[0][0] = str[0];
    J.m[0][1] = str[1];
    J.m[0][2] = str[2];
    cubeToLagrange1Hexagon(v, cubePoint, dif_XYZ3[1], str);   // производная по Y
    J.m[1][0] = str[0];
    J.m[1][1] = str[1];
    J.m[1][2] = str[2];
    cubeToLagrange1Hexagon(v, cubePoint, dif_XYZ3[2], str);   // производная по Z
    J.m[2][0] = str[0];
    J.m[2][1] = str[1];
    J.m[2][2] = str[2];
    // детерминант
    detJ = J.calcDet();
}
void lagrange2Hexagon_mappingJ(const POINT3 *v, const POINT3_CUBE &cubePoint,
                             MATR3x3 &J, double &detJ)
{
    VECTOR3 str;
    // матрица J
    cubeToLagrange2Hexagon(v, cubePoint, dif_XYZ3[0], str);    // производная по X
    J.m[0][0] = str[0];
    J.m[0][1] = str[1];
    J.m[0][2] = str[2];
    cubeToLagrange2Hexagon(v, cubePoint, dif_XYZ3[1], str);    // производная по Y
    J.m[1][0] = str[0];
    J.m[1][1] = str[1];
    J.m[1][2] = str[2];
    cubeToLagrange2Hexagon(v, cubePoint, dif_XYZ3[2], str);    // производная по Z
    J.m[2][0] = str[0];
    J.m[2][1] = str[1];
    J.m[2][2] = str[2];
    // детерминант
    detJ = J.calcDet();
}
double lagrange1Hexagon_mappingDetJ(const POINT3 *v, const POINT3_CUBE &cubePoint)
{
    MATR3x3 J;
    double detJ;
    lagrange1Hexagon_mappingJ(v, cubePoint,
                              J, detJ);
    return detJ;
}
double lagrange2Hexagon_mappingDetJ(const POINT3 *v, const POINT3_CUBE &cubePoint)
{
    MATR3x3 J;
    double detJ;
    lagrange2Hexagon_mappingJ(v, cubePoint,
                                 J, detJ);
    return detJ;
}

double lagrange1Hexagon_difLagrange1(const POINT3 *v, const POINT3_CUBE &cubePoint, const int i,
                              VECTOR3 &difBasisFunction)
{
    VECTOR3 dF;
    MATR3x3 J;
    double detJ;
    lagrange1Hexagon_mappingJ(v, cubePoint,
                              J, detJ);
    dF[0] = cubeLagrange1_3D(cubePoint, i, dif_XYZ3[0]);
    dF[1] = cubeLagrange1_3D(cubePoint, i, dif_XYZ3[1]);
    dF[2] = cubeLagrange1_3D(cubePoint, i, dif_XYZ3[2]);
    // dFi/dx
    difBasisFunction[0] = calcDet3x3(
            dF[0], J.m[0][1], J.m[0][2],
            dF[1], J.m[1][1], J.m[1][2],
            dF[2], J.m[2][1], J.m[2][2]) / detJ;
    // dFi/dy
    difBasisFunction[1] = calcDet3x3(
            J.m[0][0], dF[0], J.m[0][2],
            J.m[1][0], dF[1], J.m[1][2],
            J.m[2][0], dF[2], J.m[2][2]) / detJ;
    // dFi/dz
    difBasisFunction[2] = calcDet3x3(
            J.m[0][0], J.m[0][1], dF[0],
            J.m[1][0], J.m[1][1], dF[1],
            J.m[2][0], J.m[2][1], dF[2]) / detJ;
    return detJ;
}
double lagrange2Hexagon_difLagrange1(const POINT3 *v, const POINT3_CUBE &cubePoint, const int i,
                                 VECTOR3 &difBasisFunction)
{
    VECTOR3 dF;
    MATR3x3 J;
    double detJ;
    lagrange2Hexagon_mappingJ(v, cubePoint,
                                 J, detJ);
    dF[0] = cubeLagrange1_3D(cubePoint, i, dif_XYZ3[0]);
    dF[1] = cubeLagrange1_3D(cubePoint, i, dif_XYZ3[1]);
    dF[2] = cubeLagrange1_3D(cubePoint, i, dif_XYZ3[2]);
    // dFi/dx
    difBasisFunction[0] = calcDet3x3(
            dF[0], J.m[0][1], J.m[0][2],
            dF[1], J.m[1][1], J.m[1][2],
            dF[2], J.m[2][1], J.m[2][2]) / detJ;
    // dFi/dy
    difBasisFunction[1] = calcDet3x3(
            J.m[0][0], dF[0], J.m[0][2],
            J.m[1][0], dF[1], J.m[1][2],
            J.m[2][0], dF[2], J.m[2][2]) / detJ;
    // dFi/dz
    difBasisFunction[2] = calcDet3x3(
            J.m[0][0], J.m[0][1], dF[0],
            J.m[1][0], J.m[1][1], dF[1],
            J.m[2][0], J.m[2][1], dF[2]) / detJ;
    return detJ;
}

double lagrange1_1D(const POINT1 &p0, const double h, const POINT1 p, const int i, const DIF_STATE1 &dif)
{
    if (dif == 0)
    {
        if (i == 0) return (p0 + h - p) / h;
        if (i == 1) return (p - p0) / h;
        return 0;
    }
    if (dif == 1)	// 1-я производная
    {
        if (i == 0) return -1. / h;
        if (i == 1) return 1. / h;
        return 0;
    }
    return 0;
}
double lagrange3_1D(const POINT1 &p0, const double h, const POINT1 p, const int i, const DIF_STATE1 &dif)
{
    double ksi = (p - p0) / h;
    if (dif == 0)
    {
        if (i == 0) return -9./2. *(ksi - 1./3.)*(ksi - 2./3.)*(ksi - 1.);     //d/dx (-9./2.*(x - 1./3.)*(x - 2./3.)*(x - 1.))
        if (i == 1) return 27./2. *ksi          *(ksi - 2./3.)*(ksi - 1.);     //d/dx (27./2.*x*(x - 2./3.)*(x - 1.))
        if (i == 2) return -27./2.*ksi          *(ksi - 1./3.)*(ksi - 1.);     //d/dx (-27./2.*x*(x - 1./3.)*(x - 1.))
        if (i == 3) return 9./2.  *ksi          *(ksi - 1./3.)*(ksi - 2./3.);  //d/dx (9./2.*x*(x - 1./3.)*(x - 2./3.))
        return 0;
    }
    if (dif == 1)	// 1-я производная по a
    {
        double dksi = 1. / h;
        if (i == 0) return (-5.5 + 18*ksi - 13.5*ksi*ksi)*dksi;
        if (i == 1) return (9    - 45*ksi + 40.5*ksi*ksi)*dksi;
        if (i == 2) return (-4.5 + 36*ksi - 40.5*ksi*ksi)*dksi;
        if (i == 3) return (1    - 9*ksi  + 13.5*ksi*ksi)*dksi;
        return 0;
    }
    if (dif == 2)	// 2-я производная по a
    {
        double dksi = 1. / h;
        if (i == 0) return (18  - 27*ksi)*dksi*dksi;
        if (i == 1) return (-45 + 81*ksi)*dksi*dksi;
        if (i == 2) return (36 - 81*ksi)*dksi*dksi;
        if (i == 3) return (-9 + 27*ksi)*dksi*dksi;
        return 0;
    }
    return 0;
}
double hermite_1D(const POINT1 &p0, const double h, const POINT1 p, const int i, const DIF_STATE1 &dif)
{
    double ksi = (p - p0) / h;
    if (dif == 0)
    {
        if (i == 0) return 1 - 3 * ksi*ksi + 2 * ksi*ksi*ksi;
        if (i == 1) return (ksi - 2 * ksi*ksi + ksi*ksi*ksi)*h;	// чтобы производная в левом узле = 1 (функции 1 и 3)
        if (i == 2) return 3 * ksi*ksi - 2 * ksi*ksi*ksi;
        if (i == 3) return (-ksi*ksi + ksi*ksi*ksi)*h;
        return 0;
    }
    if (dif == 1)	// первая производная по a
    {
        //double dksi = 1. / h;
        if (i == 0) return (-6*ksi + 6 * ksi*ksi)/h;    //(-6*ksi + 6 * ksi*ksi)*dksi;
        if (i == 1) return (1 - 4 * ksi + 3*ksi*ksi);   //(1 - 4 * ksi + 3*ksi*ksi)*h*dksi;
        if (i == 2) return (6 * ksi - 6 * ksi*ksi)/h;   //(6 * ksi - 6 * ksi*ksi)*dksi;
        if (i == 3) return (-2*ksi + 3*ksi*ksi);        //(-2*ksi + 3*ksi*ksi)*h*dksi;
        return 0;
    }
    return 0;
}

double lagrange3_2D(const POINT2 &p0, const double hx, const double hy, const POINT2 &p, const int i, const DIF_STATE2 &dif)
{
#define lXI(i)	((i) % 4)
#define lYI(i)	((i) / 4)
    int k1 = lXI(i);
    int k2 = lYI(i);
    return	lagrange3_1D(p0[0], hx, p[0], k1, dif[0])*
            lagrange3_1D(p0[1], hy, p[1], k2, dif[1]);
}
double hermite_2D(const POINT2 &p0, const double hx, const double hy, const POINT2 &p, const int i, const DIF_STATE2 &dif)
{
#define hXI(i)	(2 * (((i) / 4) % 2) + ((i) % 2))
#define hYI(i)	(2 * ((i) / 8) + (((i) / 2) % 2))
    int k1 = hXI(i);
    int k2 = hYI(i);
    //return	hermite1D_ttt(p0[0], p0[0] + hx, p[0], k1, dif[0])*
    //        hermite1D_ttt(p0[1], p0[1] + hy, p[1], k2, dif[1]);
    return	hermite_1D(p0[0], hx, p[0], k1, dif[0])*
            hermite_1D(p0[1], hy, p[1], k2, dif[1]);
}















void calcCubeLinearBasisFuncValues(const Integration::IntegrationType integrationType,
                          double *integrationw, double *linearBasCube, VECTOR3 *dLinearBasCube)
{
    using namespace Integration;
    Integrator integrationCube;
    integrationCube.init3D(integrationType);
    for (int valIndex = 0; valIndex < integrationCube.size; valIndex++)
    {
        POINT3 point = integrationCube.p3[valIndex];
        integrationw[valIndex] = integrationCube.w[valIndex] * integrationCube.detJ;
        for (int m = 0; m < 8; m++)
        {
            dLinearBasCube[valIndex*8 + m][0] = Fem::cubeLagrange1_3D(point, m, Fem::dif_XYZ3[0]);
            dLinearBasCube[valIndex*8 + m][1] = Fem::cubeLagrange1_3D(point, m, Fem::dif_XYZ3[1]);
            dLinearBasCube[valIndex*8 + m][2] = Fem::cubeLagrange1_3D(point, m, Fem::dif_XYZ3[2]);
            linearBasCube[m*integrationCube.size + valIndex] = Fem::cubeLagrange1_3D(point, m, Fem::dif_NULL3);
        }
    }
}

// заполнение таблиц точек Гаусса, предварительных коэффициентов для численного интегрирования,
// значений безисных функций и их производных на шаблонном кубе
void calcCubeLinearBasisFuncValues_ttt(const Integration::IntegrationType integrationType,
                          std::vector<double> &integrationw, std::vector<double> &linearBasCube, std::vector<VECTOR3> &dLinearBasCube)
{
    using namespace Integration;
    Integrator integrationCube;
    integrationCube.init3D(integrationType);
    integrationw.resize(integrationCube.size);
    linearBasCube.resize(integrationCube.size*8);
    dLinearBasCube.resize(integrationCube.size*8);
    for (int valIndex = 0; valIndex < integrationCube.size; valIndex++)
    {
        POINT3 point = integrationCube.p3[valIndex];
        integrationw[valIndex] = integrationCube.w[valIndex] * integrationCube.detJ;
        for (int m = 0; m < 8; m++)
        {
            dLinearBasCube[valIndex*8 + m][0] = Fem::cubeLagrange1_3D(point, m, Fem::dif_XYZ3[0]);
            dLinearBasCube[valIndex*8 + m][1] = Fem::cubeLagrange1_3D(point, m, Fem::dif_XYZ3[1]);
            dLinearBasCube[valIndex*8 + m][2] = Fem::cubeLagrange1_3D(point, m, Fem::dif_XYZ3[2]);
            linearBasCube[m*integrationCube.size + valIndex] = Fem::cubeLagrange1_3D(point, m, Fem::dif_NULL3);
        }
    }
}

void calcCubeQuadraticBasisFuncValues(const Integration::IntegrationType integrationType,
                                       VECTOR3 *dQuadraticBasCube)
{
    using namespace Integration;
    Integrator integrationCube;
    integrationCube.init3D(integrationType);
    for (int valIndex = 0; valIndex < integrationCube.size; valIndex++)
    {
        POINT3 point = integrationCube.p3[valIndex];
        for (int m = 0; m < 27; m++)
        {
            dQuadraticBasCube[valIndex*27 + m][0] = Fem::cubeLagrange2_3D(point, m, Fem::dif_XYZ3[0]);
            dQuadraticBasCube[valIndex*27 + m][1] = Fem::cubeLagrange2_3D(point, m, Fem::dif_XYZ3[1]);
            dQuadraticBasCube[valIndex*27 + m][2] = Fem::cubeLagrange2_3D(point, m, Fem::dif_XYZ3[2]);
        }
    }
}


// заполнение таблиц точек Гаусса, предварительных коэффициентов для численного интегрирования,
// значений безисных функций и их производных на шаблонном кубе
void calcCubeQuadraticBasisFuncValues_ttt(const Integration::IntegrationType integrationType,
                          std::vector<VECTOR3> &dQuadraticBasCube)
{
    using namespace Integration;
    Integrator integrationCube;
    integrationCube.init3D(integrationType);
    dQuadraticBasCube.resize(integrationCube.size*27);
    for (int valIndex = 0; valIndex < integrationCube.size; valIndex++)
    {
        POINT3 point = integrationCube.p3[valIndex];
        for (int m = 0; m < 27; m++)
        {
            dQuadraticBasCube[valIndex*27 + m][0] = Fem::cubeLagrange2_3D(point, m, Fem::dif_XYZ3[0]);
            dQuadraticBasCube[valIndex*27 + m][1] = Fem::cubeLagrange2_3D(point, m, Fem::dif_XYZ3[1]);
            dQuadraticBasCube[valIndex*27 + m][2] = Fem::cubeLagrange2_3D(point, m, Fem::dif_XYZ3[2]);
        }
    }
}

// заполнение таблиц модулей детерминантов, окончательных множителей для численного интегрирования,
// и производных базисных функций в точках Гаусса шестигранника с линейным отображением
void calcLinearHexagonLinearBasisFuncValues(const POINT3 *v, const int numPoints, const double *integrationw, const VECTOR3 *dLinearBasCube,
                          double *w, double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3])
{
    // значения производных базисных функций, детерминантов и коэффициентов численного интегрирования
    for (int valIndex = 0; valIndex < numPoints; valIndex++)
    {
        VECTOR3 J[3] = {};
        for (int i = 0; i < 8; i++)
        {
            J[0] += dLinearBasCube[valIndex*8 + i][0] * v[i]; // производная по X
            J[1] += dLinearBasCube[valIndex*8 + i][1] * v[i]; // производная по Y
            J[2] += dLinearBasCube[valIndex*8 + i][2] * v[i]; // производная по Z
        }
        // детерминант
        double det = calcDet3x3(
                    J[0][0], J[0][1], J[0][2],
                    J[1][0], J[1][1], J[1][2],
                    J[2][0], J[2][1], J[2][2]);
        // множители для интегрирования
        w[valIndex] = integrationw[valIndex]*fabs(det);
        // производные по x, y, z
        for (int i = 0; i < 8; i++)
        {
            const VECTOR3 &dLinearBasCube0 = dLinearBasCube[valIndex*8 + i];
            // dbasi/dx
            dbas[i][valIndex][0] = calcDet3x3(
                    dLinearBasCube0[0], J[0][1], J[0][2],
                    dLinearBasCube0[1], J[1][1], J[1][2],
                    dLinearBasCube0[2], J[2][1], J[2][2]) / det;
            // dbasi/dy
            dbas[i][valIndex][1] = calcDet3x3(
                    J[0][0], dLinearBasCube0[0], J[0][2],
                    J[1][0], dLinearBasCube0[1], J[1][2],
                    J[2][0], dLinearBasCube0[2], J[2][2]) / det;
            // dbasi/dz
            dbas[i][valIndex][2] = calcDet3x3(
                    J[0][0], J[0][1], dLinearBasCube0[0],
                    J[1][0], J[1][1], dLinearBasCube0[1],
                    J[2][0], J[2][1], dLinearBasCube0[2]) / det;
        }
    }
}

// заполнение таблиц модулей детерминантов, окончательных множителей для численного интегрирования,
// и производных базисных функций в точках Гаусса шестигранника с квадратичным отображением
void calcQuadraticHexagonLinearBasisFuncValues(const POINT3 *v, const int numPoints, const double *integrationw, const VECTOR3 *dLinearBasCube, const VECTOR3 *dQuadraticBasCube,
                          double *w, double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3])
{
    // значения производных базисных функций, детерминантов и коэффициентов численного интегрирования
    for (int valIndex = 0; valIndex < numPoints; valIndex++)
    {
        VECTOR3 J[3] = {};
        for (int i = 0; i < 27; i++)
        {
            J[0] += dQuadraticBasCube[valIndex*27 + i][0] * v[i]; // производная по X
            J[1] += dQuadraticBasCube[valIndex*27 + i][1] * v[i]; // производная по Y
            J[2] += dQuadraticBasCube[valIndex*27 + i][2] * v[i]; // производная по Z
        }
        // детерминант
        double det = calcDet3x3(
                    J[0][0], J[0][1], J[0][2],
                    J[1][0], J[1][1], J[1][2],
                    J[2][0], J[2][1], J[2][2]);
        // множители для интегрирования
        w[valIndex] = integrationw[valIndex]*fabs(det);
        // производные по x, y, z
        for (int i = 0; i < 8; i++)
        {
            const VECTOR3 &dLinearCube0 = dLinearBasCube[valIndex*8 + i];
            // dbasi/dx
            dbas[i][valIndex][0] = calcDet3x3(
                    dLinearCube0[0], J[0][1], J[0][2],
                    dLinearCube0[1], J[1][1], J[1][2],
                    dLinearCube0[2], J[2][1], J[2][2]) / det;
            // dbasi/dy
            dbas[i][valIndex][1] = calcDet3x3(
                    J[0][0], dLinearCube0[0], J[0][2],
                    J[1][0], dLinearCube0[1], J[1][2],
                    J[2][0], dLinearCube0[2], J[2][2]) / det;
            // dbasi/dz
            dbas[i][valIndex][2] = calcDet3x3(
                    J[0][0], J[0][1], dLinearCube0[0],
                    J[1][0], J[1][1], dLinearCube0[1],
                    J[2][0], J[2][1], dLinearCube0[2]) / det;
        }
    }
}

// интегралы произведений intForG_mnjl[m][n][j][l] = integral(df_m/dx_j * df_n/dx_l)
void calcIntForG(const int numPoints, const double *w, const double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3],
                          double (&intForG_mnjl)[8][8][3][3])
{
    for (int m = 0; m < 8; m++)
        for (int n = 0; n < 8; n++) // m, n - локальные номера базисных функций (вершин)
            for (int j = 0; j < 3; j++)
                for (int l = 0; l < 3; l++)
                {
                    double E = 0;
                    for (int valIndex = 0; valIndex < numPoints; valIndex++)
                        E += w[valIndex] *
                             dbas[m][valIndex][j] * dbas[n][valIndex][l];
                    intForG_mnjl[m][n][j][l] = E;
                }
}

// интегралы произведений intForM_mn[m][n] = integral(f_m * f_n)
void calcIntForM_ttt(const int numPoints, const double *w, const std::vector<double> &basCube,
                          double (&intForM_mn)[8][8])
{
    for (int m = 0; m < 8; m++)
        for (int n = 0; n < 8; n++) // m, n - локальные номера базисных функций (вершин)
        {
            double E = 0;
            for (int valIndex = 0; valIndex < numPoints; valIndex++)
                E += w[valIndex] *
                     basCube[m*numPoints + valIndex] * basCube[n*numPoints + valIndex];
            intForM_mn[m][n] = E;
        }
}

// интегралы произведений intForb1_m[m] = integral(f_m)
void calcIntForb_f(const int numPoints, const double *w, const double *basCube,
                          double (&intForb1_m)[8])
{
    for (int m = 0; m < 8; m++)
    {
        double E = 0;
        for (int valIndex = 0; valIndex < numPoints; valIndex++)
            E += w[valIndex] *
                 basCube[m*numPoints + valIndex];
        intForb1_m[m] = E;
    }
}

// интегралы произведений intForb23_mj[m][j] = integral(df_m/dx_j)
void calcIntForb_df(const int numPoints, const double *w, const double (&dbas)[8][Integration::Integration3D_Size_Gauss3][3],
                          double (&intForb2_mj)[8][3])
{
    for (int m = 0; m < 8; m++)
        for (int j = 0; j < 3; j++)
        {
            double E = 0;
            for (int valIndex = 0; valIndex < numPoints; valIndex++)
                E += w[valIndex] *
                     dbas[m][valIndex][j];
            intForb2_mj[m][j] = E;
        }
}


/*
// одномерная линейная функция
double linear1D(const double a0, const double a1, const double a, const int i)
{
    if (i == 0) return (a1 - a) / (a1 - a0);
    if (i == 1) return (a - a0) / (a1 - a0);
    return 0;
}
// одномерная скалярная базисная функция 1-го порядка в точке a (i = 0..1)
double linear1D(const double a0, const double a1, const double a, const int i, const int dif)
{
    double h = a1 - a0;
    if (dif == 0)
    {
        if (i == 0) return (a1 - a) / h;
        if (i == 1) return (a - a0) / h;
    }
    else
        if (dif == 1)   // первая производная по a
        {
            if (i == 0) return -1. / h;
            if (i == 1) return  1. / h;
        }
    return 0;
}
// вычисление скалярной трилинейной базисной функции в точке p
double rectangleScalarLinear3D(const POINT3 &p0, const POINT3 &p, const int i, const double hx, const double hy, const double hz)
{
    int k1, k2, k3;
    k1 = i % 2;
    k2 = (i / 2) % 2;
    k3 = i / 4;
    return linear1D(p0[0], p0[0] + hx, p[0], k1)*
        linear1D(p0[1], p0[1] + hy, p[1], k2)*
        linear1D(p0[2], p0[2] + hz, p[2], k3);
}
void solveTriangleMatrix(const double x1, const double x2, const double x3, const double y1, const double y2, const double y3, MATR3x3 &inverse_D, double &det)
{
    //MATR3x3 D;
    //D.m[0][0] = 1;
    //D.m[0][1] = 1;
    //D.m[0][2] = 1;
    //D.m[1][0] = x1;
    //D.m[1][1] = x2;
    //D.m[1][2] = x3;
    //D.m[2][0] = y1;
    //D.m[2][1] = y2;
    //D.m[2][2] = y3;
    det = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
    inverse_D.m[0][0] = (x2*y3 - x3*y2) / det;
    inverse_D.m[0][1] = (y2 - y3) / det;
    inverse_D.m[0][2] = (x3 - x2) / det;
    inverse_D.m[1][0] = (x3*y1 - x1*y3) / det;
    inverse_D.m[1][1] = (y3 - y1) / det;
    inverse_D.m[1][2] = (x1 - x3) / det;
    inverse_D.m[2][0] = (x1*y2 - x2*y1) / det;
    inverse_D.m[2][1] = (y1 - y2) / det;
    inverse_D.m[2][2] = (x2 - x1) / det;
    det = fabs(det);
}
double triangleScalarLinear(const MATR3x3 &inverse_D, const double x, const double y, const int i, const int difx, const int dify)
{
    if (difx + dify >= 2)
        return 0;
    if (difx == 1)
        return inverse_D.m[i][1];
    if (dify == 1)
        return inverse_D.m[i][2];
    return inverse_D.m[i][0] + inverse_D.m[i][1] * x + inverse_D.m[i][2] * y;
}
double prism3ScalarLinear(const MATR3x3 &inverse_D, const double z1, const double z2, const POINT3 &p, const int i, const DIF_STATE3 &dif)
{
    return
        triangleScalarLinear(inverse_D, p[0], p[1], i % 3, dif[0], dif[1])*
        linear1D(z1, z2, p[2], i / 3, dif[2]);
}
*/

}   // namespace Fem
