#define _CRT_SECURE_NO_WARNINGS

#include "elementary.h"
#include "stdio.h"

namespace Elementary
{

VECTOR3 vector3Mul(const VECTOR3 &a, const VECTOR3 &b)
{
    return VECTOR3(a[1]*b[2] - a[2]*b[1],
                   a[2]*b[0] - a[0]*b[2],
                   a[0]*b[1] - a[1]*b[0]);
}


void C1mulVector1PlusC2mulVector2(const double c1, const Vector &x, const double c2, const Vector &y, Vector &z)
{
    for (size_t i = 0; i < x.size(); i++)
        z[i] = c1*x[i] + c2*y[i];
}
/*void Vector1PlusCmulVector2(const Vector &x, const double c, const Vector &y, Vector &z)
{
    for (size_t i = 0; i < x.size(); i++)
        z[i] = x[i] + c * y[i];
}
void Vector3_1PlusCmulVector3_2(const Vector3 &x, const double c, const Vector3 &y, Vector3 &z)
{
    for (size_t i = 0; i < x.size(); i++)
        z[i] = x[i] + c * y[i];
}*/

double VectorScalMul(const Vector &x, const Vector &y)
{
    double E = 0;
    for (size_t i = 0; i < x.size(); i++)
        E += x[i] * y[i];
    return E;
}
double Vector3ScalMul(const Vector3 &x, const Vector3 &y)
{
    double E = 0;
    for (size_t i = 0; i < x.size(); i++)
        E += x[i] * y[i];
    return E;
}

void VectorMulc(const Vector &x, const double c, Vector &y)
{
    for (size_t i = 0; i < x.size(); i++)
        y[i] = x[i] * c;
}

namespace Operations
{

// M1+M2 -> Mres
void c1M1plusc2M26(const double &c1, const MATR6x6 &M1, const double &c2, const MATR6x6 &M2, MATR6x6 &Mres)
{
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            Mres.m[i][j] = c1*M1.m[i][j] + c2*M2.m[i][j];
}
// M1*M2 -> Mres, Mres != M1, Mres != M2
void MmulM6(const MATR6x6 &M1, const MATR6x6 &M2, MATR6x6 &Mres)
{
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
        {
            Mres.m[i][j] = 0;
            for (int k = 0; k < 6; k++)
                Mres.m[i][j] += M1.m[i][k] * M2.m[k][j];
        }
}
// c*M1*M2 -> Mres
void cMmulM6(const double &c, const MATR6x6 &M1, const MATR6x6 &M2, MATR6x6 &Mres)
{
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
        {
            Mres.m[i][j] = 0;
            for (int k = 0; k < 6; k++)
                Mres.m[i][j] += M1.m[i][k] * M2.m[k][j];
            Mres.m[i][j] *= c;
        }
}
// M*c -> Mres
void Mmulc6(const MATR6x6 &M, const double &c, MATR6x6 &Mres)
{
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            Mres.m[i][j] = M.m[i][j] * c;
}
// M+cI -> Mres
void cIplusM6(const double &c, MATR6x6 &Mres)
{
    for (int i = 0; i < 6; i++)
        Mres.m[i][i] += c;
}

// пересчет индексов для конвертации матрицы C в матрицу D
int indexConvertion3x3to6(const int i, const int j)
{
    if (i == j) return i;
    if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 3;
    if ((i == 0 && j == 2) || (i == 2 && j == 0)) return 4;
    if ((i == 0 && j == 1) || (i == 1 && j == 0)) return 5;
    return 0;
}

// перестройка матрицы упругих констант к тензорному виду
void convert6x6to3x3x3x3(const MATR6x6 &D, MATR3x3x3x3 &C)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    C.m[i][j][k][l] = D.m[indexConvertion3x3to6(i, j)][indexConvertion3x3to6(k, l)];
}

// перестройка матрицы упругих констант к матричному виду
void convert3x3x3x3to6x6(const MATR3x3x3x3 &C, MATR6x6 &D)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    D.m[indexConvertion3x3to6(i, j)][indexConvertion3x3to6(k, l)] = C.m[i][j][k][l];
}

// нормаль к плоскости, заданной векторами p0-p1 и p0-p2
void calcNormal(const POINT3 &p0, const POINT3 &p1, const POINT3 &p2, VECTOR3 &normal)
{
    double m;
    VECTOR3 v1, v2;
    v1.x[0] = p1.x[0] - p0.x[0];
    v1.x[1] = p1.x[1] - p0.x[1];
    v1.x[2] = p1.x[2] - p0.x[2];
    v2.x[0] = p2.x[0] - p0.x[0];
    v2.x[1] = p2.x[1] - p0.x[1];
    v2.x[2] = p2.x[2] - p0.x[2];
    normal.x[0] = v1.x[1] * v2.x[2] - v1.x[2] * v2.x[1];
    normal.x[1] = v1.x[2] * v2.x[0] - v1.x[0] * v2.x[2];
    normal.x[2] = v1.x[0] * v2.x[1] - v1.x[1] * v2.x[0];
    m = sqrt(SQR(normal.x[0]) + SQR(normal.x[1]) + SQR(normal.x[2]));
    normal.x[0] /= m;
    normal.x[1] /= m;
    normal.x[2] /= m;
}

// функция строит точку (p) на дуге (a->b) сферы (0, r), такую, что угол (a,c) равен fi
// условие: дуга (a->b) не должна быть полуокружностью
void sphereArcPoint(const double R, const POINT3 &a, const POINT3 &b, const double fi, POINT3 &p)
{
    POINT3 e = a * ((a * b) / R / R);
    POINT3 n = (b - e) / (b - e).abs();
    p = a * cos(fi) + n * (R * sin(fi));
    //printf("|n| = %le, |a / a.abs()| = %le\n", n.abs(), (a / a.abs()).abs());
    //printf("(n,a) = %le\n", n * a);
    //printf("(n,a) = %le\n", p.abs());
}

void findPointOnTheLine_1d(const double p1, const double p2, const int N, const double q, const int ind, double &p)
{
    if (q != 1.)
    {
        double c = (pow(q, ind) - 1) / (pow(q, N) - 1);	// 0 <= c <= 1
        // p = p1+(p2-p1)*c = p1*(1-c)+p2*c
        p = p1*(1 - c) + p2*c;
    }
    else
    {
        p = p1 + (p2 - p1) * ind / N;
    }
}

void findPointOnTheLine_3d(const POINT3 p1, const POINT3 p2, const int N, const double q, const int ind, POINT3 &p)
{
    if (q != 1.)
    {
        double c = (pow(q, ind) - 1) / (pow(q, N) - 1);	// 0 <= c <= 1
        // p = p1+(p2-p1)*c = p1*(1-c)+p2*c
        p = p1*(1 - c) + p2*c;
    }
    else
    {
        p = p1 + (p2 - p1) * ind / N;
    }
}

void findPointOnTheLine_1d_conc(const double p1, const double p2, const int N, const double q, const int ind, double &p)
{
    if (q != 1.)
    {
        double p_center = (p1 + p2)/2;
        if(ind <= N/2)
        {
            double c = (pow(q, ind) - 1) / (pow(q, N/2) - 1);	// 0 <= c <= 1
            // p = p1+(p2-p1)*c = p1*(1-c)+p2*c
            p = p1*(1 - c) + p_center*c;
        }
        else
        {
            double c = (pow(q, N - ind) - 1) / (pow(q, N/2) - 1);	// 0 <= c <= 1
            // p = p1+(p2-p1)*c = p1*(1-c)+p2*c
            p = p2*(1 - c) + p_center*c;
        }
    }
    else
    {
        p = p1 + (p2 - p1) * ind / N;
    }
}

void solvePolynom3(double a, double b, double c, double &r1, double &r2, double &r3)
{
    //double p = c/a - b*b/3/a/a;
    //double q = 2*b*b*b/27/a/a/a - b*c/3/a/a + d/a;
    double p = b - a * a / 3.0;
    double q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
    double A = sqrt(- 4.0 * p / 3.0);
    double c3phi = - 4.0 * q / (A * A * A);
    double phi = acos(c3phi) / 3.0;
    //printf("%lf\n", - 4.0 * p / 3.0);
    r1 = A * cos(phi) - a / 3.0;
    r2 = A * cos(phi + 2 * PI / 3.0) - a / 3.0;
    r3 = A * cos(phi - 2 * PI / 3.0) - a / 3.0;
}


}

// namespace Elementary_operations
}   // namespace Elementary


