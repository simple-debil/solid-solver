/* --------------------------------------------------------- */
// ОСНОВНЫЕ СТРУКТУРЫ И ФУНКЦИИ
/* --------------------------------------------------------- */

#ifndef ELEMENTARY_H
#define ELEMENTARY_H

#include <vector>
#include "math.h"

#define MALENKOECHISLO (1.e-20)
#define MALENKOECHISLO1 (1.e-10)


template <typename T>
inline T MIN(T a, T b)
{
    return a < b ? a : b;
}

template <typename T>
inline T MAX(T a, T b)
{
    return a > b ? a : b;
}

template <typename T>
inline int SIGN(T a)
{
    return (a < 0) ? (-1) : (1);
}

template <typename T>
inline T SQR(T a)
{
    return a*a;
}
const double PI = 3.14159265358979323846264338327;

// элементарные типы данных и операции
namespace Elementary
{

// 3D точка, вектор
struct VECTOR3
{
    double x[3];
    VECTOR3(VECTOR3 const&) = default;
    VECTOR3(std::initializer_list<VECTOR3> v)
    {
        x[0] = (*v.begin()).x[0];
        x[1] = (*v.begin()).x[1];
        x[2] = (*v.begin()).x[2];
    }

    VECTOR3()
    {
        x[0] = 0;
        x[1] = 0;
        x[2] = 0;
    }
    VECTOR3(const double x0, const double x1, const double x2)
    {
        x[0] = x0;
        x[1] = x1;
        x[2] = x2;
    }
    VECTOR3(const VECTOR3 &&v0):
        x{v0.x[0], v0.x[1], v0.x[2]}
    {
    }
    inline void clear()
    {
        x[0] = 0;
        x[1] = 0;
        x[2] = 0;
    }
    inline void operator=(const VECTOR3 &v)
    {
        x[0] = v[0];
        x[1] = v[1];
        x[2] = v[2];
    }
    inline double &operator[](const int &i)
    {
        return x[i];
    }
    inline const double &operator[](const int &i)const
    {
        return x[i];
    }
    inline double abs()const
    {
        return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    }
    inline void operator*=(const double &t)
    {
        x[0] *= t;
        x[1] *= t;
        x[2] *= t;
    }
    inline void operator/=(const double &t)
    {
        x[0] /= t;
        x[1] /= t;
        x[2] /= t;
    }
    inline void operator+=(const VECTOR3 &v)
    {
        x[0] += v[0];
        x[1] += v[1];
        x[2] += v[2];
    }
    inline void operator-=(const VECTOR3 &v)
    {
        x[0] -= v[0];
        x[1] -= v[1];
        x[2] -= v[2];
    }
    inline bool operator==(const VECTOR3 &v)const
    {
        return x[0] == v[0] &&
               x[1] == v[1] &&
               x[2] == v[2];
    }

    friend inline const VECTOR3 operator+(const VECTOR3 &v)
    {
        return v;
    }
    friend inline const VECTOR3 operator-(const VECTOR3 &v)
    {
        return VECTOR3(-v[0], -v[1], -v[2]);
    }

    friend inline const VECTOR3 operator+(const VECTOR3 &l, const VECTOR3 &r)
    {
        return VECTOR3(l[0] + r[0], l[1] + r[1], l[2] + r[2]);
    }
    friend inline const VECTOR3 operator-(const VECTOR3 &l, const VECTOR3 &r)
    {
        return VECTOR3(l[0] - r[0], l[1] - r[1], l[2] - r[2]);
    }
    friend inline double operator*(const VECTOR3 &l, const VECTOR3 &r)
    {
        return l[0] * r[0] + l[1] * r[1] + l[2] * r[2];
    }
    friend inline const VECTOR3 operator*(const double &l, const VECTOR3 &r)
    {
        return VECTOR3(l * r[0], l * r[1], l * r[2]);
    }
    friend inline const VECTOR3 operator*(const VECTOR3 &l, const double &r)
    {
        return VECTOR3(l[0] * r, l[1] * r, l[2] * r);
    }
    friend inline const VECTOR3 operator/(const VECTOR3 &l, const double &r)
    {
        return VECTOR3(l[0] / r, l[1] / r, l[2] / r);
    }
    friend inline double operator^(const VECTOR3 &l, const VECTOR3 &r)
    {
        double t = (l*r)/(l.abs()*r.abs());
        if(t < -1)
            t = -1;
        if(t > 1)
            t = 1;
        return acos(t);
    }
};
// векторное произведение
VECTOR3 vector3Mul(const VECTOR3 &a, const VECTOR3 &b);



// 3D точка, вектор с целыми координатами
struct VECTOR3_int
{
    int x[3];
    VECTOR3_int(VECTOR3_int const&) = default;
    VECTOR3_int(std::initializer_list<VECTOR3_int> v)
    {
        x[0] = (*v.begin()).x[0];
        x[1] = (*v.begin()).x[1];
        x[2] = (*v.begin()).x[2];
    }
    VECTOR3_int()
    {
        x[0] = 0;
        x[1] = 0;
        x[2] = 0;
    }
    VECTOR3_int(const int x0, const int x1, const int x2)
    {
        x[0] = x0;
        x[1] = x1;
        x[2] = x2;
    }
    VECTOR3_int(const VECTOR3_int &&v0):
        x{v0.x[0], v0.x[1], v0.x[2]}
    {
    }
    inline void operator=(const VECTOR3_int &v)
    {
        x[0] = v[0];
        x[1] = v[1];
        x[2] = v[2];
    }
    inline int &operator[](const int &i)
    {
        return x[i];
    }
    inline const int &operator[](const int &i)const
    {
        return x[i];
    }
    inline void operator*=(const int &t)
    {
        x[0] *= t;
        x[1] *= t;
        x[2] *= t;
    }
    inline void operator/=(const int &t)
    {
        x[0] /= t;
        x[1] /= t;
        x[2] /= t;
    }
    inline void operator+=(const VECTOR3_int &v)
    {
        x[0] += v[0];
        x[1] += v[1];
        x[2] += v[2];
    }
    inline void operator-=(const VECTOR3_int &v)
    {
        x[0] -= v[0];
        x[1] -= v[1];
        x[2] -= v[2];
    }

    friend inline const VECTOR3_int operator+(const VECTOR3_int &v)
    {
        return v;
    }
    friend inline const VECTOR3_int operator-(const VECTOR3_int &v)
    {
        return VECTOR3_int(-v[0], -v[1], -v[2]);
    }

    friend inline const VECTOR3_int operator+(const VECTOR3_int &l, const VECTOR3_int &r)
    {
        return VECTOR3_int(l[0] + r[0], l[1] + r[1], l[2] + r[2]);
    }
    friend inline const VECTOR3_int operator-(const VECTOR3_int &l, const VECTOR3_int &r)
    {
        return VECTOR3_int(l[0] - r[0], l[1] - r[1], l[2] - r[2]);
    }
    friend inline double operator*(const VECTOR3_int &l, const VECTOR3_int &r)
    {
        return l[0] * r[0] + l[1] * r[1] + l[2] * r[2];
    }
    friend inline const VECTOR3_int operator*(const int &l, const VECTOR3_int &r)
    {
        return VECTOR3_int(l * r[0], l * r[1], l * r[2]);
    }
    friend inline const VECTOR3_int operator*(const VECTOR3_int &l, const int &r)
    {
        return VECTOR3_int(l[0] * r, l[1] * r, l[2] * r);
    }
    friend inline const VECTOR3_int operator/(const VECTOR3_int &l, const int &r)
    {
        return VECTOR3_int(l[0] / r, l[1] / r, l[2] / r);
    }
};

// 3D точка, вектор с бесзнаковыми целыми координатами
struct VECTOR3_uint
{
    unsigned int x[3];
    VECTOR3_uint(VECTOR3_uint const&) = default;
    VECTOR3_uint(std::initializer_list<VECTOR3_uint> v)
    {
        x[0] = (*v.begin()).x[0];
        x[1] = (*v.begin()).x[1];
        x[2] = (*v.begin()).x[2];
    }
    VECTOR3_uint()
    {
        x[0] = 0;
        x[1] = 0;
        x[2] = 0;
    }
    VECTOR3_uint(const unsigned int x0, const unsigned int x1, const unsigned int x2)
    {
        x[0] = x0;
        x[1] = x1;
        x[2] = x2;
    }
    VECTOR3_uint(const VECTOR3_uint &&v0):
        x{v0.x[0], v0.x[1], v0.x[2]}
    {
    }
    inline void operator=(const VECTOR3_uint &v)
    {
        x[0] = v[0];
        x[1] = v[1];
        x[2] = v[2];
    }
    inline unsigned int &operator[](const int &i)
    {
        return x[i];
    }
    inline const unsigned int &operator[](const unsigned int &i)const
    {
        return x[i];
    }
    friend inline bool operator==(const VECTOR3_uint &v1, const VECTOR3_uint &v2)
    {
        return v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2];
    }
    friend inline const VECTOR3_uint operator+(const VECTOR3_uint &l, const VECTOR3_uint &r)
    {
        return VECTOR3_uint(l[0] + r[0], l[1] + r[1], l[2] + r[2]);
    }

};
typedef VECTOR3_int POINT3_int;
typedef VECTOR3 POINT3_CUBE, POINT3;
// нулевой вектор
const VECTOR3 VECTOR3_NULL = { 0, 0, 0 };
const POINT3_CUBE CUBE_VERTEX[8] =
{
    { -1, -1, -1 },
    { +1, -1, -1 },
    { -1, +1, -1 },
    { +1, +1, -1 },
    { -1, -1, +1 },
    { +1, -1, +1 },
    { -1, +1, +1 },
    { +1, +1, +1 },
};
// 2D точка, вектор
struct VECTOR2
{
    double x[2];
    VECTOR2(VECTOR2 const&) = default;
    VECTOR2(std::initializer_list<VECTOR2> v)
    {
        x[0] = (*v.begin()).x[0];
        x[1] = (*v.begin()).x[1];
    }

    VECTOR2()
    {
        x[0] = 0;
        x[1] = 0;
    }
    VECTOR2(const double x0, const double x1)
    {
        x[0] = x0;
        x[1] = x1;
    }
    VECTOR2(const VECTOR2 &&v0):
        x{v0.x[0], v0.x[1]}
    {
    }
    inline double abs()const
    {
        return sqrt(x[0]*x[0] + x[1]*x[1]);
    }
    inline void operator*=(const double &t)
    {
        x[0] *= t;
        x[1] *= t;
    }
    inline void operator/=(const double &t)
    {
        x[0] /= t;
        x[1] /= t;
    }
    inline void operator+=(const VECTOR2 &v)
    {
        x[0] += v[0];
        x[1] += v[1];
    }
    inline void operator-=(const VECTOR2 &v)
    {
        x[0] -= v[0];
        x[1] -= v[1];
    }
    inline void operator=(const VECTOR2 &v)
    {
        x[0] = v[0];
        x[1] = v[1];
    }
    inline double &operator[](const int &i)
    {
        return x[i];
    }
    inline const double &operator[](const int &i)const
    {
        return x[i];
    }
    friend inline const VECTOR2 operator+(const VECTOR2 &v)
    {
        return v;
    }
    friend inline const VECTOR2 operator-(const VECTOR2 &v)
    {
        return VECTOR2(-v[0], -v[1]);
    }
    friend inline const VECTOR2 operator+(const VECTOR2 &l, const VECTOR2 &r)
    {
        return VECTOR2(l[0] + r[0], l[1] + r[1]);
    }
    friend inline const VECTOR2 operator-(const VECTOR2 &l, const VECTOR2 &r)
    {
        return VECTOR2(l[0] - r[0], l[1] - r[1]);
    }
    friend inline const VECTOR2 operator*(const double &l, const VECTOR2 &r)
    {
        return VECTOR2(l * r[0], l * r[1]);
    }
    friend inline const VECTOR2 operator*(const VECTOR2 &l, const double &r)
    {
        return VECTOR2(l[0] * r, l[1] * r);
    }
    friend inline const VECTOR2 operator/(const VECTOR2 &l, const double &r)
    {
        return VECTOR2(l[0] / r, l[1] / r);
    }

};
typedef VECTOR2 POINT2_CUBE, POINT2;
// 1D точка, вектор
typedef double VECTOR1, POINT1_CUBE, POINT1;

// сколько раз дифференцировать
typedef int DIF_STATE1;
// отсутствие дифференцирования
const DIF_STATE1 dif_NULL1 = 0;
// сколько раз дифференцировать по каждой переменной (0, 1)
struct DIF_STATE2
{
    DIF_STATE1 difNumber[2];
    inline int &operator[](const int &i)
    {
        return difNumber[i];
    }
    inline const int &operator[](const int &i)const
    {
        return difNumber[i];
    }
};
// отсутствие дифференцирования
const DIF_STATE2 dif_NULL2 = { 0, 0 };
// дифференцирование по одной из координат (0|1)
const DIF_STATE2 dif_XYZ2[2] =
{
    { 1, 0 },
    { 0, 1 }
};
// сколько раз дифференцировать по каждой переменной (0,1,2)
struct DIF_STATE3
{
    DIF_STATE1 difNumber[3];
    inline int &operator[](const int &i)
    {
        return difNumber[i];
    }
    inline const int &operator[](const int &i)const
    {
        return difNumber[i];
    }
};
// отсутствие дифференцирования
const DIF_STATE3 dif_NULL3 = { 0, 0, 0 };
// дифференцирование по одной из координат (0|1|2)
const DIF_STATE3 dif_XYZ3[3] =
{
    { 1, 0, 0 },
    { 0, 1, 0 },
    { 0, 0, 1 }
};

// сфера
struct Sphere
{
    POINT3 O;
    double R;
    // инициализировать шар нулевого радиуса с заданным центром
    void init0(const POINT3 set_O)
    {
        O = set_O;
        R = 0;
    }
    // расширить шар так, чтобы он стал включать заданную точку
    void expand(const POINT3 &p)
    {
        VECTOR3 v = p - O;      // вектор O1-O2
        double v_abs = v.abs(); // длина O1-O2
        if(v_abs > R)
        {
            // точка p вне шара - нужно расширить
            O = O + (v/v_abs)*(v_abs - R) / 2;
            R = (v_abs + R) / 2;
        }
    }
    // расширить шар так, чтобы он стал включать набор точек
    /*
    void expand(const POINT3 *p, const int size)
    {
        //init0(p[0]);
        for(int i = 0; i < size; i++)
            expand(p[i]);
    }
    */
    // расширить шар так, чтобы он стал включать заданный шар
    void expand(const Sphere &s2)
    {
        VECTOR3 v = s2.O - O;   // вектор O1-O2
        double v_abs = v.abs(); // длина O1-O2
        double min1 = -R;
        double max1 = R;
        double min2 = v_abs - s2.R;
        double max2 = v_abs + s2.R;
        double min = MIN(min1, min2);
        double max = MAX(max1, max2);
        O = O + (v/v_abs)*(min + max) / 2;
        R = (max - min) / 2;
        /*
        O = O + (v/v_abs)*(v_abs + s2.R - R) / 2;
        R = (v_abs + R + s2.R) / 2;
        */
    }
    // пересекается ли шар с заданным шаром?
    bool intersect(const Sphere &s2)const
    {
        return (s2.O - O)*(s2.O - O) <= (s2.R + R)*(s2.R + R);
    }
    // пересекается ли шар с заданной точкой?
    bool intersect(const POINT3 &p)const
    {
        return (p - O)*(p - O) <= R*R;
    }
    // пересекается ли шар с заданным отрезком?
    bool intersect(const POINT3 &p1, const POINT3 &p2)const
    {
        POINT3 d = p2 - p1;        // направление отрезка
        double d_abs = d.abs();    // длина отрезка
        // отрезок представляем в виде сферы
        Sphere s0;
        s0.O = (p1 + p2) / 2;
        s0.R = d_abs / 2;
        if(intersect(s0))
        {
            // шар пересекается со сферированным отрезком
            // более точные вычисления
            if(p1 == p2)
            {
                // отрезок является точкой
                //return intersect(p1);
                // в этом случае сферированный отрезок точно описывает отрезок
                return true;
            }
            else
            {
                // отрезок не является точкой
                if(intersect(p1) || intersect(p2))
                {
                    // один из концов отрезка внутри сферы
                    return true;
                }
                else
                {
                    // оба конца отрезка вне сферы
                    if(d_abs*1.e5 < R)  //###
                    {
                        // отрезок очень мал
                        // чтобы избежать высокой погрешности вычислений, принимаем что пересечение есть
                        return true;
                    }
                    // строим перпендикуляр из центра сферы на прямую, образованную отрезком
                    d = d / d_abs;                    // единичное направление отрезка
                    POINT3 O_p1 = p1 - O;
                    double H_p1 = O_p1*d;               // проекция отрезка O_p1 на прямую p1-p2
                    return O_p1*O_p1 - H_p1*H_p1 <= R*R;
                }
            }
        }
        else
        {
            // шар не пересекается со сферированным отрезком
            // значит шар не пересекается и с самим отрезком
            return false;
        }
    }
};

// интервал
struct INTERVAL
{
    double x[2];

    INTERVAL(INTERVAL const&) = default;
    INTERVAL(std::initializer_list<INTERVAL> v)
    {
        x[0] = (*v.begin()).x[0];
        x[1] = (*v.begin()).x[1];
    }

    INTERVAL()
    {
        x[0] = 0;
        x[1] = 0;
    }
    INTERVAL(const double x0, const double x1)
    {
        x[0] = x0;
        x[1] = x1;
    }
    INTERVAL(const INTERVAL &&v0):
        x{v0.x[0], v0.x[1]}
    {
    }
    inline void operator=(const INTERVAL &v)
    {
        x[0] = v[0];
        x[1] = v[1];
    }
    inline double &operator[](const int &i)
    {
        return x[i];
    }
    inline const double &operator[](const int &i)const
    {
        return x[i];
    }
    inline double abs()const
    {
        return ::abs(x[1] - x[0]);
    }
};
// квадрат
struct SQUARE
{
    INTERVAL i[2];
};
// параллелепипед
struct CUBE
{
    INTERVAL i[3];
    // инициализироветь кубом, описанным вокруг заданной сферы
    void initBySphere(const Sphere &s)
    {
        for(int k = 0; k < 3; k++)
        {
            i[k][0] = s.O[k] - s.R;
            i[k][1] = s.O[k] + s.R;
        }
    }
    void getCenter(POINT3 &O)const
    {
        for(int k = 0; k < 3; k++)
        {
            O[k] = (i[k][0] + i[k][1]) / 2;
        }
    }
    // расширить параллелепипед так, чтобы он стал включать заданную сферу
    void expand(const Sphere &s)
    {
        for(int k = 0; k < 3; k++)
        {
            double x0_new = s.O[k] - s.R;
            double x1_new = s.O[k] + s.R;
            if(x0_new < i[k][0])
                i[k][0] = x0_new;
            if(x1_new > i[k][1])
                i[k][1] = x1_new;
        }
    }
    // преобразовать параллелепипед в куб
    void toCube()
    {
        POINT3 O;       // центр куба
        double a = 0;   // половина стороны куба
        for(int k = 0; k < 3; k++)
        {
            O[k] = (i[k][0] + i[k][1]) / 2;
            double h = fabs(O[k] - i[k][0]);
            if(h > a)
                a = h;
        }
        for(int k = 0; k < 3; k++)
        {
            i[k][0] = O[k] - a;
            i[k][1] = O[k] + a;
        }
    }
    // пересекается ли параллелепипед с заданной точкой?
    bool intersect(const POINT3 &p)const
    {
        return
                p[0] >= i[0][0] && p[0] <= i[0][1] &&
                p[1] >= i[1][0] && p[1] <= i[1][1] &&
                p[2] >= i[2][0] && p[2] <= i[2][1];
    }
};

struct TRIANGLE
{
    POINT3 v[3];
    // поиск точки o пересечения треугольника с прямой p+nt
    // возвращает true если точка пересечения есть
    bool intersect(const POINT3 &p, const VECTOR3 &n,
                   POINT3 &o, double &t0)
    {
        VECTOR3 normal = vector3Mul(v[1] - v[0], v[2] - v[0]);
        double znamen = (n*normal);
        if(fabs(znamen) < MALENKOECHISLO)
        {
            // прямая параллельна плоскости треугольника - пересечения нет
            t0 = 0;
            o = VECTOR3_NULL;
            return false;
        }
        else
        {
            t0 = ((v[0] - p)*normal) / znamen;
            o = p + n*t0;   // o - точка пересечения прямой n*t с плоскостью треугольника (при t = t0)
            // проверка принадлежности точки o треугольнику или его границе
            VECTOR3 OT1 = v[0] - o;
            VECTOR3 OT2 = v[1] - o;
            VECTOR3 OT3 = v[2] - o;
            if(vector3Mul(OT1, OT2)*normal >= -MALENKOECHISLO1 &&
               vector3Mul(OT2, OT3)*normal >= -MALENKOECHISLO1 &&
               vector3Mul(OT3, OT1)*normal >= -MALENKOECHISLO1)
            {
                // точка o принадлежит треугольнику или его границе
                return true;
            }
            else
            {
                // точка пересечения находится вне треугольника
                return false;
            }
        }
    }
};




// матрица 3х3
struct MATR3x3
{
    double m[3][3];
    MATR3x3(const MATR3x3 &) = default;
    MATR3x3()
    {
        clear();
    }
    // обнуление матрицы
    void clear()
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                m[i][j] = 0;
    }
    // инициализация диагональной матрицы с одинаковыми значениями на диагонали
    void initDiag(const double d)
    {
        clear();
        m[0][0] = d;
        m[1][1] = d;
        m[2][2] = d;
    }
    // сделать полностью симметричной (усредняет симметричные внедиагональные элементы)
    void fixAsymmetry()
    {
        m[1][0] = (m[1][0] + m[0][1]) / 2;
        m[0][1] = m[1][0];
        m[2][0] = (m[2][0] + m[0][2]) / 2;
        m[0][2] = m[2][0];
        m[2][1] = (m[2][1] + m[1][2]) / 2;
        m[1][2] = m[2][1];
    }
    // детерминант матрицы
    double inline calcDet()const
    {
        return
                + m[0][0] * m[1][1] * m[2][2]
                + m[0][1] * m[1][2] * m[2][0]
                + m[0][2] * m[1][0] * m[2][1]
                - m[0][2] * m[1][1] * m[2][0]
                - m[0][0] * m[1][2] * m[2][1]
                - m[0][1] * m[1][0] * m[2][2];
    }
    friend inline const MATR3x3 operator*(const MATR3x3 &l, const MATR3x3 &r)
    {
        MATR3x3 res;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                res.m[i][j] = 0;
                for (int k = 0; k < 3; k++)
                    res.m[i][j] += l.m[i][k] * r.m[k][j];
            }
        return res;
    }
    friend inline const MATR3x3 operator+(const MATR3x3 &l, const MATR3x3 &r)
    {
        MATR3x3 res;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                res.m[i][j] = l.m[i][j] + r.m[i][j];
        return res;
    }
    friend inline const MATR3x3 operator-(const MATR3x3 &l, const MATR3x3 &r)
    {
        MATR3x3 res;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                res.m[i][j] = l.m[i][j] - r.m[i][j];
        return res;
    }
    friend inline const VECTOR3 operator*(const MATR3x3 &l, const VECTOR3 &r)
    {
        VECTOR3 t;
        for(int i = 0; i < 3; i++)
        {
            t[i] = 0;
            for(int j = 0; j < 3; j++)
            {
                t[i] += l.m[i][j]*r[j];
            }
        }
        return t;
    }
    friend inline const MATR3x3 operator*(const MATR3x3 &l, const double &r)
    {
        MATR3x3 res;
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                res.m[i][j] = l.m[i][j]*r;
            }
        }
        return res;
    }
    inline void operator*=(const double &t)
    {
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                m[i][j] *= t;
            }
        }
    }
    inline void operator+=(const MATR3x3 &t)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                m[i][j] += t.m[i][j];
    }
    inline void operator-=(const MATR3x3 &t)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                m[i][j] -= t.m[i][j];
    }
    // обратная матрица
    void inline calcInvert(MATR3x3 &inv)const
    {
        const double det = calcDet();
        inv.m[0][0] = +(m[1][1]*m[2][2] - m[1][2]*m[2][1])/det;
        inv.m[0][1] = -(m[0][1]*m[2][2] - m[0][2]*m[2][1])/det;
        inv.m[0][2] = +(m[0][1]*m[1][2] - m[0][2]*m[1][1])/det;

        inv.m[1][0] = -(m[1][0]*m[2][2] - m[1][2]*m[2][0])/det;
        inv.m[1][1] = +(m[0][0]*m[2][2] - m[0][2]*m[2][0])/det;
        inv.m[1][2] = -(m[0][0]*m[1][2] - m[0][2]*m[1][0])/det;

        inv.m[2][0] = +(m[1][0]*m[2][1] - m[1][1]*m[2][0])/det;
        inv.m[2][1] = -(m[0][0]*m[2][1] - m[0][1]*m[2][0])/det;
        inv.m[2][2] = +(m[0][0]*m[1][1] - m[0][1]*m[1][0])/det;
    }
    // обратная матрица если исходная матрица симметричная
    void inline calcInvertSym(MATR3x3 &inv)const
    {
        const double det = calcDet();
        inv.m[0][0] = +(m[1][1]*m[2][2] - m[1][2]*m[2][1])/det;
        //inv.m[0][1] = -(m[0][1]*m[2][2] - m[0][2]*m[2][1])/det;
        //inv.m[0][2] = +(m[0][1]*m[1][2] - m[0][2]*m[1][1])/det;

        inv.m[1][0] = -(m[1][0]*m[2][2] - m[1][2]*m[2][0])/det;
        inv.m[1][1] = +(m[0][0]*m[2][2] - m[0][2]*m[2][0])/det;
        //inv.m[1][2] = -(m[0][0]*m[1][2] - m[0][2]*m[1][0])/det;

        inv.m[2][0] = +(m[1][0]*m[2][1] - m[1][1]*m[2][0])/det;
        inv.m[2][1] = -(m[0][0]*m[2][1] - m[0][1]*m[2][0])/det;
        inv.m[2][2] = +(m[0][0]*m[1][1] - m[0][1]*m[1][0])/det;

        inv.m[0][1] = inv.m[1][0];
        inv.m[0][2] = inv.m[2][0];
        inv.m[1][2] = inv.m[2][1];
    }
    // матрица-матрица
    // res = a * b
    inline static void a_mul_b(const MATR3x3 &a, const MATR3x3 &b,
                               MATR3x3 &res)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                res.m[i][j] = 0;
                for (int k = 0; k < 3; k++)
                    res.m[i][j] += a.m[i][k] * b.m[k][j];
            }
        }
    }
    // E += a * b_T_transpose - умножение матрицы a на транспонированную матрицу b_T
    inline static void a_Mul_bT_add(const MATR3x3 &a, const MATR3x3 &b_T,
                                    MATR3x3 &E)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                    E.m[i][j] += a.m[i][k] * b_T.m[j][k];
            }
        }
    }
    // матрица-вектор
    // res = a * b - умножение матрицы a на вектор b
    inline static void a_mul_b(const MATR3x3 &a, const VECTOR3 &b,
                               VECTOR3 &res)
    {
        res.clear();
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                res[i] += a.m[i][j]*b[j];
            }
        }
    }
    // E += a * b - умножение матрицы a на вектор b
    inline static void a_mul_b_add(const MATR3x3 &a, const VECTOR3 &b,
                                   VECTOR3 &E)
    {
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                E[i] += a.m[i][j]*b[j];
            }
        }
    }
    // E += a_T_transpose * b - умножение транспонированной матрицы a_T на вектор b
    inline static void aT_mul_b_add(const MATR3x3 &a_T, const VECTOR3 &b,
                                    VECTOR3 &E)
    {
        for(int j = 0; j < 3; j++)
        {
            for(int i = 0; i < 3; i++)
            {
                E[i] += a_T.m[j][i]*b[j];
            }
        }
    }
    // E -= a_T_transpose * b - умножение транспонированной матрицы a_T на вектор b
    inline static void aT_mul_b_sub(const MATR3x3 &a_T, const VECTOR3 &b,
                                    VECTOR3 &E)
    {
        for(int j = 0; j < 3; j++)
        {
            for(int i = 0; i < 3; i++)
            {
                E[i] -= a_T.m[j][i]*b[j];
            }
        }
    }
};



// матрица 3х3x3x3
struct MATR3x3x3x3
{
    double m[3][3][3][3];
};
// матрица 6x6
struct MATR6x6
{
    double m[6][6];
    inline void clear()
    {
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                m[i][j] = 0;
    }
};

// вектор
typedef std::vector<double> Vector;
typedef std::vector<VECTOR3> Vector3;

// c1*x + c2*y -> z
void C1mulVector1PlusC2mulVector2(const double c1, const Vector &x, const double c2, const Vector &y, Vector &z);
// x + c * y -> z
template <typename T>
inline void Vector1PlusCmulVector2(const T &x, const double c, const T &y, T &z)
{
    for (size_t i = 0; i < x.size(); i++)
        z[i] = x[i] + c * y[i];
}
//void Vector1PlusCmulVector2(const Vector &x, const double c, const Vector &y, Vector &z);

// (x, y)
double VectorScalMul(const Vector &x, const Vector &y);
double Vector3ScalMul(const Vector3 &x, const Vector3 &y);
// x=x*c
void VectorMulc(const Vector &x, const double c, Vector &y);
// x->y
template <typename T>
void VectorCopy(const T &x, T &y)
{
    for (size_t i = 0; i < x.size(); i++)
        y[i] = x[i];
}

inline double calcDet3x3(
        const double m00, const double m01, const double m02,
        const double m10, const double m11, const double m12,
        const double m20, const double m21, const double m22)
{
    return
        + m00 * m11 * m22
        + m01 * m12 * m20
        + m02 * m10 * m21
        - m02 * m11 * m20
        - m00 * m12 * m21
        - m01 * m10 * m22;
}

// операции
namespace Operations
{
    void c1M1plusc2M26(const double &c1, const MATR6x6 &M1, const double &c2, const MATR6x6 &M2, MATR6x6 &Mres);
    void MmulM6(const MATR6x6 &M1, const MATR6x6 &M2, MATR6x6 &Mres);
    void cMmulM6(const double &c, const MATR6x6 &M1, const MATR6x6 &M2, MATR6x6 &Mres);
    void Mmulc6(const MATR6x6 &M, const double &c, MATR6x6 &Mres);
    void cIplusM6(const double &c, MATR6x6 &Mres);
    int indexConvertion3x3to6(const int i, const int j);
    void convert6x6to3x3x3x3(const MATR6x6 &D, MATR3x3x3x3 &C);
    void convert3x3x3x3to6x6(const MATR3x3x3x3 &C, MATR6x6 &D);
    void calcNormal(const POINT3 &p0, const POINT3 &p1, const POINT3 &p2, VECTOR3 &normal);
    void sphereArcPoint(const double R, const POINT3 &a, const POINT3 &b, const double fi, POINT3 &p);
    void findPointOnTheLine_1d(const double p1, const double p2, const int N, const double q, const int ind, double &p);
    void findPointOnTheLine_1d_conc(const double p1, const double p2, const int N, const double q, const int ind, double &p);
    void findPointOnTheLine_3d(const POINT3 p1, const POINT3 p2, const int N, const double q, const int ind, POINT3 &p);
    void solvePolynom3(double a, double b, double c, double &r1, double &r2, double &r3);
}   // namespace Elementary_operations
}   // namespace Elementary

/*
class Complex         // класс "Комплексное число"
{
//private:
public:
    double re, im;      // действительная и мнимая части

public:
    // конструкторы

    Complex()
    {
    }

    Complex(double r)
    {
        re = r;
        im = 0;
    }

    Complex(double r, double i)
    {
        re = r;
        im = i;
    }

    Complex(const Complex &c)   // конструктор копирования
    {
        re = c.re;
        im = c.im;
    }

    // деструктор
    ~Complex()
    {
    }

    // остальные функции

    // Модуль комплексного числа
    double abs()
    {
        return ::sqrt(re * re + im * im);
    }

    // оператор присваивания
    Complex& operator = (Complex &c)
    {
        re = c.re;
        im = c.im;

        return (*this);
    }


    // оператор +=
    Complex& operator += (Complex &c)
    {
        re += c.re;
        im += c.im;
        return *this;
    }

    // оператор сложения
    Complex operator + (const Complex &c)
    {
        return Complex(re + c.re, im + c.im);
    }

    // оператор вычитания
    Complex operator - (const Complex &c)
    {
        return Complex(re - c.re, im - c.im);
    }

    // оператор умножения
    Complex operator * (const Complex &c)
    {
        return Complex(re * c.re - im * c.im, re * c.im + im * c.re);
    }

    // оператор деления
    Complex operator / (const Complex &c)
    {
        Complex temp;

        double r = c.re * c.re + c.im * c.im;
        temp.re = (re * c.re + im * c.im) / r;
        temp.im = (im * c.re - re * c.im) / r;

        return temp;
    }

    // оператор отрецания
    Complex operator - ()
    {
        return Complex(-re, -im);
    }

    // корень комплексного числа
    Complex sqrt()
    {
        double s = ::sqrt(re*re + im*im);
        return Complex(::sqrt((s + re) / 2), SIGN(im)*::sqrt((s - re) / 2));
    }

};
*/
#endif  // ELEMENTARY_H
