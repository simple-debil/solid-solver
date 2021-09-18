/* --------------------------------------------------------- */
// ЧИСЛЕННОЕ ИНТЕГРИРОВАНИЕ
/* --------------------------------------------------------- */

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "elementary.h"

// взятие определенного интеграла численно
// интегрирование по шаблонному интервалу/квадрату/кубу
namespace Integration
{
using namespace Elementary;

const int Integration3D_Size_Gauss2 = 8;
const int Integration3D_Size_Gauss3 = 27;
const int Integration3D_Size_Gauss5 = 125;
const int Integration3D_Size_Max = 125;

const int Integration2D_Size_Gauss2 = 4;
const int Integration2D_Size_Gauss3 = 9;
const int Integration2D_Size_Gauss5 = 25;
const int Integration2D_Size_Max = 25;

const int Integration1D_Size_Gauss2 = 2;
const int Integration1D_Size_Gauss3 = 3;
const int Integration1D_Size_Gauss5 = 5;
const int Integration1D_Size_Max = 5;

// метод интегрирования
enum class IntegrationType
{
    Gauss2 = 0,
    Gauss3 = 1,
    Gauss5 = 2,
};

// интегратор
class Integrator
{
public:
    Integrator();
    ~Integrator();
    void init1D(const IntegrationType Type = IntegrationType::Gauss3, const INTERVAL &Interval = INTERVAL{-1, 1});
    void init2D(const IntegrationType Type = IntegrationType::Gauss3, const SQUARE &Square = SQUARE{{{-1, 1}, {-1, 1}}});
    void init3D(const IntegrationType Type = IntegrationType::Gauss3, const CUBE &Cube = CUBE{{{-1, 1}, {-1, 1}, {-1, 1}}});
    void release();
    inline double integrate()
    {
        double E = 0;
        for (int i = 0; i < size; i++)
            E += w[i] * value[i];
        return E * detJ;
    }
    std::vector<POINT1> p1;
    std::vector<POINT2> p2;
    std::vector<POINT3> p3;       // точки
    std::vector<double> value;    // значения
    int size = 0;               // размер value и (p1|p2|p3)
public://private:
    void resize(const int Size);
    int dim = 0;                // размерность
    IntegrationType type;       // метод
    INTERVAL interval;          // интервал
    SQUARE square;              // квадрат
    CUBE cube;                  // куб
    std::vector<double> w;
    double detJ;
};
}   // namespace Integration

#endif  // INTEGRATION_H
