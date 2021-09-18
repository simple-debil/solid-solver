#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"

#include "integration.h"

const double gauss5_w[5] = { 128. / 225.,
                             (322.+13.*sqrt(70.))/900.,
                             (322.+13.*sqrt(70.))/900.,
                             (322.-13.*sqrt(70.))/900.,
                             (322.-13.*sqrt(70.))/900.,
                           };
const double gauss5_ksi[5] = { 0,
                               -1./3.*sqrt(5.-2.*sqrt(10./7.)),
                                1./3.*sqrt(5.-2.*sqrt(10./7.)),
                               -1./3.*sqrt(5.+2.*sqrt(10./7.)),
                                1./3.*sqrt(5.+2.*sqrt(10./7.)),
                             };
const double gauss3_w[3] = { 8. / 9., 5. / 9., 5. / 9. };
const double gauss3_ksi[3] = { 0, sqrt(0.6), -sqrt(0.6) };
const double gauss2_w[2] = { 1., 1. };
const double gauss2_ksi[2] = { -1. / sqrt(3.), 1. / sqrt(3.) };

namespace Integration
{

Integrator::Integrator()
{
}

Integrator::~Integrator()
{
    release();
}

void Integrator::init1D(const IntegrationType Type, const INTERVAL &Interval)
{
    release();
    dim = 1;
    type = Type;
    interval = Interval;
    detJ = (interval[1] - interval[0]) / 2;
    switch (type)
    {
    case IntegrationType::Gauss5:
        resize(Integration1D_Size_Gauss5);
        for (int k = 0; k < size; k++)
        {
            p1[k] = interval[0] + (interval[1] - interval[0])*(gauss5_ksi[k] + 1) / 2;
            w[k] = gauss5_w[k];
        }
        break;
    case IntegrationType::Gauss3:
        resize(Integration1D_Size_Gauss3);
        for (int k = 0; k < size; k++)
        {
            p1[k] = interval[0] + (interval[1] - interval[0])*(gauss3_ksi[k] + 1) / 2;
            w[k] = gauss3_w[k];
        }
        break;
    case IntegrationType::Gauss2:
        resize(Integration1D_Size_Gauss2);
        for (int k = 0; k < size; k++)
        {
            p1[k] = interval[0] + (interval[1] - interval[0])*(gauss2_ksi[k] + 1) / 2;
            w[k] = gauss2_w[k];
        }
        break;
    }
}

void Integrator::init2D(const IntegrationType Type, const SQUARE &Square)
{
    release();
    dim = 2;
    type = Type;
    square = Square;
    detJ = (square.i[0][1] - square.i[0][0])*(square.i[1][1] - square.i[1][0]) / 4;
    switch (type)
    {
    case IntegrationType::Gauss5:
        resize(Integration2D_Size_Gauss5);
        for (int k = 0; k < size; k++)
        {
            p2[k].x[0] = square.i[0][0] + (square.i[0][1] - square.i[0][0])*(gauss5_ksi[k%5] + 1) / 2;
            p2[k].x[1] = square.i[1][0] + (square.i[1][1] - square.i[1][0])*(gauss5_ksi[k/5] + 1) / 2;
            w[k] = gauss5_w[k%5] * gauss5_w[k/5];
        }
        break;
    case IntegrationType::Gauss3:
        resize(Integration2D_Size_Gauss3);
        for (int k = 0; k < size; k++)
        {
            p2[k].x[0] = square.i[0][0] + (square.i[0][1] - square.i[0][0])*(gauss3_ksi[k%3] + 1) / 2;
            p2[k].x[1] = square.i[1][0] + (square.i[1][1] - square.i[1][0])*(gauss3_ksi[k/3] + 1) / 2;
            w[k] = gauss3_w[k%3] * gauss3_w[k/3];
        }
        break;
    case IntegrationType::Gauss2:
        resize(Integration2D_Size_Gauss2);
        for (int k = 0; k < size; k++)
        {
            p2[k].x[0] = square.i[0][0] + (square.i[0][1] - square.i[0][0])*(gauss2_ksi[k%2] + 1) / 2;
            p2[k].x[1] = square.i[1][0] + (square.i[1][1] - square.i[1][0])*(gauss2_ksi[k/2] + 1) / 2;
            w[k] = gauss2_w[k%2] * gauss2_w[k/2];
        }
        break;
    }
}

void Integrator::init3D(const IntegrationType Type, const CUBE &Cube)
{
    release();
    dim = 3;
    type = Type;
    cube = Cube;
    detJ = (cube.i[0][1] - cube.i[0][0])*(cube.i[1][1] - cube.i[1][0])*(cube.i[2][1] - cube.i[2][0]) / 8;
    switch (type)
    {
    case IntegrationType::Gauss5:
        resize(Integration3D_Size_Gauss5);
        for (int k = 0; k < size; k++)
        {
            p3[k].x[0] = cube.i[0][0] + (cube.i[0][1] - cube.i[0][0])*(gauss5_ksi[k%5] + 1) / 2;
            p3[k].x[1] = cube.i[1][0] + (cube.i[1][1] - cube.i[1][0])*(gauss5_ksi[(k/5)%5] + 1) / 2;
            p3[k].x[2] = cube.i[2][0] + (cube.i[2][1] - cube.i[2][0])*(gauss5_ksi[k/25] + 1) / 2;
            w[k] = gauss5_w[k%5]*gauss5_w[(k/5)%5]*gauss5_w[k/25];
        }
        break;
    case IntegrationType::Gauss3:
        resize(Integration3D_Size_Gauss3);
        for (int k = 0; k < size; k++)
        {
            p3[k].x[0] = cube.i[0][0] + (cube.i[0][1] - cube.i[0][0])*(gauss3_ksi[k%3] + 1) / 2;
            p3[k].x[1] = cube.i[1][0] + (cube.i[1][1] - cube.i[1][0])*(gauss3_ksi[(k/3)%3] + 1) / 2;
            p3[k].x[2] = cube.i[2][0] + (cube.i[2][1] - cube.i[2][0])*(gauss3_ksi[k/9] + 1) / 2;
            w[k] = gauss3_w[k%3]*gauss3_w[(k/3)%3]*gauss3_w[k/9];
        }
        break;
    case IntegrationType::Gauss2:
        resize(Integration3D_Size_Gauss2);
        for (int k = 0; k < size; k++)
        {
            p3[k].x[0] = cube.i[0][0] + (cube.i[0][1] - cube.i[0][0])*(gauss2_ksi[k%2] + 1) / 2;
            p3[k].x[1] = cube.i[1][0] + (cube.i[1][1] - cube.i[1][0])*(gauss2_ksi[(k/2)%2] + 1) / 2;
            p3[k].x[2] = cube.i[2][0] + (cube.i[2][1] - cube.i[2][0])*(gauss2_ksi[k/4] + 1) / 2;
            w[k] = gauss2_w[k%2]*gauss2_w[(k/2)%2]*gauss2_w[k/4];
        }
        break;
    }
}

void Integrator::release()
{
    resize(0);
}

void Integrator::resize(const int Size)
{
    if(dim == 1) p1.resize(Size);
    if(dim == 2) p2.resize(Size);
    if(dim == 3) p3.resize(Size);
    value.resize(Size);
    w.resize(Size);
    size = Size;
}

}   // namespace Integration
