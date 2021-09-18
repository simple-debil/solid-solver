#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"

#include "slausolving.h"
#include "integration.h"
#include "fem.h"

#include "interpolation.h"

using namespace SlauSolving;
using namespace Integration;
using namespace Fem;

namespace Interpolation
{
double Interpolant1D_Lagrange1::G[2][2];
double Interpolant1D_Lagrange3::G[4][4];
double Interpolant1D_Hermite3::G[4][4];
double Interpolant2D_Hermite3::G[16][16];
double Interpolant2D_Lagrange3::G[16][16];

int SIZE_X = 12;
int SIZE_Y = 14;	// размер одного символа
int LEN = 11;		// длина числа в символах
int K = 9;

#define PAINT(x_, y_, r, g, b) {if(x_ >= 0 && x_ < bmp.x && y_ >= 0 && y_ < bmp.y) bmp.c[(y_)*bmp.x + (x_)] = RGB32(0, r, g, b);}
#define NORM(a)	((a) < 0 ? (0) : ((a) > 255 ? 255 : (a)))

// вычисление цвета по значению скаляра
void calc_rgb(const int colorMode, const double value, const double valueMin, const double valueNull, const double valueMax, int &r, int &g, int &b)
{
    if(colorMode == 3)	// 3 цвета
    {
        if(value >= valueNull)
        {
            r = NORM((int)((value - valueNull) / (valueMax - valueNull) * (255)));
            b = 0;
            g = (255 - r)/2;
        }
        else
        {
            r = 0;
            b = NORM((int)((valueNull - value) / (valueNull - valueMin) * (255)));
            g = (255 - b)/2;
        }
    }else
    if(colorMode == 2)	// 2 цвета
    {
        r = NORM((int)((value - valueMin) / (valueMax - valueMin) * 255));
        g = 0;
        b = 255 - r;
    }else
    if(colorMode == 1)	// 1 цвет
    {
        double delta = ((value - valueMin) / (valueMax - valueMin));
        r = NORM((int)(pow(delta, 1) * 255));       // k_mono = 1
        r = 255 - r;
        g = r;
        b = r;
    }
}
// сохраняет 32-битное изображение в формате .bmp
void bmp_save(const char *file_name, const BMP_INF &inf)
{
    int size1 = 14, size2 = 12;
    BMP_INF1 inf1;
    BMP_INF2 inf2;
    FILE *f = fopen(file_name,"wb");
    inf1.tip = 0x4D42;
    inf1.sizef = size1 + size2 + 4 * inf.x * inf.y;
    inf1.rez1 = 0;
    inf1.rez2 = 0;
    inf1.gotobits = size1 + size2;
    inf2.lenstruct = size2;
    inf2.x = (unsigned short) inf.x;
    inf2.y = (unsigned short) inf.y;
    inf2.planes = 1;
    inf2.bits = 32;
    fseek(f, 0, SEEK_SET); fwrite(&inf1.tip, 1, 2, f);
    fseek(f, 2, SEEK_SET); fwrite(&inf1.sizef, 1, 4, f);
    fseek(f, 6, SEEK_SET); fwrite(&inf1.rez1, 1, 2, f);
    fseek(f, 8, SEEK_SET); fwrite(&inf1.rez2, 1, 2, f);
    fseek(f, 10, SEEK_SET); fwrite(&inf1.gotobits, 1, 4, f);
    fseek(f, 14, SEEK_SET); fwrite(&inf2.lenstruct, 1, 4, f);
    fseek(f, 18, SEEK_SET); fwrite(&inf2.x, 1, 2, f);
    fseek(f, 20, SEEK_SET); fwrite(&inf2.y, 1, 2, f);
    fseek(f, 22, SEEK_SET); fwrite(&inf2.planes, 1, 2, f);
    fseek(f, 24, SEEK_SET); fwrite(&inf2.bits, 1, 2, f);
    fseek(f, 26, SEEK_SET); fwrite(inf.c, 1, 4 * inf.x * inf.y, f);
    fclose(f);
};
// загружает 32-битное изображение в формате .bmp
void bmp_load(const char *file_name, BMP_INF &inf)
{
    BMP_INF2 inf2;
    FILE *f = fopen(file_name,"rb");
    fseek(f, 18, SEEK_SET); fread(&inf2.x, 1, 2, f);
    fseek(f, 20, SEEK_SET); fread(&inf2.y, 1, 2, f);
    inf.x = inf2.x;
    inf.y = inf2.y;
    inf.c = new unsigned int[inf.x * inf.y];
    fseek(f, 26, SEEK_SET); fread(inf.c, 1, 4 * inf.x * inf.y, f);
    fclose(f);
};
// выводит строку s символами symbols в картинку out с координатами x, y
void print(const char *s, int x0, int y0, const int colorMode, const BMP_INF &symbols, const BMP_INF &out)
{
    unsigned int voidc = RGB32(0, 255, 255, 255);
    y0 -= SIZE_Y/2;
    int c, color_real;
    if(colorMode != 1)
    {
        color_real = RGB32(0, 0, 0, 0);
    }
    else
    {
        //color_real = ~out.c[y*out.x + x];
    }

    for(int i = 0; (c = s[i]) != 0 && i < LEN; i++)
    {
        if(c >= '0' && c <= '9')
            c -= '0';
        else
        if(c == '-')
            c = 10;
        else
        if(c == '.')
            c = 11;
        else
        if(c == 'e')
            c = 12;
        else
        if(c == '+')
            c = 13;
        else
        if(c == ' ')
            goto next;
        else
        {
            //printf("<%c>\n", c);
            return;
        }
        for(int y = y0; y < y0 + SIZE_Y; y++)
        for(int x = x0; x < x0 + SIZE_X; x++)
        {
            unsigned int color = symbols.c[(y-y0)*symbols.x + (x-x0) + c*SIZE_X];
            if(color != voidc)
            {
                if(colorMode != 1)
                {
                    out.c[y*out.x + x] = color_real;
                }
                else
                {
                    int c = (out.c[y*out.x + x]%256 + 128)%256;
                    //out.c[y*out.x + x] = out.c[y*out.x + x]^color_real;
                    if(c >= 128)
                        out.c[y*out.x + x] = RGB32(0, 255, 255, 255);// RGB32(0, c, c, c);
                    else
                        out.c[y*out.x + x] = RGB32(0, 0, 0, 0);
                }
            }
        }
next:;
        x0 += SIZE_X;
    }
}



void Interpolant1D_base::init(const Grid::GridRectangleRegular1D &set_grid, const double set_alpha)
{
    gr = set_grid;
    alpha = set_alpha;
    l = gr.rect[1] - gr.rect[0];
    h = l / gr.N;
    m.clear();
    res.clear();
}
void Interpolant1D_base::release()
{
    m.clear();
    res.clear();
}
void Interpolant1D_base::addPoint(const POINT1 &p, const double value)
{
    POINT1_VALUE p1vEl;
    int i;
    findFe(p,
           i);
    p1vEl.p = p;
    p1vEl.value = value;
    p1vEl.feIndex = i;
    m.push_back(p1vEl);
}
void Interpolant1D_base::findFe(const POINT1 &p, int &i) const
{
    i = (int)((p - gr.rect[0]) / l * gr.N);  // 0..gr.N
}

int POINT1_VALUE_cmp(const void* x1, const void* x2)
{
    POINT1_VALUE *px1 = (POINT1_VALUE *)x1;
    POINT1_VALUE *px2 = (POINT1_VALUE *)x2;
    if (px1->feIndex < px2->feIndex) return -1;
    if (px1->feIndex > px2->feIndex) return  1;
    return 0;
}
int POINT1_VALUE_UNREGULAR_cmp(const void* x1, const void* x2)
{
    POINT1_VALUE_UNREGULAR *px1 = (POINT1_VALUE_UNREGULAR *)x1;
    POINT1_VALUE_UNREGULAR *px2 = (POINT1_VALUE_UNREGULAR *)x2;
    if (px1->p < px2->p) return -1;
    if (px1->p > px2->p) return  1;
    return 0;
}


void Interpolant1D_Lagrange1::buildInterpolant()
{
    buildLagrange1DLocalMatrix();
    int matrixSize = (gr.N + 1); // размер матрицы СЛАУ
    // сортировка точек КЭ в которые они попали
    qsort(m.data(), m.size(), sizeof(POINT1_VALUE), POINT1_VALUE_cmp);
    SlauSolving::SSCM matrix;                           // матрица
    SlauSolving::SSCMBulder mBulder;                    // сборщик
    SlauSolving::SSCMPreconditioner_base *mPreconditioner;    // предобусловливатель
    SlauSolving::SSCMSolver_base *mSolver;                    // решатель СЛАУ
    Vector b;                                           // правая часть
    SlauSolving::SolverParameters ssp;
    ssp.type = SlauSolving::SolverType::Direct;
    ssp.preconditioning = SlauSolving::Preconditioning::Profile_LLT;

    matrix.init();
    b.resize(matrixSize);
    res.resize(matrixSize);
    for (int i = 0; i < matrixSize; i++)
    {
        b[i] = 0;
        res[i] = 0;
    }
    mBulder.setMatrixSize(matrixSize);
    mBulder.start();

    // занесение элементов локальных матриц в глобальную
    size_t gcount = 0;
    size_t count;
    for(int i = 0; i < gr.N; i++)
    {                   // выбран конечный элемент i
        int gind[2];	// глобальные номера локальных функций конечного элемента i
        int feIndex = i;
        POINT1 p0;
        p0 = gr.rect[0] + i * h;// одномерный конечный элемент задается одной вершиной с локальным номером 0
        for(int lvi = 0; lvi < 2; lvi++)
        {
            gind[lvi] = getBasFunGlobaIndex(feIndex, lvi);
        }
        // подсчет суммы базисных функций от точек на этом конечном элементе
        // и сборка локальной матрицы конечного элемента i
        for(int lvi = 0; lvi < 2; lvi++)
        for(int lvj = 0; lvj < 2; lvj++)	// lvi,lvj - локальные номера базисных функций (координаты в G)
        {
            // локальные номера вершин lvi, lvj соответствуют глобальным номерам вершин gind[lvi], gind[lvj]
            double sum = 0;
            count = gcount;
            for (;;)
            {
                if(count >= m.size() || m[count].feIndex > feIndex) break;
                sum +=  Fem::lagrange1_1D(p0, h, m[count].p, lvi, Fem::dif_NULL1) *
                        Fem::lagrange1_1D(p0, h, m[count].p, lvj, Fem::dif_NULL1);
                count++;
            }
            if(gind[lvi] >= gind[lvj])
                mBulder.addElement_not_null(sum + alpha*G[lvi][lvj], gind[lvi], gind[lvj]);
        }
        // сборка правой части
        for(int lvi = 0; lvi < 2; lvi++)	// lvi - локальные номера базисных функций
        {
            // локальные номера вершин lvi соответствуют глобальным номерам вершин gind[lvi]
            double sum = 0;
            count = gcount;
            for(;;)
            {
                if (count >= m.size() || m[count].feIndex > feIndex) break;
                sum += Fem::lagrange1_1D(p0, h, m[count].p, lvi, Fem::dif_NULL1) * m[count].value;
                count++;
            }
            b[gind[lvi]] += sum;
        }
        gcount = count;	// переходим к точкам следующего КЭ
    }
    mBulder.complete(matrix);
    double residual;
    double relativeResidual;
    int iterations;
    double time;
    double timePreconditioner;
    mPreconditioner = SlauSolving::SSCMPreconditioner_base::gen(ssp.preconditioning);
    mSolver = SlauSolving::SSCMSolver_base::gen(ssp.type);
    mPreconditioner->bulid(matrix, 0, timePreconditioner);
    mSolver->init(matrix.getMatrixSize());
    mSolver->solve(matrix, b,
                  mPreconditioner, res, ssp,
                  res, residual, relativeResidual, iterations, time);
    matrix.release();
    mBulder.release();
    delete mPreconditioner;
    delete mSolver;
    b.clear();
    valueMin = +1.e+200;
    valueMax = -1.e+200;
    for (size_t k = 0; k < m.size(); k++)
    {
        double value = m[k].value;
        if (value < valueMin) valueMin = value;
        if (value > valueMax) valueMax = value;
    }
    valueNull = (valueMin + valueMax) / 2;
}
void Interpolant1D_Lagrange1::getNodesCoordinates(std::vector<POINT1> &coordinates) const
{
    coordinates.clear();
    for(int i = 0; i <= gr.N; i++)
    {
        POINT1 p0;      // координата узла i
        p0 = gr.rect[0] + i * h;
        coordinates.push_back(p0);
    }
}
void Interpolant1D_Lagrange1::buildInterpolantByAllNodes(const Vector &nodeValue)
{
    int matrixSize = (gr.N + 1); // размер матрицы СЛАУ
    res.resize(matrixSize);
    int count = 0;
    valueMin = +1.e+200;
    valueMax = -1.e+200;
    for(int i = 0; i <= gr.N; i++)
    {
        double value = nodeValue[count];
        if (value < valueMin) valueMin = value;
        if (value > valueMax) valueMax = value;
        res[count] = value;
        count++;
    }
    valueNull = (valueMin + valueMax) / 2;
}
double Interpolant1D_Lagrange1::fun(const POINT1 &p) const
{
    return difFun(p, Fem::dif_NULL1);
}
double Interpolant1D_Lagrange1::difFun(const POINT1 &p, const DIF_STATE1 &dif) const
{
    int i;
    findFe(p,
           i);
    if(i < 0)
        i = 0;
    if(i >= gr.N)
        i = gr.N - 1;
    int feIndex = i;
    double E = 0;
    POINT1 p0;
    p0 = gr.rect[0] + i * h;// конечный элемент задается одной вершиной с локальным номером 0
    for(int lvi = 0; lvi < 2; lvi++)
    {
        E += res[getBasFunGlobaIndex(feIndex, lvi)] * Fem::lagrange1_1D(p0, h, p, lvi, dif);
    }
    return E;
}
int Interpolant1D_Lagrange1::getBasFunGlobaIndex(const int feIndex, const int localBusFunIndex) const
{
    return feIndex + localBusFunIndex;
}
void Interpolant1D_Lagrange1::buildLagrange1DLocalMatrix()
{
    G[0][0] = 1. / h;
    G[0][1] = -1. / h;
    G[1][0] = -1. / h;
    G[1][1] = 1. / h;
}

void Interpolant1D_Lagrange1Unregular::addPoint(const POINT1 &p, const double value)
{
    POINT1_VALUE_UNREGULAR p1vEl;
    p1vEl.p = p;
    p1vEl.value = value;
    mu.push_back(p1vEl);
}
void Interpolant1D_Lagrange1Unregular::findFe(const POINT1 &p,
                                              int &i)const
{
    int a, b, c;// индексы задают КЭ
    a = 0;
    b = mu.size() - 2;
    for(;;)
    {
        nextiter:;
        c = (a + b) / 2;
        // КЭ: [mu[c].p, mu[c + 1].p]
        if(p < mu[c].p)
        {
            b = c - 1;
            if(a > b)
            {
                i = a;
                return;   // возвращается ближайший КЭ
            }
            goto nextiter;
        }
        if(p > mu[c + 1].p)
        {
            a = c + 1;
            if(a > b)
            {
                i = b;
                return;   // возвращается ближайший КЭ
            }
            goto nextiter;
        }
        i = c;
        return;
    }
    i = 0;
    return;
}
void Interpolant1D_Lagrange1Unregular::buildInterpolant()
{
    // сортировка точек
    qsort(mu.data(), mu.size(), sizeof(POINT1_VALUE_UNREGULAR), POINT1_VALUE_UNREGULAR_cmp);
}
void Interpolant1D_Lagrange1Unregular::getNodesCoordinates(std::vector<POINT1> &) const
{
}
void Interpolant1D_Lagrange1Unregular::buildInterpolantByAllNodes(const Vector &Value)
{
}
double Interpolant1D_Lagrange1Unregular::fun(const POINT1 &p) const
{
    return difFun(p, Fem::dif_NULL1);
}
double Interpolant1D_Lagrange1Unregular::difFun(const POINT1 &p, const DIF_STATE1 &dif) const
{
    int feIndex;
    findFe(p,
           feIndex);
    double p1 = mu[feIndex].p;
    double value1 = mu[feIndex].value;
    double p2 = mu[feIndex + 1].p;
    double value2 = mu[feIndex + 1].value;
    // f(x) = a*x + b
    // f(p1) = value1
    // f(p2) = value2
    // a*p1 + b = value1 => b = value1 - a*p1
    // a*p2 + b = value2 => a*p2 + value1 - a*p1 = value2 => a*(p2 - p1) = value2 - value1

    // a = (value2 - value1) / (p2 - p1)
    // b = (value1*p2 - value2*p1) / (p2 - p1)
    double a = (value2 - value1) / (p2 - p1);
    double b = (value1*p2 - value2*p1) / (p2 - p1);
    if(dif == 0)
        return a*p + b;
    if(dif == 1)
        return a;
    return 0;
}


void Interpolant1D_Hermite3::buildInterpolant()
{
    buildHermite1DLocalMatrix();
    int matrixSize = (gr.N + 1)*2; // размер матрицы СЛАУ
    // сортировка точек по КЭ в которые они попали
    qsort(m.data(), m.size(), sizeof(POINT1_VALUE), POINT1_VALUE_cmp);
    SlauSolving::SSCM matrix;                           // матрица
    SlauSolving::SSCMBulder mBulder;                    // сборщик
    SlauSolving::SSCMPreconditioner_base *mPreconditioner;    // предобусловливатель
    SlauSolving::SSCMSolver_base *mSolver;                    // решатель СЛАУ
    Vector b;                                           // правая часть
    SlauSolving::SolverParameters ssp;
    ssp.type = SlauSolving::SolverType::Direct;
    ssp.preconditioning = SlauSolving::Preconditioning::Profile_LLT;//SlauPreconditioning_none;

    matrix.init();
    b.resize(matrixSize);
    res.resize(matrixSize);
    for (int i = 0; i < matrixSize; i++)
    {
        b[i] = 0;
        res[i] = 0;
    }
    mBulder.setMatrixSize(matrixSize);
    mBulder.start();

    // занесение элементов локальных матриц в глобальную
    size_t gcount = 0;
    size_t count;
    for(int i = 0; i < gr.N; i++)
    {					// выбран конечный элемент i
        int gind[4];	// глобальные номера локальных функций конечного элемента i
        int feIndex = i;
        POINT1 p0;
        p0 = gr.rect[0] + i * h;// одномерный конечный элемент задается одной вершиной с локальным номером 0

        for(int lvi = 0; lvi < 4; lvi++)
        {
            gind[lvi] = getBasFunGlobaIndex(feIndex, lvi);
        }
        // подсчет суммы базисных функций от точек на этом конечном элементе
        // и сборка локальной матрицы конечного элемента i
        for(int lvi = 0; lvi < 4; lvi++)
        for(int lvj = 0; lvj < 4; lvj++)	// lvi,lvj - локальные номера базисных функций (координаты в G)
        {
            // локальные номера вершин lvi, lvj соответствуют глобальным номерам вершин gind[lvi], gind[lvj]
            double sum = 0;
            count = gcount;
            for (;;)
            {
                if(count >= m.size() || m[count].feIndex > feIndex) break;
                sum +=  Fem::hermite_1D(p0, h, m[count].p, lvi, Fem::dif_NULL1) *
                        Fem::hermite_1D(p0, h, m[count].p, lvj, Fem::dif_NULL1);
                count++;
            }
            if(gind[lvi] >= gind[lvj])
                mBulder.addElement_not_null(sum + alpha*G[lvi][lvj], gind[lvi], gind[lvj]);
        }
        // сборка правой части
        for(int lvi = 0; lvi < 4; lvi++)	// lvi - локальные номера базисных функций
        {
            // локальные номера вершин lvi соответствуют глобальным номерам вершин gind[lvi]
            double sum = 0;
            count = gcount;
            for(;;)
            {
                if (count >= m.size() || m[count].feIndex > feIndex) break;
                sum += Fem::hermite_1D(p0, h, m[count].p, lvi, Fem::dif_NULL1) * m[count].value;
                count++;
            }
            b[gind[lvi]] += sum;
        }
        gcount = count;	// переходим к точкам следующего КЭ
    }
    mBulder.complete(matrix);
    double residual;
    double relativeResidual;
    int iterations;
    double time;
    double timePreconditioner;
    mPreconditioner = SlauSolving::SSCMPreconditioner_base::gen(ssp.preconditioning);
    mSolver = SlauSolving::SSCMSolver_base::gen(ssp.type);
    mPreconditioner->bulid(matrix, 0, timePreconditioner);
    mSolver->init(matrix.getMatrixSize());
    mSolver->solve(matrix, b,
                  mPreconditioner, res, ssp,
                  res, residual, relativeResidual, iterations, time);
    matrix.release();
    mBulder.release();
    delete mPreconditioner;
    delete mSolver;
    b.clear();
    valueMin = +1.e+200;
    valueMax = -1.e+200;
    for (size_t k = 0; k < m.size(); k++)
    {
        double value = m[k].value;
        if(value < valueMin) valueMin = value;
        if(value > valueMax) valueMax = value;
    }
    valueNull = (valueMin + valueMax) / 2;
}
void Interpolant1D_Hermite3::getNodesCoordinates(std::vector<POINT1> &coordinates) const
{
    coordinates.clear();
    for(int i = 0; i <= gr.N; i++)
    {
        POINT1 p0;      // координата узла i
        p0 = gr.rect[0] + i * h;
        coordinates.push_back(p0);
    }
}
void Interpolant1D_Hermite3::buildInterpolantByAllNodes(const Vector &nodeValue)
{
    buildHermite1DLocalMatrix();
    int matrixSize = (gr.N + 1)*2; // размер матрицы СЛАУ
    SlauSolving::SSCM matrix;                           // матрица
    SlauSolving::SSCMBulder mBulder;                    // сборщик
    SlauSolving::SSCMPreconditioner_base *mPreconditioner;    // предобусловливатель
    SlauSolving::SSCMSolver_base *mSolver;                    // решатель СЛАУ
    Vector b;                                           // правая часть
    SlauSolving::SolverParameters ssp;
    ssp.type = SlauSolving::SolverType::Direct;
    ssp.preconditioning = SlauSolving::Preconditioning::Profile_LLT;

    matrix.init();
    b.resize(matrixSize);
    res.resize(matrixSize);
    for (int i = 0; i < matrixSize; i++)
    {
        b[i] = 0;
        res[i] = 0;
    }
    mBulder.setMatrixSize(matrixSize);
    mBulder.start();
    // занесение элементов локальных матриц в глобальную
    for(int i = 0; i < gr.N; i++)
    {					// выбран конечный элемент i
        int gind[4];	// глобальные номера локальных функций конечного элемента i
        int feIndex = i;
        POINT1 p0;
        p0 = gr.rect[0] + i * h;// КЭ задается одной вершиной с локальным номером 0
        for(int lvi = 0; lvi < 4; lvi++)
        {
            gind[lvi] = getBasFunGlobaIndex(feIndex, lvi);
        }
        // узлы конечного элемента i
        POINT1 feVertex[2];
        feVertex[0] = p0;
        feVertex[1] = p0 + h;
        // значение интерполируемой ф-и в этих узлах
        double value[2];
        int nodeIndex0 = i;
        value[0] = nodeValue[nodeIndex0];
        value[1] = nodeValue[nodeIndex0 + 1];

        // подсчет суммы базисных функций от точек на этом конечном элементе
        // и сборка локальной матрицы конечного элемента i
        for(int lvi = 0; lvi < 4; lvi++)
        for(int lvj = 0; lvj < 4; lvj++)// lvi,lvj - локальные номера базисных функций (координаты в G)
        {
            // локальные номера вершин li, lj соответствуют глобальным номерам вершин gind[li], gind[lj]
            double sum = 0;
            sum +=  Fem::hermite_1D(p0, h, feVertex[0], lvi, Fem::dif_NULL1) *
                    Fem::hermite_1D(p0, h, feVertex[0], lvj, Fem::dif_NULL1);
            if(i == gr.N - 1)
            {
                sum +=  Fem::hermite_1D(p0, h, feVertex[1], lvi, Fem::dif_NULL1) *
                        Fem::hermite_1D(p0, h, feVertex[1], lvj, Fem::dif_NULL1);
            }
            if(gind[lvi] >= gind[lvj])
                mBulder.addElement_not_null(sum + alpha*G[lvi][lvj], gind[lvi], gind[lvj]);
        }
        // сборка правой части
        for(int lvi = 0; lvi < 4; lvi++)	// lvi - локальные номера базисных функций
        {
            // локальные номера вершин lvi соответствуют глобальным номерам вершин gind[lvi]
            double sum = 0;
            sum += Fem::hermite_1D(p0, h, feVertex[0], lvi, Fem::dif_NULL1) * value[0];
            if(i == gr.N - 1)
            {
                sum += Fem::hermite_1D(p0, h, feVertex[1], lvi, Fem::dif_NULL1) * value[1];
            }
            b[gind[lvi]] += sum;
        }
    }
    mBulder.complete(matrix);
    double residual;
    double relativeResidual;
    int iterations;
    double time;
    double timePreconditioner;
    mPreconditioner = SlauSolving::SSCMPreconditioner_base::gen(ssp.preconditioning);
    mSolver = SlauSolving::SSCMSolver_base::gen(ssp.type);
    mPreconditioner->bulid(matrix, 0, timePreconditioner);
    mSolver->init(matrix.getMatrixSize());
    mSolver->solve(matrix, b,
                  mPreconditioner, res, ssp,
                  res, residual, relativeResidual, iterations, time);
    matrix.release();
    mBulder.release();
    delete mPreconditioner;
    delete mSolver;
    b.clear();
    valueMin = +1.e+200;
    valueMax = -1.e+200;
    for (size_t k = 0; k < nodeValue.size(); k++)
    {
        double value = nodeValue[k];
        if(value < valueMin) valueMin = value;
        if(value > valueMax) valueMax = value;
    }
    valueNull = (valueMin + valueMax) / 2;
}
double Interpolant1D_Hermite3::fun(const POINT1 &p) const
{
    return difFun(p, Fem::dif_NULL1);
}
double Interpolant1D_Hermite3::difFun(const POINT1 &p, const DIF_STATE1 &dif) const
{
    int i;
    findFe(p,
           i);
    if(i < 0)
        i = 0;
    if(i >= gr.N)
        i = gr.N - 1;
    int feIndex = i;
    double E = 0;
    POINT1 p0;
    p0 = gr.rect[0] + i * h;// КЭ задается одной вершиной с локальным номером 0
    for(int lvi = 0; lvi < 4; lvi++)
    {
        E += res[getBasFunGlobaIndex(feIndex, lvi)] * Fem::hermite_1D(p0, h, p, lvi, dif);
    }
    return E;
}
int Interpolant1D_Hermite3::getBasFunGlobaIndex(const int feIndex, const int localBusFunIndex) const
{
    int i = feIndex;
    int basFun0Ind = i * 2; // глобальный индекс базисной ф-и с локальным индексом 0
    return basFun0Ind + localBusFunIndex;
}
void Interpolant1D_Hermite3::buildHermite1DLocalMatrix()
{
    // локальные матрицы для одномерных эрмитовых базисных функций 3-го порядка
     G[0][0] = 36;
     G[0][1] = 3 * h;
     G[0][2] = -36;
     G[0][3] = 3 * h;
      G[1][0] = 3 * h;
      G[1][1] = 4 * h*h;
      G[1][2] = -3 * h;
      G[1][3] = -h*h;
     G[2][0] = -36;
     G[2][1] = -3 * h;
     G[2][2] = 36;
     G[2][3] = -3 * h;
      G[3][0] = 3 * h;
      G[3][1] = -h*h;
      G[3][2] = -3 * h;
      G[3][3] = 4 * h*h;
      for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
          G[i][j] *= 1.L / (30 * h);

}

void Interpolant1D_Lagrange3::buildInterpolant()
{
    buildLagrange1DLocalMatrix();
    int matrixSize = (gr.N * 3 + 1); // размер матрицы СЛАУ
    // сортировка точек КЭ в которые они попали
    qsort(m.data(), m.size(), sizeof(POINT1_VALUE), POINT1_VALUE_cmp);
    SlauSolving::SSCM matrix;                           // матрица
    SlauSolving::SSCMBulder mBulder;                    // сборщик
    SlauSolving::SSCMPreconditioner_base *mPreconditioner;    // предобусловливатель
    SlauSolving::SSCMSolver_base *mSolver;                    // решатель СЛАУ
    Vector b;                                           // правая часть
    SlauSolving::SolverParameters ssp;
    ssp.type = SlauSolving::SolverType::Direct;
    ssp.preconditioning = SlauSolving::Preconditioning::Profile_LLT;

    matrix.init();
    b.resize(matrixSize);
    res.resize(matrixSize);
    for (int i = 0; i < matrixSize; i++)
    {
        b[i] = 0;
        res[i] = 0;
    }
    mBulder.setMatrixSize(matrixSize);
    mBulder.start();

    // занесение элементов локальных матриц в глобальную
    size_t gcount = 0;
    size_t count;
    for(int i = 0; i < gr.N; i++)
    {                   // выбран конечный элемент i
        int gind[4];	// глобальные номера локальных функций конечного элемента i
        int feIndex = i;
        POINT1 p0;
        p0 = gr.rect[0] + i * h;// одномерный конечный элемент задается одной вершиной с локальным номером 0
        for(int lvi = 0; lvi < 4; lvi++)
        {
            gind[lvi] = getBasFunGlobaIndex(feIndex, lvi);
        }
        // подсчет суммы базисных функций от точек на этом конечном элементе
        // и сборка локальной матрицы конечного элемента i
        for(int lvi = 0; lvi < 4; lvi++)
        for(int lvj = 0; lvj < 4; lvj++)	// lvi,lvj - локальные номера базисных функций (координаты в G)
        {
            // локальные номера вершин lvi, lvj соответствуют глобальным номерам вершин gind[lvi], gind[lvj]
            double sum = 0;
            count = gcount;
            for (;;)
            {
                if(count >= m.size() || m[count].feIndex > feIndex) break;
                sum +=  Fem::lagrange3_1D(p0, h, m[count].p, lvi, Fem::dif_NULL1) *
                        Fem::lagrange3_1D(p0, h, m[count].p, lvj, Fem::dif_NULL1);
                count++;
            }
            if(gind[lvi] >= gind[lvj])
                mBulder.addElement_not_null(sum + alpha*G[lvi][lvj], gind[lvi], gind[lvj]);
        }
        // сборка правой части
        for(int lvi = 0; lvi < 4; lvi++)	// lvi - локальные номера базисных функций
        {
            // локальные номера вершин lvi соответствуют глобальным номерам вершин gind[lvi]
            double sum = 0;
            count = gcount;
            for(;;)
            {
                if (count >= m.size() || m[count].feIndex > feIndex) break;
                sum += Fem::lagrange3_1D(p0, h, m[count].p, lvi, Fem::dif_NULL1) * m[count].value;
                count++;
            }
            b[gind[lvi]] += sum;
        }
        gcount = count;	// переходим к точкам следующего КЭ
    }
    mBulder.complete(matrix);
    double residual;
    double relativeResidual;
    int iterations;
    double time;
    double timePreconditioner;
    mPreconditioner = SlauSolving::SSCMPreconditioner_base::gen(ssp.preconditioning);
    mSolver = SlauSolving::SSCMSolver_base::gen(ssp.type);
    mPreconditioner->bulid(matrix, 0, timePreconditioner);
    mSolver->init(matrix.getMatrixSize());
    mSolver->solve(matrix, b,
                  mPreconditioner, res, ssp,
                  res, residual, relativeResidual, iterations, time);
    matrix.release();
    mBulder.release();
    delete mPreconditioner;
    delete mSolver;
    b.clear();
    valueMin = +1.e+200;
    valueMax = -1.e+200;
    for (size_t k = 0; k < m.size(); k++)
    {
        double value = m[k].value;
        if (value < valueMin) valueMin = value;
        if (value > valueMax) valueMax = value;
    }
    valueNull = (valueMin + valueMax) / 2;
}
void Interpolant1D_Lagrange3::getNodesCoordinates(std::vector<POINT1> &coordinates) const
{
    coordinates.clear();
    for(int i = 0; i <= gr.N*3; i++)
    {
        POINT1 p0;      // координата узла i
        p0 = gr.rect[0] + i * h / 3;
        coordinates.push_back(p0);
    }
}
void Interpolant1D_Lagrange3::buildInterpolantByAllNodes(const Vector &nodeValue)
{
    int matrixSize = (gr.N * 3 + 1); // размер матрицы СЛАУ
    res.resize(matrixSize);
    int count = 0;
    valueMin = +1.e+200;
    valueMax = -1.e+200;
    for(int i = 0; i <= gr.N*3; i++)
    {
        double value = nodeValue[count];
        if (value < valueMin) valueMin = value;
        if (value > valueMax) valueMax = value;
        res[count] = value;
        count++;
    }
    valueNull = (valueMin + valueMax) / 2;
}
double Interpolant1D_Lagrange3::fun(const POINT1 &p) const
{
    return difFun(p, Fem::dif_NULL1);
}
double Interpolant1D_Lagrange3::difFun(const POINT1 &p, const DIF_STATE1 &dif) const
{
    int i;
    findFe(p,
           i);
    /*if(i < 0 || i >= gr.N)
    {
        // выход за пределы сетки
        return 0;
    }*/
    if(i < 0)
        i = 0;
    if(i >= gr.N)
        i = gr.N - 1;
    int feIndex = i;
    double E = 0;
    POINT1 p0;
    p0 = gr.rect[0] + i * h;// конечный элемент задается одной вершиной с локальным номером 0
    for(int lvi = 0; lvi < 4; lvi++)
    {
        E += res[getBasFunGlobaIndex(feIndex, lvi)] * Fem::lagrange3_1D(p0, h, p, lvi, dif);
    }
    return E;
}
int Interpolant1D_Lagrange3::getBasFunGlobaIndex(const int feIndex, const int localBusFunIndex) const
{
    int basFun0Ind = feIndex * 3;
    // basFun0Ind - глобальный индекс базисной ф-и с локальным индексом 0,
    // он же - глобальный индекс вершины с локальным индексом 0
    return basFun0Ind + localBusFunIndex;
}
void Interpolant1D_Lagrange3::buildLagrange1DLocalMatrix()
{
    // построение локальной матрицы G для лагранжевых базисных функций 3-го порядка
     G[0][0] = 148;
     G[0][1] = -189;
     G[0][2] = 54;
     G[0][3] = -13;
      G[1][0] = -189;
      G[1][1] = 432;
      G[1][2] = -297;
      G[1][3] = 54;
     G[2][0] = 54;
     G[2][1] = -297;
     G[2][2] = 432;
     G[2][3] = -189;
      G[3][0] = -13;
      G[3][1] = 54;
      G[3][2] = -189;
      G[3][3] = 148;
      for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
          G[i][j] *= 1.L / (40 * h);
}

void Interpolant2D_base::init(const Grid::GridRectangleRegular2D &set_grid, const double set_alpha)
{
    gr = set_grid;
    alpha = set_alpha;
    l[0] = gr.rect[0][1] - gr.rect[0][0];
    l[1] = gr.rect[1][1] - gr.rect[1][0];
    h[0] = l[0] / gr.N[0];
    h[1] = l[1] / gr.N[1];
    m.clear();
    res.clear();
}
void Interpolant2D_base::release()
{
    m.clear();
    res.clear();
}
void Interpolant2D_base::addPoint(const POINT2 &p, const double value)
{
    POINT2_VALUE p2vEl;
    int i, j;
    findFe_ij(p,
              i, j);
    p2vEl.p = p;
    p2vEl.value = value;
    p2vEl.feIndex = getFeIndex(i, j);
    m.push_back(p2vEl);
}
void Interpolant2D_base::save(const char *file_name, const int bmp_x, const int bmp_y, const int colorMode, const int pictureMode)
{
    int r,g,b;
    BMP_INF bmp;		//
    // инициализация картинки
    bmp.x = bmp_x;
    bmp.y = bmp_y;
    bmp.c = new unsigned int[bmp.x * bmp.y];
    // рисование картинки
    for (int y = 0; y < bmp.y; y++)
        for (int x = 0; x < bmp.x; x++)
        {
            POINT2 p;
            p[0] = gr.rect[0][0] + l[0] * x / bmp.x;
            p[1] = gr.rect[1][0] + l[1] * y / bmp.y;
            calc_rgb(colorMode, fun(p), valueMin, valueNull, valueMax, r, g, b);
            PAINT(x, y, r, g, b);
        }
    // вывод картинки без градаций
    if(pictureMode == 0)
    {
        bmp_save(file_name, bmp);
        delete[] bmp.c;
        return;
    }
    // вывод картинки с градациями
    if(pictureMode == 1)
    {
        BMP_INF symbols;	// изображения символов
        BMP_INF bmp2;		// интерполянт + градации
        bmp_load("_symbols.bmp", symbols);
        // инициализация второй картинки
        bmp2.x = bmp_x+LEN*SIZE_X;
        bmp2.y = bmp_y;
        bmp2.c = new unsigned int[bmp2.x * bmp2.y];
        for(int y = 0; y < bmp.y; y++)
        for(int x = 0; x < bmp.x; x++)
        {
            bmp2.c[y*bmp2.x + x] = bmp.c[y*bmp.x + x];
        }
        for(int y = 0; y < bmp.y; y++)
        for(int x = bmp.x; x < bmp2.x; x++)
        {
            bmp2.c[y*bmp2.x + x] = RGB32(0, 255, 255, 255);;
        }

        for(int y = 0; y < SIZE_Y/2; y++)
        for(int x = bmp.x; x < bmp2.x; x++)
        {
            double value = valueMin;
            calc_rgb(colorMode, value, valueMin, valueNull, valueMax, r, g, b);
            bmp2.c[y*bmp2.x + x] = RGB32(0, r, g, b);;
        }
        for(int y = bmp.y-SIZE_Y/2; y < bmp.y; y++)
        for(int x = bmp.x; x < bmp2.x; x++)
        {
            double value = valueMax;
            calc_rgb(colorMode, value, valueMin, valueNull, valueMax, r, g, b);
            bmp2.c[y*bmp2.x + x] = RGB32(0, r, g, b);;
        }

        for(int y = SIZE_Y/2; y < bmp.y-SIZE_Y/2; y++)
        for(int x = bmp.x; x < bmp2.x; x++)
        {
            double value = valueMin + (valueMax - valueMin)*(y-SIZE_Y / 2)/(bmp.y - SIZE_Y);
            calc_rgb(colorMode, value, valueMin, valueNull, valueMax, r, g, b);
            bmp2.c[y*bmp2.x + x] = RGB32(0, r, g, b);;
        }
        double val = valueMin;
        double dval = (valueMax - valueMin)/K;
        for(int i = 0; i <= K; i++)
        {
            char str[100];
            str[0] = ' ';
            if(val >=0 )
                sprintf(str+1, "%.3le", val);
            else
                sprintf(str, "%.3le", val);
            print(str, bmp.x, SIZE_Y/2 + i*(bmp.y-SIZE_Y)/K, colorMode, symbols, bmp2);
            val += dval;
        }
        //print(0, 100, "-1.234e980", bmp2);
        // сохранение изображение
        bmp_save(file_name, bmp2);
        delete[] symbols.c;
        delete[] bmp2.c;
        delete[] bmp.c;
        return;
    }
}
void Interpolant2D_base::findFe_ij(const POINT2 &p,
                                  int &i, int &j)const
{
    i = (int)((p[1] - gr.rect[1][0]) / l[1] * gr.N[1]);  // 0..gr.N[1]
    j = (int)((p[0] - gr.rect[0][0]) / l[0] * gr.N[0]);  // 0..gr.N[0]
}
int Interpolant2D_base::getFeIndex(const int i, const int j)const
{
    return i * gr.N[0] + j;
}

int POINT2_VALUE_cmp(const void* x1, const void* x2)
{
    POINT2_VALUE *px1 = (POINT2_VALUE *)x1;
    POINT2_VALUE *px2 = (POINT2_VALUE *)x2;
    if (px1->feIndex < px2->feIndex) return -1;
    if (px1->feIndex > px2->feIndex) return  1;
    return 0;
}
void Interpolant2D_Hermite3::buildInterpolant()
{
    buildHermite2DLocalMatrix();
    int matrixSize = (gr.N[0] + 1)*(gr.N[1] + 1)*4; // размер матрицы СЛАУ
    // сортировка точек по параллелепипедам в которые они попали
    qsort(m.data(), m.size(), sizeof(POINT2_VALUE), POINT2_VALUE_cmp);
    SlauSolving::SSCM matrix;                           // матрица
    SlauSolving::SSCMBulder mBulder;                    // сборщик
    SlauSolving::SSCMPreconditioner_base *mPreconditioner;    // предобусловливатель
    SlauSolving::SSCMSolver_base *mSolver;                    // решатель СЛАУ
    Vector b;                                           // правая часть
    SlauSolving::SolverParameters ssp;
    /*
    ssp.type = SlauMetod_LOS;
    ssp.preconditioning = SlauPreconditioning_none;//SlauPreconditioning_none;
    ssp.eps = 1.e-14;
    ssp.maxIter = 500;
    */
    ssp.type = SlauSolving::SolverType::Direct;
    ssp.preconditioning = SlauSolving::Preconditioning::Profile_LLT;//SlauPreconditioning_none;

    matrix.init();
    b.resize(matrixSize);
    res.resize(matrixSize);
    for (int i = 0; i < matrixSize; i++)
    {
        b[i] = 0;
        res[i] = 0;
    }
    mBulder.setMatrixSize(matrixSize);
    mBulder.start();

    // занесение элементов локальных матриц в глобальную
    size_t gcount = 0;
    size_t count;
    for(int i = 0; i < gr.N[1]; i++)
    for(int j = 0; j < gr.N[0]; j++)
    {					// выбран конечный элемент (i,j)
        int gind[16];	// глобальные номера локальных функций конечного элемента (i,j)
        int feIndex = getFeIndex(i, j);
        POINT2 p0;
        p0[0] = gr.rect[0][0] + j * h[0];
        p0[1] = gr.rect[1][0] + i * h[1];		// прямоугольный конечный элемент задается одной вершиной с локальным номером 0
        for(int lvi = 0; lvi < 16; lvi++)
        {
            gind[lvi] = getBasFunGlobaIndex(feIndex, lvi);
        }
        // подсчет суммы базисных функций от точек на этом конечном элементе
        // и сборка локальной матрицы конечного элемента (i,j)
        for(int lvi = 0; lvi < 16; lvi++)
        for(int lvj = 0; lvj < 16; lvj++)	// lvi,lvj - локальные номера базисных функций (координаты в G)
        {
            // локальные номера вершин li, lj соответствуют глобальным номерам вершин gind[li], gind[lj]
            double sum = 0;
            count = gcount;
            for (;;)
            {
                if(count >= m.size() || m[count].feIndex > feIndex) break;
                sum +=  Fem::hermite_2D(p0, h[0], h[1], m[count].p, lvi, Fem::dif_NULL2) *
                        Fem::hermite_2D(p0, h[0], h[1], m[count].p, lvj, Fem::dif_NULL2);
                count++;
            }
            if(gind[lvi] >= gind[lvj])
                mBulder.addElement_not_null(sum + alpha*G[lvi][lvj], gind[lvi], gind[lvj]);
        }
        // сборка правой части
        for(int lvi = 0; lvi < 16; lvi++)	// lvi - локальные номера базисных функций
        {
            // локальные номера вершин li соответствуют глобальным номерам вершин gind[li]
            double sum = 0;
            count = gcount;
            for(;;)
            {
                if (count >= m.size() || m[count].feIndex > feIndex) break;
                sum += Fem::hermite_2D(p0, h[0], h[1], m[count].p, lvi, Fem::dif_NULL2) * m[count].value;
                count++;
            }
            b[gind[lvi]] += sum;
        }
        gcount = count;	// переходим к точкам следующего КЭ
    }
    mBulder.complete(matrix);
    double residual;
    double relativeResidual;
    int iterations;
    double time;
    double timePreconditioner;
    mPreconditioner = SlauSolving::SSCMPreconditioner_base::gen(ssp.preconditioning);
    mSolver = SlauSolving::SSCMSolver_base::gen(ssp.type);
    mPreconditioner->bulid(matrix, 0, timePreconditioner);
    mSolver->init(matrix.getMatrixSize());
    mSolver->solve(matrix, b,
                  mPreconditioner, res, ssp,
                  res, residual, relativeResidual, iterations, time);
    matrix.release();
    mBulder.release();
    delete mPreconditioner;
    delete mSolver;
    b.clear();
    valueMin = +1.e+200;
    valueMax = -1.e+200;
    for (size_t k = 0; k < m.size(); k++)
    {
        double value = m[k].value;
        if(value < valueMin) valueMin = value;
        if(value > valueMax) valueMax = value;
    }
    valueNull = (valueMin + valueMax) / 2;
}
void Interpolant2D_Hermite3::getNodesCoordinates(std::vector<POINT2> &coordinates)const
{
    coordinates.clear();
    for(int i = 0; i <= gr.N[1]; i++)
    for(int j = 0; j <= gr.N[0]; j++)
    {
        POINT2 p0;      // координаты узла (i,j)
        p0[0] = gr.rect[0][0] + j * h[0];
        p0[1] = gr.rect[1][0] + i * h[1];
        coordinates.push_back(p0);
    }
}
void Interpolant2D_Hermite3::buildInterpolantByAllNodes(const Vector &nodeValue)
{
    buildHermite2DLocalMatrix();
    int matrixSize = (gr.N[0] + 1)*(gr.N[1] + 1)*4; // размер матрицы СЛАУ
    SlauSolving::SSCM matrix;                           // матрица
    SlauSolving::SSCMBulder mBulder;                    // сборщик
    SlauSolving::SSCMPreconditioner_base *mPreconditioner;    // предобусловливатель
    SlauSolving::SSCMSolver_base *mSolver;                    // решатель СЛАУ
    Vector b;                                           // правая часть
    SlauSolving::SolverParameters ssp;
    ssp.type = SlauSolving::SolverType::Direct;
    ssp.preconditioning = SlauSolving::Preconditioning::Profile_LLT;

    matrix.init();
    b.resize(matrixSize);
    res.resize(matrixSize);
    for (int i = 0; i < matrixSize; i++)
    {
        b[i] = 0;
        res[i] = 0;
    }
    mBulder.setMatrixSize(matrixSize);
    mBulder.start();
    // занесение элементов локальных матриц в глобальную
    for(int i = 0; i < gr.N[1]; i++)
    for(int j = 0; j < gr.N[0]; j++)
    {					// выбран конечный элемент (i,j)
        int gind[16];	// глобальные номера локальных функций конечного элемента (i,j)
        int feIndex = getFeIndex(i, j);
        POINT2 p0;
        p0[0] = gr.rect[0][0] + j * h[0];
        p0[1] = gr.rect[1][0] + i * h[1];		// прямоугольный конечный элемент задается одной вершиной с локальным номером 0
        for(int lvi = 0; lvi < 16; lvi++)
        {
            gind[lvi] = getBasFunGlobaIndex(feIndex, lvi);
        }
        // узлы конечного элемента (i,j)
        POINT2 feVertex[4];
        feVertex[0][0] = p0[0];
        feVertex[0][1] = p0[1];
         feVertex[1][0] = p0[0] + h[0];
         feVertex[1][1] = p0[1];
        feVertex[2][0] = p0[0];
        feVertex[2][1] = p0[1] + h[1];
         feVertex[3][0] = p0[0] + h[0];
         feVertex[3][1] = p0[1] + h[1];
        // значение интерполируемой ф-и в этих узлах
        double value[4];
        int nodeIndex0 = i*(gr.N[0] + 1) + j;
        value[0] = nodeValue[nodeIndex0];
        value[1] = nodeValue[nodeIndex0 + 1];
        value[2] = nodeValue[nodeIndex0 + (gr.N[0] + 1)];
        value[3] = nodeValue[nodeIndex0 + (gr.N[0] + 1) + 1];

        // подсчет суммы базисных функций от точек на этом конечном элементе
        // и сборка локальной матрицы конечного элемента (i,j)
        for(int lvi = 0; lvi < 16; lvi++)
        for(int lvj = 0; lvj < 16; lvj++)	// lvi,lvj - локальные номера базисных функций (координаты в G)
        {
            // локальные номера вершин li, lj соответствуют глобальным номерам вершин gind[li], gind[lj]
            double sum = 0;
            sum +=  Fem::hermite_2D(p0, h[0], h[1], feVertex[0], lvi, Fem::dif_NULL2) *
                    Fem::hermite_2D(p0, h[0], h[1], feVertex[0], lvj, Fem::dif_NULL2);
            if(j == gr.N[0] - 1)
            {
                sum +=  Fem::hermite_2D(p0, h[0], h[1], feVertex[1], lvi, Fem::dif_NULL2) *
                        Fem::hermite_2D(p0, h[0], h[1], feVertex[1], lvj, Fem::dif_NULL2);
            }
            if(i == gr.N[1] - 1)
            {
                sum +=  Fem::hermite_2D(p0, h[0], h[1], feVertex[2], lvi, Fem::dif_NULL2) *
                        Fem::hermite_2D(p0, h[0], h[1], feVertex[2], lvj, Fem::dif_NULL2);
            }
            if(j == gr.N[0] - 1 && i == gr.N[1] - 1)
            {
                sum +=  Fem::hermite_2D(p0, h[0], h[1], feVertex[3], lvi, Fem::dif_NULL2) *
                        Fem::hermite_2D(p0, h[0], h[1], feVertex[3], lvj, Fem::dif_NULL2);
            }
            if(gind[lvi] >= gind[lvj])
                mBulder.addElement_not_null(sum + alpha*G[lvi][lvj], gind[lvi], gind[lvj]);
        }
        // сборка правой части
        for(int lvi = 0; lvi < 16; lvi++)	// lvi - локальные номера базисных функций
        {
            // локальные номера вершин li соответствуют глобальным номерам вершин gind[li]
            double sum = 0;
            sum += Fem::hermite_2D(p0, h[0], h[1], feVertex[0], lvi, Fem::dif_NULL2) * value[0];
            if(j == gr.N[0] - 1)
            {
                sum += Fem::hermite_2D(p0, h[0], h[1], feVertex[1], lvi, Fem::dif_NULL2) * value[1];
            }
            if(i == gr.N[1] - 1)
            {
                sum += Fem::hermite_2D(p0, h[0], h[1], feVertex[2], lvi, Fem::dif_NULL2) * value[2];
            }
            if(j == gr.N[0] - 1 && i == gr.N[1] - 1)
            {
                sum += Fem::hermite_2D(p0, h[0], h[1], feVertex[3], lvi, Fem::dif_NULL2) * value[3];
            }
            b[gind[lvi]] += sum;
        }
    }
    mBulder.complete(matrix);
    double residual;
    double relativeResidual;
    int iterations;
    double time;
    double timePreconditioner;
    mPreconditioner = SlauSolving::SSCMPreconditioner_base::gen(ssp.preconditioning);
    mSolver = SlauSolving::SSCMSolver_base::gen(ssp.type);
    mPreconditioner->bulid(matrix, 0, timePreconditioner);
    mSolver->init(matrix.getMatrixSize());
    mSolver->solve(matrix, b,
                  mPreconditioner, res, ssp,
                  res, residual, relativeResidual, iterations, time);
    matrix.release();
    mBulder.release();
    delete mPreconditioner;
    delete mSolver;
    b.clear();
    valueMin = +1.e+200;
    valueMax = -1.e+200;
    for (size_t k = 0; k < nodeValue.size(); k++)
    {
        double value = nodeValue[k];
        if(value < valueMin) valueMin = value;
        if(value > valueMax) valueMax = value;
    }
    valueNull = (valueMin + valueMax) / 2;
}
double Interpolant2D_Hermite3::fun(const POINT2 &p)const
{
    return difFun(p, Fem::dif_NULL2);
}
double Interpolant2D_Hermite3::difFun(const POINT2 &p, const DIF_STATE2 &dif)const
{
    int i, j;
    findFe_ij(p,
              i, j);
    /*if(i < 0 || i >= gr.N[1] || j < 0 || j >= gr.N[0])
    {
        // выход за пределы сетки
        return 0;
    }*/
    if(i < 0)
        i = 0;
    if(j < 0)
        j = 0;
    if(i >= gr.N[1])
        i = gr.N[1] - 1;
    if(j >= gr.N[0])
        j = gr.N[0] - 1;

    int feIndex = getFeIndex(i, j);
    double E = 0;
    POINT2 p0;
    p0[0] = gr.rect[0][0] + j * h[0];
    p0[1] = gr.rect[1][0] + i * h[1];		// прямоугольный конечный элемент задается одной вершиной с локальным номером 0
    for(int lvi = 0; lvi < 16; lvi++)
    {
        E += res[getBasFunGlobaIndex(feIndex, lvi)] * Fem::hermite_2D(p0, h[0], h[1], p, lvi, dif);
    }
    return E;
}
int Interpolant2D_Hermite3::getBasFunGlobaIndex(const int feIndex, const int localBusFunIndex)const
{
    int i = feIndex / gr.N[0];
    int j = feIndex % gr.N[0];
    int basFun0Ind = i * (gr.N[0] + 1) * 4 + j * 4; // глобальный индекс базисной ф-и с локальным индексом 0
    if(localBusFunIndex < 8)
        return basFun0Ind + localBusFunIndex;
    else
        return basFun0Ind + (gr.N[0] + 1) * 4 + (localBusFunIndex - 8);
}
void Interpolant2D_Hermite3::buildHermite2DLocalMatrix()
{
    double G1x[4][4];
    double M1x[4][4];
    double G1y[4][4];
    double M1y[4][4];
    // локальные матрицы для одномерных эрмитовых базисных функций 3-го порядка
     G1x[0][0] = 36;
     G1x[0][1] = 3 * h[0];
     G1x[0][2] = -36;
     G1x[0][3] = 3 * h[0];
      G1x[1][0] = 3 * h[0];
      G1x[1][1] = 4 * h[0]*h[0];
      G1x[1][2] = -3 * h[0];
      G1x[1][3] = -h[0]*h[0];
     G1x[2][0] = -36;
     G1x[2][1] = -3 * h[0];
     G1x[2][2] = 36;
     G1x[2][3] = -3 * h[0];
      G1x[3][0] = 3 * h[0];
      G1x[3][1] = -h[0]*h[0];
      G1x[3][2] = -3 * h[0];
      G1x[3][3] = 4 * h[0]*h[0];

     G1y[0][0] = 36;
     G1y[0][1] = 3 * h[1];
     G1y[0][2] = -36;
     G1y[0][3] = 3 * h[1];
      G1y[1][0] = 3 * h[1];
      G1y[1][1] = 4 * h[1]*h[1];
      G1y[1][2] = -3 * h[1];
      G1y[1][3] = -h[1]*h[1];
     G1y[2][0] = -36;
     G1y[2][1] = -3 * h[1];
     G1y[2][2] = 36;
     G1y[2][3] = -3 * h[1];
      G1y[3][0] = 3 * h[1];
      G1y[3][1] = -h[1]*h[1];
      G1y[3][2] = -3 * h[1];
      G1y[3][3] = 4 * h[1]*h[1];

    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
        G1x[i][j] *= 1.L / (30 * h[0]);
        G1y[i][j] *= 1.L / (30 * h[1]);
    }

     M1x[0][0] = 156;
     M1x[0][1] = 22 * h[0];
     M1x[0][2] = 54;
     M1x[0][3] = -13 * h[0];
      M1x[1][0] = 22 * h[0];
      M1x[1][1] = 4 * h[0]*h[0];
      M1x[1][2] = 13 * h[0];
      M1x[1][3] = -3*h[0]*h[0];
     M1x[2][0] = 54;
     M1x[2][1] = 13 * h[0];
     M1x[2][2] = 156;
     M1x[2][3] = -22 * h[0];
      M1x[3][0] = -13 * h[0];
      M1x[3][1] = -3*h[0]*h[0];
      M1x[3][2] = -22 * h[0];
      M1x[3][3] = 4 * h[0]*h[0];

     M1y[0][0] = 156;
     M1y[0][1] = 22 * h[1];
     M1y[0][2] = 54;
     M1y[0][3] = -13 * h[1];
      M1y[1][0] = 22 * h[1];
      M1y[1][1] = 4 * h[1]*h[1];
      M1y[1][2] = 13 * h[1];
      M1y[1][3] = -3*h[1]*h[1];
     M1y[2][0] = 54;
     M1y[2][1] = 13 * h[1];
     M1y[2][2] = 156;
     M1y[2][3] = -22 * h[1];
      M1y[3][0] = -13 * h[1];
      M1y[3][1] = -3*h[1]*h[1];
      M1y[3][2] = -22 * h[1];
      M1y[3][3] = 4 * h[1]*h[1];


    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
        M1x[i][j] *= h[0] / 420;
        M1y[i][j] *= h[1] / 420;
    }

    // построение локальной матрицы G для эрмитовых базисных функций 3-го порядка
    for(int i = 0; i < 16; i++)
    for(int j = 0; j < 16; j++)
    {
#define XI(i)	(2 * (((i) / 4) % 2) + ((i) % 2))
#define YI(i)	(2 * ((i) / 8) + (((i) / 2) % 2))
        int xi = XI(i);
        int yi = YI(i);
        int xj = XI(j);
        int yj = YI(j);
        G[i][j] =
            G1x[xi][xj]*M1y[yi][yj] +
            M1x[xi][xj]*G1y[yi][yj];
    }
}

void Interpolant2D_Lagrange3::buildInterpolant()
{

    buildLagrange2DLocalMatrix();
    int matrixSize = (gr.N[0] * 3 + 1)*(gr.N[1] * 3 + 1); // размер матрицы СЛАУ
    // сортировка точек по параллелепипедам в которые они попали
    qsort(m.data(), m.size(), sizeof(POINT2_VALUE), POINT2_VALUE_cmp);
    SlauSolving::SSCM matrix;                           // матрица
    SlauSolving::SSCMBulder mBulder;                    // сборщик
    SlauSolving::SSCMPreconditioner_base *mPreconditioner;    // предобусловливатель
    SlauSolving::SSCMSolver_base *mSolver;                    // решатель СЛАУ
    Vector b;                                           // правая часть
    SlauSolving::SolverParameters ssp;
    ssp.type = SlauSolving::SolverType::Direct;
    ssp.preconditioning = SlauSolving::Preconditioning::Profile_LLT;

    matrix.init();
    b.resize(matrixSize);
    res.resize(matrixSize);
    for (int i = 0; i < matrixSize; i++)
    {
        b[i] = 0;
        res[i] = 0;
    }
    mBulder.setMatrixSize(matrixSize);
    mBulder.start();

    // занесение элементов локальных матриц в глобальную
    size_t gcount = 0;
    size_t count;
    for(int i = 0; i < gr.N[1]; i++)
    for(int j = 0; j < gr.N[0]; j++)
    {					// выбран конечный элемент (i,j)
        int gind[16];	// глобальные номера локальных функций конечного элемента (i,j)
        int feIndex = getFeIndex(i, j);
        POINT2 p0;
        p0[0] = gr.rect[0][0] + j * h[0];
        p0[1] = gr.rect[1][0] + i * h[1];		// прямоугольный конечный элемент задается одной вершиной с локальным номером 0
        for(int lvi = 0; lvi < 16; lvi++)
        {
            gind[lvi] = getBasFunGlobaIndex(feIndex, lvi);
        }
        // подсчет суммы базисных функций от точек на этом конечном элементе
        // и сборка локальной матрицы конечного элемента (i,j)
        for(int lvi = 0; lvi < 16; lvi++)
        for(int lvj = 0; lvj < 16; lvj++)	// lvi,lvj - локальные номера базисных функций (координаты в G)
        {
            // локальные номера вершин li, lj соответствуют глобальным номерам вершин gind[li], gind[lj]
            double sum = 0;
            count = gcount;
            for (;;)
            {
                if(count >= m.size() || m[count].feIndex > feIndex) break;
                sum +=  Fem::lagrange3_2D(p0, h[0], h[1], m[count].p, lvi, Fem::dif_NULL2) *
                        Fem::lagrange3_2D(p0, h[0], h[1], m[count].p, lvj, Fem::dif_NULL2);
                count++;
            }
            if(gind[lvi] >= gind[lvj])
                mBulder.addElement_not_null(sum + alpha*G[lvi][lvj], gind[lvi], gind[lvj]);
        }
        // сборка правой части
        for(int lvi = 0; lvi < 16; lvi++)	// lvi - локальные номера базисных функций
        {
            // локальные номера вершин li соответствуют глобальным номерам вершин gind[li]
            double sum = 0;
            count = gcount;
            for(;;)
            {
                if (count >= m.size() || m[count].feIndex > feIndex) break;
                sum += Fem::lagrange3_2D(p0, h[0], h[1], m[count].p, lvi, Fem::dif_NULL2) * m[count].value;
                count++;
            }
            b[gind[lvi]] += sum;
        }
        gcount = count;	// переходим к точкам следующего КЭ
    }
    mBulder.complete(matrix);
    double residual;
    double relativeResidual;
    int iterations;
    double time;
    double timePreconditioner;
    mPreconditioner = SlauSolving::SSCMPreconditioner_base::gen(ssp.preconditioning);
    mSolver = SlauSolving::SSCMSolver_base::gen(ssp.type);
    mPreconditioner->bulid(matrix, 0, timePreconditioner);
    mSolver->init(matrix.getMatrixSize());
    mSolver->solve(matrix, b,
                  mPreconditioner, res, ssp,
                  res, residual, relativeResidual, iterations, time);
    matrix.release();
    mBulder.release();
    delete mPreconditioner;
    delete mSolver;
    b.clear();
    valueMin = +1.e+200;
    valueMax = -1.e+200;
    for (size_t k = 0; k < m.size(); k++)
    {
        double value = m[k].value;
        if (value < valueMin) valueMin = value;
        if (value > valueMax) valueMax = value;
    }
    valueNull = (valueMin + valueMax) / 2;
}
void Interpolant2D_Lagrange3::getNodesCoordinates(std::vector<POINT2> &coordinates)const
{
    coordinates.clear();
    for(int i = 0; i <= gr.N[1]*3; i++)
    for(int j = 0; j <= gr.N[0]*3; j++)
    {
        POINT2 p0;      // координаты узла (i,j)
        p0[0] = gr.rect[0][0] + j * h[0] / 3;
        p0[1] = gr.rect[1][0] + i * h[1] / 3;
        coordinates.push_back(p0);
    }
}
void Interpolant2D_Lagrange3::buildInterpolantByAllNodes(const Vector &nodeValue)
{
    int matrixSize = (gr.N[0] * 3 + 1)*(gr.N[1] * 3 + 1); // размер матрицы СЛАУ
    res.resize(matrixSize);
    int count = 0;
    valueMin = +1.e+200;
    valueMax = -1.e+200;
    for(int i = 0; i <= gr.N[1]*3; i++)
    for(int j = 0; j <= gr.N[0]*3; j++)
    {
        double value = nodeValue[count];
        if (value < valueMin) valueMin = value;
        if (value > valueMax) valueMax = value;
        res[count] = value;
        count++;
    }
    valueNull = (valueMin + valueMax) / 2;
}
double Interpolant2D_Lagrange3::fun(const POINT2 &p)const
{
    return difFun(p, Fem::dif_NULL2);
}
double Interpolant2D_Lagrange3::difFun(const POINT2 &p, const DIF_STATE2 &dif)const
{
    int i, j;
    findFe_ij(p,
              i, j);
    /*if(i < 0 || i >= gr.N[1] || j < 0 || j >= gr.N[0])
    {
        // выход за пределы сетки
        return 0;
    }*/
    if(i < 0)
        i = 0;
    if(j < 0)
        j = 0;
    if(i >= gr.N[1])
        i = gr.N[1] - 1;
    if(j >= gr.N[0])
        j = gr.N[0] - 1;
    int feIndex = getFeIndex(i, j);
    double E = 0;
    POINT2 p0;
    p0[0] = gr.rect[0][0] + j * h[0];
    p0[1] = gr.rect[1][0] + i * h[1];		// прямоугольный конечный элемент задается одной вершиной с локальным номером 0
    for(int lvi = 0; lvi < 16; lvi++)
    {
        E += res[getBasFunGlobaIndex(feIndex, lvi)] * Fem::lagrange3_2D(p0, h[0], h[1], p, lvi, dif);
    }
    return E;

}
int Interpolant2D_Lagrange3::getBasFunGlobaIndex(const int feIndex, const int localBusFunIndex)const
{
    int i = feIndex / gr.N[0];
    int j = feIndex % gr.N[0];
    int basFun0Ind = i * 3 * (gr.N[0] * 3 + 1) + j * 3;
    // basFun0Ind - глобальный индекс базисной ф-и с локальным индексом 0,
    // он же - глобальный индекс вершины с локальным индексом 0
    int dj = localBusFunIndex % 4;
    int di = localBusFunIndex / 4;
    return basFun0Ind + di * (gr.N[0] * 3 + 1) + dj;
}
void Interpolant2D_Lagrange3::buildLagrange2DLocalMatrix()
{
    double G1x[4][4];
    double G1y[4][4];
    double M1x[4][4];
    double M1y[4][4];

    // локальные матрицы для одномерных лагранжевых базисных функций 3-го порядка
     G1x[0][0] = 148;
     G1x[0][1] = -189;
     G1x[0][2] = 54;
     G1x[0][3] = -13;
      G1x[1][0] = -189;
      G1x[1][1] = 432;
      G1x[1][2] = -297;
      G1x[1][3] = 54;
     G1x[2][0] = 54;
     G1x[2][1] = -297;
     G1x[2][2] = 432;
     G1x[2][3] = -189;
      G1x[3][0] = -13;
      G1x[3][1] = 54;
      G1x[3][2] = -189;
      G1x[3][3] = 148;

    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
        G1y[i][j] = G1x[i][j];

        G1x[i][j] *= 1.L / (40 * h[0]);
        G1y[i][j] *= 1.L / (40 * h[1]);
    }

     M1x[0][0] = 128;
     M1x[0][1] = 99;
     M1x[0][2] = -36;
     M1x[0][3] = 19;
      M1x[1][0] = 99;
      M1x[1][1] = 618;
      M1x[1][2] = -81;
      M1x[1][3] = -36;
     M1x[2][0] = -36;
     M1x[2][1] = -81;
     M1x[2][2] = 648;
     M1x[2][3] = 99;
      M1x[3][0] = 19;
      M1x[3][1] = -36;
      M1x[3][2] = 99;
      M1x[3][3] = 128;

    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
        M1y[i][j] = M1x[i][j];

        M1x[i][j] *= h[0] / 1680;
        M1y[i][j] *= h[1] / 1680;
    }

    // построение локальной матрицы G для эрмитовых базисных функций 3-го порядка
    for(int i = 0; i < 16; i++)
    for(int j = 0; j < 16; j++)
    {
#define lXI(i)	((i) % 4)
#define lYI(i)	((i) / 4)
        int xi = lXI(i);
        int yi = lYI(i);
        int xj = lXI(j);
        int yj = lYI(j);
        G[i][j] =
            G1x[xi][xj]*M1y[yi][yj] +
            M1x[xi][xj]*G1y[yi][yj];
    }
}



InterpolantSurfaceDecomposition_base *InterpolantSurfaceDecomposition_base::newEl(const Grid::RegionDecompositionType elType)
{
    switch (elType)
    {
    case Grid::RegionDecompositionType::none:
        return new InterpolantSurfaceDecomposition_none;
        break;
    case Grid::RegionDecompositionType::cubeVertexIndexation:
        return new InterpolantSurfaceDecomposition_cubeVertexIndexation;
        break;
    case Grid::RegionDecompositionType::surfaceAsSpheres:
        return new InterpolantSurfaceDecomposition_surfaceAsSpheres;
        break;
    }
}

void InterpolantSurfaceDecomposition_none::buildRegions(const Grid::RegionDecompositionParameters *, const InterpolantSurface_base *)
{
}
void InterpolantSurfaceDecomposition_none::excludeIntersectionCase(const POINT3 &, const POINT3 &, const Grid::SurfacePositionData &, const Grid::SurfacePositionData &,
                                                                   Grid::SurfacePositionData &, Grid::SurfacePositionData &, bool &mayBeFound)const
{
    mayBeFound = true;
}

void InterpolantSurfaceDecomposition_cubeVertexIndexation::buildRegions(const Grid::RegionDecompositionParameters *set_param, const InterpolantSurface_base *surface)
{
    param = *set_param;
    Grid::RegionIndexation_vertexIndex v0Index;
    v0Index[0] = 0;
    v0Index[1] = 0;
    v0Index[2] = 0;
    Grid::SurfacePositionData ssd0;
    ssd0.index = -1;
    ri.init(param, v0Index, 0, surface, ssd0,
            &map);
}
void InterpolantSurfaceDecomposition_cubeVertexIndexation::excludeIntersectionCase(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &prevSolutionData1, const Grid::SurfacePositionData &prevSolutionData2,
                                                                                   Grid::SurfacePositionData &solutionData1, Grid::SurfacePositionData &solutionData2, bool &mayBeFound)const
{
    bool p1_defined;
    int p1_side;
    bool p1_onBorder;
    Grid::SurfacePositionData p1_ssd;
    bool p2_defined;
    int p2_side;
    bool p2_onBorder;
    Grid::SurfacePositionData p2_ssd;
    ri.getInformation(param, &map, p1, p1_defined, p1_side, p1_onBorder, p1_ssd);
    ri.getInformation(param, &map, p2, p2_defined, p2_side, p2_onBorder, p2_ssd);
    if(p1_defined && p2_defined)
    {
        // получена точная информация для обоих точек
        if(p1_side == p2_side)
        {
            // точки с одной и той же стороны от поверхности - пересечения с поверхностью нет
            mayBeFound = false;
            return;
        }
    }
    mayBeFound = true;
    // не известно есть ли пересечение, значит необходим полноценный поиск
    // и начальные приближения
    // при выборе начальных приближений предпочтение отдаём входным данным
    if(prevSolutionData1.index == -1)
    {
        solutionData1 = p1_ssd;
    }
    else
    {
        solutionData1 = prevSolutionData1;
    }
    if(prevSolutionData2.index == -1)
    {
        solutionData2 = p2_ssd;
    }
    else
    {
        solutionData2 = prevSolutionData2;
    }
}

void InterpolantSurfaceDecomposition_surfaceAsSpheres::buildRegions(const Grid::RegionDecompositionParameters *set_param, const InterpolantSurface_base *surface)
{
    param = *set_param;
    int div = param.div;
    std::vector<POINT3> vertex;
    vertex.clear();
    // 0) составление массива узлов
    for(int i = 0; i <= surface->gr.N[1]*div; i++)
    for(int j = 0; j <= surface->gr.N[0]*div; j++)
    {
        POINT2 uv;  // узел
        uv[0] = surface->gr.rect[0][0] + j * surface->h[0]/div;
        uv[1] = surface->gr.rect[1][0] + i * surface->h[1]/div;
        POINT3 v0 = surface->s(uv);
        vertex.push_back(v0);
    }
    // 1) составление массива сферированых тел
    std::vector<Grid::SpheredBody> sb;
    sb.clear();
    for(int i = 0; i < surface->gr.N[1]*div; i++)
    for(int j = 0; j < surface->gr.N[0]*div; j++)
    {
        Grid::SpheredBody sbEl;
        int vInd0 = i*(surface->gr.N[0]*div + 1) + j;
        sbEl.bodyInd = vInd0;           // индекс сферированного тела = индекс нулевой вершины
        POINT3 localVertex[4];          // тело - это набор из четырёх узлов (которые составляют участок поверхности)
        int localVertexIndex = 0;
        for(int di = 0; di <= 1; di++)
        for(int dj = 0; dj <= 1; dj++)
        {
            int vInd = vInd0 + (surface->gr.N[0]*div + 1)*di + dj;       // индекс узла
            localVertex[localVertexIndex] = vertex[vInd];   // координаты узла
            localVertexIndex++;
        }
        sbEl.s.init0(localVertex[0]);
        sbEl.s.expand(localVertex[1]);  // узлы расширяют сферу
        sbEl.s.expand(localVertex[2]);
        sbEl.s.expand(localVertex[3]);
        // сферированное тело sbEl добавляется в массив
        sb.push_back(sbEl);
    }
    // 2) поиск параллелепипеда, который содержит все сферированные тела и инициализация корневой области
    CUBE q0;
    q0.initBySphere(sb[0].s);
    for(size_t i = 1; i < sb.size(); i++)
    {
        // каждое новое сферированное тело может расширить параллелепипед
        q0.expand(sb[i].s);
    }
    // нашли параллелепипед, который содержит все тела
    // преобразуем его в куб
    q0.toCube();
    // 3) инициализация корневой области
    ri.init(q0);
    // 4) добавление сферированных тел (индекс + шар) в корневую область ri
    for(size_t i = 0; i < sb.size(); i++)
    {
        ri.insert(sb[i]);
    }
}
void InterpolantSurfaceDecomposition_surfaceAsSpheres::excludeIntersectionCase(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &, const Grid::SurfacePositionData &,
                                                                               Grid::SurfacePositionData &, Grid::SurfacePositionData &, bool &mayBeFound)const
{
    // быстрая проверка на отсутствие пересечения
    Grid::SpheredBodyArray bodyInd;
    bodyInd.clear();
    ri.findIntersectingBodies(p1, p2,
                              bodyInd);
    //Sphere s0;  // отрезок представляем в виде сферы
    //s0.init0(p1);
    //s0.expand(p2);
    //ri.findIntersectingBodies(s0,
    //                          bodyInd);
    if(bodyInd.size() == 0)
    {
        mayBeFound = false;
    }
    else
    {
        mayBeFound = true;
    }
}



InterpolantSurface_base::InterpolantSurface_base(const Grid::RegionDecompositionType set_decompositionType, const Grid::GridRectangleRegular2D &set_gr, const int set_maxIter, const double set_eps_uv_2d, const double set_eps_uv_1d, const double set_eps_grad, const double set_eps_inerception)
{
    maxIter = set_maxIter;
    eps_uv_2d = set_eps_uv_2d;
    eps_uv_1d = set_eps_uv_1d;
    eps_grad = set_eps_grad;
    eps_inerception = set_eps_inerception;
    gr = set_gr;
    l[0] = gr.rect[0][1] - gr.rect[0][0];
    l[1] = gr.rect[1][1] - gr.rect[1][0];
    h[0] = l[0] / gr.N[0];
    h[1] = l[1] / gr.N[1];
    decomposition = InterpolantSurfaceDecomposition_base::newEl(set_decompositionType);
    decompositionType = set_decompositionType;
    regionIndexationInited = false;
}
void InterpolantSurface_base::findIntersection(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &prevSolutionData1, const Grid::SurfacePositionData &prevSolutionData2,
                                               Grid::SurfacePositionData &solutionData1, Grid::SurfacePositionData &solutionData2, POINT3 &intersectionPoint, POINT3 &normal, int &side, bool &found)const
{
    // быстрая проверка на отсутствие пересечения
    bool mayBeFound;
    decomposition->excludeIntersectionCase(p1, p2, prevSolutionData1, prevSolutionData2,
                                           solutionData1, solutionData2, mayBeFound);
    if(!mayBeFound)
    {
        // пересечений нет
        found = false;
        return;
    }
    // пересечение может быть
    // полноценный поиск
    double a = 0;
    double b = (p2 - p1).abs();
    if(b == 0)
    {
        found = false;
        return;
    }
    /*
    {
        // информация о ближайшей к x2 точке поверхности
        POINT3 x2_nearestPoint;
        POINT3 x2_normal;
        int x2_side;
        bool x2_onBorder;
        Grid::SurfaceSolutionData ssd;
        ssd.index = -1;
        findNearestPoint(p2, p2 - p1, ssd,
                         ssd, x2_nearestPoint, x2_normal, x2_side, x2_onBorder);
        if(x2_side == 1)
        {
            found = false;
            return;
        }
    }
    */
    POINT3 d = (p2 - p1) / b; // p1 + d*t - прямая,  a <= t <= b
    POINT3 a_nearestPoint;
    POINT3 b_nearestPoint;
    POINT3 c_nearestPoint;
    POINT3 a_normal;
    POINT3 b_normal;
    POINT3 c_normal;
    int a_side;
    int b_side;
    int c_side;
    bool a_onBorder;
    bool b_onBorder;
    bool c_onBorder;
    findNearestPoint(p1 + d*a, {0,0,0}, prevSolutionData1,
                     solutionData1, a_nearestPoint, a_normal, a_side, a_onBorder);
    findNearestPoint(p1 + d*b, {0,0,0}, prevSolutionData2,
                     solutionData2, b_nearestPoint, b_normal, b_side, b_onBorder);
    if(a_side == b_side)
    {
        // концы отрезка находятся с одной стороны от поверхности => пересечения нет
        found = false;
        return;
    }

    if(a_onBorder == true && b_onBorder == true)
    {
        // концы отрезка находятся за пределами поверхности => пересечения нет
        found = false;
        return;
    }

    // поиск пересечения методом дихотомии
    double c;
    Grid::SurfacePositionData ssd_c;
    ssd_c.uv0 = (solutionData1.uv0 + solutionData2.uv0) / 2;
    ssd_c.index = 0;    // данные есть
    for(;;)
    {
        c = (a + b) / 2;
        findNearestPoint(p1 + d*c, {0,0,0}, ssd_c,
                         ssd_c, c_nearestPoint, c_normal, c_side, c_onBorder);
        if(fabs(a - b) < eps_inerception)
            break;
        if(c_side == a_side)
        {
            // границу a смещаем в положение c
            a = c;
        }
        else
        {
            // c_side == b_side
            // границу b смещаем в положение c
            b = c;
        }
    }

    if(c_onBorder)
    {
        // пересечение за пределами поверхности
        found = false;
        return;
    }

    intersectionPoint = p1 + d*c;
    normal = c_normal;
    side = b_side;  // +1 если точка p2 находится в внешней стороны поверхности
    found = true;
}
void InterpolantSurface_base::update(const double, const Grid::RegionDecompositionParameters *set_param)
{
    if(!regionIndexationInited)
    {
        decomposition->buildRegions(set_param, this);
        regionIndexationInited = true;
    }
}
void InterpolantSurface_base::findFe(const POINT2 &uv, int &i, int &j, POINT2 &uv0) const
{
    i = (int)((uv[1] - gr.rect[1][0]) / l[1] * gr.N[1]);  // 0..gr.N[1]
    j = (int)((uv[0] - gr.rect[0][0]) / l[0] * gr.N[0]);  // 0..gr.N[0]
    if(i < 0)
        i = 0;
    if(j < 0)
        j = 0;
    if(i >= gr.N[1])
        i = gr.N[1] - 1;
    if(j >= gr.N[0])
        j = gr.N[0] - 1;
    uv0[0] = gr.rect[0][0] + j * h[0];
    uv0[1] = gr.rect[1][0] + i * h[1];		// прямоугольный конечный элемент задается одной вершиной с локальным номером 0
}
double InterpolantSurface_base::findNearestPoint_find1d_zolot_sechenie(const POINT3 &p, const POINT2 &uv0, const int coordinateIndex, double a, double b, double eps) const
{
    double x1, x2, f1, f2;
    POINT2 uv = uv0;
    x1 = a + 0.381966011250105*(b - a);
    x2 = b - 0.381966011250105*(b - a);
    uv[coordinateIndex] = x1; f1 = (p - s(uv)).abs();
    uv[coordinateIndex] = x2; f2 = (p - s(uv)).abs();
    for(;;)
    {
        if(fabs(a - b) < eps)break;
        if(f1 > f2)
        {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = b - 0.381966011250105*(b - a);
            uv[coordinateIndex] = x2; f2 = (p - s(uv)).abs();
        }
        else
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + 0.381966011250105*(b - a);
            uv[coordinateIndex] = x1; f1 = (p - s(uv)).abs();
        }
    }
    return x1;
}
double InterpolantSurface_base::findNearestPoint_find1d_zolot_sechenie_vector(const POINT3 &p, const POINT2 &uv0, const VECTOR2 &duv, double a, double b, double eps) const
{
    double x1, x2, f1, f2;
    x1 = a + 0.381966011250105*(b - a);
    x2 = b - 0.381966011250105*(b - a);
    f1 = (p - s(uv0 + duv*x1)).abs();
    f2 = (p - s(uv0 + duv*x2)).abs();
    for(;;)
    {
        if(fabs(a - b) < eps)break;
        if(f1 > f2)
        {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = b - 0.381966011250105*(b - a);
            f2 = (p - s(uv0 + duv*x2)).abs();
        }
        else
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + 0.381966011250105*(b - a);
            f1 = (p - s(uv0 + duv*x1)).abs();
        }
    }
    return x1;
}
bool InterpolantSurface_base::pointOnBorder(const POINT2 &uv, const double eps_uv)const
{
    if(fabs(uv[0] - gr.rect[0][0]) < eps_uv ||
       fabs(uv[0] - gr.rect[0][1]) < eps_uv ||
       fabs(uv[1] - gr.rect[1][0]) < eps_uv ||
       fabs(uv[1] - gr.rect[1][1]) < eps_uv)
    {
        return true;    // точка uv находится на границе поверхности
    }
    else
    {
        return false;   // точка uv не находится на границе поверхности
    }
}
void InterpolantSurface_base::calcInterval(const POINT2 &uv0, const POINT2 duv_1,
                                            double &a, double &b)const
{
    // gr.rect[0][0] <= uv0[0] + duv_1[0]*t <= gr.rect[0][1]
    // gr.rect[1][0] <= uv0[1] + duv_1[1]*t <= gr.rect[1][1]

    // uv0[0] + duv_1[0]*t = gr.rect[0][0] => tu_1 = (gr.rect[0][0] - uv0[0]) / duv_1[0]
    // uv0[0] + duv_1[0]*t = gr.rect[0][1] => tu_2 = (gr.rect[0][1] - uv0[0]) / duv_1[0]
    // duv_1[0] == 0 => tu_1 = -1.e200, tu_2 = +1.e200
    // MIN(tu_1, tu_2) <= t <= MAX(tu_1, tu_2)

    // uv0[1] + duv_1[1]*t = gr.rect[1][0] => tv_1 = (gr.rect[1][0] - uv0[1]) / duv_1[1]
    // uv0[1] + duv_1[1]*t = gr.rect[1][1] => tv_2 = (gr.rect[1][1] - uv0[1]) / duv_1[1]
    // duv_1[1] == 0 => tu_1 = -1.e200, tu_2 = +1.e200
    // MIN(tv_1, tv_2) <= t <= MAX(tv_1, tv_2)

    double tu_1, tu_2, tv_1, tv_2;
    if(duv_1[0] == 0)
    {
        tu_1 = -1.e200;
        tu_2 = +1.e200;
    }
    else
    {
        tu_1 = (gr.rect[0][0] - uv0[0]) / duv_1[0];
        tu_2 = (gr.rect[0][1] - uv0[0]) / duv_1[0];
        if(tu_1 > tu_2)
        {
            double temp = tu_1;
            tu_1 = tu_2;
            tu_2 = temp;
        }
    }
    if(duv_1[1] == 0)
    {
        tv_1 = -1.e200;
        tv_2 = +1.e200;
    }
    else
    {
        tv_1 = (gr.rect[1][0] - uv0[1]) / duv_1[1];
        tv_2 = (gr.rect[1][1] - uv0[1]) / duv_1[1];
        if(tv_1 > tv_2)
        {
            double temp = tv_1;
            tv_1 = tv_2;
            tv_2 = temp;
        }
    }
    a = MAX(tu_1, tv_1);
    b = MIN(tu_2, tv_2);
}
void InterpolantSurface_base::loadSolutionData(const Grid::SurfacePositionData &solutionData,
                                               POINT2 &uv0)const
{
    if(solutionData.index == -1)
    {
        // данные о начальном приближении отсутствуют
        // начальное приближение - центр области
        uv0[0] = (gr.rect[0][0] + gr.rect[0][1]) / 2;
        uv0[1] = (gr.rect[1][0] + gr.rect[1][1]) / 2;
    }
    else
    {
        // начальное приближение - предыдущее решение
        uv0 = solutionData.uv0;
    }
}
void InterpolantSurface_base::saveSolutionData(const POINT2 &uv0,
                                               Grid::SurfacePositionData &solutionData)const
{
    solutionData.index = 0; // данные есть
    solutionData.uv0 = uv0;
}

AnaliticalSurface_Cylinder::AnaliticalSurface_Cylinder(const Grid::RegionDecompositionType set_decompositionType, const Grid::GridRectangleRegular2D &set_gr, Interpolant1D_base *set_R, Interpolant1D_base *set_x, Interpolant1D_base *set_y, const int set_maxIter, const double set_eps_uv_2d, const double set_eps_uv_1d, const double set_eps_grad, const double set_eps_inerception):
    InterpolantSurface_base(set_decompositionType, set_gr, set_maxIter, set_eps_uv_2d, set_eps_uv_1d, set_eps_grad, set_eps_inerception)
{
    R_fun = set_R;
    x_fun = set_x;
    y_fun = set_y;
    surfaceType = Grid::SurfaceType::AnaliticalSurface_Cylinder;
}
void AnaliticalSurface_Cylinder::findNearestPoint(const POINT3 &p, const VECTOR3 &, const Grid::SurfacePositionData &,
                                                  Grid::SurfacePositionData &, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const
{
    POINT3 C_v = {C[0], C[1], p[2]};         // центр в плоскости вершины по Z
    normal = (p - C_v) / (p - C_v).abs();    // нормаль наружу
    nearestPoint = C_v + R*normal;           // ближайшая к точке p точка поверхности
    if(sqrt(SQR(p[0] - C[0]) + SQR(p[1] - C[1])) - R >= 0)
        side = +1;
    else
        side = -1;
    onBorder = false;
    //normal = VECTOR3(0, -1, 0);//#######
}
void AnaliticalSurface_Cylinder::findIntersection_analitical(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &, const Grid::SurfacePositionData &,
                                                             Grid::SurfacePositionData &, Grid::SurfacePositionData &, POINT3 &intersectionPoint, POINT3 &normal, int &side, bool &found)const
{
    VECTOR3 n = p2 - p1;    // направление p1->p2
    n = n / n.abs();        // ищем пересечение прямой p1 + n*t с цилиндром (координата z любая)
    POINT3 d = p1 - C;
    double D = SQR(d[0]*n[0] + d[1]*n[1]) - (SQR(n[0]) + SQR(n[1]))*(SQR(d[0]) + SQR(d[1]) - SQR(R));
    if(D < 0)
    {
        // пересечения нет
        found = false;
        return;
    }
    else
    {
        // пересечение есть
        double t1 = (-(d[0]*n[0] + d[1]*n[1]) + sqrt(D)) / (SQR(n[0]) + SQR(n[1]));
        double t2 = (-(d[0]*n[0] + d[1]*n[1]) - sqrt(D)) / (SQR(n[0]) + SQR(n[1]));
        double t0;
        // интересует наименьшее неотрицательное решение
        if(t1 >= 0 && t2 >= 0)
        {
            t0 = MIN(t1, t2);
        }
        else
        {
            if(t1 < 0 && t2 < 0)
            {
                // оба решения < 0
                // точка пересечения вне отрезка
                found = false;   // пересечения нет
                return;
            }
            else
            {
                // ровно 1 решение < 0
                if(t1 < 0)
                {
                    // t1 < 0
                    t0 = t2;
                }
                else
                {
                    // t2 < 0
                    t0 = t1;
                }
            }
        }
        // найдено решение t0 > 0
        // проверим, находится ли точка p1 + n*t0 на отрезке p1-p2
        if((n*t0).abs() <= (p2 - p1).abs())
        {
            // точка p1 + n*t0 находится на отрезке p1-p2
            intersectionPoint = p1 + n*t0;   // точка пересечения
            POINT3 C_v = {C[0], C[1], intersectionPoint[2]};         // центр в плоскости вершины по Z
            normal = (intersectionPoint - C_v);                      // направление нормали
            normal = normal / normal.abs();                     // единичная нормаль к поверхности в точке пересечения
            if(normal*(p2 - p1) >= 0)
            {
                side = +1;     // p2 находится с внешней стороны поверхности
            }
            else
            {
                side = -1;     // p2 находится с внутренней стороны поверхности
            }
            found = true;   // пересечение найдено
            return;

        }
        else
        {
            found = false;   // пересечения нет
            return;
        }
    }
}
void AnaliticalSurface_Cylinder::findIntersection(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &, const Grid::SurfacePositionData &prevSolutionData2,
                                                             Grid::SurfacePositionData &, Grid::SurfacePositionData &solutionData2, POINT3 &intersectionPoint, POINT3 &normal, int &side, bool &found)const
{
    // информация о ближайшей к x2 точке поверхности
    POINT3 x2_nearestPoint;
    POINT3 x2_normal;
    int x2_side;
    bool x2_onBorder;
    findNearestPoint(p2, p2 - p1, prevSolutionData2,
                     solutionData2, x2_nearestPoint, x2_normal, x2_side, x2_onBorder);
    if(x2_side == -1)
    {
        // точка x2 находится с внутренней стороны поверхности
        // есть контакт
        intersectionPoint = x2_nearestPoint;
        normal = x2_normal;
        side = x2_side;
        found = true;
    }
    else
    {
        found = false;
    }
}
void AnaliticalSurface_Cylinder::update(const double time, const Grid::RegionDecompositionParameters *set_param)
{
    C[0] = x_fun->fun(time);
    C[1] = y_fun->fun(time);
    R = R_fun->fun(time);
}

POINT3 AnaliticalSurface_Cylinder::s(const POINT2 &) const
{
    return VECTOR3_NULL;
}
POINT3 AnaliticalSurface_Cylinder::difs(const POINT2 &, const DIF_STATE2 &) const
{
    return VECTOR3_NULL;
}

AnaliticalSurface_Sphere::AnaliticalSurface_Sphere(const Grid::RegionDecompositionType set_decompositionType, const Grid::GridRectangleRegular2D &set_gr, Interpolant1D_base *set_R, Interpolant1D_base *set_x, Interpolant1D_base *set_y, Interpolant1D_base *set_z, const int set_maxIter, const double set_eps_uv_2d, const double set_eps_uv_1d, const double set_eps_grad, const double set_eps_inerception):
    InterpolantSurface_base(set_decompositionType, set_gr, set_maxIter, set_eps_uv_2d, set_eps_uv_1d, set_eps_grad, set_eps_inerception)
{
    R_fun = set_R;
    x_fun = set_x;
    y_fun = set_y;
    z_fun = set_z;
    surfaceType = Grid::SurfaceType::AnaliticalSurface_Sphere;
}
void AnaliticalSurface_Sphere::findNearestPoint(const POINT3 &p, const VECTOR3 &, const Grid::SurfacePositionData &,
                                                  Grid::SurfacePositionData &, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const
{
    normal = (p - C) / (p - C).abs();    // нормаль наружу
    nearestPoint = C + R*normal;         // ближайшая к точке p точка поверхности
    if((p - C).abs() - R >= 0)
        side = +1;
    else
        side = -1;
    onBorder = false;
    //normal = VECTOR3(0, 0, -1);//#######cheat
}
void AnaliticalSurface_Sphere::findIntersection(const POINT3 &p1, const POINT3 &p2, const Grid::SurfacePositionData &, const Grid::SurfacePositionData &prevSolutionData2,
                                                             Grid::SurfacePositionData &, Grid::SurfacePositionData &solutionData2, POINT3 &intersectionPoint, POINT3 &normal, int &side, bool &found)const
{
    // информация о ближайшей к x2 точке поверхности
    POINT3 x2_nearestPoint;
    POINT3 x2_normal;
    int x2_side;
    bool x2_onBorder;
    findNearestPoint(p2, p2 - p1, prevSolutionData2,
                     solutionData2, x2_nearestPoint, x2_normal, x2_side, x2_onBorder);
    if(x2_side == -1)
    {
        // точка x2 находится с внутренней стороны поверхности
        // есть контакт
        intersectionPoint = x2_nearestPoint;
        normal = x2_normal;
        side = x2_side;
        found = true;
    }
    else
    {
        found = false;
    }
}
void AnaliticalSurface_Sphere::update(const double time, const Grid::RegionDecompositionParameters *set_param)
{
    C[0] = x_fun->fun(time);
    C[1] = y_fun->fun(time);
    C[2] = z_fun->fun(time);
    R = R_fun->fun(time);
}
POINT3 AnaliticalSurface_Sphere::s(const POINT2 &) const
{
    return VECTOR3_NULL;
}
POINT3 AnaliticalSurface_Sphere::difs(const POINT2 &, const DIF_STATE2 &) const
{
    return VECTOR3_NULL;
}



// двумерные функции эрмита
// вычисление локальных индексов одномерных базисных ф-й по локальному индексу двумерной базисной ф-и
#define hXI(i)	(2 * (((i) / 4) % 2) + ((i) % 2))
#define hYI(i)	(2 * ((i) / 8) + (((i) / 2) % 2))
// значение одномерных базисных функций
#define hermite1D_0_dif0(ksi, h) (1 - 3 * (ksi)*(ksi) + 2 * (ksi)*(ksi)*(ksi))
#define hermite1D_1_dif0(ksi, h) (((ksi) - 2 * (ksi)*(ksi) + (ksi)*(ksi)*(ksi))*(h))
#define hermite1D_2_dif0(ksi, h) (3 * (ksi)*(ksi) - 2 * (ksi)*(ksi)*(ksi))
#define hermite1D_3_dif0(ksi, h) ((-(ksi)*(ksi) + (ksi)*(ksi)*(ksi))*(h))
// значение первых производных одномерных базисных функций
#define hermite1D_0_dif1(ksi, h) ((-6*(ksi) + 6*(ksi)*(ksi))/(h))
#define hermite1D_1_dif1(ksi, h) (1 - 4*(ksi) + 3*(ksi)*(ksi))
#define hermite1D_2_dif1(ksi, h) ((6*(ksi) - 6*(ksi)*(ksi))/(h))
#define hermite1D_3_dif1(ksi, h) (-2*(ksi) + 3*(ksi)*(ksi))

InterpolantSurface_Hermite3::InterpolantSurface_Hermite3(const Grid::RegionDecompositionType set_decompositionType, const Interpolant2D_Hermite3 *set_it_x, const Interpolant2D_Hermite3 *set_it_y, const Interpolant2D_Hermite3 *set_it_z, const int set_maxIter, const double set_eps_uv_2d, const double set_eps_uv_1d, const double set_eps_grad, const double set_eps_inerception):
    InterpolantSurface_base(set_decompositionType, set_it_x->gr, set_maxIter, set_eps_uv_2d, set_eps_uv_1d, set_eps_grad, set_eps_inerception)
{
    interpolantRes.clear();
    interpolantRes.resize(set_it_x->res.size()*3);
    for(size_t i = 0; i < set_it_x->res.size(); i++)
    {
        interpolantRes[i*3 + 0] = set_it_x->res[i];
        interpolantRes[i*3 + 1] = set_it_y->res[i];
        interpolantRes[i*3 + 2] = set_it_z->res[i];
    }
    surfaceType = Grid::SurfaceType::InterpolantSurface_Hermite3;
}
POINT3 InterpolantSurface_Hermite3::s(const POINT2 &uv) const
{
    return difs(uv, dif_NULL2);
}
POINT3 InterpolantSurface_Hermite3::difs(const POINT2 &uv, const DIF_STATE2 &dif) const
{
    int i, j;
    POINT2 uv0;     // прямоугольный конечный элемент задается одной вершиной с локальным номером 0
    findFe(uv,
           i, j, uv0);
    int globalBasFunIndex_0 = i * (gr.N[0] + 1) * 4 + j * 4; // глобальный индекс базисной ф-и с локальным индексом 0
    // вычисление значений 4-х одномерных функций по парметрам u и v
    double fun0[4];
    double fun1[4];
    double ksi[2];
    ksi[0] = (uv[0] - uv0[0]) / h[0];
    ksi[1] = (uv[1] - uv0[1]) / h[1];
    if(dif[0] == 0)
    {
        fun0[0] = hermite1D_0_dif0(ksi[0], h[0]);
        fun0[1] = hermite1D_1_dif0(ksi[0], h[0]);
        fun0[2] = hermite1D_2_dif0(ksi[0], h[0]);
        fun0[3] = hermite1D_3_dif0(ksi[0], h[0]);
    }
    else
    {
        fun0[0] = hermite1D_0_dif1(ksi[0], h[0]);
        fun0[1] = hermite1D_1_dif1(ksi[0], h[0]);
        fun0[2] = hermite1D_2_dif1(ksi[0], h[0]);
        fun0[3] = hermite1D_3_dif1(ksi[0], h[0]);
    }

    if(dif[1] == 0)
    {
        fun1[0] = hermite1D_0_dif0(ksi[1], h[1]);
        fun1[1] = hermite1D_1_dif0(ksi[1], h[1]);
        fun1[2] = hermite1D_2_dif0(ksi[1], h[1]);
        fun1[3] = hermite1D_3_dif0(ksi[1], h[1]);
    }
    else
    {
        fun1[0] = hermite1D_0_dif1(ksi[1], h[1]);
        fun1[1] = hermite1D_1_dif1(ksi[1], h[1]);
        fun1[2] = hermite1D_2_dif1(ksi[1], h[1]);
        fun1[3] = hermite1D_3_dif1(ksi[1], h[1]);
    }
    // вычисление значений двумерных функций и интерполянта
    POINT3 result = VECTOR3_NULL;
    /*for(int localBusFunIndex = 0; localBusFunIndex < 16; localBusFunIndex++)
    {
        int globalBasFunIndex;                                    // глобальный индекс базисной ф-и с локальным индексом localBusFunIndex
        if(localBusFunIndex < 8)
            globalBasFunIndex = globalBasFunIndex_0 + localBusFunIndex;
        else
            globalBasFunIndex = globalBasFunIndex_0 + (gr.N[0] + 1) * 4 + (localBusFunIndex - 8);
        double f = fun0[hXI(localBusFunIndex)] *
                   fun1[hYI(localBusFunIndex)];
        result[0] += interpolantRes[globalBasFunIndex*3 + 0] * f;
        result[1] += interpolantRes[globalBasFunIndex*3 + 1] * f;
        result[2] += interpolantRes[globalBasFunIndex*3 + 2] * f;
    }*/
    for(int localBusFunIndex = 0; localBusFunIndex < 8; localBusFunIndex++)
    {
        int globalBasFunIndex = 3*(globalBasFunIndex_0 + localBusFunIndex);// утроенный глобальный индекс базисной ф-и с локальным индексом localBusFunIndex
        double f = fun0[hXI(localBusFunIndex)] *
                   fun1[hYI(localBusFunIndex)];
        result[0] += interpolantRes[globalBasFunIndex + 0] * f;
        result[1] += interpolantRes[globalBasFunIndex + 1] * f;
        result[2] += interpolantRes[globalBasFunIndex + 2] * f;
    }
    for(int localBusFunIndex = 8; localBusFunIndex < 16; localBusFunIndex++)
    {
        int globalBasFunIndex = 3*(globalBasFunIndex_0 + (gr.N[0] + 1) * 4 + (localBusFunIndex - 8));// утроенный глобальный индекс базисной ф-и с локальным индексом localBusFunIndex
        double f = fun0[hXI(localBusFunIndex)] *
                   fun1[hYI(localBusFunIndex)];
        result[0] += interpolantRes[globalBasFunIndex + 0] * f;
        result[1] += interpolantRes[globalBasFunIndex + 1] * f;
        result[2] += interpolantRes[globalBasFunIndex + 2] * f;
    }
    return result;
}
// метод сопряжённых градиентов
void InterpolantSurface_Hermite3::findNearestPoint(const POINT3 &p, const VECTOR3 &, const Grid::SurfacePositionData &prevSolutionData,
                                                   Grid::SurfacePositionData &solutionData, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const
{
    POINT2 uv0;
    // учёт данных о предыдущем решении
    loadSolutionData(prevSolutionData,
                     uv0);
    POINT2 duv;
    POINT2 grad;
    double h = 1.e200;
    for(int iter = 0; iter < maxIter; iter++)
    {
        // сохранение информации с предыдущего шага
        POINT2 uv0_prev = uv0;      // приближение на предыдущем шаге
        POINT2 duv_prev = duv;      // направление на предыдущем шаге
        POINT2 grad_prev = grad;    // градиент на предыдущем шаге
        // вычисление градиента
        grad[0] = (s(uv0) - p) * difs(uv0, Fem::dif_XYZ2[0]) / (s(uv0) - p).abs();
        grad[1] = (s(uv0) - p) * difs(uv0, Fem::dif_XYZ2[1]) / (s(uv0) - p).abs();
        // определение направления одномерного поиска
        if(grad.abs() < eps_grad)
            break;
        if(iter == 0)
        {
            duv = -grad;
        }
        else
        {
            double w = SQR(grad.abs() / grad_prev.abs());
            //if(w > 0.5)
            //    w = 0.5;
            duv = -grad + w*duv_prev;
        }
        POINT2 duv_1 = duv / duv.abs();    // единичное направление
        // определение промежутка одномерного поиска
        double a;
        double b;
        calcInterval(uv0, duv_1,
                      a, b);
        a = 0;
        if(b > h*2)
            b = h*2;
        // поиск минимума по направлению duv_1
        double t = findNearestPoint_find1d_zolot_sechenie_vector(p, uv0, duv_1, a, b, eps_uv_1d);
        uv0 += duv_1*t;
        // проверка достижения точности
        double h_new = (uv0_prev - uv0).abs();
        if(h_new < eps_uv_2d)
            break;
        if(iter == 0)
            h = h_new;
        else
            h = h*0.5 + h_new*0.5;
    }
    // координаты ближайшей точки поверхности
    nearestPoint = s(uv0);
    // проверка нахождения на границе
    onBorder = pointOnBorder(uv0, eps_uv_2d);
    // нормаль в к поверхности
    // ## точно данная поверхность позволяет вычислить производные?
    normal = vector3Mul(difs(uv0, {1, 0}), difs(uv0, {0, 1}));
    normal = normal / normal.abs();
    // с какой стороны от поверхности находится точка p
    if((p - nearestPoint)*normal < 0)
        side = -1;  // точка p находится с внутренней стороны поверхности
    else
        side = +1;  // точка p находится с внешней стороны поверхности
    //normal = (p - nearestPoint)/(p - nearestPoint).abs()*side;
    // сохранение данных о решении
    saveSolutionData(uv0,
                     solutionData);
}
// наискорейший спуск
/*
void InterpolantSurface_Hermite3::findNearestPoint(const POINT3 &p, const VECTOR3 &, const Grid::SurfaceSolutionData prevSolutionData,
                                                   Grid::SurfaceSolutionData &solutionData, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const
{
    POINT2 uv0;
    // учёт данных о предыдущем решении
    loadSolutionData(prevSolutionData,
                     uv0);
    POINT2 duv;
    POINT2 grad;
    double h = 1.e200;
    for(int iter = 0; iter < 20; iter++)
    {
        // сохранение предыдущего приближения
        POINT2 uv0_prev = uv0;
        // вычисление градиента
        grad[0] = (s(uv0) - p) * difs(uv0, Fem::dif_XYZ2[0]) / (s(uv0) - p).abs();
        grad[1] = (s(uv0) - p) * difs(uv0, Fem::dif_XYZ2[1]) / (s(uv0) - p).abs();
        if(grad.abs() < eps_grad)
            break;
        // определение направления одномерного поиска
        POINT2 duv_1 = -grad / grad.abs();    // единичное направление в сторону против градиената
        // определение промежутка одномерного поиска
        double a;
        double b;
        calcInterval(uv0, duv_1,
                      a, b);
        a = 0;
        if(b > h*1.5)
            b = h*1.5;
        // поиск минимума по направлению duv_1
        double t = findNearestPoint_find1d_zolot_sechenie_vector(p, uv0, duv_1, a, b, eps_uv_1d);
        // перемещение точки к минимому
        uv0 += duv_1*t;
        // проверка достижения точности
        h = (uv0_prev - uv0).abs();
        if(h < eps_uv_2d)
            break;
    }
    // координаты ближайшей точки поверхности
    nearestPoint = s(uv0);
    // проверка нахождения на границе
    onBorder = pointOnBorder(uv0, eps_uv_2d);
    // нормаль в к поверхности
    // ## точно данная поверхность позволяет вычислить производные?
    normal = vector3Mul(difs(uv0, {1, 0}), difs(uv0, {0, 1}));
    normal = normal / normal.abs();
    // с какой стороны от поверхности находится точка p
    if((p - nearestPoint)*normal < 0)
        side = -1;  // точка p находится с внутренней стороны поверхности
    else
        side = +1;  // точка p находится с внешней стороны поверхности
    // сохранение данных о решении
    saveSolutionData(uv0,
                     solutionData);
}
*/
// покоординатный поиск
/*void InterpolantSurface_Hermite3::findNearestPoint(const POINT3 &p, const VECTOR3 &, const Grid::SurfaceSolutionData prevSolutionData,
                                                   Grid::SurfaceSolutionData &solutionData, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const
{
    POINT2 uv0;
    // учёт данных о предыдущем решении
    loadSolutionData(prevSolutionData,
                     uv0);
    double hu = 1.e200;
    double hv = 1.e200;
    for(int iter = 0; iter < 20; iter++)
    {
        // сохранение предыдущего приближения
        POINT2 uv0_prev = uv0;
        // поиск по координате u
        {
            double a = MAX(uv0[0] - hu/2*1.5, gr.rect[0][0]);
            double b = MIN(uv0[0] + hu/2*1.5, gr.rect[0][1]);
            double new_u = findNearestPoint_find1d_zolot_sechenie(p, uv0, 0, a, b, eps_uv_1d);
            hu = fabs(new_u - uv0[0]);
            uv0[0] = new_u;
            //POINT2 duv = {1, 0};
            //double a = gr.rect[0][0] - uv0[0];// uv0[0] + 1*a = gr.rect[0][0]
            //double b = gr.rect[0][1] - uv0[0];// uv0[0] + 1*b = gr.rect[0][1]
            //double t = findNearestPoint_find1d_zolot_sechenie_vector(p, uv0, duv, a, b, eps_uv_1d);
            //uv0 += duv*t;
        }
        // поиск по координате v
        {
            double a = MAX(uv0[1] - hv/2*1.5, gr.rect[1][0]);
            double b = MIN(uv0[1] + hv/2*1.5, gr.rect[1][1]);
            double new_v = findNearestPoint_find1d_zolot_sechenie(p, uv0, 1, a, b, eps_uv_1d);
            hv = fabs(new_v - uv0[1]);
            uv0[1] = new_v;
            //POINT2 duv = {0, 1};
            //double a = gr.rect[1][0] - uv0[1];// uv0[1] + 1*a = gr.rect[1][0]
            //double b = gr.rect[1][1] - uv0[1];// uv0[1] + 1*b = gr.rect[1][1]
            //double t = findNearestPoint_find1d_zolot_sechenie_vector(p, uv0, duv, a, b, eps_uv_1d);
            //uv0 += duv*t;
        }
        // проверка достижения точности
        if((uv0_prev - uv0).abs() < eps_uv_2d)
            break;
    }
    // координаты ближайшей точки поверхности
    nearestPoint = s(uv0);
    // проверка нахождения на границе
    onBorder = pointOnBorder(uv0, eps_uv_2d);
    // нормаль в к поверхности
    // ## точно данная поверхность позволяет вычислить производные?
    normal = vector3Mul(difs(uv0, {1, 0}), difs(uv0, {0, 1}));
    normal = normal / normal.abs();
    // с какой стороны от поверхности находится точка p
    if((p - nearestPoint)*normal < 0)
        side = -1;  // точка p находится с внутренней стороны поверхности
    else
        side = +1;  // точка p находится с внешней стороны поверхности
    // сохранение данных о решении
    saveSolutionData(uv0,
                     solutionData);
}*/
// метод ньютона
/*
bool InterpolantSurface::project_newton(const POINT3 &p,
                                 POINT3 &nearestPoint, POINT3 &normal, double &h)const
{
    POINT2 uv;
    POINT2 duv;
    const Grid::GridRectangleRegular2D &gr = it[0]->gr;
    // начальное приближение - центр области
    uv[0] = (gr.rect[0][0] + gr.rect[0][1]) / 2;
    uv[1] = (gr.rect[1][0] + gr.rect[1][1]) / 2;
    // метод Ньютона
    for(int iter = 0; iter < 1000; iter++)
    {
        VECTOR3 s11 = difs(uv, {2, 0});
        VECTOR3 s22 = difs(uv, {0, 2});
        VECTOR3 s12 = difs(uv, {1, 1});
        VECTOR3 s21 = s12;
        VECTOR3 s1 = difs(uv, {1, 0});
        VECTOR3 s2 = difs(uv, {0, 1});
        VECTOR3 s0 = difs(uv, {0, 0});
        // a1*x + b1*y = c1
        // a2*x + b2*y = c2
        double a1 = (p - s0)*s11 - s1*s1;
        double b1 = (p - s0)*s12 - s2*s1;
        double c1 = (p - s0)*s1;
        double a2 = (p - s0)*s21 - s1*s2;
        double b2 = (p - s0)*s22 - s2*s2;
        double c2 = (p - s0)*s2;
        duv[0] = (b1*c2 - b2*c1) / (b1*a2 - b2*a1);
        duv[1] = (a1*c2 - a2*c1) / (a1*b2 - a2*b1);
        duv[0] = duv[0] / 100;
        duv[1] = duv[1] / 100;
        POINT2 uv_new = uv + duv;
        // запрет выхода за область определения сплайна
        if(uv_new[0] < gr.rect[0][0])
            uv_new[0] = gr.rect[0][0];
        if(uv_new[0] > gr.rect[0][1])
            uv_new[0] = gr.rect[0][1];
        if(uv_new[1] < gr.rect[1][0])
            uv_new[1] = gr.rect[1][0];
        if(uv_new[1] > gr.rect[1][1])
            uv_new[1] = gr.rect[1][1];
        uv = uv_new;
        //if(fabs(uv_new[0] - uv[0]) < 0.000001 &&
        //   fabs(uv_new[1] - uv[1]) < 0.000001)
        //    break;
    }
    nearestPoint = s(uv);
    normal = vector3Mul(difs(uv, {1, 0}), difs(uv, {0, 1}));
    normal = normal / normal.abs();
    // предполагаем что p - nearestPoint || normal
    h = (p - nearestPoint)*normal;
    return false;
}
*/

// двумерные функции лагранжа
// вычисление локальных индексов одномерных базисных ф-й по локальному индексу двумерной базисной ф-и
#define lXI(i)	((i) % 4)
#define lYI(i)	((i) / 4)
// значение одномерных базисных функций
#define lagrange1D_0_dif0(ksi, h) (-9./2. *((ksi) - 1./3.)*((ksi) - 2./3.)*((ksi) - 1.))
#define lagrange1D_1_dif0(ksi, h) (27./2. *(ksi)          *((ksi) - 2./3.)*((ksi) - 1.))
#define lagrange1D_2_dif0(ksi, h) (-27./2.*(ksi)          *((ksi) - 1./3.)*((ksi) - 1.))
#define lagrange1D_3_dif0(ksi, h) (9./2.  *(ksi)          *((ksi) - 1./3.)*((ksi) - 2./3.))
// значение первых производных одномерных базисных функций
#define lagrange1D_0_dif1(ksi, h) ((-5.5 + 18*(ksi) - 13.5*(ksi)*(ksi))/(h))
#define lagrange1D_1_dif1(ksi, h) ((9    - 45*(ksi) + 40.5*(ksi)*(ksi))/(h))
#define lagrange1D_2_dif1(ksi, h) ((-4.5 + 36*(ksi) - 40.5*(ksi)*(ksi))/(h))
#define lagrange1D_3_dif1(ksi, h) ((1    - 9*(ksi)  + 13.5*(ksi)*(ksi))/(h))

InterpolantSurface_Lagrange3::InterpolantSurface_Lagrange3(const Grid::RegionDecompositionType set_decompositionType, const Interpolant2D_Lagrange3 *set_it_x, const Interpolant2D_Lagrange3 *set_it_y, const Interpolant2D_Lagrange3 *set_it_z, const int set_maxIter, const double set_eps_uv_2d, const double set_eps_uv_1d, const double set_eps_grad, const double set_eps_inerception):
    InterpolantSurface_base(set_decompositionType, set_it_x->gr, set_maxIter, set_eps_uv_2d, set_eps_uv_1d, set_eps_grad, set_eps_inerception)
{
    interpolantRes.clear();
    interpolantRes.resize(set_it_x->res.size()*3);
    for(size_t i = 0; i < set_it_x->res.size(); i++)
    {
        interpolantRes[i*3 + 0] = set_it_x->res[i];
        interpolantRes[i*3 + 1] = set_it_y->res[i];
        interpolantRes[i*3 + 2] = set_it_z->res[i];
    }
    surfaceType = Grid::SurfaceType::InterpolantSurface_Lagrange3;
}
POINT3 InterpolantSurface_Lagrange3::s(const POINT2 &uv) const
{
    return difs(uv, Fem::dif_NULL2);
}
POINT3 InterpolantSurface_Lagrange3::difs(const POINT2 &uv, const DIF_STATE2 &dif) const
{
    int i, j;
    POINT2 uv0;     // прямоугольный конечный элемент задается одной вершиной с локальным номером 0
    findFe(uv,
           i, j, uv0);
    int globalBasFunIndex_0 = i * 3 * (gr.N[0] * 3 + 1) + j * 3; // глобальный индекс базисной ф-и с локальным индексом 0
    // вычисление значений 4-х одномерных функций по парметрам u и v
    double fun0[4];
    double fun1[4];
    double ksi[2];
    ksi[0] = (uv[0] - uv0[0]) / h[0];
    ksi[1] = (uv[1] - uv0[1]) / h[1];
    if(dif[0] == 0)
    {
        fun0[0] = lagrange1D_0_dif0(ksi[0], h[0]);
        fun0[1] = lagrange1D_1_dif0(ksi[0], h[0]);
        fun0[2] = lagrange1D_2_dif0(ksi[0], h[0]);
        fun0[3] = lagrange1D_3_dif0(ksi[0], h[0]);
    }
    else
    {
        fun0[0] = lagrange1D_0_dif1(ksi[0], h[0]);
        fun0[1] = lagrange1D_1_dif1(ksi[0], h[0]);
        fun0[2] = lagrange1D_2_dif1(ksi[0], h[0]);
        fun0[3] = lagrange1D_3_dif1(ksi[0], h[0]);
    }

    if(dif[1] == 0)
    {
        fun1[0] = lagrange1D_0_dif0(ksi[1], h[1]);
        fun1[1] = lagrange1D_1_dif0(ksi[1], h[1]);
        fun1[2] = lagrange1D_2_dif0(ksi[1], h[1]);
        fun1[3] = lagrange1D_3_dif0(ksi[1], h[1]);
    }
    else
    {
        fun1[0] = lagrange1D_0_dif1(ksi[1], h[1]);
        fun1[1] = lagrange1D_1_dif1(ksi[1], h[1]);
        fun1[2] = lagrange1D_2_dif1(ksi[1], h[1]);
        fun1[3] = lagrange1D_3_dif1(ksi[1], h[1]);
    }
    // вычисление значений двумерных функций и интерполянта
    POINT3 result = VECTOR3_NULL;
    for(int localBusFunIndex = 0; localBusFunIndex < 16; localBusFunIndex++)
    {
        int dj = lXI(localBusFunIndex);
        int di = lYI(localBusFunIndex);
        int globalBasFunIndex = 3*(globalBasFunIndex_0 + di * (gr.N[0] * 3 + 1) + dj); // глобальный индекс базисной ф-и с локальным индексом localBusFunIndex
        double f = fun0[dj] *
                   fun1[di];
        result[0] += interpolantRes[globalBasFunIndex + 0] * f;
        result[1] += interpolantRes[globalBasFunIndex + 1] * f;
        result[2] += interpolantRes[globalBasFunIndex + 2] * f;
    }
    return result;
}
// покоординатный поиск
void InterpolantSurface_Lagrange3::findNearestPoint(const POINT3 &p, const VECTOR3 &, const Grid::SurfacePositionData &prevSolutionData,
                                                    Grid::SurfacePositionData &solutionData, POINT3 &nearestPoint, POINT3 &normal, int &side, bool &onBorder)const
{
    POINT2 uv0;
    // учёт данных о предыдущем решении
    loadSolutionData(prevSolutionData,
                     uv0);
    double h = 1.e200;
    for(int iter = 0; iter < maxIter; iter++)
    {
        POINT2 uv0_prev = uv0;
        double hu_new;
        double hv_new;
        // поиск по координате u
        {
            double a = MAX(uv0[0] - h*2, gr.rect[0][0]);
            double b = MIN(uv0[0] + h*2, gr.rect[0][1]);
            double new_u = findNearestPoint_find1d_zolot_sechenie(p, uv0, 0, a, b, eps_uv_1d);
            hu_new = fabs(new_u - uv0[0]);
            uv0[0] = new_u;
        }
        // поиск по координате v
        {
            double a = MAX(uv0[1] - h*2, gr.rect[1][0]);
            double b = MIN(uv0[1] + h*2, gr.rect[1][1]);
            double new_v = findNearestPoint_find1d_zolot_sechenie(p, uv0, 1, a, b, eps_uv_1d);
            hv_new = fabs(new_v - uv0[1]);
            uv0[1] = new_v;
        }
        if((uv0_prev - uv0).abs() < eps_uv_2d)
            break;
        double h_new = sqrt(hu_new*hu_new + hv_new*hv_new);
        if(iter == 0)
            h = h_new;
        else
            h = h*0.5 + h_new*0.5;
        //hu = fabs(new_u - uv0[0]);
        //hv = fabs(new_v - uv0[1]);
    }
    // координаты ближайшей точки поверхности
    nearestPoint = s(uv0);
    // проверка нахождения на границе
    onBorder = pointOnBorder(uv0, eps_uv_2d);
    // нормаль в к поверхности
    // ## точно данная поверхность позволяет вычислить производные?
    normal = vector3Mul(difs(uv0, {1, 0}), difs(uv0, {0, 1}));
    normal = normal / normal.abs();
    // с какой стороны от поверхности находится точка p
    if((p - nearestPoint)*normal < 0)
        side = -1;  // точка p находится с внутренней стороны поверхности
    else
        side = +1;  // точка p находится с внешней стороны поверхности
    // сохранение данных о решении
    saveSolutionData(uv0,
                     solutionData);
}

}   // namespace Interpolation
