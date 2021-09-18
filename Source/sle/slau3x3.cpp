#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"

#include "timeinterval.h"
#include "slau3x3.h"

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl.h"

namespace SlauSolving
{

// пустой предобусловливатель (отсутствие)
struct SSCM3x3Preconditioner_none: public SSCM3x3Preconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCM3x3Elements &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector3 &x)const final;
    virtual void InversedQmulVector(Vector3 &x)const final;
};
// диагональный предобусловливатель (не полный)
struct SSCM3x3Preconditioner_diag: public SSCM3x3Preconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCM3x3Elements &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    // saveBMP ничего не делает
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector3 &x)const final;
    virtual void InversedQmulVector(Vector3 &x)const final;
private:
    std::vector<MATR3x3> d; // диагональные элементы
};
/*
// LLT предобусловливатель, симметричный строчно-столбцовый формат, портрет не меняется (не полный)
struct SSCMPreconditioner3x3_SSCM_LLT_0: public SSCMPreconditioner3x3_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCMElements3x3 &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector3 &x)const final;
    virtual void InversedQmulVector(Vector3 &x)const final;
private:
    SSCMElements3x3 m; // элементы матрицы L, партрет такой же как и у исходной
};
// LDLT предобусловливатель, симметричный строчно-столбцовый формат, портрет не меняется (не полный)
struct SSCMPreconditioner3x3_SSCM_LDLT_0: public SSCMPreconditioner3x3_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCMElements3x3 &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector3 &x)const final;
    virtual void InversedQmulVector(Vector3 &x)const final;
private:
    SSCMElements3x3 m; // элементы матрицы L, партрет такой же как и у исходной
};
// LLT предобусловливатель, симметричный профильный формат, портрет расширяется (полный)
struct SSCMPreconditioner3x3_Profile_LLT: public SSCMPreconditioner3x3_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCMElements3x3 &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector3 &x)const final;
    virtual void InversedQmulVector(Vector3 &x)const final;
private:
    ProfilePortrait p_profile;  // профильный портрет матрицы полного предобусловливателя
    SSCMElements3x3 m; // элементы матрицы L
};
// LLT предобусловливатель, симметричный строчно-столбцовый формат, портрет расширяется (полный)
struct SSCMPreconditioner3x3_SSCM_LLT: public SSCMPreconditioner3x3_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCMElements3x3 &set_e, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector3 &x)const final;
    virtual void InversedQmulVector(Vector3 &x)const final;
private:
    SSCMPortrait p_SSCM; // строчно-столбцовый портрет матрицы полного предобусловливателя(может получиться расширенным)
    SSCMElements3x3 m_SSCM;
    std::vector<MATR3x3> s; // строка - временный массив
};
*/
// LDLT предобусловливатель, симметричный строчно-столбцовый формат, портрет расширяется (полный)
struct SSCM3x3Preconditioner_SSCM_LDLT: public SSCM3x3Preconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCM3x3Elements &set_e0, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector3 &x)const final;
    virtual void InversedQmulVector(Vector3 &x)const final;
private:
    SSCMPortrait p; // строчно-столбцовый портрет матрицы полного предобусловливателя(может получиться расширенным)
    //SSCMElements3x3 e;
    std::vector<MATR3x3> a; // элементы нижнего треугольника (доступ по индексу)
    std::vector<MATR3x3> amuld; // элементы нижнего треугольника, домноженные на элементы диагонали (доступ по индексу)
    std::vector<MATR3x3_Pair> d; // диагональные элементы + обратные
    std::vector<size_t> s; // строка - временный массив ссылок на блоки или -1 если блок нулевой
};
// LDLT предобусловливатель, симметричный профильный формат, портрет расширяется (полный)
struct SSCM3x3Preconditioner_Profile_LDLT: public SSCM3x3Preconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCM3x3Elements &set_e0, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector3 &x)const final;
    virtual void InversedQmulVector(Vector3 &x)const final;
private:
    ProfilePortrait p; // строчно-столбцовый портрет матрицы полного предобусловливателя(может получиться расширенным)
    std::vector<MATR3x3> a; // элементы нижнего треугольника (доступ по индексу)
    std::vector<MATR3x3> amuld; // элементы нижнего треугольника, домноженные на элементы диагонали (доступ по индексу)
    std::vector<MATR3x3_Pair> d; // диагональные элементы + обратные
};

// прямой решатель (предобусловие должно быть полное)
struct SSCM3x3Solver_direct: public SSCM3x3Solver_base
{
    virtual void init(const size_t matrixSize) final;
    virtual void release() final;
    virtual void solve(const SSCM3x3 &matrix, const Vector3 &b, const SSCM3x3Preconditioner_base *p, const Vector3 &x0, const SolverParameters &ssp,
                       Vector3 &x, double &residual, double &relativeResidual, int &iterations, double &time) final;
private:
    Vector3 r;
    Vector3 r0;
};
// итерационный решатель (предобусловие может быть не полное)
struct SSCM3x3Solver_iterative_base: public SSCM3x3Solver_base
{
    virtual void solve(const SSCM3x3 &matrix, const Vector3 &b, const SSCM3x3Preconditioner_base *p, const Vector3 &x0, const SolverParameters &ssp,
                       Vector3 &x, double &residual, double &relativeResidual, int &iterations, double &time) final;
protected:
    virtual void start(const SSCM3x3 &matrix, const Vector3 &b, const SSCM3x3Preconditioner_base *p, Vector3 &x) = 0;
    virtual void iter(const SSCM3x3 &matrix, const Vector3 &, const SSCM3x3Preconditioner_base *p, Vector3 &x) = 0;
    Vector3 r;
    Vector3 r0;
};
// итерационный решатель LOS
struct SSCM3x3Solver_LOS: public SSCM3x3Solver_iterative_base
{
    virtual void init(const size_t matrixSize) final;
    virtual void release() final;
private:
    virtual void start(const SSCM3x3 &matrix, const Vector3 &b, const SSCM3x3Preconditioner_base *p, Vector3 &x) final;
    virtual void iter(const SSCM3x3 &matrix, const Vector3 &, const SSCM3x3Preconditioner_base *p, Vector3 &x) final;
    Vector3 ss;
    Vector3 z;
    Vector3 w;
    Vector3 aa;
    Vector3 pp;
};
// итерационный решатель CGM
struct SSCM3x3Solver_CGM: public SSCM3x3Solver_iterative_base
{
    virtual void init(const size_t matrixSize) final;
    virtual void release() final;
private:
    virtual void start(const SSCM3x3 &matrix, const Vector3 &b, const SSCM3x3Preconditioner_base *p, Vector3 &x) final;
    virtual void iter(const SSCM3x3 &matrix, const Vector3 &, const SSCM3x3Preconditioner_base *p, Vector3 &x) final;
    Vector3 ss;
    Vector3 z;
    Vector3 w;
    Vector3 aa;
    Vector3 pp;
};
// предобусловливатель PARDISO, BSR 3x3
struct SSCM3x3Preconditioner_PARDISO: public SSCM3x3Preconditioner_base
{
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) final;
    virtual void updatePreconditioner(const SSCM3x3Elements &set_e0, const size_t firstStr, double &time) final;
    virtual void release() final;
    virtual void saveBMP(const char *fileName, size_t size)const final;
    virtual void saveProperties(const char *fileName)const final;
    virtual size_t getElementsNumber()const final;
    virtual void InversedSmulVector(Vector3 &x) const final;
    virtual void InversedQmulVector(Vector3 &x) const final;
private:
    size_t matrixSize;
    size_t elementsNumber;
    MKL_INT *ia;
    MKL_INT *ja;
    double *a;
    double *b;
    double *x;
    int *k_convert;
    // переменные для PARDISO
    MKL_INT n;
    const MKL_INT block_size = 3;
    MKL_INT mtype = -2;   // Real symmetric matrix
    // RHS and solution vectors
    MKL_INT nrhs = 1;     // Number of right hand sides
    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    mutable void *pt[64];
    // Pardiso control parameters
    mutable MKL_INT iparm[64];
    mutable MKL_INT maxfct, mnum, phase, error, msglvl;
    // Auxiliary variables
    double ddum;          // Double dummy
    mutable MKL_INT idum;         // Integer dummy
};


SSCM3x3Elements::SSCM3x3Elements(const size_t elementsNumber, const size_t matrixSize)
{
    init(elementsNumber, matrixSize);
}
void SSCM3x3Elements::init(const size_t elementsNumber, const size_t matrixSize)
{
    a.resize(elementsNumber);
    d.resize(matrixSize);
}
void SSCM3x3Elements::fill(const MATR3x3 value)
{
    for(size_t i = 0; i < a.size(); i++)
        a[i] = value;
    for(size_t i = 0; i < d.size(); i++)
        d[i] = value;
}
void SSCM3x3Elements::release()
{
    a.clear();
    d.clear();
}

void SSCM3x3::init()
{
    p = new SSCMPortrait;
    e = new SSCM3x3Elements;
}
void SSCM3x3::initPortrait3x3by1x1AndAllocateMemory(const SSCMPortrait &p1x1)
{
    const size_t matrixSize1x1 = p1x1.getMatrixSize();
    const size_t matrixSize = matrixSize1x1 / 3;
    SSCMPortrait &p = *this->p;
    SSCM3x3Elements &e = *this->e;
    // сборка портрета
    SSCMPortraitBulder pb;
    pb.init(matrixSize);
    for (size_t i1x1 = 0; i1x1 < matrixSize1x1; i1x1++) // i1x1 - номер строки
    {
        const size_t k1x1_min = p1x1.ind[i1x1];
        const size_t k1x1_max = p1x1.ind[i1x1 + 1];  // строго меньше
        size_t k1x1 = k1x1_min;
        while(k1x1 < k1x1_max)
        {
            size_t j1x1 = p1x1.ai[k1x1];// i1x1, j1x1 - ненулевая ячейка исходной матрицы 1х1
            size_t i = i1x1/3;
            size_t j = j1x1/3;
            pb.addElement(i, j);// i, j - ненулевая ячейка матрицы 3х3
            k1x1++;
        }
    }
    pb.completePortrait(p);
    e.init(p.getElementsNumber(), matrixSize);
}
void SSCM3x3::initMatrixByMatrix1x1(const SSCMPortrait &p1x1, const SSCMElements &e1x1)
{
    const size_t matrixSize1x1 = p1x1.getMatrixSize();
    //const size_t matrixSize = matrixSize1x1 / 3;
    SSCMPortrait &p = *this->p;
    SSCM3x3Elements &e = *this->e;
    // очистка
    {
        MATR3x3 t;
        t.clear();
        e.fill(t);
    }
    // копирование матрицы
//fprintf(stderr, "3x3 e.a.size() = %d e.d.size() = %d\n", (int)e.a.size(), (int)e.d.size());
    for (size_t i1x1 = 0; i1x1 < matrixSize1x1; i1x1++) // i1x1 - номер строки
    {
        // вне диагонали
        const size_t k1x1_min = p1x1.ind[i1x1];
        const size_t k1x1_max = p1x1.ind[i1x1 + 1];  // строго меньше
        size_t k1x1 = k1x1_min;
        while(k1x1 < k1x1_max)
        {
            size_t j1x1 = p1x1.ai[k1x1];
            double a1x1 = e1x1.a[k1x1];    // a = m0[i0][j0]
            size_t i = i1x1/3;
            size_t j = j1x1/3;
            if(i == j)
            {
                // диагональ
                e.d[i].m[i1x1%3][j1x1%3] = a1x1;
                e.d[i].m[j1x1%3][i1x1%3] = a1x1;
            }
            else
            {
                // матрица
                size_t k = p.findSorted(i, j);
                e.a[k].m[i1x1%3][j1x1%3] = a1x1;
                //fprintf(stderr, "k =%d\n", (int)k);
            }
            k1x1++;
        }
        // диагональ
        {
            double d1x1 = e1x1.d[i1x1];
            size_t i = i1x1/3;
            e.d[i].m[i1x1%3][i1x1%3] = d1x1;
        }
    }
}
/*void SSCM3x3::init3x3by1x1(const SSCMPortrait &p1x1, const SSCMElements &e1x1)
{
    const size_t matrixSize1x1 = p1x1.getMatrixSize();
    const size_t matrixSize = matrixSize1x1 / 3;
    SSCMPortrait &p = *this->p;
    SSCMElements3x3 &e = *this->e;
    // 1) сборка портрета
    SSCMPortraitBulder pb;
    pb.init(matrixSize);
    for (size_t i1x1 = 0; i1x1 < matrixSize1x1; i1x1++) // i1x1 - номер строки
    {
        const size_t k1x1_min = p1x1.ind[i1x1];
        const size_t k1x1_max = p1x1.ind[i1x1 + 1];  // строго меньше
        size_t k1x1 = k1x1_min;
        while(k1x1 < k1x1_max)
        {
            size_t j1x1 = p1x1.ai[k1x1];// i1x1, j1x1 - ненулевая ячейка исходной матрицы 1х1
            size_t i = i1x1/3;
            size_t j = j1x1/3;
            pb.addElement(i, j);// i, j - ненулевая ячейка матрицы 3х3
            k1x1++;
        }
    }
    pb.completePortrait(p);
    // 2) копирование матрицы
    e.init(p.getElementsNumber(), matrixSize);
//fprintf(stderr, "3x3 e.a.size() = %d e.d.size() = %d\n", (int)e.a.size(), (int)e.d.size());
    for (size_t i1x1 = 0; i1x1 < matrixSize1x1; i1x1++) // i1x1 - номер строки
    {
        // вне диагонали
        const size_t k1x1_min = p1x1.ind[i1x1];
        const size_t k1x1_max = p1x1.ind[i1x1 + 1];  // строго меньше
        size_t k1x1 = k1x1_min;
        while(k1x1 < k1x1_max)
        {
            size_t j1x1 = p1x1.ai[k1x1];
            double a1x1 = e1x1.a[k1x1];    // a = m0[i0][j0]
            size_t i = i1x1/3;
            size_t j = j1x1/3;
            if(i == j)
            {
                // диагональ
                e.d[i].m[i1x1%3][j1x1%3] = a1x1;
                e.d[i].m[j1x1%3][i1x1%3] = a1x1;
            }
            else
            {
                // матрица
                size_t k = p.findSorted(i, j);
                e.a[k].m[i1x1%3][j1x1%3] = a1x1;
                //fprintf(stderr, "k =%d\n", (int)k);
            }
            k1x1++;
        }
        // диагональ
        {
            double d1x1 = e1x1.d[i1x1];
            size_t i = i1x1/3;
            e.d[i].m[i1x1%3][i1x1%3] = d1x1;
        }
    }
}*/
void SSCM3x3::release()
{
    p->release();
    delete p;
    e->release();
    delete e;
}
size_t SSCM3x3::getMatrixSize()const
{
    return p->matrixSize;
}
size_t SSCM3x3::getElementsNumber()const
{
    return p->elementsNumber;
}

void Vector1x1CopyToVector3x3(const Vector &x, Vector3 &y)
{
    for (size_t i = 0; i < x.size(); i++)
        y[i/3][i%3] = x[i];
}
void Vector3x3CopyToVector1x1(const Vector3 &x, Vector &y)
{
    for (size_t i = 0; i < y.size(); i++)
        y[i] = x[i/3][i%3];
}
void SSCM3x3copy(const SSCM3x3Elements &e1, SSCM3x3Elements &e2)
{
    if(e2.a.size() < e1.a.size())
        e2.a.resize(e1.a.size());
    if(e2.d.size() < e1.d.size())
        e2.d.resize(e1.d.size());
    for (size_t i = 0; i < e1.a.size(); i++)
        e2.a[i] = e1.a[i];
    for (size_t i = 0; i < e1.d.size(); i++)
        e2.d[i] = e1.d[i];
}
void SSCM3x3mulScalar(SSCM3x3 &matrix, const double c)
{
    SSCM3x3Elements &e = *matrix.e;
    for (size_t i = 0; i < e.a.size(); i++)
        e.a[i] *= c;
    for (size_t i = 0; i < e.d.size(); i++)
        e.d[i] *= c;
}
void SSCM3x3mulVector3(const SSCM3x3 &matrix, const Vector3 &x, Vector3 &y)
{
    SSCMPortrait &p = *matrix.p;
    SSCM3x3Elements &e = *matrix.e;
    for (size_t i = 0; i < p.matrixSize; i++)
        y[i].clear();
    // вне диагонали
    for (size_t i = 0; i < p.matrixSize; i++)           // i - номер строки
        for (size_t k = p.ind[i]; k < p.ind[i + 1]; k++)
        {                             // k - индекс элемента i-й строки
            size_t j = p.ai[k];          // j - номер столбца
                                      // A(i, j) = e.a[k]
            MATR3x3::a_mul_b_add(e.a[k], x[j], y[i]);
            //y[i] += e.a[k] * x[j];    // нижний треугольник
            MATR3x3::aT_mul_b_add(e.a[k], x[i],
                                          y[j]);
            //y[j] += e.a[k].transpose() * x[i];    // верхний треугольник
        }
    // диагональ
    for (size_t i = 0; i < p.matrixSize; i++)
        y[i] += e.d[i] * x[i];
}
void SSCM3x3calcResidual(const SSCM3x3 &matrix, const Vector3 &b, const Vector3 &x, Vector3 &r)
{
    SSCM3x3mulVector3(matrix, x, r);            // A*x -> r
    Vector1PlusCmulVector2(b, -1, r, r);    // b - r -> r
}
 void SSCM3x3addBoundaryCondition1(SSCM3x3 &matrix, Vector3 &b, const std::vector<double>&u0, const std::vector<bool> &state)
{
    SSCMPortrait &p = *matrix.p;
    SSCM3x3Elements &e = *matrix.e;
    // вне диагонали
    for (size_t i_3x3 = 0; i_3x3 < p.matrixSize; i_3x3++) // i_3x3 - номер строки
        for (size_t k = p.ind[i_3x3]; k < p.ind[i_3x3 + 1]; k++)
        {                             // k - индекс элемента i_3x3-й строки
            size_t j_3x3 = p.ai[k];   // j_3x3 - номер столбца
                                      // A(i_3x3, j_3x3) = e.a[k]
            for(size_t bi = 0; bi < 3; bi++)
            {
                for(size_t bj = 0; bj < 3; bj++)
                {
                    size_t i = i_3x3*3 + bi;
                    size_t j = j_3x3*3 + bj;    // (i,j) - координаты исходной матрицы (не блочной)
                    // b[i_3x3].x[bi] = b[i] - выделение элемента блочного вектора
                    // b[j_3x3].x[bj] = b[j] - выделение элемента блочного вектора
                    // e.a[k].m[bi][bj] = m[i][j] - выделение элемента блочной матрицы
                    double &m_el_ij = e.a[k].m[bi][bj]; // m[i][j] = e.a[k]
                    double &b_el_i = b[i_3x3].x[bi];    // b[i]
                    double &b_el_j = b[j_3x3].x[bj];    // b[j]
                    if(state[i])
                    {

                        b_el_j -= u0[i]*m_el_ij;
                        m_el_ij = 0;
                    }
                    if(state[j])
                    {
                        b_el_i -= u0[j]*m_el_ij;
                        m_el_ij = 0;
                    }
                }
            }
            /*if(state[i])
            {
                b[j] -= u0[i]*e.a[k];   //j - номер строки зеркального элемента к e.a[k]
                e.a[k] = 0;
            }
            if(state[j])
            {
                b[i] -= u0[j]*e.a[k];
                e.a[k] = 0;
            }*/
        }
    // диагональ
    for (size_t i_3x3 = 0; i_3x3 < p.matrixSize; i_3x3++)
    {
        for(size_t bi = 0; bi < 3; bi++)
        {
            size_t i = i_3x3*3 + bi;    // координата исходного вектора (не блочного)
            if(state[i])
            {
                double &d_el_ii = e.d[i_3x3].m[bi][bi]; // d[i] = e.d[i]
                double &b_el_i = b[i_3x3].x[bi]; // b[i]
                d_el_ii = 1;
                b_el_i = u0[i];
                // плюс нужно обнулить внедиагональные элементы строки bi и столбца bi блока e.d[i_3x3]
                for(size_t bj = 0; bj < 3; bj++)
                {
                    if(bj != bi)
                    {
                        // вне диагонали
                        e.d[i_3x3].m[bi][bj] = 0;   // строка
                        if(bj > bi)// при обнулении элементов столбца наже диагонали меняется вектор правой части
                        {
                            // столбец, ниже диагонали
                            // для неоднородных условий
                            //size_t j = i_3x3*3 + bj;
                            double &b_el_j = b[i_3x3].x[bj]; // b[j]
                            b_el_j -= u0[i]*e.d[i_3x3].m[bj][bi];
                        }
                        e.d[i_3x3].m[bj][bi] = 0;   // столбец
                    }
                }
            }
        }
        /*if(state[i])
        {
            e.d[i] = 1;
            b[i] = u0[i];
        }*/
    }
}
void SSCM3x3_1addSimilarM3x3_2(SSCM3x3 &matrix1, SSCM3x3 &matrix2)
{
    SSCM3x3Elements &e1 = *matrix1.e;
    SSCM3x3Elements &e2 = *matrix2.e;
    for (size_t i = 0; i < e1.a.size(); i++)
        e1.a[i] += e2.a[i];
    for (size_t i = 0; i < e1.d.size(); i++)
        e1.d[i] += e2.d[i];
}
void SSCM3x3_1addEnclosedM3x3_2mulScalar(SSCM3x3 &matrix1, SSCM3x3 &matrix2, const double c)
{
    SSCMPortrait &p1 = *matrix1.p;
    SSCM3x3Elements &e1 = *matrix1.e;
    SSCMPortrait &p2 = *matrix2.p;
    SSCM3x3Elements &e2 = *matrix2.e;
    //FILE *out = fopen("_", "w");
    size_t matrixSize = p1.matrixSize;         // размеры матриц должны быть одинаковыми
    // вне диагонали
    for (size_t i = 0; i < matrixSize; i++)    // i - номер строки обоих матриц
    {
        size_t k1 = p1.ind[i];                 // k1 - индекс элемента i-й строки матрицы 1
        size_t k2 = p2.ind[i];                 // k2 - индекс элемента i-й строки матрицы 2
        for (; k2 < p2.ind[i + 1]; k2++)    // обходим элементы матрицы 2
        {
            size_t j2 = p2.ai[k2];             // j2 - номер столбца элемента e2.a[k2]
                                            // M2(i, j2) = e2.a[k2]
            // поиск элемента (i, j2) в матрице M1
            // он должен существовать ввиду вложенности портретов поэтому проверка k1 < p1.ind[i + 1] не нужна
            for(;;)
            {
                size_t j1 = p1.ai[k1];   // j1 - номер столбца элемента e1.a[k1]
                                      // M1(i, j1) = e1.a[k1]
                if(j1 == j2)
                {
                    // M1(i, j1) += M2(i, j2)
                    e1.a[k1] += e2.a[k2] * c;
                    k1++;
                    break;
                }
                k1++;
            }
        }
    }
    //fclose(out);
    // диагональ
    for (size_t i = 0; i < matrixSize; i++)
        e1.d[i] += e2.d[i] * c;
}

// копирование строки i из матрицы p2, e2 в матрицу p1, e1
// матрицы должны быть одного размера, p2 вложен в p1
// элементы которых нет в портрете p2 задаются нулями в матрице p1
void SSCM3x3StrCopyToSSCM3x3Str(const size_t i, const SSCMPortrait &p1, std::vector<MATR3x3> &a1, std::vector<MATR3x3_Pair> &d1,
                                                const SSCMPortrait &p2, const std::vector<MATR3x3> &a2, const std::vector<MATR3x3> &d2)
{
    // вне диагонали
    size_t k1 = p1.ind[i];          // k1 - индекс элемента i-й строки матрицы 1
    size_t k1_max = p1.ind[i + 1];
    size_t k2 = p2.ind[i];          // k2 - индекс элемента i-й строки матрицы 2
    size_t k2_max = p2.ind[i + 1];
    while(k2 < k2_max) // обход элементов строки матрицы 2
    {
        size_t j2 = p2.ai[k2];      // j2 - номер столбца элемента e2.a[k2]
                                    // M2(i, j2) = e2.a[k2]
        // поиск элемента (i, j2) в матрице M1
        // он должен существовать ввиду вложенности портретов поэтому проверка k1 < p1.ind[i + 1] не нужна
        for(;;) // обход элементов строки матрицы 1
        {
            size_t j1 = p1.ai[k1];  // j1 - номер столбца элемента e1.a[k1]
                                    // M1(i, j1) = e1.a[k1]
            if(j1 == j2)
            {
                // M1(i, j1) = M2(i, j2)
                a1[k1] = a2[k2];
                k1++;
                k2++;
                break;
            }
            else
            {
                a1[k1].clear();
                k1++;
            }
        }
    }
    // заполнение нулями остатка строки матрицы M1
    while(k1 < k1_max)  // обход элементов матрицы 1
    {
        a1[k1].clear();
        k1++;
    }
    // диагональ
    d1[i].value = d2[i];
}
// копирование строки i из SSCM3x3 матрицы p2, e2 в профильную 3x3 матрицу p1, e1
// матрицы должны быть одного размера, p2 вложен в p1
// элементы которых нет в портрете p2 задаются нулями в матрице p1
void SSCM3x3StrCopyToProfileStr(const size_t i, const ProfilePortrait &p1, std::vector<MATR3x3> &a1, std::vector<MATR3x3_Pair> &d1,
                                             const SSCMPortrait &p2, const std::vector<MATR3x3> &a2, const std::vector<MATR3x3> &d2)
{
    // копия матрицы
    if(p2.ind[i + 1] - p2.ind[i] != 0)  // строка не пустая
    {
        size_t k1_min = p1.ind[i];
        size_t k1_max = p1.ind[i + 1];
        size_t k2_min = p2.ind[i];
        size_t k2_max = p2.ind[i + 1];
        size_t j_min = p1.getFirstCol(i);// j_min - номер столбца первого ненулевого элемента строки i
        //size_t stringSize = i - j_min;     // в строке i плотный профиль занимает stringSize элементов
        for(size_t k1 = k1_min; k1 < k1_max; k1++)
            a1[k1].clear();     // a1[k1] = L(i, i-stringSize + k1-p1.ind[i])
        for(size_t k2 = k2_min; k2 < k2_max; k2++)
        {
            size_t j = p2.ai[k2];           // a2[k2] = A(i, j)
            size_t k1 = k1_min + (j - j_min); // a1[k1] = L(i, j)
            a1[k1] = a2[k2];
        }
    }
    d1[i].value = d2[i];
}
// вычисление скалярного произведения строки i матрицы (p, a, d) на плотный вектор s(ссылки на блоки или -1), домножается дополнительно на d[]
// s_j_min - первый ненулевой столбец строки, соответствующей вектору s
MATR3x3 SSCM3x3StrScalarMulVecto3MulD(const size_t i, const SSCMPortrait &p, const std::vector<MATR3x3> &a, const std::vector<MATR3x3> &amuld, const std::vector<size_t> &s, const size_t s_j_min)
{
    // справа налево
    MATR3x3 E;
    E.clear();
    size_t k_min = p.ind[i];
    size_t k_max = p.ind[i + 1];
    if(k_min == k_max)
        return E;
    size_t m_j_min = p.ai[k_min];
    size_t j_min = MAX(s_j_min, m_j_min); // столбцы j < j_min нулевые для строки матрицы и строки вектора s
    size_t k = k_max - 1;
    for(;;) // обход элементов строки
    {
        size_t j = p.ai[k]; // j - номер столбца элемента e.a[k]
        if(s[j] != (size_t)-1)
        {
            MATR3x3::a_Mul_bT_add(amuld[s[j]], a[k], E);
            //E += amuld[s[j]] * a[k].transpose();// j = p.ai[k] - номер столбца элемента e.a[k],  M2(i, j) = e.a[k]
            //E += a[s[j]] * d[j].value * a[k].transpose();// j = p.ai[k] - номер столбца элемента e.a[k],  M2(i, j) = e.a[k]
        }
        if(j <= j_min)
            break;
        k--;
    }
    return E;
}

// вычисление скалярного произведения строк i и j матрицы (p, a, d), домножается дополнительно на d[]
MATR3x3 Profile3x3StrScalarMulProfile3x3StrMulD(const size_t i, const size_t j, const ProfilePortrait &p, const std::vector<MATR3x3> &a, const std::vector<MATR3x3> &amuld)
{
    const size_t k1_min = p.ind[i];
    const size_t k1_max = p.ind[i + 1]; // строго меньше
    const size_t j1_min = i - (k1_max - k1_min);//p.getFirstCol(i);
    const size_t k2_min = p.ind[j];
    const size_t k2_max = p.ind[j + 1];
    const size_t j2_min = j - (k2_max - k2_min);
    size_t k1 = k1_min;
    size_t k2 = k2_min;
    size_t j_count; // j_count = max(j1_min, j2_min) - номер столбца, одинаковый для перемножаемых строк i и j
    if(j2_min >= j1_min)
    {
        k1 += j2_min - j1_min;
        j_count = j2_min;
    }
    else
    {
        k2 -= j2_min - j1_min;
        j_count = j1_min;
    }
    const size_t k1_max_for_cycle = k1 + (k2_max - k2);// всегда (k2_max - k2) < (k1_max - k1) т.к. начинаем идти с одного столбца
    MATR3x3 E;
    E.clear();
    while(k1 < k1_max_for_cycle)
    {
        MATR3x3::a_Mul_bT_add(amuld[k1], a[k2],
                              E);
        //E += a[k1] * d[j_count].value * a[k2].transpose();
        k1++; // L(i, j_count) = a[k1]
        k2++; // L(j, j_count) = a[k2]
        j_count++;
    }
    return E;
}

SSCM3x3Preconditioner_base *SSCM3x3Preconditioner_base::gen(const Preconditioning preconditioningType)
{
    switch (preconditioningType)
    {
    case Preconditioning::none:
    {
        return new SSCM3x3Preconditioner_none;
    }break;
    case Preconditioning::diag:
    {
        return new SSCM3x3Preconditioner_diag;
    }break;
    case Preconditioning::SSCM_LLT_0:
    {
        return new SSCM3x3Preconditioner_none;//new SSCMPreconditioner3x3_SSCM_LLT_0;
    }break;
    case Preconditioning::SSCM_LDLT_0:
    {
        return new SSCM3x3Preconditioner_none;//new SSCMPreconditioner3x3_SSCM_LDLT_0;
    }break;
    case Preconditioning::Profile_LLT:
    {
        return new SSCM3x3Preconditioner_none;//new SSCMPreconditioner3x3_Profile_LLT;
    }break;
    case Preconditioning::Profile_LDLT:
    {
        return new SSCM3x3Preconditioner_Profile_LDLT;
    }break;
    case Preconditioning::SSCM_LLT:
    {
        return new SSCM3x3Preconditioner_none;//new SSCMPreconditioner3x3_SSCM_LLT;
    }break;
    case Preconditioning::SSCM_LDLT:
    {
        return new SSCM3x3Preconditioner_SSCM_LDLT;
    }break;
    case Preconditioning::PARDISO:
    {
        return new SSCM3x3Preconditioner_PARDISO;
    }break;

    }
}
void SSCM3x3Preconditioner_base::bulid(const SSCM3x3 &matrix, const size_t firstStr, double &time)
{
    double time1;
    initPortraitAndAllocateMemory(*matrix.p, time1);
    double time2;
    updatePreconditioner(*matrix.e, firstStr, time2);
    time = time1 + time2;
}
void SSCM3x3Preconditioner_base::InversedMmulVector(Vector3 &x) const
{
    InversedSmulVector(x);
    InversedQmulVector(x);
}

void SSCM3x3Preconditioner_none::initPortraitAndAllocateMemory(const SSCMPortrait &, double &time)
{
    fprintf(stderr, "SSCMPreconditioner3x3_none initPortraitAndAllocateMemory..\n");
    time = 0;
}
void SSCM3x3Preconditioner_none::updatePreconditioner(const SSCM3x3Elements &, const size_t, double &time)
{
    time = 0;
}
void SSCM3x3Preconditioner_none::release()
{
}
void SSCM3x3Preconditioner_none::saveBMP(const char *, size_t) const
{
}
void SSCM3x3Preconditioner_none::saveProperties(const char *) const
{
}
size_t SSCM3x3Preconditioner_none::getElementsNumber() const
{
    return 0;
}
void SSCM3x3Preconditioner_none::InversedSmulVector(Vector3 &) const
{
}
void SSCM3x3Preconditioner_none::InversedQmulVector(Vector3 &) const
{
}

void SSCM3x3Preconditioner_diag::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
fprintf(stderr, "SSCMPreconditioner3x3_diag initPortraitAndAllocateMemory..\n");
    TimeIntervals::timeInterval t;
    t.begin();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    d.resize(p0.matrixSize);
    time = t.getCurrentDuration();
}
void SSCM3x3Preconditioner_diag::updatePreconditioner(const SSCM3x3Elements &set_e, const size_t, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    const SSCM3x3Elements &e0 = set_e;         // элементы исходной матрицы
    const size_t matrixSize = d.size();
    for (size_t i = 0; i < matrixSize; i++)
        e0.d[i].calcInvertSym(d[i]);//1./e0.d[i];
        //d[i] = e0.d[i].calcInvert();//1./e0.d[i];
    time = t.getCurrentDuration();
}
void SSCM3x3Preconditioner_diag::release()
{
    d.clear();
}
void SSCM3x3Preconditioner_diag::saveBMP(const char *, size_t) const
{
}
void SSCM3x3Preconditioner_diag::saveProperties(const char *fileName) const
{
    FILE *f = fopen(fileName, "w");
    fprintf(f, "Size = %dx\n", (int)d.size());
    fclose(f);
}
size_t SSCM3x3Preconditioner_diag::getElementsNumber() const
{
    return d.size();
}
void SSCM3x3Preconditioner_diag::InversedSmulVector(Vector3 &x) const
{
    // (d^-1)*x -> x
    const size_t matrixSize = d.size();
    for (size_t i = 0; i < matrixSize; i++)
        x[i] = d[i]*x[i];//x[i] *= d[i];
    // (создаётся копия т.к. в процессе умножения нельзя изменять множители)
}
void SSCM3x3Preconditioner_diag::InversedQmulVector(Vector3 &) const
{
}
/*
// не отлажено
void SSCMPreconditioner3x3_SSCM_LLT_0::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
    clock_t time1 = clock();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    m.init(p0.elementsNumber, p0.matrixSize);
    time = (clock() - time1) / (double)CLOCKS_PER_SEC;
}
void SSCMPreconditioner3x3_SSCM_LLT_0::updatePreconditioner(const SSCMElements3x3 &set_e, const size_t firstStr, double &time)
{
    clock_t time1 = clock();
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const SSCMElements3x3 &e0 = set_e;         // элементы исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    MATR3x3 E;
    for (size_t j = 0; j < matrixSize; j++)               // j - номер столбца
    {                                       // i - номер строки
        size_t i;
        for (size_t k = p0.ind[j]; k < p0.ind[j + 1] && (i = p0.ai[k]) < j; k++)
        {
            size_t k1 = p0.ind[j];
            size_t k2 = p0.ind[i];
            E.clear();
            size_t i1;
            size_t i2;
            for (; k1 < p0.ind[j + 1] && k2 < p0.ind[i + 1] && (i1 = p0.ai[k1]) < i && (i2 = p0.ai[k2]) < i;)
            {
                if (i1 == i2)               // L(j, i1) = L[k1]
                    E -= m.a[k1++] * m.a[k2++]; // L(i, i2) = L[k2]
                else
                    if (i1 < i2)
                        k1++;
                    else
                        k2++;
            }
            m.a[k] = (E + e0.a[k]) / m.d[i];       // L(j, i) = L[k]
        }
        E.clear();
        for (size_t k = p0.ind[j]; k < p0.ind[j + 1] && p0.ai[k] < j; k++)
            E -= m.a[k] * m.a[k];           // i = ai[k] - номер строки
        m.d[j] = (E + e0.d[j]).sqrt();          // A(j, j) = d[j]
        //m.d[j] = sqrt(E + e0.d[j]);          // A(j, j) = d[j]
                                            // L(j, j) = D[j]
    }
    time = (clock() - time1) / (double)CLOCKS_PER_SEC;
}
void SSCMPreconditioner3x3_SSCM_LLT_0::release()
{
    m.release();
}
void SSCMPreconditioner3x3_SSCM_LLT_0::saveBMP(const char *fileName, size_t size) const
{
    p_ptr0_SSCM->saveBMP(fileName, size);
}
void SSCMPreconditioner3x3_SSCM_LLT_0::saveProperties(const char *fileName) const
{
    p_ptr0_SSCM->saveProperties(fileName);
}
size_t SSCMPreconditioner3x3_SSCM_LLT_0::getElementsNumber() const
{
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    return p0.getElementsNumber();
}
void SSCMPreconditioner3x3_SSCM_LLT_0::InversedSmulVector(Vector3 &x) const
{
    // (L^-1)*x -> x
    const SSCMPortrait &p = *p_ptr0_SSCM;
    size_t matrixSize = p.matrixSize;
    VECTOR3 E;
    for (size_t j = 0; j < matrixSize; j++) // j - номер столбца
    {
        E = x[j];
        size_t i; // i - номер строки
        for (size_t k = p.ind[j]; k < p.ind[j + 1] && (i = p.ai[k]) < j; k++)
            E -= m.a[k] * x[i];               // L(j, i) = L[k]
        x[j] = E / m.d[j];
    }
}
void SSCMPreconditioner3x3_SSCM_LLT_0::InversedQmulVector(Vector3 &x) const
{
    // (U^-1)*x -> x
    // U*x_new = x
    const SSCMPortrait &p = *p_ptr0_SSCM;
    size_t j = p.matrixSize - 1;
    for (;;)            // j - номер столбца
    {
        x[j] = x[j] / m.d[j];
        //x[j] /= m.d[j];                       // i = ai[k] - номер строки
        for(size_t k = p.ind[j]; k < p.ind[j + 1]; k++)
            x[p.ai[k]] -=  m.a[k]*x[j];        // U(i, j) = L[k]
        if(j == 0)
            break;
        j--;
    }
}

// не отлажено
void SSCMPreconditioner3x3_SSCM_LDLT_0::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
    clock_t time1 = clock();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    m.init(p0.elementsNumber, p0.matrixSize);
    time = (clock() - time1) / (double)CLOCKS_PER_SEC;
}
void SSCMPreconditioner3x3_SSCM_LDLT_0::updatePreconditioner(const SSCMElements3x3 &set_e, const size_t firstStr, double &time)
{
    clock_t time1 = clock();
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const SSCMElements3x3 &e0 = set_e;         // элементы исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    MATR3x3 E;
    for (size_t j = 0; j < matrixSize; j++)               // j - номер столбца
    {                                       // i - номер строки
        size_t i;
        for (size_t k = p0.ind[j]; k < p0.ind[j + 1] && (i = p0.ai[k]) < j; k++)
        {
            size_t k1 = p0.ind[j];
            size_t k2 = p0.ind[i];
            E.clear();
            size_t i1;
            size_t i2;
            for (; k1 < p0.ind[j + 1] && k2 < p0.ind[i + 1] && (i1 = p0.ai[k1]) < i && (i2 = p0.ai[k2]) < i;)
            {                               // L(j, i1) = L[k1]
                if (i1 == i2)               // L(i, i2) = L[k2]
                    E -= m.a[k1++] * m.a[k2++] * m.d[i1];
                else
                    if (i1 < i2)
                        k1++;
                    else
                        k2++;
            }
            m.a[k] = (E + e0.a[k]) / m.d[i];       // L(j, i) = L[k]
        }
        E.clear();
        for (size_t k = p0.ind[j]; k < p0.ind[j + 1] && (i = p0.ai[k]) < j; k++)
            E -= m.a[k] * m.a[k] * m.d[i];  // L(j, i) = L[k]
                                            // L(i, i) = D[i]
        m.d[j] = E + e0.d[j];                    // A(i, i) = d[i]
    }
    time = (clock() - time1) / (double)CLOCKS_PER_SEC;
}
void SSCMPreconditioner3x3_SSCM_LDLT_0::release()
{
    m.release();
}
void SSCMPreconditioner3x3_SSCM_LDLT_0::saveBMP(const char *fileName, size_t size) const
{
    p_ptr0_SSCM->saveBMP(fileName, size);
}
void SSCMPreconditioner3x3_SSCM_LDLT_0::saveProperties(const char *fileName) const
{
    p_ptr0_SSCM->saveProperties(fileName);
}
size_t SSCMPreconditioner3x3_SSCM_LDLT_0::getElementsNumber() const
{
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    return p0.getElementsNumber();
}
void SSCMPreconditioner3x3_SSCM_LDLT_0::InversedSmulVector(Vector3 &x) const
{
    // (D^-1*L^-1)*x -> x
    // LD*x_new = x
    // LD(i, j) = L(i, j)*D(j, j)
    // LD(i, i) = D(i, i)
    const SSCMPortrait &p = *p_ptr0_SSCM;
    size_t matrixSize = p.matrixSize;
    VECTOR3 E;
    for (size_t j = 0; j < matrixSize; j++)                 // j - номер столбца
    {
        E = x[j];                           // i - номер строки
        size_t i;
        for (size_t k = p.ind[j]; k < p.ind[j + 1] && (i = p.ai[k]) < j; k++)
            E -= m.a[k] * m.d[i] * x[i];        // L(j, i) = L[k]
        x[j] = E / m.d[j];
    }
}
void SSCMPreconditioner3x3_SSCM_LDLT_0::InversedQmulVector(Vector3 &x) const
{
    // (U^-1)*x -> x
    // U*x_new = x
    const SSCMPortrait &p = *p_ptr0_SSCM;
    size_t j = p.matrixSize - 1;
    for(;;)                    // j - номер столбца
    {
        for(size_t k = p.ind[j]; k < p.ind[j + 1]; k++)  // i = ai[k] - номер строки
            x[p.ai[k]] -= m.a[k] * x[j];                // U(i, j) = L[k]
        if(j == 0)
            break;
        j--;
    }

    //for (j = (int)p.matrixSize - 1; j >= 0; j--)                    // j - номер столбца
    //    for (k = p.ind[j + 1] - 1; k >= p.ind[j]; k--)  // i = ai[k] - номер строки
    //        x[p.ai[k]] -= x[j] * m.a[k];                // U(i, j) = L[k]
}

// не отлажено
void SSCMPreconditioner3x3_Profile_LLT::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
    clock_t time1 = clock();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    ProfilePortrait &p = p_profile;
    // инициализация портрета и матрицы
    size_t elementsNumber = 0; // общее количество элементов матрицы с плотным профилем
    for (size_t i = 0; i < matrixSize; i++) // i - номер строки
    {
        if(p0.ind[i + 1] - p0.ind[i] != 0)  // строка не пустая
        {
            size_t j0 = p0.ai[p0.ind[i]];   // j0 - номер столбца первого элемента строки i
            size_t stringSize = i - j0;     // в строке i плотный профиль занимает (i - j0) элементов
            elementsNumber += stringSize;
        }
    }
    p.init(matrixSize);
    m.init(elementsNumber, matrixSize);
    // заполнение профильного портрета
    p.ind[0] = 0;
    for (size_t i = 0; i < matrixSize; i++) // i - номер строки
    {
        if(p0.ind[i + 1] - p0.ind[i] != 0)  // строка не пустая
        {
            size_t j0 = p0.ai[p0.ind[i]];   // j0 - номер столбца первого элемента строки i
            size_t stringSize = i - j0;     // в строке i плотный профиль занимает (i - j0) элементов
            p.ind[i + 1] = p.ind[i] + stringSize;
        }
        else
            p.ind[i + 1] = p.ind[i];    // строка пустая
    }
    time = (clock() - time1) / (double)CLOCKS_PER_SEC;
}
void SSCMPreconditioner3x3_Profile_LLT::updatePreconditioner(const SSCMElements3x3 &set_e, const size_t firstStr, double &time)
{
    clock_t time1 = clock();
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const SSCMElements3x3 &e0 = set_e;         // элементы исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    ProfilePortrait &p = p_profile;
    // копирование исходной матрицы в профильный формат
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        if(p0.ind[i + 1] - p0.ind[i] != 0)  // строка не пустая
        {
            size_t k_min = p.ind[i];
            size_t k_max = p.ind[i + 1];
            size_t k0_min = p0.ind[i];
            size_t k0_max = p0.ind[i + 1];
            size_t j_min = p.getFirstCol(i);//p0.ai[p0.ind[i]];   // j0 - номер столбца первого элемента строки i
            //size_t stringSize = i - j_min;     // в строке i плотный профиль занимает (i - j0) элементов
            for(size_t k = k_min; k < k_max; k++)
                m.a[k].clear();     // m.a[k] = L(i, i-stringSize + k-p.ind[i])
            for(size_t k0 = k0_min; k0 < k0_max; k0++)
            {
                size_t j = p0.ai[k0];           // e.a[k] = A(i, j)
                size_t k = k_min + (j - j_min); // m.a[kk] = L(i, j)
                m.a[k] = e0.a[k0];
            }
        }
        m.d[i] = e0.d[i];
    }
fprintf(stderr, "1) copy time %le\n", (clock() - time1) / (double)CLOCKS_PER_SEC);
int lastPerc = 0;
double lastTime = 0;
    MATR3x3 *__restrict a = m.a.data();
    MATR3x3 *__restrict d = m.d.data();
    size_t *__restrict ind = p.ind.data();
    // построение матрицы предобусловливателя
    for (size_t i = firstStr; i < matrixSize; i++)           // i - номер строки
    {                                                // j - номер стоблца
        const size_t k1_min = ind[i];
        const size_t k1_max = ind[i + 1];  // строго меньше
        const size_t j1_min = i - (ind[i + 1] - ind[i]);//p.getFirstCol(i);
        MATR3x3 E_sqr_a;
        E_sqr_a.clear();
        size_t k = k1_min;
        size_t j = j1_min;
        for (;k < k1_max;)
        {
            const size_t k2_min = ind[j];
            const size_t k2_max = ind[j + 1];
            size_t k1 = k1_min;
            size_t k2 = k2_min;
            size_t firstCol_j = j - (k2_max - k2_min);
            if(firstCol_j >= j1_min)
            {
                k1 += firstCol_j - j1_min;
            }
            else
            {
                k2 -= firstCol_j - j1_min;
            }
            MATR3x3 E;
            E.clear();
            size_t k1_max_for_cycle = k1 + (k2_max - k2);// всегда (k2_max - k2) < (k1_max - k1) т.к. начинаем идти с одного столбца
            while(k1 < k1_max_for_cycle)
            {
                E += a[k1] * a[k2];
                k1++;                   // L(i, j1) = m.a[k1]
                k2++;                   // L(j, j2) = m.a[k2]
            }
            a[k] = (a[k] - E) / d[j];       // L(i, j) = m.a[k]
            E_sqr_a += a[k] * a[k];
            k++;
            j++;
        }
        d[i] = (d[i] - E_sqr_a).sqrt();
        //d[i] = sqrt(d[i] - E_sqr_a);
// вывод сообщения о проценте проделанной работы не чаще 1 раза в 1024 циклов и не чаще 1 раза в 5 секунд
if(i%1024 == 1023)
{
double newTime = (clock() - time1) / (double)CLOCKS_PER_SEC;
if(newTime - lastTime > 5)
{
lastTime = newTime;
int newPerc = round((double)k/p.ind[matrixSize]*100);
if(newPerc > lastPerc)
{
lastPerc = newPerc;
fprintf(stderr, "     %d\n", lastPerc);
}
}
}
    }
fprintf(stderr, "2) prec time = %le\n", (clock() - time1) / (double)CLOCKS_PER_SEC);
    time = (clock() - time1) / (double)CLOCKS_PER_SEC;
}
void SSCMPreconditioner3x3_Profile_LLT::release()
{
    p_profile.release();
    m.release();
}
void SSCMPreconditioner3x3_Profile_LLT::saveBMP(const char *fileName, size_t size) const
{
    p_profile.saveBMP(fileName, size);
}
void SSCMPreconditioner3x3_Profile_LLT::saveProperties(const char *fileName) const
{
    p_profile.saveProperties(fileName);
}
size_t SSCMPreconditioner3x3_Profile_LLT::getElementsNumber() const
{
    return p_profile.getElementsNumber();
}
void SSCMPreconditioner3x3_Profile_LLT::InversedSmulVector(Vector3 &x) const
{
    const ProfilePortrait &p = p_profile;
    size_t matrixSize = p.matrixSize;
    // (L^-1)*x -> x
    for (size_t i = 0; i < matrixSize; i++)    // i - номер строки
    {
        VECTOR3 E;
        E.clear();
        size_t k = p.ind[i];
        size_t j = p.getFirstCol(i);
        size_t k_max = p.ind[i + 1];
        while(k < k_max)
        {
            E += m.a[k] * x[j];               // L(i, j) = m.a[k]
            k++;
            j++;
        }
        x[i] = (x[i] - E) / m.d[i];
    }
}
void SSCMPreconditioner3x3_Profile_LLT::InversedQmulVector(Vector3 &x) const
{
    const ProfilePortrait &p = p_profile;
    // (U^-1)*x -> x
    // U*x_new = x
    size_t j = p.matrixSize - 1;
    for (;;)  // j - номер столбца
    {
        x[j] = x[j] / m.d[j];
       // x[j] /= m.d[j];
        const size_t k_min = p.ind[j];
        const size_t k_max = p.ind[j + 1];  // строго меньше
        size_t i = j - (k_max - k_min);
        size_t k = k_min;
        while(k < k_max)
        {
            x[i] -= m.a[k] * x[j];         // U(i, j) = L[k]
            i++;    // i < j
            k++;
        }
        if(j == 0)
            break;
        j--;
    }
}

// не отлажено
void SSCMPreconditioner3x3_SSCM_LLT::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
    clock_t time1 = clock();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    // p0, e0 - портрет и элементы исходной матрицы в строчно-столбцовом формате
    // требуется построить портрет p и элементы m матрицы L
    SSCMPortrait &p = p_SSCM;
    SSCMElements3x3 &m = m_SSCM;
    p.init_LLT_SSCM(p0);
    m.init(p.elementsNumber, matrixSize);
    s.resize(matrixSize);
    for (size_t i = 0; i < matrixSize; i++)
        s[i].clear();
fprintf(stderr, "Portrait time = %le\n", (clock() - time1) / (double)CLOCKS_PER_SEC);
    time = (clock() - time1) / (double)CLOCKS_PER_SEC;
}
void SSCMPreconditioner3x3_SSCM_LLT::updatePreconditioner(const SSCMElements3x3 &set_e, const size_t firstStr, double &time)
{
    clock_t time1 = clock();
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const SSCMElements3x3 &e0 = set_e;         // элементы исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    SSCMPortrait &p = p_SSCM;
    SSCMElements3x3 &m = m_SSCM;
    // 2) копирование исходной матрицы в L
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        SSCM3x3StrCopyToSSCM3x3Str(i, p, m, p0, e0);    // p0,e -> p,m
    }
fprintf(stderr, "1) copy time = %le\n", (clock() - time1) / (double)CLOCKS_PER_SEC);
int lastPerc = 0;
double lastTime = 0;
    time1 = clock();
    // 3) построение разложения - заполнение элементов m
    for (size_t i = firstStr; i < matrixSize; i++)           // i - номер строки
    {
        const size_t k_min = p.ind[i];
        const size_t k_max = p.ind[i + 1];  // строго меньше
        MATR3x3 E_sqr_a;
        E_sqr_a.clear();
        size_t k = k_min;
        while(k < k_max)
        {
            size_t j = p.ai[k]; // номер столбца в строке i
            // вычисление скалярного произведения
            size_t s_j_min = p.ai[k_min];
            MATR3x3 E = SSCM3x3StrScalarMulVector3(j, p, m, s, s_j_min);
            MATR3x3 newValue = (m.a[k] - E) / m.d[j];
            m.a[k] = newValue; // L(i, j) = m.a[k]
            E_sqr_a += newValue * newValue;
            s[j] = newValue;
            k++;

            //double E = SSCMStrScalarMulSSCMStr(i, j, p, m);
            //m.a[k] = (m.a[k] - E) / m.d[j]; // L(i, j) = m.a[k]
            //E_sqr_a += m.a[k] * m.a[k];
            //k++;
        }
        // очистка всех ненулевых элементов временного вектора
        k = k_min;
        while(k < k_max)
        {
            s[p.ai[k]].clear();
            k++;
        }
        m.d[i] = (m.d[i] - E_sqr_a).sqrt();
        //m.d[i] = sqrt(m.d[i] - E_sqr_a);
// вывод сообщения о проценте проделанной работы не чаще 1 раза в 1024 циклов и не чаще 1 раза в 5 секунд
if(i%1024 == 1023)
{
double newTime = (clock() - time1) / (double)CLOCKS_PER_SEC;
if(newTime - lastTime > 5)
{
lastTime = newTime;
int newPerc = round((double)k/p.ind[matrixSize]*100);
if(newPerc > lastPerc)
{
lastPerc = newPerc;
fprintf(stderr, "     %d\n", lastPerc);
}
}
}
    }
fprintf(stderr, "2) prec time = %le\n", (clock() - time1) / (double)CLOCKS_PER_SEC);
    time = (clock() - time1) / (double)CLOCKS_PER_SEC;
}
void SSCMPreconditioner3x3_SSCM_LLT::release()
{
    p_SSCM.release();
    m_SSCM.release();
    s.clear();
}
void SSCMPreconditioner3x3_SSCM_LLT::saveBMP(const char *fileName, size_t size) const
{
    p_SSCM.saveBMP(fileName, size);
}
void SSCMPreconditioner3x3_SSCM_LLT::saveProperties(const char *fileName) const
{
    p_SSCM.saveProperties(fileName);
}
size_t SSCMPreconditioner3x3_SSCM_LLT::getElementsNumber() const
{
    return p_SSCM.getElementsNumber();
}
void SSCMPreconditioner3x3_SSCM_LLT::InversedSmulVector(Vector3 &x) const
{
    // (L^-1)*x -> x
    const SSCMPortrait &p = p_SSCM;
    const SSCMElements3x3 &m = m_SSCM;
    const size_t matrixSize = p.matrixSize;
    for (size_t i = 0; i < matrixSize; i++) // i - номер строки
    {
        const size_t k_min = p.ind[i];
        const size_t k_max = p.ind[i + 1];  // строго меньше
        VECTOR3 E;
        E.clear();
        size_t k = k_min;
        while(k < k_max)
        {
            size_t j = p.ai[k]; // номер столбца в строке i
            E += m.a[k] * x[j]; // L(i, j) = m.a[k]
            k++;
        }
        x[i] = (x[i] - E) / m.d[i];
    }
}
void SSCMPreconditioner3x3_SSCM_LLT::InversedQmulVector(Vector3 &x) const
{
    // (U^-1)*x -> x
    const SSCMPortrait &p = p_SSCM;
    const SSCMElements3x3 &m = m_SSCM;
    size_t j = p.matrixSize - 1;
    for (;;)  // j - номер столбца
    {
        x[j] = x[j] / m.d[j];
        //x[j] /= m.d[j];
        const size_t k_min = p.ind[j];
        const size_t k_max = p.ind[j + 1];  // строго меньше
        size_t k = k_min;
        while(k < k_max)
        {
            size_t i = p.ai[k]; // номер строки
            x[i] -= m.a[k] * x[j]; // U(i, j) = L[k]
            k++;
        }
        if(j == 0)
            break;
        j--;
    }
}
*/

void SSCM3x3Preconditioner_Profile_LDLT::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
fprintf(stderr, "SSCMPreconditioner3x3_Profile_LDLT initPortraitAndAllocateMemory..\n");
    TimeIntervals::timeInterval t;
    t.begin();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    // инициализация портрета
    p.init_LLT(p0);
fprintf(stderr, "Portrait time = %le\n", t.getCurrentDuration());
    // матрицы
    a.resize(p.getElementsNumber());
    amuld.resize(p.getElementsNumber());
    d.resize(p.getMatrixSize());
    time = t.getCurrentDuration();
}
void SSCM3x3Preconditioner_Profile_LDLT::updatePreconditioner(const SSCM3x3Elements &set_e0, const size_t firstStr, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const SSCM3x3Elements &e0 = set_e0;     // элементы исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    // копирование исходной матрицы в профильный формат
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        SSCM3x3StrCopyToProfileStr(i, p, a, d, p0, e0.a, e0.d); // копирование строки
    }
fprintf(stderr, "1) copy time %le\n", t.getCurrentDuration());
TimeIntervals::timeInterval debug_t;
debug_t.begin();
int lastPerc = 0;
double lastTime = 0;
    // построение матрицы предобусловливателя
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        const size_t k1_min = p.ind[i];
        const size_t k1_max = p.ind[i + 1]; // строго меньше
        const size_t j1_min = i - (k1_max - k1_min);//p.getFirstCol(i);
        MATR3x3 E_lld;
        E_lld.clear();
        size_t k = k1_min;
        size_t j = j1_min; // j - номер стоблца
        for (;k < k1_max;)
        {
            MATR3x3 E = Profile3x3StrScalarMulProfile3x3StrMulD(i, j, p, a, amuld);
            MATR3x3 newValue;
            MATR3x3::a_mul_b(a[k] - E, d[j].invertValue,
                             newValue);
            a[k] = newValue; // L(i, j) = m.a[k]
            //a[k] = (a[k] - E) * d[j].invertValue; // L(i, j) = a[k]
            MATR3x3::a_mul_b(newValue, d[j].value,
                             amuld[k]);
            //amuld[k] = newValue * d[j].value;
            MATR3x3::a_Mul_bT_add(amuld[k], newValue,
                                  E_lld);
            //E_lld += amuld[k] * a[k];
            //E_lld += a[k] * d[j].value * a[k];
            k++;
            j++;
        }
        d[i].value -= E_lld;
        d[i].value.fixAsymmetry();
        d[i].value.calcInvertSym(d[i].invertValue);// вычисление обратной матрицы
        //d[i].invertValue = d[i].value.calcInvert();// вычисление обратной матрицы
// вывод сообщения о проценте проделанной работы не чаще 1 раза в 1024 циклов и не чаще 1 раза в 5 секунд
if(i%1024 == 1023)
{
double newTime = debug_t.getCurrentDuration();
if(newTime - lastTime > 5)
{
lastTime = newTime;
int newPerc = round((double)k/p.ind[matrixSize]*100);
if(newPerc > lastPerc)
{
lastPerc = newPerc;
fprintf(stderr, "     %d\n", lastPerc);
}
}
}
    }
time = t.getCurrentDuration();
fprintf(stderr, "2) prec time = %le\n", debug_t.getCurrentDuration());
}
void SSCM3x3Preconditioner_Profile_LDLT::release()
{
    p.release();
    a.clear();
    amuld.clear();
    d.clear();
}
void SSCM3x3Preconditioner_Profile_LDLT::saveBMP(const char *fileName, size_t size) const
{
    p.saveBMP(fileName, size);
}
void SSCM3x3Preconditioner_Profile_LDLT::saveProperties(const char *fileName) const
{
    p.saveProperties(fileName);
}
size_t SSCM3x3Preconditioner_Profile_LDLT::getElementsNumber() const
{
    return p.getElementsNumber();
}
void SSCM3x3Preconditioner_Profile_LDLT::InversedSmulVector(Vector3 &x) const
{
    size_t matrixSize = p.matrixSize;
    // (D^-1*L^-1)*x -> x
    // LD*x_new = x
    // LD(i, j) = L(i, j)*D(j, j)
    // LD(i, i) = D(i, i)
    for (size_t i = 0; i < matrixSize; i++)
    {
        VECTOR3 E;
        E.clear();
        size_t k = p.ind[i];
        size_t j = p.getFirstCol(i);
        size_t k_max = p.ind[i + 1];
        while(k < k_max)
        {
            MATR3x3::a_mul_b_add(amuld[k], x[j], E);
            //E += amuld[k] * x[j]; // L(i, j) = m.a[k]
            k++;
            j++;
        }
        MATR3x3::a_mul_b(d[i].invertValue, x[i] - E,
                         x[i]);
        //x[i] = d[i].invertValue * (x[i] - E);
    }
}
void SSCM3x3Preconditioner_Profile_LDLT::InversedQmulVector(Vector3 &x) const
{
    // (U^-1)*x -> x
    // U*x_new = x
    size_t j = p.matrixSize - 1;
    for (;;) // j - номер столбца
    {
        const size_t k_min = p.ind[j];
        const size_t k_max = p.ind[j + 1]; // строго меньше
        size_t i = j - (k_max - k_min);
        size_t k = k_min;
        while(k < k_max)
        {
            MATR3x3::aT_mul_b_sub(a[k], x[j],
                                  x[i]);
            //x[i] -= a[k].transpose() * x[j]; // U(i, j) = L[k]
            i++;    // i < j
            k++;
        }
        if(j == 0)
            break;
        j--;
    }
}

void SSCM3x3Preconditioner_SSCM_LDLT::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
fprintf(stderr, "SSCMPreconditioner3x3_SSCM_LDLT initPortraitAndAllocateMemory..\n");
TimeIntervals::timeInterval t;
t.begin();
    p_ptr0_SSCM = &set_p;
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    // p0, e0 - портрет и элементы исходной матрицы в строчно-столбцовом формате
    // требуется построить портрет p и элементы e матрицы L
    p.init_LLT(p0);
fprintf(stderr, "Portrait time = %le\n", t.getCurrentDuration());
    a.resize(p.elementsNumber);
    amuld.resize(p.elementsNumber);
    d.resize(matrixSize);
    s.resize(matrixSize);
    for (size_t i = 0; i < matrixSize; i++)
        s[i] = -1;
    time = t.getCurrentDuration();
}
void SSCM3x3Preconditioner_SSCM_LDLT::updatePreconditioner(const SSCM3x3Elements &set_e0, const size_t firstStr, double &time)
{
    if(firstStr >= p_ptr0_SSCM->matrixSize)
    {
        time = 0;
        return;
    }
    fprintf(stderr, "...firstStr/matrixSize = %ld / %ld\n", firstStr, p_ptr0_SSCM->matrixSize);
TimeIntervals::timeInterval t;
t.begin();
    const SSCMPortrait &p0 = *p_ptr0_SSCM;  // портрет исходной матрицы
    const SSCM3x3Elements &e0 = set_e0;         // элементы исходной матрицы
    const size_t matrixSize = p0.matrixSize;
    // 2) копирование исходной матрицы в L
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        SSCM3x3StrCopyToSSCM3x3Str(i, p, a, d, p0, e0.a, e0.d);    // p0,e -> p,m
    }
fprintf(stderr, "1) copy time %le\n", t.getCurrentDuration());
TimeIntervals::timeInterval debug_t;
debug_t.begin();
int lastPerc = 0;
double lastTime = 0;
    // 3) построение разложения - заполнение элементов m
    for (size_t i = firstStr; i < matrixSize; i++) // i - номер строки
    {
        const size_t k_min = p.ind[i];
        const size_t k_max = p.ind[i + 1];  // строго меньше
        MATR3x3 E_lld;
        E_lld.clear();
        size_t k = k_min;
        while(k < k_max)
        {
            size_t j = p.ai[k]; // номер столбца в строке i
            // вычисление скалярного произведения
            size_t s_j_min = p.ai[k_min];
            MATR3x3 E = SSCM3x3StrScalarMulVecto3MulD(j, p, a, amuld, s, s_j_min);
            MATR3x3 newValue;
            MATR3x3::a_mul_b(a[k] - E, d[j].invertValue,
                             newValue);
            //MATR3x3 newValue = (a[k] - E) * d[j].invertValue;
            //MATR3x3 newValue = (a[k] - E) * d[j].value.calcInvert();
            a[k] = newValue; // L(i, j) = m.a[k]
            MATR3x3::a_mul_b(newValue, d[j].value,
                             amuld[k]);
            //amuld[k] = newValue * d[j].value;
            MATR3x3::a_Mul_bT_add(amuld[k], newValue,
                                  E_lld);
            //E_lld += amuld[k] * newValue.transpose();
            //E_lld += newValue * d[j].value * newValue.transpose();
            s[j] = k;
            k++;
        }
        // очистка всех ненулевых элементов временного вектора
        k = k_min;
        while(k < k_max)
        {
            s[p.ai[k]] = -1;
            k++;
        }
        d[i].value -= E_lld;
        d[i].value.fixAsymmetry();
        d[i].value.calcInvertSym(d[i].invertValue);// вычисление обратной матрицы
        //d[i].invertValue = d[i].value.calcInvert();// вычисление обратной матрицы
// вывод сообщения о проценте проделанной работы не чаще 1 раза в 1024 циклов и не чаще 1 раза в 5 секунд
if(i%1024 == 1023)
{
double newTime = debug_t.getCurrentDuration();
if(newTime - lastTime > 5)
{
lastTime = newTime;
int newPerc = round((double)k/p.ind[matrixSize]*100);
if(newPerc > lastPerc)
{
lastPerc = newPerc;
fprintf(stderr, "     %d\n", lastPerc);
}
}
}
    }
    time = t.getCurrentDuration();
fprintf(stderr, "2) prec time = %le\n", debug_t.getCurrentDuration());
}
void SSCM3x3Preconditioner_SSCM_LDLT::release()
{
    p.release();
    a.clear();
    amuld.clear();
    d.clear();
    s.clear();
}
void SSCM3x3Preconditioner_SSCM_LDLT::saveBMP(const char *fileName, size_t size) const
{
    p.saveBMP(fileName, size);
}
void SSCM3x3Preconditioner_SSCM_LDLT::saveProperties(const char *fileName) const
{
    p.saveProperties(fileName);
}
size_t SSCM3x3Preconditioner_SSCM_LDLT::getElementsNumber() const
{
    return p.getElementsNumber();
}
void SSCM3x3Preconditioner_SSCM_LDLT::InversedSmulVector(Vector3 &x) const
{
    // (D^-1*L^-1)*x -> x
    // LD*x_new = x
    // LD(i, j) = L(i, j)*D(j, j)
    // LD(i, i) = D(i, i)
    size_t matrixSize = p.matrixSize;
    for (size_t j = 0; j < matrixSize; j++) // j - номер столбца
    {
        VECTOR3 E;
        E.clear();
        for (size_t k = p.ind[j]; k < p.ind[j + 1]; k++)
        {
            size_t i = p.ai[k]; // i - номер строки
            MATR3x3::a_mul_b_add(amuld[k], x[i], E);
            //E += amuld[k] * x[i]; // L(j, i) = L[k]
            //E += a[k] * d[i].value * x[i]; // L(j, i) = L[k]
        }
        MATR3x3::a_mul_b(d[j].invertValue, x[j] - E,
                         x[j]);
        //x[j] = d[j].invertValue * (x[j] - E);
        //x[j] = (x[j] - E) / d[j].value;
    }
}
void SSCM3x3Preconditioner_SSCM_LDLT::InversedQmulVector(Vector3 &x) const
{
    // (U^-1)*x -> x
    // U*x_new = x
    size_t j = p.matrixSize - 1;
    for(;;)                    // j - номер столбца
    {
        for(size_t k = p.ind[j]; k < p.ind[j + 1]; k++)  // i = ai[k] - номер строки
        {
            MATR3x3::aT_mul_b_sub(a[k], x[j],
                                  x[p.ai[k]]);
            //x[p.ai[k]] -= a[k].transpose() * x[j];     // U(i, j) = L[k]
        }
        if(j == 0)
            break;
        j--;
    }
}


struct convert_to_pardiso_el
{
    int j;   // номер столбца
    int k;   // номер элемента
    friend inline bool operator<(const convert_to_pardiso_el &l, const convert_to_pardiso_el &r)
    {
        return l.j < r.j;
    }
};

void SSCM3x3Preconditioner_PARDISO::initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time)
{
    mkl_set_num_threads(4);
fprintf(stderr, "SSCMPreconditioner3x3_SSCM_LDLT initPortraitAndAllocateMemory..\n");
TimeIntervals::timeInterval t;
t.begin();
    matrixSize = set_p.getMatrixSize();
    std::vector<std::set<convert_to_pardiso_el>> m; // m[i] - упорядоченное множество номеров столбцов с ненулевыми значениями элементов для строки i
    m.resize(matrixSize);
    // вне диагонали
    for(size_t i = 0; i < matrixSize; i++)
    {
        size_t k_min = set_p.ind[i];
        size_t k_max = set_p.ind[i + 1]; // строго меньше
        for(size_t k = k_min; k < k_max; k++)
        {
            size_t j = set_p.ai[k];
            convert_to_pardiso_el el;
            el.j = i;
            el.k = k;
            m[j].insert(el);
        }
    }
    // диагональ
    for(size_t i = 0; i < matrixSize; i++)
    {
        convert_to_pardiso_el el;
        el.j = i;
        el.k = -1 - i;  // диагональ
        m[i].insert(el);
    }
    elementsNumber = 0;
    for(size_t i = 0; i < matrixSize; i++)
        elementsNumber += m[i].size();
    ia = new MKL_INT[matrixSize + 1];
    ja = new MKL_INT[elementsNumber];
    a = new double[elementsNumber*9];
    b = new double[matrixSize*3];
    x = new double[matrixSize*3];
    k_convert = new int[elementsNumber];
    size_t count_a = 0;
    for(size_t i = 0; i < matrixSize; i++)
    {
        ia[i] = count_a;
        for(std::set<convert_to_pardiso_el>::iterator it = m[i].begin(); it != m[i].end(); ++it)
        {
            ja[count_a] = (*it).j;
            k_convert[count_a] = (*it).k;
            count_a++;
        }
    }
    ia[matrixSize] = count_a;
    m.clear();
    // --------------------------------------------------------------------
    // .. Setup Pardiso control parameters
    // --------------------------------------------------------------------
    n = (MKL_INT)matrixSize;
    MKL_INT i;
    for ( i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;         // No solver default
    iparm[1] = 2;         // Fill-in reordering from METIS
    iparm[3] = 0;         // No iterative-direct algorithm
    iparm[4] = 0;         // No user fill-in reducing permutation
    iparm[5] = 0;         // Write solution into x
    iparm[6] = 0;         // Not in use
    iparm[7] = 2;         // Max numbers of iterative refinement steps
    iparm[8] = 0;         // Not in use
    iparm[9] = 8;        // Perturb the pivot elements with 1E-13
    iparm[10] = 0;        // Use nonsymmetric permutation and scaling MPS
    iparm[11] = 0;        // Not in use
    iparm[12] = 0;        // Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy
    iparm[13] = 0;        // Output: Number of perturbed pivots
    iparm[14] = 0;        // Not in use
    iparm[15] = 0;        // Not in use
    iparm[16] = 0;        // Not in use
    iparm[17] = -1;       // Output: Number of nonzeros in the factor LU
    iparm[18] = -1;       // Output: Mflops for LU factorization
    iparm[19] = 0;        // Output: Numbers of CG Iterations
    iparm[34] = 1;    // нумерация с 0
    iparm[36] = block_size;        // Size of block
    maxfct = 1;           // Maximum number of numerical factorizations
    mnum = 1;         // Which factorization to use
    msglvl = 0;           // Print statistical information in file
    error = 0;            // Initialize error flag
    // --------------------------------------------------------------------
    // .. Initialize the internal solver memory pointer. This is only
    // necessary for the FIRST call of the PARDISO solver
    // --------------------------------------------------------------------
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
    // --------------------------------------------------------------------
    // .. Reordering and Symbolic Factorization. This step also allocates
    // all memory that is necessary for the factorization.
    // --------------------------------------------------------------------
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }
    printf ("\nReordering completed ... ");
    printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
    printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
fprintf(stderr, "Portrait time = %le\n", t.getCurrentDuration());
    time = t.getCurrentDuration();
}
void SSCM3x3Preconditioner_PARDISO::updatePreconditioner(const SSCM3x3Elements &set_e0, const size_t firstStr, double &time)
{
    if(firstStr >= matrixSize)
    {
        time = 0;
        return;
    }
    fprintf(stderr, "...firstStr/matrixSize = %ld / %ld\n", firstStr, matrixSize);
TimeIntervals::timeInterval t;
t.begin();
    // копирование матрицы в формат PARDISO
    for(size_t k = 0; k < elementsNumber; k++)
    {
        if(k_convert[k] >= 0)
        {
            int k0 = k_convert[k];
            for(size_t ii = 0; ii < 3; ii++)
                for(size_t jj = 0; jj < 3; jj++)
                    a[k*9 + ii*3 + jj] = set_e0.a[k0].m[jj][ii];//транспонировано
        }
        else
        {
            int k0 = -1 - k_convert[k];
            for(size_t ii = 0; ii < 3; ii++)
                for(size_t jj = 0; jj < 3; jj++)
                    a[k*9 + ii*3 + jj] = set_e0.d[k0].m[ii][jj];//не транспонировано
        }
    }
fprintf(stderr, "1) copy time %le\n", t.getCurrentDuration());
TimeIntervals::timeInterval debug_t;
debug_t.begin();
// --------------------------------------------------------------------
// .. Numerical factorization
// --------------------------------------------------------------------
    const MKL_INT n = (MKL_INT)matrixSize;
    msglvl = 0;           // Print statistical information in file
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        fprintf (stderr, "\nERROR during numerical factorization: %d", error);
        exit (2);
    }
    //printf ("\nFactorization completed ... ");
    time = t.getCurrentDuration();
fprintf(stderr, "2) prec time = %le\n", debug_t.getCurrentDuration());
}
void SSCM3x3Preconditioner_PARDISO::release()
{
    delete ia;
    delete ja;
    delete a;
    delete b;
    delete x;
    delete k_convert;
    phase = -1;           // Release internal memory
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);
}
void SSCM3x3Preconditioner_PARDISO::saveBMP(const char * /*fileName*/, size_t /*size*/) const
{
}
void SSCM3x3Preconditioner_PARDISO::saveProperties(const char * /*fileName*/) const
{
}
size_t SSCM3x3Preconditioner_PARDISO::getElementsNumber() const
{
    return 0;
    //return elementsNumber;
}
void SSCM3x3Preconditioner_PARDISO::InversedSmulVector(Vector3 &x3x3) const
{
    for(size_t i = 0; i < matrixSize; i++)
    {
        b[i*3 + 0] = x3x3[i][0];
        b[i*3 + 1] = x3x3[i][1];
        b[i*3 + 2] = x3x3[i][2];
    }
    // --------------------------------------------------------------------
    // .. Back substitution and iterative refinement
    // --------------------------------------------------------------------
    phase = 33;
    iparm[7] = 2;         // Max numbers of iterative refinement steps
    // Set right hand side to one
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during solution: %d", error);
        exit (3);
    }
    //printf ("\nSolve completed ... ");
    for(size_t i = 0; i < matrixSize; i++)
    {
        x3x3[i][0] = x[i*3 + 0];
        x3x3[i][1] = x[i*3 + 1];
        x3x3[i][2] = x[i*3 + 2];
    }
}
void SSCM3x3Preconditioner_PARDISO::InversedQmulVector(Vector3 &) const
{
}



SSCM3x3Solver_base *SSCM3x3Solver_base::gen(const SolverType solverType)
{
    switch (solverType)
    {
    case SolverType::Direct:
    {
        return new SSCM3x3Solver_direct;
    }break;
    case SolverType::LOS:
    {
        return new SSCM3x3Solver_LOS;
    }break;
    case SolverType::CGM:
    {
        return new SSCM3x3Solver_CGM;
    }break;
    }
}
void SSCM3x3Solver_direct::init(const size_t matrixSize)
{
    r.resize(matrixSize);
    r0.resize(matrixSize);
}
void SSCM3x3Solver_direct::release()
{
    r.clear();
    r0.clear();
}
void SSCM3x3Solver_direct::solve(const SSCM3x3 &matrix, const Vector3 &b, const SSCM3x3Preconditioner_base *p, const Vector3 &x0, const SolverParameters &, Vector3 &x, double &residual, double &relativeResidual, int &iterations, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    double norma_b = sqrt(Vector3ScalMul(b, b));
    if(norma_b < 1.e-20)
    {
        iterations = 0;
        relativeResidual = -1;
        residual = 0;
        for(size_t i = 0; i < x.size(); i++)
            x[i].clear();
        goto ext;
    }
    // прямой решатель
    iterations = 0;
    // невязка первого приближения
    // r0 = A-bx0
    SSCM3x3calcResidual(matrix, b, x0, r0);
    // решение СЛАУ
    // x = Mm1 * b
    VectorCopy(b, x);  // b -> x
    p->InversedMmulVector(x);
    // невязка решения
    // r = A-bx
    SSCM3x3calcResidual(matrix, b, x, r);
    // residual = |b-Ax|/|b|
    residual = sqrt(Vector3ScalMul(r, r)) / norma_b;

    if(Vector3ScalMul(r0, r0) == 0)
        relativeResidual = 100000;
    else
        relativeResidual = sqrt(Vector3ScalMul(r, r) / Vector3ScalMul(r0, r0));
ext:;
    time = t.getCurrentDuration();
}

void SSCM3x3Solver_iterative_base::solve(const SSCM3x3 &matrix, const Vector3 &b, const SSCM3x3Preconditioner_base *p, const Vector3 &x0, const SolverParameters &ssp, Vector3 &x, double &residual, double &relativeResidual, int &iterations, double &time)
{
    TimeIntervals::timeInterval t;
    t.begin();
    TimeIntervals::timeInterval debug_t;
    debug_t.begin();
    double norma_b = sqrt(Vector3ScalMul(b, b));

    if(norma_b < 1.e-20)
    {
        iterations = 0;
        relativeResidual = 1;
        residual = 0;
        for(size_t i = 0; i < x.size(); i++)
            x[i].clear();
        goto ext;
    }
    VectorCopy(x0, x);
    // r0 = A-bx0
    SSCM3x3calcResidual(matrix, b, x, r0);
    residual = sqrt(Vector3ScalMul(r0, r0)) / norma_b;
    iterations = 0;
    // итерации
    for(;;)
    {
        //printf("residual = %le\n", residual);
        if (iterations >= ssp.maxIter)
            break;
        if (residual <= ssp.eps)
            break;
        /*
        if (residual <= ssp.eps &&
            residual >= lastResidual)     // дожимаем
            break;*/
        iterations++;
        //if(iterations % 16 == 0)    // периодические обновления процесса
        //    start(portrait, e, b, p, x);
        start(matrix, b, p, x);
        iter(matrix, b, p, x);
        //lastResidual = residual;
        residual = sqrt(Vector3ScalMul(r, r)) / norma_b;
        //residual = sqrt(r_mul_r) / norma_F;
        if (debug_t.getCurrentDuration() >= 60)   // сообщаем каждую минуту
        {
            //printf("\niter = %6d	residual = %le\n", iterations, residual);
            debug_t.begin();
        }
        //printf("\niter = %6d	residual = %le\n", iterations, residual);
    }
    // r = A-bx
    SSCM3x3calcResidual(matrix, b, x, r);
    // residual = |b-Ax|/|b|
    residual = sqrt(Vector3ScalMul(r, r)) / norma_b;
    // relativeResidual = |b-Ax|/|b-Ax0| = |r|/|r0|
    if(Vector3ScalMul(r0, r0) == 0)
        relativeResidual = 100000;
    else
        relativeResidual = sqrt(Vector3ScalMul(r, r) / Vector3ScalMul(r0, r0));
ext:;
    time = t.getCurrentDuration();
}

void SSCM3x3Solver_LOS::init(const size_t matrixSize)
{
    r.resize(matrixSize);
    r0.resize(matrixSize);
    ss.resize(matrixSize);
    z.resize(matrixSize);
    w.resize(matrixSize);
    aa.resize(matrixSize);
    pp.resize(matrixSize);
}
void SSCM3x3Solver_LOS::release()
{
    r.clear();
    r0.clear();
    ss.clear();
    z.clear();
    w.clear();
    aa.clear();
    pp.clear();
}
void SSCM3x3Solver_LOS::start(const SSCM3x3 &matrix, const Vector3 &b, const SSCM3x3Preconditioner_base *p, Vector3 &x)
{
    // r0 = (b - A*x0)
    SSCM3x3mulVector3(matrix, x, r);        //A*x -> r
    Vector1PlusCmulVector2(b, -1, r, r);//b - r -> r
    // p0 = s0 = Mm1*r0
    VectorCopy(r, pp);  // r -> p
    p->InversedMmulVector(pp);
    VectorCopy(pp, ss); // p -> ss
    // z0 = aa0 = A*p0
    SSCM3x3mulVector3(matrix, pp, z);   // A*p -> z
    VectorCopy(z, aa);  // z -> aa
    // w0 = Mm1*z0
    VectorCopy(z, w);   // z -> w
    p->InversedMmulVector(w);
}
void SSCM3x3Solver_LOS::iter(const SSCM3x3 &matrix, const Vector3 &, const SSCM3x3Preconditioner_base *p, Vector3 &x)
{
    double ca = Vector3ScalMul(w, r) / Vector3ScalMul(w, z);
    Vector1PlusCmulVector2(x, ca, pp, x);   // x + ca * p -> x
    Vector1PlusCmulVector2(r, -ca, z, r);   // r - ca * z -> r
    Vector1PlusCmulVector2(ss, -ca, w, ss); // ss - ca * w -> ss
    SSCM3x3mulVector3(matrix, ss, aa);          // A*ss -> aa
    double cb = -(Vector3ScalMul(w, aa) / Vector3ScalMul(w, z));
    Vector1PlusCmulVector2(ss, cb, pp, pp); // s + cb * p -> p
    Vector1PlusCmulVector2(aa, cb, z, z);   // aa + cb * z -> z
    // w = Mm1*z
    VectorCopy(z, w);   // z -> w
    p->InversedMmulVector(w);
}

void SSCM3x3Solver_CGM::init(const size_t matrixSize)
{
    r.resize(matrixSize);
    r0.resize(matrixSize);
    ss.resize(matrixSize);
    z.resize(matrixSize);
    w.resize(matrixSize);
    aa.resize(matrixSize);
    pp.resize(matrixSize);
}
void SSCM3x3Solver_CGM::release()
{
    r.clear();
    r0.clear();
    ss.clear();
    z.clear();
    w.clear();
    aa.clear();
    pp.clear();
}
void SSCM3x3Solver_CGM::start(const SSCM3x3 &matrix, const Vector3 &b, const SSCM3x3Preconditioner_base *p, Vector3 &x)
{
    // r0 = (b - A*x0)
    SSCM3x3mulVector3(matrix, x, r);        //A*x -> r
    Vector1PlusCmulVector2(b, -1, r, r);//b - r -> r
    // z0 = Mm1*r0
    VectorCopy(r, z);   // r -> z
    p->InversedMmulVector(z);
}
void SSCM3x3Solver_CGM::iter(const SSCM3x3 &matrix, const Vector3 &, const SSCM3x3Preconditioner_base *p, Vector3 &x)
{
    // aa = Mm1*r
    VectorCopy(r, aa);  // r -> aa
    p->InversedMmulVector(aa);
    // A*z -> pp
    SSCM3x3mulVector3(matrix, z, pp);
    double ca1 = Vector3ScalMul(aa, r);
    double ca = ca1 / Vector3ScalMul(pp, z);
    // x + ca * z -> x
    Vector1PlusCmulVector2(x, ca, z, x);
    // r - ca * (A*z) -> r
    Vector1PlusCmulVector2(r, -ca, pp, r);
    // ss = Mm1*r
    VectorCopy(r, ss);  // r -> ss
    p->InversedMmulVector(ss);
    double cb = Vector3ScalMul(ss, r) / ca1;
    // ss + cb*z -> z
    Vector1PlusCmulVector2(ss, cb, z, z);
}


}   // namespace SlauSolving
