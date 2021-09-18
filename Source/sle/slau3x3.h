#ifndef SLAU3X3_H
#define SLAU3X3_H

#include "slausolving.h"

// блочное представление матрицы (3x3)
// размер указывается в блоках
namespace SlauSolving
{

using namespace Elementary;

// матрица+обратная матрица для диагональных элементов
struct MATR3x3_Pair
{
    MATR3x3 value;
    MATR3x3 invertValue;
};

// элементы матрицы в симметричном строчно-столбцовом формате
// или в симметричном профильном формате
struct SSCM3x3Elements
{
    SSCM3x3Elements(const size_t elementsNumber = 0, const size_t matrixSize = 0);
    // задать портрет матрицы, для которого будут храниться элементы матрицы
    void init(const size_t elementsNumber, const size_t matrixSize);
    // заполнить заданным значением
    void fill(const MATR3x3 value);
    // освободить память
    void release();
    std::vector<MATR3x3> a; // элементы нижнего треугольника (доступ по индексу)
    std::vector<MATR3x3> d; // диагональные элементы
};

// матрица в симметричном строчно-столбцовом формате (портрет + элементы)
struct SSCM3x3
{
    SSCMPortrait *p;
    SSCM3x3Elements *e;
    void init();
    // ##не эффективная сборка портрета
    void initPortrait3x3by1x1AndAllocateMemory(const SSCMPortrait &p1x1);
    void initMatrixByMatrix1x1(const SSCMPortrait &p1x1, const SSCMElements &e1x1);
    //void init3x3by1x1(const SSCMPortrait &p1x1, const SSCMElements &e1x1);
    void release();
    size_t getMatrixSize()const;
    size_t getElementsNumber()const;
    ~SSCM3x3()
    {
        //release();
    }
};

// x -> y
void Vector1x1CopyToVector3x3(const Vector &x, Vector3 &y);
// x -> y
void Vector3x3CopyToVector1x1(const Vector3 &x, Vector &y);
// скопировать элементы матрицы
void SSCM3x3copy(const SSCM3x3Elements &e1, SSCM3x3Elements &e2);
// умножение элементов матрицы на число
void SSCM3x3mulScalar(SSCM3x3 &matrix, const double c);
// умножение матрицы на вектор
void SSCM3x3mulVector3(const SSCM3x3 &matrix, const Vector3 &x, Vector3 &y);
// вычислить невязку
void SSCM3x3calcResidual(const SSCM3x3 &matrix, const Vector3 &b, const Vector3 &x, Vector3 &r);
// учесть сразу все краевые условия первого рода
void SSCM3x3addBoundaryCondition1(SSCM3x3 &matrix, Vector3 &b, const std::vector<double>&u0, const std::vector<bool> &state);
// прибавить к одной матрице вторую матрицу со вложенным портретом
void SSCM3x3_1addEnclosedM3x3_2mulScalar(SSCM3x3 &matrix1, SSCM3x3 &matrix2, const double c);
// прибавить к одной матрице вторую матрицу с тем же портретом (сем портрет тут не нужен)
void SSCM3x3_1addSimilarM3x3_2(SSCM3x3 &matrix1, SSCM3x3 &matrix2);

// предобусловливатель
struct SSCM3x3Preconditioner_base
{
    static SSCM3x3Preconditioner_base *gen(const Preconditioning preconditioningType);
    // построение портрета разложения и выделение памяти
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) = 0;
    // обновление предобусловливателя начиная со строки firstStr
    virtual void updatePreconditioner(const SSCM3x3Elements &set_e, const size_t firstStr, double &time) = 0;
    // инициализировать (вычислить предобусловливатель к матрице M=SQ, неполное предобусловливание с заданным портретом)
     virtual void bulid(const SSCM3x3 &matrix, const size_t firstStr, double &time);
    // освободление памяти
    virtual void release() = 0;
    virtual void saveBMP(const char *fileName, size_t size)const = 0;
    virtual void saveProperties(const char *fileName)const = 0;
    virtual size_t getElementsNumber()const = 0;
    // (S^-1)*x -> x
    virtual void InversedSmulVector(Vector3 &x)const = 0;
    // (Q^-1)*x -> x
    virtual void InversedQmulVector(Vector3 &x)const = 0;
    virtual ~SSCM3x3Preconditioner_base() = default;
    // (Q^-1)*((S^-1)*x) -> x
    void InversedMmulVector(Vector3 &x)const;
    const SSCMPortrait *p_ptr0_SSCM = nullptr;   // строчно-столбцовый портрет исходной матрицы
};

// решатель
struct SSCM3x3Solver_base
{
    static SSCM3x3Solver_base *gen(const SolverType solverType);
    virtual void init(const size_t matrixSize) = 0;
    virtual void release() = 0;
    virtual void solve(const SSCM3x3 &matrix, const Vector3 &b, const SSCM3x3Preconditioner_base *p, const Vector3 &x0, const SolverParameters &ssp,
                       Vector3 &x, double &residual, double &relativeResidual, int &iterations, double &time) = 0;
    virtual ~SSCM3x3Solver_base() = default;
};

}   // namespace SlauSolving

#endif // SLAU3X3_H
