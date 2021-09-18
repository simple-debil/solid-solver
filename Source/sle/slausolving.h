/* --------------------------------------------------------- */
// МАТРИЦЫ В РАЗРЕЖЕННОМ СИММЕТРИЧНОМ СТРОЧНО-СТОЛБЦОВОМ ФОРМАТЕ
/* --------------------------------------------------------- */

#ifndef SLAUSOLVING_H
#define SLAUSOLVING_H

#include <set>
#include "elementary.h"


namespace SlauSolving
{
using namespace Elementary;
// итерационный метод решения СЛАУ
enum class SolverType
{
    LOS = 0,  // ЛОС
    CGM = 1,  // Метод сопряжонных градиентов
    Direct = 2,			// Прямой метод (зависит от предобусловливателя)
};
// метод неполного предобуславливания матрицы
enum class Preconditioning
{
    none = 0,
    diag = 1,
    SSCM_LLT_0 = 2,      // строчно-столбцовый, портрет не меняется
    SSCM_LDLT_0 = 3,     // строчно-столбцовый, портрет не меняется
    Profile_LLT = 4,     // профильный, полное предобусловливание
    Profile_LDLT = 5,     // профильный, полное предобусловливание
    SSCM_LLT = 6,        // строчно-столбцовый, полное предобусловливание, портрет может немного расшириться
    SSCM_LDLT = 7,       //строчно-столбцовый, полное предобусловливание, портрет может немного расшириться
    PARDISO = 8,       //PARDISO    ## реализовано только для блочного формата (BSR)
};
// параметры решателя СЛАУ
struct SolverParameters
{
    SolverType type;                    // id решателя
    bool blocks;                        // использовать блочный формат хранения матрицы
    Preconditioning preconditioning;    // предобусловливатель
    double eps;                         // желаемая точность
    int maxIter;                        // максимальное число итераций
};
/*
// параметры решения СЛАУ
struct SolverResultInformation
{
    double slauResidual;                        // невязка решения
    double slauRelativeResidual;                // невязка решения относительно начального приближения
    size_t slauIterations;                         // количество итераций
    double slauTime;                            // время решения в секундах
    size_t slauSolvedWithoutPreconditioner;        // количество раз решали СЛАУ без предобусловливания
    double slauTimeSumWithoutPreconditioner;    // общее время, затраченное на решение СЛАУ (+ предобусловливание), начиная с последнего построения предобусловливателя
};*/

// симметричная строчно-столбцовая матрица (Symmetric String-Column Matrix)

// портрет матрицы в симметричном строчно-столбцовом формате
struct SSCMPortrait
{
    SSCMPortrait(const size_t newElementsNumber = 0, const size_t newMatrixSize = 0);
    // задать размер матрицы и количество внедиагональных элементов матрицы
    void init(const size_t newElementsNumber, const size_t newMatrixSize);
    // построение портрета для LLT разложения
    void init_LLT(const SSCMPortrait &p0);
    // освободить память
    void release();
    void saveBMP(const char *fileName, size_t size)const;
    void saveProperties(const char *fileName)const;
    size_t getMatrixSize()const;
    size_t getElementsNumber()const;
    size_t findUnsorted(const size_t i, const size_t j)const;
    size_t findSorted(const size_t i, const size_t j)const;
    size_t matrixSize;            // размер матрицы
    size_t elementsNumber;        // количество ненулевых элементов в нижнем треугольнике матрицы
    //#### может привести к ошибкам!
    std::vector<size_t> ind;   // a[ind[i]] - первый элемент строки i. (k = ind[i]...(ind[i+1]-1) для строки i)
    std::vector<size_t> ai;    // ai[k] - номер столбца элемента a[k]
};

// симметричная матрица в симметричном профильном формате

// портрет матрицы в симметричном профильном формате
struct ProfilePortrait
{
    ProfilePortrait(const size_t newMatrixSize = 0);
    // задать размер матрицы
    void init(const size_t newMatrixSize);
    // построение портрета для LLT разложения
    void init_LLT(const SSCMPortrait &p0);
    // освободить память
    void release();
    void saveBMP(const char *fileName, size_t size)const;
    void saveProperties(const char *fileName)const;
    size_t getMatrixSize()const;
    size_t getElementsNumber()const;
    inline size_t getFirstCol(const size_t i) const
    {
        return i - (ind[i + 1] - ind[i]);
    }
    size_t matrixSize;         // размер матрицы
    std::vector<size_t> ind;   // a[ind[i]] - первый элемент строки i. (k = ind[i]...(ind[i+1]-1) для строки i)
                            // (ind[i + 1] - ind[i]) - количество элементов в строке i
};

// элементы матрицы в симметричном строчно-столбцовом формате
// или в симметричном профильном формате
struct SSCMElements
{
    SSCMElements(const size_t elementsNumber = 0, const size_t matrixSize = 0);
    // задать портрет матрицы, для которого будут храниться элементы матрицы
    void init(const size_t elementsNumber, const size_t matrixSize);
    // заполнить заданным значением
    void fill(const double value);
    // освободить память
    void release();
    std::vector<double> a; // элементы нижнего треугольника (доступ по индексу)
    std::vector<double> d; // диагональные элементы
};

// матрица в симметричном строчно-столбцовом формате (портрет + элементы)
struct SSCM
{
    SSCMPortrait *p;
    SSCMElements *e;
    void init();
    void release();
    size_t getMatrixSize()const;
    size_t getElementsNumber()const;
    ~SSCM()
    {
        //release();
    }
};

// возвращает true если портрет строки i матрицы (с портретом p) пересекает строку с портретом s
//bool SSCMintersection_fast(const size_t i, const SSCMPortrait &p, const std::vector<bool> &s, const size_t s_j_min, const size_t s_j_max);
bool SSCMintersection(const SSCMPortrait &p, const size_t i1, const size_t i2, const size_t k1_max);

// скопировать элементы матрицы
void SSCMcopy(const SSCMElements &e1, SSCMElements &e2);
// умножение элементов матрицы на число
void SSCMmulScalar(SSCM &matrix, const double c);
// умножение матрицы на вектор
void SSCMmulVector(const SSCM &matrix, const Vector &x, Vector &y);
// вычислить невязку
void SSCMsolveResidual(const SSCM &matrix, const Vector &b, const Vector &x, Vector &r);
// учесть сразу все краевые условия первого рода
void SSCMaddBoundaryCondition1(SSCM &matrix, Vector &b, const std::vector<double>&u0, const std::vector<bool> &state);
// преобразовать строку i к виду (00..010...0|t), 1 на диагонали и обнулить её весь столбец (ниже диагонали) матрицы
void SSCMaddBoundaryCondition1(SSCM &matrix, Vector &b, const size_t i0, const double t);
/*
// вычислить сумму 2-х матриц
void SSCM1plusSSCM2(const double c1, const SSCM &matrix1,
                    const double c2, const SSCM &matrix2,
                    SSCM &matrixRes);
*/
// прибавить к одной матрице вторую матрицу со вложенным портретом
void SSCM1addEnclosedM2mulScalar(SSCM &matrix1, SSCM &matrix2, const double c);
// прибавить к одной матрице вторую матрицу с тем же портретом (сем портрет тут не нужен)
void SSCM1addSimilarM2(SSCM &matrix1, SSCM &matrix2);

// построитель портрета, строчно-столбцовый формат
struct SSCMPortraitBulder
{
    SSCMPortraitBulder(const size_t newMatrixSize = 0);
    void init(const size_t newMatrixSize);
    inline void addElement(const size_t i, const size_t j)
    {
        if(i > j)
            m[i].insert(j);
    }
    void completePortrait(SlauSolving::SSCMPortrait &p);
    std::vector<std::set<size_t>> m; // m[i] - упорядоченное множество номеров столбцов с ненулевыми значениями элементов для строки i
};

// предобусловливатель
struct SSCMPreconditioner_base
{
    static SSCMPreconditioner_base *gen(const Preconditioning preconditioningType);
    // построение портрета разложения и выделение памяти
    virtual void initPortraitAndAllocateMemory(const SSCMPortrait &set_p, double &time) = 0;
    // обновление предобусловливателя начиная со строки firstStr
    virtual void updatePreconditioner(const SSCMElements &set_e, const size_t firstStr, double &time) = 0;
    // инициализировать (вычислить предобусловливатель к матрице M=SQ, неполное предобусловливание с заданным портретом)
     virtual void bulid(const SSCM &matrix, const size_t firstStr, double &time);
    // освободление памяти
    virtual void release() = 0;
    virtual void saveBMP(const char *fileName, size_t size)const = 0;
    virtual void saveProperties(const char *fileName)const = 0;
    virtual size_t getElementsNumber()const = 0;
    // (S^-1)*x -> x
    virtual void InversedSmulVector(Vector &x)const = 0;
    // (Q^-1)*x -> x
    virtual void InversedQmulVector(Vector &x)const = 0;
    virtual ~SSCMPreconditioner_base() = default;
    // (Q^-1)*((S^-1)*x) -> x
    void InversedMmulVector(Vector &x)const;
    const SSCMPortrait *p_ptr0_SSCM = nullptr;   // строчно-столбцовый портрет исходной матрицы
};

// решатель
struct SSCMSolver_base
{
    static SSCMSolver_base *gen(const SolverType solverType);
    virtual void init(const size_t matrixSize) = 0;
    virtual void release() = 0;
    virtual void solve(const SSCM &matrix, const Vector &b, const SSCMPreconditioner_base *p, const Vector &x0, const SolverParameters &ssp,
                       Vector &x, double &residual, double &relativeResidual, int &iterations, double &time) = 0;
    virtual ~SSCMSolver_base() = default;
};

// структура для хранения множества элементов матрицы без определенного формата
struct SlauElement
{
    double a;       // значение элемента матрицы
    size_t j;          // координата элемента в матрице
};

// сборщик матрицы ## устарел?
class SSCMBulder
{
public:
    SSCMBulder(const size_t matrixSize = 0);
    ~SSCMBulder();
    // освободить память
    void release();
    // задать размер матрицы
    void setMatrixSize(const size_t matrixSize);
    // вызвать перед началом сборки
    void start();
    // добавление элемента el в ячейку (i, j) во временный массив, если el не нулевой
    inline void addElement_not_null(const double el, const size_t i, const size_t j)
    {
        if(el != 0 && i >= j)   //####
            str[i].push_back(SlauElement{el, j});
    }
    // добавление элемента el в ячейку (i, j) во временный массив
    inline void addElement(const double el, const size_t i, const size_t j)
    {
        if(i >= j)   //####
            str[i].push_back(SlauElement{el, j});
    }
    // сборка (суммирование) всех добавленных элементов матрицы в разреженую матрицу
    void complete(SSCM &matrix);
    void completeWithPortrait(SSCM &matrix);
    void fixReservedMemory();
private:
    // сравнение 2-х элементов типа SlauElement для сортировки
    static int elementsCmp(const void* x1, const void* x2);
    std::vector<std::vector<SlauElement>> str;  // массив строк из элементов
};

}   // namespace SlauSolving
#endif  // SLAUSOLVING_H
