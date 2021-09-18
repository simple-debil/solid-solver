/* --------------------------------------------------------- */
// ВХОДНЫЕ И ВЫХОДНЫЕ ПАРАМЕТРЫ
/* --------------------------------------------------------- */

#ifndef SOLID_BASE_H
#define SOLID_BASE_H

#include "elementary.h"
#include "funParser.h"
#include "slausolving.h"
#include "grid.h"


namespace Solid
{
using namespace Elementary;

// типы материала
enum class MechPlasticityMethodType
{
    Elasticity = 0,                 // упругий
    D_pl = 1,                       // только C_pl
    D_pl_Solov = 2,                 // только C_pl, алгоритм Соловейчика
    InitialSigma = 4,               // метод начальных напряжений
    Combo_D_pl_InitialSigma = 5,    // метод начальных напряжений и C_pl
};
// от каких параметров зависит кривая деформирования
enum class MechPlasticityCurveDependenceType
{
    None,   // напряжения зависят от деформаций / деформации зависят от напряжений
    Time,   // напряжения зависят от деформаций и времени / деформации зависят от напряжений и времени
};
// поведение при разгрузке
enum class MechPlasticityCurveUnloadingType
{
    Curve,   // разгрузка по кривой (нелинейная упругость, ползучесть)
    Elastic, // упругая разгрузка (пластичность)
};
// вариант кривой деформирования (sigma(eps) или eps(sigma))
enum class MechPlasticityCurveType
{
    Sigma_eps,  // напряжения зависят от деформаций
    Eps_sigma,  // ####деформации зависят от напряжений
};

// общий результат итерации
enum class MechIterationsGeneralResult
{
    Interrupted_Ok,                     // итерации успешно завершились после достижения требуемой точности
    Interrupted_OutOfIterationsLimit,   // точность не достигнута, итерации прекратились из-за ограничения на количество итераций
    Interrupted_Terminated,             // итерационный процесс прерван пользователем
    Continue_NeedMoreIterations,        // нужно продолжать итерации
};

enum class SwitchIterationsMode
{
    Serial = 0,     // итерации по контакту и пластичности запускаются по-очереди
    Parallel = 1,   // итерации по контакту и пластичности происходят параллельно
};

// режим итераций для учёта различных нелинейностей
enum class IterationMode
{
    ReadOnly = 0,   // в конце итерации вычисляется следующее приближение, но состояния материала/контактов не меняются
    Free = 1,       // полноценные итерации
};

// режим учёта сил в правой части СЛАУ
enum class IncForsesMode
{
    MinusIntegral = 1,  // Я: b = integral(sumP) - integral(sumSigma) - вычисляется для текущей сетки в начале шага
    IncrementP = 2,     // Без коррекции: R = 0
    bPlusR = 3,        // Темис: R = integral(sumSigma+deltaSigma) - integral(sumP+deltaP) - вычисляется для предыдущей сетки в конце шага
    //bMinusR = 3,        // учёт невязки в правой части: b = sumP - Sum(R1 + R2 + ...)
};

// невязки напряжений и деформаций
struct MaxSigmaEpsResidual
{
    double sigma;
    double eps;
    void set_null();
    void set_undefined();
    void set_abs(const MaxSigmaEpsResidual &r);
    void compareAndSaveIfLarger(const MaxSigmaEpsResidual &r);
};
const MaxSigmaEpsResidual MaxSigmaEpsResidual_NULL = {0, 0};

// упругие параметры материала
struct ElasticParameters
{
    double Talpha;      // коэффициент линейного расширения
    double ro;          // плотность
    double K;           // модуль объёмной упругости (модуль объёмного сжатия)
    double E;           // модуль Юнга (модуль упругости)
    double LAMBDA;      // параметр Ламе 1
    double G;           // модуль сдвига (модуль жесткости, параметр Ламе 2) = MU
    double NU;          // коэффициент Пуассона
    void clear()
    {
        Talpha = 0;
        ro = 0;
    }
};
struct TanAngle
{
    bool is_elastic;    // true - упругость
    bool is_infinum;    // true - угол равен pi/2 или -pi/2
    double value;       // угол наклона, по которому строится D_pl (игнорируется если tan_betta_is_elastic = true или tan_betta_is_infinum = true)
};

struct TanAngle_Yeld
{
    bool is_elastic;    // true - упругость
    double value;       // значение, если не упругость = dsigma/deps_pl = dsigma/dq
    void setValue(const double tan_elastic, const double dsigma_want, const double deps_want);
};


struct TENSOR2
{
    double m[3][3];
    TENSOR2(const TENSOR2 &) = default;
    TENSOR2()
    {
        clear();
    }
    // обнуление
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
    // эквивалентное напряжение
    double deviatorEuclideanNorm() const
    {
        TENSOR2 dev = this->deviator();
        double E = 0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                E += dev.m[i][j]*dev.m[i][j];
        return sqrt(E);
    }
    // эквивалентное напряжение
    double eqv_SIGMA() const
    {
        TENSOR2 dev = this->deviator();
        double E = 0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                E += dev.m[i][j]*dev.m[i][j];
        return sqrt(3./2.*E);
    }
    // эквивалентная деформация
    double eqv_EPS() const
    {
        TENSOR2 dev = this->deviator();
        double E = 0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                E += dev.m[i][j]*dev.m[i][j];
        return sqrt(2./3.*E);
    }
    // объёмная часть девиатора
    double hydrostatic() const
    {
        return (m[0][0] + m[1][1] + m[2][2]) / 3.;
    }
    // единичный девиатор напряжений
    TENSOR2 deviator() const
    {
        TENSOR2 deviator = *this;
        double hydrostatic = deviator.hydrostatic();
        deviator.m[0][0] -= hydrostatic;
        deviator.m[1][1] -= hydrostatic;
        deviator.m[2][2] -= hydrostatic;
        return deviator;
    }
    // детерминант
    double calcDet()const
    {
        return
                + m[0][0] * m[1][1] * m[2][2]
                + m[0][1] * m[1][2] * m[2][0]
                + m[0][2] * m[1][0] * m[2][1]
                - m[0][2] * m[1][1] * m[2][0]
                - m[0][0] * m[1][2] * m[2][1]
                - m[0][1] * m[1][0] * m[2][2];
    }
    void calcInvariants(VECTOR3 &I)const
    {
        I[0] = m[0][0] + m[1][1] + m[2][2];
        I[1] = m[0][0]*m[1][1] + m[1][1]*m[2][2] + m[2][2]*m[0][0] - SQR(m[1][2]) - SQR(m[2][0]) - SQR(m[0][1]);
        I[2] = calcDet();
    }
    void ToMATR3x3_sigma(MATR3x3 &sigma3x3)const
    {
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
            {
                sigma3x3.m[i][j] = m[i][j];
            }
    }

    // сделать полностью симметричной (усредняет симметричные внедиагональные элементы)
    /*
    void fixAsymmetry()
    {
        m[1][0] = (m[1][0] + m[0][1]) / 2;
        m[0][1] = m[1][0];
        m[2][0] = (m[2][0] + m[0][2]) / 2;
        m[0][2] = m[2][0];
        m[2][1] = (m[2][1] + m[1][2]) / 2;
        m[1][2] = m[2][1];
    }
    */
    // : (двойная свёртка, аналог скалярного произведения)
    friend inline double operator*(const TENSOR2 &l, const TENSOR2 &r)
    {
        double E = 0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                E += l.m[i][j] * r.m[i][j];
            }
        return E;
    }
    friend inline const TENSOR2 operator+(const TENSOR2 &l, const TENSOR2 &r)
    {
        TENSOR2 res;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                res.m[i][j] = l.m[i][j] + r.m[i][j];
        return res;
    }
    friend inline const TENSOR2 operator-(const TENSOR2 &l, const TENSOR2 &r)
    {
        TENSOR2 res;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                res.m[i][j] = l.m[i][j] - r.m[i][j];
        return res;
    }
    friend inline const TENSOR2 operator*(const TENSOR2 &l, const double &r)
    {
        TENSOR2 res;
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                res.m[i][j] = l.m[i][j]*r;
            }
        }
        return res;
    }
    friend inline const TENSOR2 operator/(const TENSOR2 &l, const double &r)
    {
        TENSOR2 res;
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                res.m[i][j] = l.m[i][j]/r;
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
    inline void operator/=(const double &t)
    {
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                m[i][j] /= t;
            }
        }
    }
    inline void operator+=(const TENSOR2 &t)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                m[i][j] += t.m[i][j];
    }
    inline void operator-=(const TENSOR2 &t)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                m[i][j] -= t.m[i][j];
    }
};

struct TENSOR4
{
    double m[3][3][3][3];
    TENSOR4(const TENSOR4 &) = default;
    TENSOR4()
    {
        clear();
    }
    // обнуление матрицы
    void clear()
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                        m[i][j][k][l] = 0;
    }
    // : (двойная свёртка, аналог скалярного произведения)
    friend inline const TENSOR2 operator*(const TENSOR4 &left, const TENSOR2 &right)
    {
        TENSOR2 res;
        res.clear();
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                    {
                        res.m[i][j] += left.m[i][j][k][l] * right.m[k][l];
                    }
        return res;
    }
    static void setCY_pl(const ElasticParameters &ep_T, const TanAngle_Yeld &A, const TENSOR2 &z,
                         TENSOR4 &C_el, TENSOR4 &C_pl, TENSOR4 &Y_pl)
    {
        double MU = ep_T.G;
        double LAMBDA = ep_T.LAMBDA;
        if(A.is_elastic)
        {
            // упругий случай (A.value = бесконечность)
            for(int i = 0; i < 3; i++)
                for(int j = 0; j < 3; j++)
                    for(int k = 0; k < 3; k++)
                        for(int l = 0; l < 3; l++)
                        {
                            double el_value = LAMBDA*(i==j)*(k==l) + MU*((i==k)*(j==l) + (i==l)*(j==k));
                            Y_pl.m[i][j][k][l] = 0;
                            C_pl.m[i][j][k][l] = el_value;
                            C_el.m[i][j][k][l] = el_value;
                        }
        }
        else
        {
            double coef = 2*MU / (A.value + 3*MU);
            for(int i = 0; i < 3; i++)
                for(int j = 0; j < 3; j++)
                    for(int k = 0; k < 3; k++)
                        for(int l = 0; l < 3; l++)
                        {
                            double el_value = LAMBDA*(i==j)*(k==l) + MU*((i==k)*(j==l) + (i==l)*(j==k));
                            double Y_value = coef * z.m[i][j]*z.m[k][l];
                            Y_pl.m[i][j][k][l] = Y_value;
                            C_pl.m[i][j][k][l] = el_value - Y_value*2*MU;
                            C_el.m[i][j][k][l] = el_value;
                        }
        }
    }
    static void setC_el(const ElasticParameters &ep_T,
                         TENSOR4 &C_el)
    {
        double MU = ep_T.G;
        double LAMBDA = ep_T.LAMBDA;
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                for(int k = 0; k < 3; k++)
                    for(int l = 0; l < 3; l++)
                    {
                        double el_value = LAMBDA*(i==j)*(k==l) + MU*((i==k)*(j==l) + (i==l)*(j==k));
                        C_el.m[i][j][k][l] = el_value;
                    }
    }
    // приращение обратного тензора упругости
    static void setInvertDeltaIsotropicC(const ElasticParameters &ep_T1, const ElasticParameters &ep_T2,
                                         TENSOR4 &deltaA)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                    {
                        double a1 = 1./ep_T1.E*((1. + ep_T1.NU)/2.*((i==k)*(j==l) + (i==l)*(j==k)) - ep_T1.NU*(i==j)*(k==l));
                        double a2 = 1./ep_T2.E*((1. + ep_T2.NU)/2.*((i==k)*(j==l) + (i==l)*(j==k)) - ep_T2.NU*(i==j)*(k==l));
                        deltaA.m[i][j][k][l] = a2 - a1;
                    }
    }
};


// параметры материала
struct MechMaterialSource
{
    MechPlasticityMethodType plasticityMethodType;       // тип материала
    MechPlasticityCurveDependenceType PCDependenceType;  // от каких параметров зависит кривая деформирования
    MechPlasticityCurveUnloadingType PCUnloadingType;    // разгрузка по кривой/упругая разгрузка
    MechPlasticityCurveType PCType;// eps(sigma) или sigma(eps) - только для MechPlasticityMethodType::Plasticity
    double w_midPoint;          // для направляющего тензора
    double w_project;      // параметр сглаживания прибавки к начальным напряжениям (1 - нет сглаживания, ~0 - сильно сглаживать)
    double cosTettaMin = 0;     // минимальный косинус угла для разгрузки заранее
    bool temperatureDependence; // = true если D зависит от температуры

    ElasticParameters elasticParameters0;    // упругие параметры материала (начальные)
     VECTOR3 F;         // объемные силы (например z=Ro*g)  (##)тут можно задать вектор-функцию F
    MATR6x6 M_sigma;    // матрица для вычисления эквивалентных деформаций (пластичная изотропия / #ортотропия)
    double Csigma;
    double Ceps;
    FunParser::Function epsFun;         // eps(sigma) - кривая пластичности
    FunParser::Function difEpsFun;      // deps/dsigma - производная кривой пластичности
    FunParser::Function sigmaFun;       // sigma(eps) - кривая пластичности
    FunParser::Function difSigmaFun;    // dsigma/deps - производная кривой пластичности

    FunParser::Function E_Fun;        // E(T) - зависимость модуля Юнга от температуры
    FunParser::Function NU_Fun;       // NU(T) - зависимость коэффициента Пуассона от температуры

     double elasticSigmaLimit;      // предел упругости (для идеально пластичного материала)
     double elasticEpsLimit;        // предел упругости (для идеально пластичного материала)
     bool is2Dxy; // плоскодеформированное НДС, eps33 = 0

    MechMaterialSource();
    // деформации в зависимости от напряжений (вспомогательная)
    double eps_dep(const double sigmaEqv, const double t);
    // deps/dsigma - тангенс угла наклона кривой деформирования (вспомогательная)
    double difEps_dep(const double sigmaEqv, const double t);
    // напряжения в зависимости от деформаций (вспомогательная)
    double sigma_dep(const double epsEqv, const double t);
    // deps/dsigma - тангенс угла наклона кривой деформирования (вспомогательная)
    double difSigma_dep(const double epsEqv, const double t);

    // деформации в зависимости от напряжений
    double eps(const ElasticParameters &ep_T, const double sigmaT, const double sigmaEqv, const double t);
    // тангенс угла наклона кривой деформирования (deps/dsigma)
    double difEps(const ElasticParameters &ep_T, const double sigmaT, const double sigmaEqv, const double t);
    // напряжения в зависимости от деформаций
    double sigma(const ElasticParameters &ep_T, const double epsT, const double epsEqv, const double t);
    // тангенс угла наклона кривой деформирования (deps/dsigma)
    double difSigma(const ElasticParameters &ep_T, const double epsT, const double epsEqv, const double t);

    // напряжения в зависимости от параметра Одквиста
    double sigma_Yeld(const ElasticParameters &ep_T, const double q, const double t);
    // тангенс угла наклона кривой sigma(q) т.е. (dsigma/dq)
    double difSigma_Yeld(const ElasticParameters &ep_T, const double q, const double t);
    // параметр Одквиста в зависимости от напряжения
    double q_Yeld(const ElasticParameters &ep_T, const double sigmaEqv, const double t);
    // 1/тангенс угла наклона кривой q(sigma) т.е. (dsigma/dq)
    double difq_Yeld(const ElasticParameters &ep_T, const double sigmaEqv, const double t);

    // расчёт тангенс угла упругого наклона
    double tan_el(const ElasticParameters &ep_T);


    // расчёт упругих параметров материала при температуре T
    void calc_elastic_parameters(const double T, ElasticParameters &ep_T);
    // расчёт матрицы определяющих соотношений D~ в зависимости от температуры
    void calcD(const ElasticParameters &ep_T, MATR6x6 &D);
    // расчёт обратного тензора определяющих соотношений
    void calcInvertIsotropicD(const ElasticParameters &ep_T, MATR3x3x3x3 &A);
    // расчёт приращения обратного тензора определяющих соотношений вследствие изменения температуры
    void calcInvertDeltaIsotropicD(const ElasticParameters &ep_T1, const ElasticParameters &ep_T2, MATR3x3x3x3 &deltaA);
    // задание матрицы определяющих соотношений для изотропного упругого материала
    // она же будет начальной матрицей для пластичного изотропного материала
    void set_K_E(const double K, const double E);
    void set_K_LAMBDA(const double K, const double LAMBDA);
    void set_K_G(const double K, const double G);
    void set_K_NU(const double K, const double NU);
    void set_E_LAMBDA(const double E, const double LAMBDA);
    void set_E_G(const double E, const double G);
    void set_E_NU(const double E, const double NU);
    void set_LAMBDA_G(const double LAMBDA, const double G);
    void set_LAMBDA_NU(const double LAMBDA, const double NU);
    // задание матрицы Мизеса, по котоорй будут расчитываться эквивалентные деформации и напряжения
    void set_M_sigma();
    // задание плоскодеформированного НДС
    // ##(только простейшая упругость без температуры, т.к. не сделано плоское обращение матрицы)
    void set_2Dxy();
protected:
    // удаление строк и столбцов 1,2,3 из матрицы определяющих соотношений (чтобы сделать 2d)
    void set_matrix_to_2D(MATR6x6 &D);
};

// ресурс 1-го краевого условия
struct MechBoundaryCondition1Source
{
    VECTOR3_int mode;   // -1 - отсутствует
    VECTOR3 u0;         // вектор перемещения
};

// ресурс 2-го краевого краевого условия
// базовый класс
struct MechBoundaryCondition2Source_base
{
    // обновление значения поверхностных сил в моменты времени t_prev и t
    virtual void updateValue(const double t_prev1, const double t0) = 0;
    virtual VECTOR3 calcVectorValueInPoint1(const VECTOR3 &normal) = 0;
    virtual VECTOR3 calcVectorValueInPoint2(const VECTOR3 &normal) = 0;
};
// ресурс 2-го краевого условие
// отсутствие краевого условия
struct MechBoundaryCondition2Source_None:
    public MechBoundaryCondition2Source_base
{
    virtual void updateValue(const double, const double) override;
    virtual VECTOR3 calcVectorValueInPoint1(const VECTOR3 &) override;
    virtual VECTOR3 calcVectorValueInPoint2(const VECTOR3 &) override;
};
// векторная константа
struct MechBoundaryCondition2Source_VectorConstant:
    public MechBoundaryCondition2Source_base
{
    VECTOR3 P1;          // значение поверхностных сил
    VECTOR3 P2;          // значение поверхностных сил
    virtual void updateValue(const double, const double) override;
    virtual VECTOR3 calcVectorValueInPoint1(const VECTOR3 &) override;
    virtual VECTOR3 calcVectorValueInPoint2(const VECTOR3 &) override;
    void init(const VECTOR3 set_P);
};
// скалярная константа (поверхностные силы перпендикулярны поверхности)
struct MechBoundaryCondition2Source_ScalarConstant:
    public MechBoundaryCondition2Source_base
{
    double P1;          // значение поверхностных сил
    double P2;          // значение поверхностных сил
    virtual void updateValue(const double, const double) override;
    virtual VECTOR3 calcVectorValueInPoint1(const VECTOR3 &normal) override;
    virtual VECTOR3 calcVectorValueInPoint2(const VECTOR3 &normal) override;
    void init(const double set_P);
};
// векторная функция
struct MechBoundaryCondition2Source_VectorFunction:
    public MechBoundaryCondition2Source_VectorConstant
{
    FunParser::Function *Pf[3];
    virtual void updateValue(const double t_prev1, const double t0) override;
    VECTOR3 value(const double t);
    void init(FunParser::Function *set_Pf1, FunParser::Function *set_Pf2, FunParser::Function *set_Pf3);
};
// скалярная функция (поверхностные силы перпендикулярны поверхности)
struct MechBoundaryCondition2Source_ScalarFunction:
    public MechBoundaryCondition2Source_ScalarConstant
{
    FunParser::Function *Pf;
    virtual void updateValue(const double t_prev1, const double t0) override;
    double value(const double t);
    void init(FunParser::Function *set_Pf);
};

// результаты в узлах сетки
struct MechOutVertexData
{
    VECTOR3 p;              // положение точки
    VECTOR3 sum_du;         // суммарные перемещения за все время движения
    bool contact;           // наличие контакта
    VECTOR3 F_sum;          // сила реакции опоры (вектор)
    double h;
    double stiffness;
};
// результаты в точке кэ
struct MechOutFePointData
{
    POINT3 p;           // координаты центра КЭ
    double q;              // параметр Одквиста - сумма эквивалентных приращений пластических деформаций
    TENSOR2 sumEpsElastic; // упругая деформация
    TENSOR2 sumEpsPlastic; // пластическая деформация
    TENSOR2 sumEpsTerm;    // температурная деформация
    TENSOR2 sumSigma;      // напряжение
    bool isYelding;  // пластичное состояние, т.е. предел упругости был превышен и после последнего активного нагружения не было разгрузки
    MaxSigmaEpsResidual residual;       // невязки
    void mainStressesXY(VECTOR3 &ms);
    void mainStresses(VECTOR3 &ms) const;
};
// результаты в КЭ
struct MechOutFeData
{
    std::vector<MechOutFePointData> pd;
    std::vector<VECTOR3> q;
    //VECTOR3 q[Grid::FE_LinearHexagon_XFEM::DOF_max];      // ###значение коэффициентов перед базисными ф-ми (из решения СЛАУ)
    void init(const int setSize);
};

// параметры глобального шага
struct MechGlobalStepParameters
{
    // решение СЛАУ
    SlauSolving::SolverParameters slausolverParameters;        // параметры решателя СЛАУ
    // метод дискретизации по времени
    int timeMode;   // 0 - квазистатический, 1 - ..
    IncForsesMode incForsesMode; // 1 - dSigma = integral(...), 2 - приращения сил, dSigma = 0
    // фиксация сетки или способ смещения сетки
    int fixGrid;                // 0 - смещать сетку после каждого шага по времени, 1 - фиксировать положение сетки
    int movingGridMode;         // 0 - простые приращения, 1## - приращения Яумана
    // параметры завершение итераций на шаге по времени
    SwitchIterationsMode switchIterationsMode;  // последовательные или параллельные итерации
    int controlMode;            // 0 - игнорировать не достижение невязок, 1## - изменять шаг по времени для достижения невязок, 2## - завершать процесс в случае не достижения невязки
    bool terminateIfAccuracyIsNotAchieving; // true - завершать итерационный процесс, если он начал расходиться
    int iterLimit;              // максимальное количество итераций на 1 шаг (решение нелинейной задачи)
    double slauResidualLimit;   // условие завершения по изменению решения СЛАУ
    MaxSigmaEpsResidual plasticResidualLimit;
    double contactEndPointResidualLimit;// условие завершения по изменению отклонения от положения контакта
    double contactDeltaFResidualLimit;  // условие завершения по результирующей силы реакции
    std::vector<MechMaterialSource> *material;                      // материалы
    std::vector<MechBoundaryCondition1Source> *bc1Source;           // 1-e краевые условия
    std::vector<MechBoundaryCondition2Source_base *> *bc2Source;    // 2-е краевые
    std::vector<Grid::MechBoundaryCondition2> *bc2;                 // 2-е краевые: индекс КЭ-поверхности + индекс MechBoundaryCondition2Source_base
};

// данные для вершин
struct MechVertexData
{
    VECTOR3 sum_du;      // суммарное перемещение за все время движения
    VECTOR3 du;          // приращения перемещений за 1 шаг по времени
};

}   // namespace Solid
#endif  // SOLID_BASE_H
