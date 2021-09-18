/* --------------------------------------------------------- */
// ПАРСЕР ФУНКЦИЙ
/* --------------------------------------------------------- */

#ifndef FUN_PARSER_H
#define FUN_PARSER_H

#include <string>

#include "elementary.h"
#include "interpolation.h"

// парсер функций
namespace FunParser
{

// тип разделителя
enum class SeparatorType
{
    Comma,   // ,
    Space,   // _
};
// тип скобки
enum BracketType
{
    Bracket1,   // ()
    Bracket2,   // []
    Bracket3,   // {}
    Bracket4,   // <>
};
// тип арифметического оператора
enum OperatorType
{
    OperatorPow,   // ^
    OperatorMul,   // *
    OperatorDiv,   // /
    OperatorAdd,   // +
    OperatorSub,   // -
};
// тип элементарной функции
enum ElementaryFunctionType
{
    EFunPow,        // 0
    EFunLog,        // 1
    EFunExp,        // 2
    EFunLn,         // 3
    EFunLog10,      // 4
    EFunSin,        // 5
    EFunCos,        // 6
    EFunTg,         // 7
    EFunArcsin,     // 8
    EFunArccos,     // 9
    EFunArctg,      // 10
    EFunAbs,        // 11
    EFunSign,        // 12
    _EFunType_Size, //
};
// тип токена
enum TokenType
{
    TokenOpen,      // 0 открывающая скобка
    TokenClose,     // 1 закрывающая скобка
    TokenSeparator, // 2 разделитель
    TokenOperator,  // 3 арифметический оператор
     TokenElementaryFunction,  // 4 элементарная функция
     TokenFunction, // 5 функция
     //TokenSpline, //(сплайн-функция)
    TokenNumber,    // 6 действительное число
    TokenVariable,  // 7 аргумент

};
// токен
struct Token
{
    TokenType type;     // тип
    size_t id;          // идентификатор
};
// аргумент функции
struct Argument
{
    std::string name;   // название аргумента
    double value;       // значение аргумента
};
/*
// строковое представление числа
class NumberString
{
public:
    NumberString(const char *expression = "");
    ~NumberString();
    void setExpression(const char *expression);
    double getValue()const;
    std::string getExpression();
protected:
    std::string s;
};
*/


// Сплайн

////////////////////////////////////////////////////////////////////////////////////////////////
class cubic_spline
{
private:
    // Структура, описывающая сплайн на каждом сегменте сетки
    struct spline_tuple
    {
        double a, b, c, d, x;
    };

    spline_tuple *splines; // Сплайн
    std::size_t n; // Количество узлов сетки

    void free_mem(); // Освобождение памяти

public:
    cubic_spline(); //конструктор
    ~cubic_spline(); //деструктор

    // Построение сплайна
    // x - узлы сетки, должны быть упорядочены по возрастанию, кратные узлы запрещены
    // y - значения функции в узлах сетки
    // n - количество узлов сетки
    void build_spline(const std::vector<Elementary::VECTOR2> &points);

    // Вычисление значения интерполированной функции в произвольной точке
    double f(double x) const;
};
////////////////////////////////////////////////////////////////////////////////////////////////
enum class myFunType
{
    Expression,     // 0
    Spline,         // 1
    Bezie,          // 2 - пластичность, eps(sigma), Безье
    difBezie,       //
    PlasticSigma,   // 3 - пластичность, sigma(eps), 2 линейных участка
    difPlasticSigma,//
     PlasticSigma_yeld,
     diffPlasticSigma_yeld,
     PlasticSigma_yeld_hardening,
     diffPlasticSigma_yeld_hardening,
    PlasticEps,
    difPlasticEps,
    CreepSigma,     // 4 - ползучесть, sigma(eps), начинается с линейного участка
    difCreepSigma,  //
    CreepEps,       // 5 - ползучесть, eps(sigma)
    difCreepEps,    //
};

// Класс парсера функций
class Function
{
public:
    Function(const std::string expression = "",
             const std::string functionName = "");
    ~Function();
    void setName(const std::string functionName);
    void setExpression(const std::string expression);
    const std::string getName()const;
    const std::string getExpression()const;
    int addArgument(const std::string name, const double value = 0);
    void addFunction(Function *fun);
    void setArgumentValue(const int id, const double value);
    void parse();
    double calcValue();
///////////////
    void setBezie(double E1, Elementary::POINT3 newp1, Elementary::POINT3 newp2, Elementary::POINT3 newp3, double E2);
    void setDifBezie(double E1, Elementary::POINT3 newp1, Elementary::POINT3 newp2, Elementary::POINT3 newp3, double E2);
    void setPlasticSigma(Elementary::POINT3 newp1, double E1, double E2);
    void setDifPlasticSigma(Elementary::POINT3 newp1, double E1, double E2);
     void setPlasticSigma_Yeld(double elasticSigmaLimit, double E1, double E2);
     void setDiffPlasticSigma_Yeld(double elasticSigmaLimit, double E1, double E2);
     void setPlasticSigma_Yeld_hardening(double elasticSigmaLimit, double E, double n);
     void setDiffPlasticSigma_Yeld_hardening(double elasticSigmaLimit, double E, double n);
    void setPlasticEps(Elementary::POINT3 newp1, double E1, double E2);
    void setDifPlasticEps(Elementary::POINT3 newp1, double E1, double E2);
    void __setCreep(const double A, const double n, const double m, const double tan0);
    void setCreepSigma(const double A, const double n, const double m, const double tan0);
    void setDifCreepSigma(const double A, const double n, const double m, const double tan0);
    void setCreepEps(const double A, const double n, const double m, const double tan0);
    void setDifCreepEps(const double A, const double n, const double m, const double tan0);
///////////////
    void setSpline(const std::vector<Elementary::VECTOR2> &points);
    double calcValue(const double x);
    double calcValue(const double x, const double y);
    double calcDif(const double x);
    double calcInverse(const double y);
//private:
    //Interpolation::Interpolant1D_Lagrange1 it;
    Interpolation::Interpolant1D_Lagrange1Unregular it;
    //Interpolation::Interpolant1D_Hermite3 it;
    double E1, E2;
     double elasticSigmaLimit_pl;
    Elementary::POINT3 p1;
    Elementary::POINT3 p2;
    Elementary::POINT3 p3;
    cubic_spline spline;
    double creep_A;
    double creep_n;
    double creep_m;
    double creep_tan0;
    myFunType mode;
    std::string name;      // имя функции
    std::string s;         // текстовое выражение функция
    std::vector<double> nums;       // действительные числа
    std::vector<Argument> args;    // аргументы
    std::vector<Function *> funs;  // функции
     std::vector<Token> tokens;    // #в нем нет необходимости#функция в виде массива токенов, исходная(инфиксная) запись
    std::vector<Token> ts;         // стек токенов
    std::vector<double> ns;        // стек чисел
    std::vector<Token> RPN;        // функция в виде массива токенов, постфиксная запись
     void printToken(const Token &t);
    //bool FirstTokenOfOperand(const Token &t);
    //bool LastTokenOfOperand(const Token &t);
    Token *pushWhileNotOpen();
    void translate();
    bool isLetter(const char c);
    bool isNumber(const char c);
    int charToNumber(const char c);
    bool firstWordEqual(const std::string str, const std::string word);
};

double strToDouble(const std::string expression);
bool strToDoubleTest(const std::string expression);

}   //namespace Fun_parser

#endif  // FUN_PARSER_H
