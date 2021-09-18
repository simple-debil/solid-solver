#include "stdio.h"
#include "string.h"

#include "funParser.h"
#include "interpolation.h"

namespace FunParser
{


cubic_spline::cubic_spline() : splines(nullptr)
{

}
cubic_spline::~cubic_spline()
{
    free_mem();
}
void cubic_spline::build_spline(const std::vector<Elementary::VECTOR2> &points)
{
    free_mem();

    this->n = points.size();

    // Инициализация массива сплайнов
    splines = new spline_tuple[n];
    for (std::size_t i = 0; i < n; ++i)
    {
        splines[i].x = points[i].x[0];
        splines[i].a = points[i].x[1];
    }
    splines[0].c = 0.;

    // Решение СЛАУ относительно коэффициентов сплайнов c[i] методом прогонки для трехдиагональных матриц
    // Вычисление прогоночных коэффициентов - прямой ход метода прогонки
    double *alpha = new double[n - 1];
    double *beta = new double[n - 1];
    double A, B, C, F, h_i, h_i1, z;
    alpha[0] = beta[0] = 0.;
    for (std::size_t i = 1; i < n - 1; ++i)
    {
        h_i = points[i].x[0] - points[i - 1].x[0], h_i1 = points[i + 1].x[0] - points[i].x[0];
        A = h_i;
        C = 2. * (h_i + h_i1);
        B = h_i1;
        F = 6. * ((points[i + 1].x[1] - points[i].x[1]) / h_i1 - (points[i].x[1] - points[i - 1].x[1]) / h_i);
        z = (A * alpha[i - 1] + C);
        alpha[i] = -B / z;
        beta[i] = (F - A * beta[i - 1]) / z;
    }

    splines[n - 1].c = (F - A * beta[n - 2]) / (C + A * alpha[n - 2]);

    // Нахождение решения - обратный ход метода прогонки
    for (std::size_t i = n - 2; i > 0; --i)
        splines[i].c = alpha[i] * splines[i + 1].c + beta[i];

    // Освобождение памяти, занимаемой прогоночными коэффициентами
    delete[] beta;
    delete[] alpha;

    // По известным коэффициентам c[i] находим значения b[i] и d[i]
    for (std::size_t i = n - 1; i > 0; --i)
    {
        double h_i = points[i].x[0] - points[i - 1].x[0];
        splines[i].d = (splines[i].c - splines[i - 1].c) / h_i;
        splines[i].b = h_i * (2. * splines[i].c + splines[i - 1].c) / 6. + (points[i].x[1] - points[i - 1].x[1]) / h_i;
    }
}
double cubic_spline::f(double x) const
{
    if (!splines)
        //return std::numeric_limits<double>::quiet_NaN(); // Если сплайны ещё не построены - возвращаем NaN
        return nan("");

    spline_tuple *s;
    if (x <= splines[0].x) // Если x меньше точки сетки x[0] - пользуемся первым эл-тов массива
        s = splines + 1;
    else if (x >= splines[n - 1].x) // Если x больше точки сетки x[n - 1] - пользуемся последним эл-том массива
        s = splines + n - 1;
    else // Иначе x лежит между граничными точками сетки - производим бинарный поиск нужного эл-та массива
    {
        std::size_t i = 0, j = n - 1;
        while (i + 1 < j)
        {
            std::size_t k = i + (j - i) / 2;
            if (x <= splines[k].x)
                j = k;
            else
                i = k;
        }
        s = splines + j;
    }

    double dx = (x - s->x);
    return s->a + (s->b + (s->c / 2. + s->d * dx / 6.) * dx) * dx; // Вычисляем значение сплайна в заданной точке.
}
void cubic_spline::free_mem()
{
    if(splines != nullptr)
        delete[] splines;
    splines = nullptr;
}

void Function::setBezie(double E1, Elementary::POINT3 newp1, Elementary::POINT3 newp2, Elementary::POINT3 newp3, double E2)
{
    this->E1 = E1;
    this->E2 = E2;
    p1 = newp1;
    p2 = newp2;
    p3 = newp3;
    mode = myFunType::Bezie;
}
void Function::setDifBezie(double E1, Elementary::POINT3 newp1, Elementary::POINT3 newp2, Elementary::POINT3 newp3, double E2)
{
    setBezie(E1, newp1, newp2, newp3, E2);
    mode = myFunType::difBezie;
}
void Function::setPlasticSigma(Elementary::POINT3 newp1, double E1, double E2)
{
    this->E1 = E1;
    this->E2 = E2;
    p1 = newp1;
    mode = myFunType::PlasticSigma;
}
void Function::setDifPlasticSigma(Elementary::POINT3 newp1, double E1, double E2)
{
    setPlasticSigma(newp1, E1, E2);
    mode = myFunType::difPlasticSigma;
}
void Function::setPlasticSigma_Yeld(double elasticSigmaLimit, double E1, double E2)
{
    this->E1 = E1;
    this->E2 = E2;
    this->elasticSigmaLimit_pl = elasticSigmaLimit;
    mode = myFunType::PlasticSigma_yeld;
}
void Function::setDiffPlasticSigma_Yeld(double elasticSigmaLimit, double E1, double E2)
{
    this->E1 = E1;
    this->E2 = E2;
    this->elasticSigmaLimit_pl = elasticSigmaLimit;
    mode = myFunType::diffPlasticSigma_yeld;
}

void Function::setPlasticSigma_Yeld_hardening(double elasticSigmaLimit, double E, double n)
{
    double q_max = 1;//0.1; // сплайн строится для аргумента 0 <= q <= q_max
    double N = 10000;//1000; // количество узлов сплайна минус 1
    Grid::GridRectangleRegular1D it_grid;
    //it_grid.init(0, q_max, N/2);// для эрмита
    it_grid.init(0, q_max, N);// для лагранжа
    double CoefAlpha = 0.0000000000001;// для Лагранжа не важно
    it.init(it_grid, CoefAlpha);    // аргументы игнорируются для нерегулярной сетки

    // построение интерполянта по значениям ф-и в узлах
    for(size_t j = 0; j < N; j++)
    {
        //double q = (j) / ((N-1)) * q_max;
        double q = (j*j) / ((N-1)*(N-1)) * q_max;
        double sigma_y = elasticSigmaLimit;
        double eps_y = sigma_y/E;
        // поиск eps из уравнения [eps - 1./E*F(eps) = q], F(eps) = sigma_y*(eps/eps_y)^n, eps >= eps_y
        // решается уравнение f(eps) = q
        // функция f(eps) = eps - 1./E*F(eps) возрастает
        double eps;
        double eps1 = eps_y;
        double eps2 = eps_y + 10;
        for(int i = 0; i < 100; i++)
        {
            eps = (eps1 + eps2) / 2.;
            double f_trial = eps - 1./E*sigma_y*pow(eps/eps_y, n);
            if(f_trial < q)
            {
                eps1 = eps;
            }
            else
            {
                eps2 = eps;
            }
        }
        double sigma = sigma_y*pow(eps/eps_y, n);
        it.addPoint(q, sigma);
    }
    it.buildInterpolant();
    mode = myFunType::PlasticSigma_yeld_hardening;

    /*
    double q_max = 1;//0.1; // сплайн строится для аргумента 0 <= q <= q_max
    double N = 10000;//1000; // количество узлов сплайна минус 1
    Grid::GridRectangleRegular1D it_grid;
    it_grid.init(0, q_max, N);
    double CoefAlpha = 0.0000000001;
    it.init(it_grid, CoefAlpha);
    // получение набора узлов интерполянта
    std::vector<Elementary::POINT1> coordinates;
    it.getNodesCoordinates(coordinates);
    // построение интерполянта по значениям ф-и в узлах
    Elementary::Vector nodeValue;
    nodeValue.resize(coordinates.size());
    for(size_t j = 0; j < coordinates.size(); j++)
    {
        double q = coordinates[j];
        double sigma_y = elasticSigmaLimit;
        double eps_y = sigma_y/E;
        // поиск eps из уравнения [eps - 1./E*F(eps) = q], F(eps) = sigma_y*(eps/eps_y)^n, eps >= eps_y
        // решается уравнение f(eps) = q
        // функция f(eps) = eps - 1./E*F(eps) возрастает
        double eps;
        double eps1 = eps_y;
        double eps2 = eps_y + 10;
        for(int i = 0; i < 100; i++)
        {
            eps = (eps1 + eps2) / 2.;
            double f_trial = eps - 1./E*sigma_y*pow(eps/eps_y, n);
            if(f_trial < q)
            {
                eps1 = eps;
            }
            else
            {
                eps2 = eps;
            }
        }
        double sigma = sigma_y*pow(eps/eps_y, n);
        nodeValue[j] = sigma;
    }
    it.buildInterpolantByAllNodes(nodeValue);
    mode = myFunType::PlasticSigma_yeld_hardening;
    */
}

void Function::setDiffPlasticSigma_Yeld_hardening(double elasticSigmaLimit, double E, double n)
{
    setPlasticSigma_Yeld_hardening(elasticSigmaLimit, E, n);
    mode = myFunType::diffPlasticSigma_yeld_hardening;
}
void Function::setPlasticEps(Elementary::POINT3 newp1, double E1, double E2)
{
    setPlasticSigma(newp1, E1, E2);
    mode = myFunType::PlasticEps;
}
void Function::setDifPlasticEps(Elementary::POINT3 newp1, double E1, double E2)
{
    setPlasticSigma(newp1, E1, E2);
    mode = myFunType::difPlasticEps;
}


void Function::__setCreep(const double A, const double n, const double m, const double tan0)
{
    creep_A = A;
    creep_n = n;
    creep_m = m;
    creep_tan0 = tan0;
}
void Function::setCreepSigma(const double A, const double n, const double m, const double tan0)
{
    __setCreep(A, n, m, tan0);
    mode = myFunType::CreepSigma;
}
void Function::setDifCreepSigma(const double A, const double n, const double m, const double tan0)
{
    __setCreep(A, n, m, tan0);
    mode = myFunType::difCreepSigma;
}
void Function::setCreepEps(const double A, const double n, const double m, const double tan0)
{
    __setCreep(A, n, m, tan0);
    mode = myFunType::CreepEps;
}
void Function::setDifCreepEps(const double A, const double n, const double m, const double tan0)
{
    __setCreep(A, n, m, tan0);
    mode = myFunType::difCreepEps;
}
void Function::setSpline(const std::vector<Elementary::VECTOR2> &points)
{
    spline.build_spline(points);
    mode = myFunType::Spline;
}
double Function::calcValue(const double x)
{
    if(mode == myFunType::Expression)
    {
        setArgumentValue(0, x);
        return calcValue();
    }
    if(mode == myFunType::Spline)
    {
        //return calcInverse(x);
        return spline.f(x);
    }
    if(mode == myFunType::Bezie)
    {
        const double &sigmaEqv = x;
        const double &elasticSigmaLimit = p2[0];
        // чуть до текучести
        if(sigmaEqv <= p1[0])
            return sigmaEqv*E1;
        // чуть после текучести
        if(sigmaEqv >= p3[0])
            return elasticSigmaLimit*E1 + (sigmaEqv-elasticSigmaLimit)*E2;
        // начало текучести
        //if(sigmaEqv > p1[0] && sigmaEqv < p3[0])
        // функция растёт с ростом t
        double t1 = 0;
        double t2 = 1;
        double t;
        double tLast = 1.e100;
        Elementary::POINT3 b;
        for(;;)
        {
            t = (t1 + t2) / 2;
            //b = SQR(1 - t)*p1 + 2*t*(1 - t)*p2 + SQR(t)*p3;
            b = SQR(1 - t)*p1 + 2*t*(1 - t)*p2 + SQR(t)*p3;
            double d = b[0] - x;
            if(d > 0)
                t2 = t;
            else
                t1 = t;
            if(t == tLast)
                break;
            else
                tLast = t;
        }
        return b[1];
        /*
        if(x <= p1[0])
        {
            double x1 = p1[0];
            double y1 = p1[1];
            double x2 = p2[0];
            double y2 = p2[1];
            return y1 + (x - x1)/(x2-x1)*(y2-y1);
        }else
        if(x >= p3[0])
        {
            double x1 = p2[0];
            double y1 = p2[1];
            double x2 = p3[0];
            double y2 = p3[1];
            return y1 + (x - x1)/(x2-x1)*(y2-y1);
        }
        else
        {
            // функция растёт с ростом t
            double t1 = 0;
            double t2 = 1;
            double t;
            double tLast = 1.e100;
            Elementary::POINT3 b;
            for(;;)
            {
                t = (t1 + t2) / 2;
                b = SQR(1 - t)*p1 + 2*t*(1 - t)*p2 + SQR(t)*p3;
                double d = b[0] - x;
                if(d > 0)
                    t2 = t;
                else
                    t1 = t;
                //if(abs(t2 - t1) < 1.e-15) break;
                if(t == tLast)
                    break;
                else
                    tLast = t;
            }
            return b[1];
        }*/
    }
    if(mode == myFunType::difBezie)
    {
        const double &sigmaEqv = x;
        // чуть до текучести
        if(sigmaEqv <= p1[0])
            return E1;
        // чуть после текучести
        if(sigmaEqv >= p3[0])
            return E2;
        return calcDif(x);// на переходном участке производная находится численно
        // начало текучести
        //if(sigmaEqv > p1[0] && sigmaEqv < p3[0])
        //return calcDif(sigmaEqv);
    }
    if(mode == myFunType::PlasticSigma || mode == myFunType::PlasticEps)
    {
        const double &epsEqv = x;
        const double &elasticEpsLimit = p1[0];
        // до текучести
        if(epsEqv <= elasticEpsLimit)
            return epsEqv*E1;
        // чуть после текучести
        if(epsEqv > elasticEpsLimit)
            return elasticEpsLimit*E1 + (epsEqv-elasticEpsLimit)*E2;
        return 0;
    }
    if(mode == myFunType::difPlasticSigma || mode == myFunType::difPlasticEps)
    {
        const double &epsEqv = x;
        const double &elasticEpsLimit = p1[0];
        // до текучести
        if(epsEqv <= elasticEpsLimit)
        {
            if(mode == myFunType::difPlasticSigma)
                return E1;
            if(mode == myFunType::difPlasticEps)
                return 1/E1;
        }
        // чуть после текучести
        if(epsEqv > elasticEpsLimit)
        {
            if(mode == myFunType::difPlasticSigma)
                return E2;
            if(mode == myFunType::difPlasticEps)
                return 1/E2;
        }
        return E1;
    }
    if(mode == myFunType::PlasticSigma_yeld)
    {
        return elasticSigmaLimit_pl + E1*E2/(E1 + E2)*x;
    }
    if(mode == myFunType::diffPlasticSigma_yeld)
    {
        return E1*E2/(E1 + E2);
    }
    if(mode == myFunType::PlasticSigma_yeld_hardening)
    {
        return it.fun(x);
    }
    if(mode == myFunType::diffPlasticSigma_yeld_hardening)
    {
        //return it.difFun(x, 1);
        return it.difFun(x, 1);
    }

    return 0;
}
double Function::calcValue(const double x, const double y)
{
    if(mode == myFunType::Expression)
    {
        setArgumentValue(0, x);
        setArgumentValue(1, y);
        return calcValue();
    }
    //x = epsEqv, y = t
    if(mode == myFunType::CreepSigma)
    {
        double f = pow(x/(creep_A*y), 1./creep_n);
        double lin_f = creep_tan0*x;
        if(lin_f < f)
            return lin_f;
        else
            return f;
    }
    if(mode == myFunType::difCreepSigma)
    {
        double df = pow(x/(creep_A*y), 1./creep_n) / (creep_n * x);
        double dlin_f = creep_tan0;

        double f = pow(x/(creep_A*y), 1./creep_n);
        double lin_f = creep_tan0*x;
        if(lin_f < f || y == 0 || x == 0)
            return dlin_f;
        else
            return df;
    }
    //x = sigmaEqv, y = t
    if(mode == myFunType::CreepEps)
    {
        if(y < 0)
            return 0;
        return creep_A*pow(x, creep_n)*pow(y, creep_m);
    }
    if(mode == myFunType::difCreepEps)
    {
        if(y < 0)
            return creep_tan0;
        return creep_n*creep_A*pow(x, creep_n - 1)*pow(y, creep_m);
    }
}
double Function::calcDif(const double x)
{
    // приблизительное вычисление производной
    double x1;
    double x2;
    if(x == 0)
    {
        x1 = -0.00001;
        x2 = 0.00001;
    }
    else
    {
        x1 = 0.99999*x;
        x2 = 1.00001*x;
    }
    double df = (calcValue(x2) - calcValue(x1)) / (x2 - x1);
    return df;
}
double Function::calcInverse(const double y)
{
    // поиск x:f(x)=y
    // для возрастающей ф-и с аргументом от 0 до 1
    double x1 = 0;
    double x2 = 0.01;
    double x;
    for(;;)
    {
        x = (x1 + x2) / 2;
        double y0;
        /*if(mode == 0)
        {
            setArgumentValue(0, x);
            y0 = calcValue();
        }
        if(mode == 1)
        {
            y0 = spline.f(x);
        }*/
        y0 = spline.f(x);
        double d = y0 - y;
        if(fabs(d / (y  + (double)(y==0))) < 1.e-15) break;
        if(d > 0)
            x2 = x;
        else
            x1 = x;
        //printf("x = %le y0 = %le y = %le\n", x, y0, y);
    }
    return (x1 + x2) / 2;
}













// приоритеты арифметических операций
// порядок соответствует OperatorType
const int OperatorPriority[5] =
{
    100,            // ^
    99,             // *
    99,             // /
    98,             // +
    98,             // -
};

// символы арифметических операций
// порядок соответствует OperatorType
const char OperatorSymbols[5] =
{
    '^',
    '*',
    '/',
    '+',
    '-',
};

// названия элементарных функций
// порядок соответствует ElementaryFunctionType
const char *ElementaryFunctionNames[_EFunType_Size] =
{
    "pow",        // 0
    "log",        // 1
    "exp",        // 2
    "ln",         // 3
    "log10",      // 4
    "sin",        // 5
    "cos",        // 6
    "tg",         // 7
    "arcsin",     // 8
    "arccos",     // 9
    "arctg",      // 10
    "abs",        // 11
    "sign",       // 12
};

// количества аргументов элементарных функций
// порядок соответствует ElementaryFunctionType
const int ElementaryFunctionNumArgs[_EFunType_Size] =
{
    2,
    2,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
};

// инициализация
Function::Function(const std::string expression,
                   const std::string functionName)
{
    setName(functionName);
    setExpression(expression);
}

// освобождение
Function::~Function()
{
}

// изменение имени функции
void Function::setName(const std::string functionName)
{
    name = functionName;
}

// изменение текстового выражения функции
void Function::setExpression(const std::string expression)
{
    s = expression;
}

// текущее имя функции
const std::string Function::getName() const
{
    return name;
}

// текущее текстового выражения функции
const std::string Function::getExpression() const
{
    return s;
}

// добление аргумента. возвращает id добавленного элемента.
int Function::addArgument(const std::string name, const double value)
{
    Argument arg;
    arg.name = name;
    arg.value = value;
    args.push_back(arg);
    return (int)args.size() - 1;
}

// добавление функции
void Function::addFunction(Function *fun)
{
    funs.push_back(fun);
}

// изменение значения аргумента по идентификатору
void Function::setArgumentValue(const int id, const double value)
{
    args[id].value = value;
}

// построение постфиксной записи функции, предоставленной в виде строки
void Function::parse()
{
    mode = myFunType::Expression;
    Token *tt;
    Token t;
    translate();
    ts.clear();
    RPN.clear();
    //printf("%s\n", s.c_str());
    for(size_t i = 0; i < tokens.size(); i++)
    {
        t = tokens[i];
        //printf("token = %d, %d\n", t.type, t.id);
        //printToken(t);
        switch (t.type)
        {
        case TokenSeparator:
            if(t.id == (size_t)SeparatorType::Comma)
                pushWhileNotOpen();
            if(t.id == (size_t)SeparatorType::Space){}
            break;
        case TokenNumber:
        case TokenVariable:
            RPN.push_back(t);
            break;
        case TokenFunction:
        case TokenElementaryFunction:
            ts.push_back(t);
            if(i + 1 == tokens.size() || tokens[i + 1].type != TokenOpen)
                throw "После названия функции необходима открывающая скобка";
            break;
        case TokenOpen:
            ts.push_back(t);
            break;
        case TokenClose:
            tt = pushWhileNotOpen();   // tt - токен открывающей скобки
            if(tt->id != t.id)
                throw "Скобки не согласованы(2)";
            ts.pop_back();
            if(!ts.empty())
            {
                tt = &ts.back();
                if(tt->type == TokenFunction || tt->type == TokenElementaryFunction)
                {
                    RPN.push_back(*tt);
                    ts.pop_back();
                }
            }
            break;
        case TokenOperator:
            for(;;)
            {
                if(ts.empty())
                    break;
                tt = &ts.back();
                if(tt->type == TokenOperator && OperatorPriority[t.id] <= OperatorPriority[tt->id])
                {
                    RPN.push_back(*tt);
                    ts.pop_back();
                }
                else
                    break;
            }
            ts.push_back(t);
            break;
        default:
            throw "Не допустимый тип токена";
            break;
        }
    }   // for(int i = 0; i < tokens.getCount(); i++)
    for(;;)
    {
        if(ts.empty())
            break;
        tt = &ts.back();
        if(!(tt->type == TokenOperator ||
             tt->type == TokenElementaryFunction ||
             tt->type == TokenFunction
             ))
            throw "Скобки не согласованы(3)";
        RPN.push_back(*tt);
        ts.pop_back();
    }
    // отладочный вывод RPN
    for(size_t i = 0; i < RPN.size(); i++)
    {
        t = RPN[i];
        //printToken(t);
        //printf("token = %d, %d\n", t.type, t.id);
    }
}

// вычисление значения функции
double Function::calcValue()
{
    //if(RPN.getCount() == 0) return 0;
    Token t;
    double value;
    double val[10];
    ns.clear();
    for(size_t i = 0; i < RPN.size(); i++)
    {
        t = RPN[i];
        switch (t.type)
        {
        case TokenNumber:
            ns.push_back(nums[t.id]);
            break;
        case TokenVariable:
            ns.push_back(args[t.id].value);
            break;
        case TokenFunction:
            for(size_t i = 0; i < 1; i++)   //##funs[t.id]->args.size()
            {
                if(ns.empty())
                    throw "Не достаточно аргументов для функции";
                value = ns.back();
                funs[t.id]->setArgumentValue((int)i, value);
                ns.pop_back();
            }
            value = funs[t.id]->calcValue();
            ns.push_back(value);
            break;
        case TokenElementaryFunction:
            for(int i = 0; i < ElementaryFunctionNumArgs[t.id]; i++)
            {
                if(ns.empty())
                    throw "Не достаточно аргументов для элементарной функции";
                val[i] = ns.back();
                ns.pop_back();
            }
            switch (t.id)
            {
            case EFunPow:
                value = pow(val[1], val[0]);
                break;
            case EFunLog:
                value = log(val[0]) / log(val[1]);
                break;
            case EFunExp:
                value = exp(val[0]);
                break;
            case EFunLn:
                value = log(val[0]);
                break;
            case EFunLog10:
                value = log10(val[0]);
                break;
            case EFunSin:
                value = sin(val[0]);
                break;
            case EFunCos:
                value = cos(val[0]);
                break;
            case EFunTg:
                value = tan(val[0]);
                break;
            case EFunArcsin:
                value = asin(val[0]);
                break;
            case EFunArccos:
                value = acos(val[0]);
                break;
            case EFunArctg:
                value = atan(val[0]);
                break;
            case EFunAbs:
                value = fabs(val[0]);
                break;
            case EFunSign:
                value = val[0]>0?1:-1;
                break;
            default:
                break;
            }
            ns.push_back(value);
            break;
        case TokenOperator:
            for(int i = 0; i < 2; i++)
            {
                if(ns.empty())
                    throw "Не достаточно аргументов для оператора";
                val[i] = ns.back();
                ns.pop_back();
            }
            switch (t.id)
            {
            case OperatorPow:
                value = pow(val[1], val[0]);
                break;
            case OperatorMul:
                value = val[1] * val[0];
                break;
            case OperatorDiv:
                value = val[1] / val[0];
                break;
            case OperatorAdd:
                value = val[1] + val[0];
                break;
            case OperatorSub:
                value = val[1] - val[0];
                break;
            }
            ns.push_back(value);
            break;
        default:
            throw "Не допустимый тип токена";
            break;
        }
    }
    if(ns.empty())
        throw "Пустая функция";
    value = ns.back();
    ns.pop_back();
    if(!ns.empty())
        throw "Возможно, слишком много аргументов для функции";
    return value;
}

//#вывод токена на экран
void Function::printToken(const Token &t)
{
    switch (t.type)
    {
    case TokenNumber:
        printf("%.1le\n", nums[t.id]);
        break;
    case TokenVariable:
        printf("%.1le\n", args[t.id].value);
        break;
    case TokenFunction:
        printf("%s\n", funs[t.id]->name.c_str());
        break;
    case TokenElementaryFunction:
        printf("%s\n", ElementaryFunctionNames[t.id]);
        break;
    case TokenOperator:
        printf("%c\n", OperatorSymbols[t.id]);
        break;
    default:
        //throw "Не допустимый тип токена";
        break;
    }
}

/*
// операнд может начинаться этим токеном?
bool Function::FirstTokenOfOperand(const Token &t)
{
    return
 char *Function::GetNamese ||
            t.type == TokenNumber ||
            t.type == TokenVariable;
}

// операнд может заканчиваться этим токеном?
bool Function::LastTokenOfOperand(const Token &t)
{
    return
            t.type == TokenOpen ||
            t.type == TokenNumber ||
            t.type == TokenVariable ||
            t.type == TokenElementaryFunction ||
            t.type == TokenFunction;
}*/

// выводить все токены из стека ts пока не встретится открывающая скобка
Token* Function::pushWhileNotOpen()
{
    Token *tt;
    for(;;)
    {
        if(ts.empty())
            throw "Скобки не согласованы(1)";
        tt = &ts.back();
        if(tt->type == TokenOpen)
            return tt;
        else
        {
            RPN.push_back(*tt);
            ts.pop_back();
        }
    }
}

// преобразование функции из строки в последовательность токенов
void Function::translate()
{
    tokens.clear();

    char c;
    for(size_t i = 0; i < s.length();)
    {
        Token t;
        c = s[i];
        switch (c)
        {
        //#лучше табличку для каждого символа
        // скобки
        case '(': t = Token{TokenOpen, Bracket1}; break;
        case '[': t = Token{TokenOpen, Bracket2}; break;
        case '{': t = Token{TokenOpen, Bracket3}; break;
        case '<': t = Token{TokenOpen, Bracket4}; break;
        case ')': t = Token{TokenClose, Bracket1}; break;
        case ']': t = Token{TokenClose, Bracket2}; break;
        case '}': t = Token{TokenClose, Bracket3}; break;
        case '>': t = Token{TokenClose, Bracket4}; break;
        // арифметические операторы
        case '^': t = Token{TokenOperator, OperatorPow}; break;
        case '*': t = Token{TokenOperator, OperatorMul}; break;
        case '/': t = Token{TokenOperator, OperatorDiv}; break;
        case '+': t = Token{TokenOperator, OperatorAdd}; break;
        case '-': t = Token{TokenOperator, OperatorSub}; break;
        case ',': t = Token{TokenSeparator, (int)SeparatorType::Comma}; break;
        case ' ': t = Token{TokenSeparator, (int)SeparatorType::Space}; break;
        // число, переменная или функция
        default:
        {
            // число
            if (isNumber(c) || c == '.')
            {
                int j = (int)i;
                int la = 0, lb = 0, lc = 0;
                int sign_multiplier = 1;
                double a = 0;
                double b = 0;
                int c = 0;
                // (num^l1)[.(num^l2)][(e|E)[+|-](num^l3)]
                for(; isNumber(s[j]); j++)
                {
                    la++;
                    a *= 10;
                    a += charToNumber(s[j]);
                }
                if(s[j] == '.')
                {
                    j++;  // .
                }
                for(double factor = 1; isNumber(s[j]); j++)
                {
                    lb++;
                    factor /= 10;
                    b += factor * charToNumber(s[j]);
                }
                if(la + lb == 0)
                {
                    throw "Перед точкой или после точки должна быть хотя бы одна цифра";
                }
                if(s[j] == 'e' || s[j] == 'E')
                {
                    j++;    // e|E
                    if(s[j] == '+')
                    {
                        sign_multiplier = +1;
                        j++;    // +
                    }else
                    if(s[j] == '-')
                    {
                        sign_multiplier = -1;
                        j++;    // -
                    }
                    for(; isNumber(s[j]); j++)
                    {
                        lc++;
                        c *= 10;
                        c += charToNumber(s[j]);
                    }
                    if(lc == 0)
                    {
                        throw "После символа e или E в записи числа должна быть хотя бы одна цифра";
                    }
                    c *= sign_multiplier;
                }
                double x = (a + b) * pow(10, c);
                nums.push_back(x);
                t = Token{TokenNumber, nums.size() - 1};
                i = j; // число прочитано, смещаемся дальше
                //printf("%.16le\n", x);
                goto switch_exit_without_inc;   // токен готов
            }
            // переменная или функция
            if(isLetter(c))
            {
                size_t len;
                // сравнение с названиями элементарных функций
                for(size_t j = 0; j < _EFunType_Size; j++)
                {
                    len = strlen(ElementaryFunctionNames[j]);
                    if(firstWordEqual(&s[i], ElementaryFunctionNames[j]))
                    {
                        t = Token{TokenElementaryFunction, j};
                        i += len;  // имя функции прочитано, смещаемся дальше
                        goto switch_exit_without_inc;   // токен готов
                    }
                }
                // сравнение с названиями переменных
                for(size_t j = 0; j < args.size(); j++)
                {
                    len = args[j].name.length();
                    if(firstWordEqual(&s[i], args[j].name))
                    {
                        t = Token{TokenVariable, j};
                        i += len;   // имя переменной прочитано, смещаемся дальше
                        goto switch_exit_without_inc;   // токен готов
                    }
                }
                // сравнение с названиями функций
                for(size_t j = 0; j < funs.size(); j++)
                {
                    len = funs[j]->name.length();
                    if(firstWordEqual(&s[i], funs[j]->name) && !isLetter(s[i+len]))
                    {
                        t = Token{TokenFunction, j};
                        i += len;  // имя функции прочитано, смещаемся дальше
                        goto switch_exit_without_inc;   // токен готов
                    }
                }
                //printf("s1=%s, s2=%s\n", &s[i], ElementaryFunctionNames[1]);
                throw "Переменная или функция не определена";
            }
            throw "Не допустимый символ";
        }
            break;
        }   // switch (c)
        i++;    // одиночный символ прочитан, смещаемся дальше
switch_exit_without_inc:
        //printf("token[%d] = %d, %d\n", tokens.getCount(), t.type, t.id);
        tokens.push_back(t);     // добавление токена
    }
}

// c является буквой?
bool Function::isLetter(const char c)
{
    return (c >= 'a' && c <= 'z') || ((c >= 'A' && c <= 'Z') || c == '_');
}

// c является цифрой?
bool Function::isNumber(const char c)
{
    return (c >= '0' && c <= '9');
}

// преобразование символа-цифры в целое число
int Function::charToNumber(const char c)
{
    return c - '0';
}

// строка str начинается со слова word?
bool Function::firstWordEqual(const std::string str, const std::string word)
{
    for(size_t i = 0; ; i++)
    {
        if(i >= word.length())
            return (!isLetter(str[i])) && (!isNumber(str[i]));
        if(str[i] != word[i] || i >= str.length()) return false;
    }
}

double strToDouble(const std::string expression)
{
    double value;
    Function tfun(expression, "");
    try
    {
        tfun.parse();
        value = tfun.calcValue();
    }
    catch(const char *err)
    {
        printf("%s\n", err);
        value = 0;
    }
    return value;
}

bool strToDoubleTest(const std::string expression)
{
    Function tfun(expression, "");
    try
    {
        tfun.parse();
        tfun.calcValue();
        return true;
    }
    catch(const char *err)
    {
        printf("strToDoubleTest: %s\n", err);
        return false;
    }
}


}   // namespace Fun_parser
