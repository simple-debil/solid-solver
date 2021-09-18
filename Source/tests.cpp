#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"
#include <fstream>

#include "tests.h"
#include "console.h"
#include "interpolation.h"

using namespace Elementary;
using namespace Solid;
using namespace Grid;

namespace Tests
{
void GraphInf::makeGraph() const
{
    fprintf(stderr, "%s...\n", fn_pic_short.c_str());
    OS::Gnuplot gnuplot;
    gnuplot.exec(fn_graph.c_str());
}

void GraphInfs::add_all_in_subdir(const std::string dir, const std::string fn_gnuplot, const std::string fn_picture)
{
    push_back({});
    back().fn_graph = dir + fn_gnuplot;
    back().fn_pic = dir + fn_picture;
    back().fn_pic_short = fn_picture;
}
void GraphInfs::add_gnuplot_in_root(const std::string dir, const std::string fn_gnuplot, const std::string fn_picture)
{
    push_back({});
    back().fn_graph = fn_gnuplot;
    back().fn_pic = dir + fn_picture;
    back().fn_pic_short = fn_picture;
}
void GraphInfs::read_gnuplot_from_file(const std::string fn_source)
{

    clear();
    std::ifstream f(fn_source);
    std::string s;
    std::getline(f, s);
    if(s.length() == 0)
        return;
    generalGraphsNomber = std::stoi(s);
    std::getline(f, s);
    graphsPerGlobalStepNomber = std::stoi(s);
    for(int i = 0; i < generalGraphsNomber; i++)
    {
        push_back({});
        back().globalStepNumber = -1;
        std::getline(f, back().fn_graph);
        std::getline(f, back().fn_pic);
        std::getline(f, back().fn_pic_short);
    }
    for(;;)
    {
        s.clear();
        std::getline(f, s);
        if(s.length() == 0)
            break;
        //fprintf(stderr, "%s\n", s.c_str());
        int globalStepNumber = std::stoi(s);
        for(int i = 0; i < graphsPerGlobalStepNomber; i++)
        {
            push_back({});
            back().globalStepNumber = globalStepNumber;
            std::getline(f, back().fn_graph);
            std::getline(f, back().fn_pic);
            std::getline(f, back().fn_pic_short);
        }
    }
    f.close();
}
void GraphInfs::makeAllGraphs() const
{
    OS::Gnuplot gnuplot;
    for(size_t ind = 0; ind < size(); ind++)
        (*this)[ind].makeGraph();
}
void GraphInfs::makeglobalStepNumberGraphs(const int globalStepNumber) const
{
    int i1 = generalGraphsNomber + graphsPerGlobalStepNomber*globalStepNumber;
    int i2 = generalGraphsNomber + graphsPerGlobalStepNomber*(globalStepNumber + 1);
    for(int i = i1; i < i2; i++)
        (*this)[i].makeGraph();
}
void GraphInfs::makeGeneralGraphs() const
{
    for(int i = 0; i < generalGraphsNomber; i++)
        (*this)[i].makeGraph();
}

Test_base *Test_base::gen(const Test_base::Type type)
{
    switch (type)
    {
    case Type::Sphere_T:
        return new Test_sphere_T;
        break;
    case Type::Sphere_ep:
        return new Test_sphere_ep;
        break;
    case Type::Sphere_creep:
        return new Test_sphere_creep;
        break;
    case Type::Sphere_Te:
        return new Test_sphere_Te;
        break;
    case Type::Formovka:
        return new Test_formovka;
        break;
    case Type::ContactCylinder:
        return new Test_contactCylinder;
        break;
    case Type::ContactSphere:
        return new Test_contactSphere;
        break;
    case Type::Creep:
        return new Test_creep;
        break;
    case Type::Kirsch:
        return new Test_Kirsch;
        break;
    case Type::Crack_plate_gen:
        return new Test_crack_plate_gen;
        break;
    //default:
    //    return nullptr;
    }
}
void Test_base::initStepsGraphs()
{
    fn_steps_all_iterations = dir + "steps.txt";
    fn_visible_steps_all_iterations = dir + "visible_steps.txt";
    fn_steps_results_of_iterations = dir + "steps_res.txt";
    fn_log = dir + "log.txt";
    // полные данные о всех итерациях
    stepsGraph_all_iterations.clear();
    stepsGraph_all_iterations.add_gnuplot_in_root(dir, "steps_slau.gnu", "steps_slau.png");
    stepsGraph_all_iterations.add_gnuplot_in_root(dir, "steps_ep.gnu", "steps_ep.png");
    stepsGraph_all_iterations.add_gnuplot_in_root(dir, "steps_contact.gnu", "steps_contact.png");
    //stepsGraph_all_iterations.add_gnuplot_in_root(dir, "steps_creep_IE.gnu", "steps_creep_IE.png");
    //stepsGraph_all_iterations.add_gnuplot_in_root(dir, "steps_creep_IS.gnu", "steps_creep_IS.png");
    // данные о результатах итераций
    stepsGraph_results_of_iterations.clear();
    stepsGraph_results_of_iterations.add_gnuplot_in_root(dir, "steps_res_iterations.gnu", "steps_res_iterations.png");
    stepsGraph_results_of_iterations.add_gnuplot_in_root(dir, "steps_res_residuals.gnu", "steps_res_residuals.png");
}
void Test_base::makeGraphsAfterGlobalStep(const int globalStepNumber, const Task &task)
{
    // построение графиков для текущего глобального шага
    resGraph.makeglobalStepNumberGraphs(globalStepNumber);
    // если шаг последний, то построение графиков, общих для всех глобальных шагов
    if(globalStepNumber == (int)task.step->size() - 1)
        resGraph.makeGeneralGraphs();
}
void Test_base::makeResultGraphs()
{
    resGraph.makeAllGraphs();
}
const GraphInfs *Test_base::get_oldResultGraphs() const
{
    return &resGraph;
}
const std::string *Test_base::get_logFileName()
{
    return &fn_log;
}
void Test_base::saveLastStepInf(const OutDataForLogger &stepsInf)
{
    Solid::MechIterInf el = (*stepsInf.stepInf)[stepsInf.stepNumber].iterInf[stepsInf.iterNumber];
    if(stepsInf.iterNumber == 0)
    {
        el.contact.max_deltaF_residual = -1;
        el.contact.max_endPoint_residual = -1;
        el.plastic.maxPlasticResidual.set_undefined();
    }
    if(stepsGraph_all_iterations.size() != 0)
    {
        FILE *f_steps;
        if(stepsInf.stepNumber == 0 && stepsInf.iterNumber == 0)
            f_steps = fopen(fn_steps_all_iterations.c_str(), "w");
        else
            f_steps = fopen(fn_steps_all_iterations.c_str(), "a");
        if(f_steps != nullptr)
        {
            fprintf(f_steps, "%d\t%d\t",
                    stepsInf.stepNumber,
                    stepsInf.iterNumber);
            el.save(f_steps);
            fprintf(f_steps, "\n");
            fclose(f_steps);
        }
    }
}
void Test_base::loadStepsInf(std::vector<MechStepInf> &stepInf)
{
    stepInf.clear();
    FILE *f_steps = fopen(fn_steps_all_iterations.c_str(), "r");
    if(f_steps == nullptr)
        return;
    int maxIterNumber = -1;
    int prevStepNumber = -1;
    for(;;)
    {
        int stepNumber;
        int iterNumber;
        fscanf(f_steps, "%d%d",
               &stepNumber,
               &iterNumber);
        if(feof(f_steps))
            break;
        if(stepNumber > prevStepNumber)
        {
            // начался новый шаг
            stepInf.push_back({});
            maxIterNumber = -1;
            prevStepNumber = stepNumber;
        }
        stepInf.back().iterInf.push_back({});
        Solid::MechIterInf &el = stepInf.back().iterInf.back();
        el.load(f_steps);
    }
    fclose(f_steps);
}
int Test_base::get_stepsInfLastSavedStepNomber()
{
    if(stepsGraph_all_iterations.size() == 0)
        return -1;
    std::vector<MechStepInf> stepInf;
    // чтение всех сохранённых невязок в stepInf
    loadStepsInf(stepInf);
    return stepInf.size() - 1;
}
void Test_base::makeStepsGraphs(const int stepIndex1, const int stepIndex2)
{
    std::vector<MechStepInf> stepInf;
    // чтение всех сохранённых невязок в stepInf
    loadStepsInf(stepInf);
    // запись данных для графиков (выводятся данные о всех итерациях)
    if(stepsGraph_all_iterations.size() != 0)
    {
        FILE *f_visible_steps = fopen(fn_visible_steps_all_iterations.c_str(), "w");
        if(f_visible_steps != nullptr)
        {
            for(int stepIndex = stepIndex1; stepIndex <= stepIndex2; stepIndex++)
            {
                int iterIndex_max = (int)stepInf[stepIndex].iterInf.size();
                for(int iterIndex = 0; iterIndex < iterIndex_max; iterIndex++)
                {
                    Solid::MechIterInf &el = stepInf[stepIndex].iterInf[iterIndex];
                    fprintf(f_visible_steps, "%le\t", (double)stepIndex + (double)iterIndex/iterIndex_max);
                    el.save(f_visible_steps);
                    fprintf(f_visible_steps, "\n");
                }
            }
            fclose(f_visible_steps);
            // построение графиков
            // пути к файлам передаются через командную строку в gnuplot
            OS::Console console;
            for(size_t ind = 0; ind < stepsGraph_all_iterations.size(); ind++)
            {
                std::string command;
                command = "gnuplot -e \"";
                // начало кода
                command += "f_in='";
                command += fn_visible_steps_all_iterations;
                command += "'; ";
                command += "f_out='";
                command += stepsGraph_all_iterations[ind].fn_pic;
                command += "';";
                // конец кода
                command +=  "\" ";
                command += stepsGraph_all_iterations[ind].fn_graph;
                console.exec(command.c_str());
            }
        }
    }
    // запись данных для графиков (выводятся данные о результатах итераций)
    if(stepsGraph_results_of_iterations.size() != 0)
    {
        FILE *f_steps_res = fopen(fn_steps_results_of_iterations.c_str(), "w");
        if(f_steps_res != nullptr)
        {
            for(int stepIndex = 0; stepIndex < (int)stepInf.size(); stepIndex++)
            {
                int iterIndex_max = (int)stepInf[stepIndex].iterInf.size();
                // подсчёт количества итераций, при которых менялась глобальная матрица  из-за разгрузок
                int GGlobalChangedCount_plastic = 1; // на 0-й итерации 1 матрица менялась 1 раз полностью
                int GGlobalChangedCount_contact = 0;
                for(int iterIndex = 0; iterIndex < iterIndex_max; iterIndex++)
                {
                    Solid::MechIterInf &el = stepInf[stepIndex].iterInf[iterIndex];
                    if(el.plastic.GLocalChangedNumber != 0)
                        GGlobalChangedCount_plastic++;
                    if(el.contact.contactChangedNumber != 0)
                        GGlobalChangedCount_contact++;
                }
                int iterIndex_last = iterIndex_max - 1;
                Solid::MechIterInf &el = stepInf[stepIndex].iterInf[iterIndex_last];
                fprintf(f_steps_res, "%d\t%d\t%d\t%d\t",
                        stepIndex,
                        iterIndex_last + 1,
                        GGlobalChangedCount_contact,
                        GGlobalChangedCount_plastic);
                el.save(f_steps_res);
                fprintf(f_steps_res, "\n");
            }
            fclose(f_steps_res);
        }
        // построение графиков
        // пути к файлам передаются через командную строку в gnuplot
        OS::Console console;
        for(size_t ind = 0; ind < stepsGraph_results_of_iterations.size(); ind++)
        {
            std::string command;
            command = "gnuplot -e \"";
            // начало кода
            command += "f_in='";
            command += fn_steps_results_of_iterations;
            command += "'; ";
            command += "f_out='";
            command += stepsGraph_results_of_iterations[ind].fn_pic;
            command += "';";
            // конец кода
            command +=  "\" ";
            command += stepsGraph_results_of_iterations[ind].fn_graph;
            console.exec(command.c_str());
        }
    }
}
void Test_base::get_oldStepsGraphs(GraphInfs const *&oldStepsGraph_all_iterations, GraphInfs const *&oldStepsGraph_results_of_iterations)
{
    oldStepsGraph_all_iterations = &stepsGraph_all_iterations;
    oldStepsGraph_results_of_iterations = &stepsGraph_results_of_iterations;
}
bool Test_base::possibleToShow2d() const
{
    return false;
}
bool Test_base::getContactSurfaceCircle(const Task &, const int, POINT2 &, double &) const
{
    return false;
}
void Test_base::notFixedCoordinates(int &cn1, int &cn2)
{
    cn1 = 0;
    cn2 = 1;
}
void Test_base::needToDrawFe(const Task &, const OutData &, const size_t, const int, std::array<bool, 6> &faceState) const
{
    faceState = {0, 0, 0, 0, 0, 0};
}

void Test_base::needToDrawSubFe(const Task &, const OutData &, const size_t, const size_t, const int, std::array<bool, 6> &faceState) const
{
    faceState = {0, 0, 0, 0, 0, 0};
}
bool Test_base::needToPaintFiniteElementSurface() const
{
    return false;
}



#define MOVE(v_ind, d0, d1, d2)	((v_ind + offset0*d0 + offset1*d1 + offset2*d2))

// формирует 1 из строковых параметров для помещения его в коммандную строку gnuplot
std::string genGnuplotCommandLineParameter(const std::string &var_name, const std::string &value)
{
    return var_name + "='" + value +"';";
}

// f = (t-A)*B + C: f(t1)=a1, f(t2)=a2
// A, B, C = ?
void calcCoef(const double t1, const double t2, const double a1, const double a2, double &A, double &B, double &C)
{
    if(a1 == a2)
    {
        A = 0;
        B = 0;
        C = a1;
    }
    else
    {
        A = (a1*t2 - a2*t1) / (a1 - a2);
        if(A != t1)
            B = a1 / (t1 - A);
        else
            B = a2 / (t2 - A);
        C = 0;
    }
}

// инициализация кривых для упруго-пластичного материала
void setLinearFunction(double P1, double P2, double t1, double t2, FunParser::Function &fun)
{
    // линейная функция
    char Px[1000];
    double PA, PB, PC;
    calcCoef(t1, t2, P1, P2, PA, PB, PC);
    if(PA >= 0)
    {
        if(PB >= 0)
        {
            if(PC >= 0)
                sprintf(Px, "(t-%.16le)*%.16le+(0+%.16le)\0", fabs(PA), fabs(PB), fabs(PC));
            else
                sprintf(Px, "(t-%.16le)*%.16le+(0-%.16le)\0", fabs(PA), fabs(PB), fabs(PC));
        }
        else
        {
            if(PC >= 0)
                sprintf(Px, "(t-%.16le)*(0-%.16le)+(0+%.16le)\0", fabs(PA), fabs(PB), fabs(PC));
            else
                sprintf(Px, "(t-%.16le)*(0-%.16le)+(0-%.16le)\0", fabs(PA), fabs(PB), fabs(PC));
        }
    }
    else
    {
        if(PB >= 0)
        {
            if(PC >= 0)
                sprintf(Px, "(t+%.16le)*%.16le+(0+%.16le)\0", fabs(PA), fabs(PB), fabs(PC));
            else
                sprintf(Px, "(t+%.16le)*%.16le+(0-%.16le)\0", fabs(PA), fabs(PB), fabs(PC));
        }
        else
        {
            if(PC >= 0)
                sprintf(Px, "(t+%.16le)*(0-%.16le)+(0+%.16le)\0", fabs(PA), fabs(PB), fabs(PC));
            else
                sprintf(Px, "(t+%.16le)*(0-%.16le)+(0-%.16le)\0", fabs(PA), fabs(PB), fabs(PC));
        }
    }
    (fun).setExpression(std::string(Px));
    (fun).args.clear();
    (fun).addArgument("t\0");
    (fun).parse();
    FILE *tt = fopen("______", "a");
    fprintf(tt, "%s\n", Px);
    //for(;;);
    //fprintf(tt, "%le _ %le\n", (*fun).solve(t1), (*fun).solve(t2));
    fclose(tt);
}
void setLinearBc2_el(double P1, double P2, double t1, double t2, Solid::MechBoundaryCondition2Source_base *&bc2Source_el)
{
    // линейная функция
    FunParser::Function *fun;
    fun = new FunParser::Function;
    setLinearFunction(P1, P2, t1, t2, *fun);
    MechBoundaryCondition2Source_ScalarFunction *mbc2Source_sf = new MechBoundaryCondition2Source_ScalarFunction;
    mbc2Source_sf->init(fun);
    bc2Source_el = mbc2Source_sf;
}
void setPlasticCurve(Solid::MechMaterialSource &mechMat, const double k_sigma_eps, const double k_eps_sigma)
{
    using namespace Solid;
    // скорректированный предел текучести. А вообще, лучше задавать кривые и обработку результата сразу подогнанные под коэффициенты Csigma, Ceps
    double elasticSigmaLimit = mechMat.elasticSigmaLimit*sqrt(mechMat.Csigma);//###
    mechMat.elasticSigmaLimit = elasticSigmaLimit;
    mechMat.elasticEpsLimit = elasticSigmaLimit/mechMat.tan_el(mechMat.elasticParameters0);
    // построение сплайнов (и вывод в файлы для отладки)
    // eps(sigma)
    {
        POINT3 p1, p2, p3;
        //double k = 10;//0.5;//100000;
        //double k = 10;//100000;
        double E1 = 1/mechMat.tan_el(mechMat.elasticParameters0);//1/(3*mechMat.elasticParameters0.G);
        double E2 = E1*k_eps_sigma;   // неидеальная пластичность
        double percent1 = 0.01;//0.01;
        double percent2 = percent1/k_eps_sigma*2;
        // чуть до текучести
        p1[0] = elasticSigmaLimit*(1 - percent1);
        p1[1] = elasticSigmaLimit*(1 - percent1)*E1;
        // начало текучести
        p2[0] = elasticSigmaLimit;
        p2[1] = elasticSigmaLimit*E1;
        // чуть после текучести
        p3[0] = elasticSigmaLimit + elasticSigmaLimit*percent2;
        p3[1] = elasticSigmaLimit*E1 + elasticSigmaLimit*percent2*E2;

        // Безье
        //mechMat.epsFun.setBezie(E1, p1, p2, p3, E2);
        //mechMat.difEpsFun.setDifBezie(E1, p1, p2, p3, E2);

        // Кусочно-линейная
        {
            POINT3 p1;
            // начало текучести
            p1[0] = elasticSigmaLimit;
            p1[1] = elasticSigmaLimit*E1;
            mechMat.epsFun.setPlasticEps(p1, E1, E2);
            mechMat.difEpsFun.setDifPlasticEps(p1, E1, E2);
        }

        // вывод кривой в файл
        FILE *f = fopen("___plasticEpsSigma.txt", "w");
        FILE *df = fopen("___dif_plasticEpsSigma.txt", "w");
        double sigma1 = 0;
        double sigma2 = elasticSigmaLimit + 20000; // strToDouble(m0.vals[1]) - предел упругости
        int N = 10;//1000000;
        for(int i = 0; i < N; i++)
        {
            double sigma = sigma1+(sigma2-sigma1)*i/N;
            double eps = mechMat.eps(mechMat.elasticParameters0, 0, sigma, 0);
            double deps = mechMat.difEps(mechMat.elasticParameters0, 0, sigma, 0);
            if(eps > 0.003) break;
            //if(sigma >= p1[0] && sigma <= p3[0])
            {
                fprintf(f, "%le %le\n", eps, sigma);
                fprintf(df, "%le %le\n", eps, deps);
            }
            int i_new = i + N/1000;
            double sigma_new = sigma1+(sigma2-sigma1)*i_new/N;
            if(sigma_new <= p3[0])
                i = i_new;
        }
        fclose(df);
        fclose(f);
    }
    // sigma(eps)
    {
        POINT3 p1;
        //double k = 1000;//0.5;//100000;//0.5;//100000;
        double E1 = mechMat.tan_el(mechMat.elasticParameters0);//3*mechMat.elasticParameters0.G;
        //double E2 = E1/10; // неидеальная пластичность
        //double E2 = -E1/10; // отрицательный наклон
        //double E2 = 0;    // идеальная пластичность
        //double E2 = E1/k;
        double E2 = E1/k_sigma_eps;
        //double E2 = 0;
        double elasticEpsLimit = elasticSigmaLimit/E1;
        // начало текучести
        p1[0] = elasticEpsLimit;
        p1[1] = elasticEpsLimit*E1;
        mechMat.sigmaFun.setPlasticSigma(p1, E1, E2);
        mechMat.difSigmaFun.setDifPlasticSigma(p1, E1, E2);
        // вывод кривой в файл
        FILE *f = fopen("___plasticSigmaEps.txt", "w");
        FILE *df = fopen("___dif_plasticSigmaEps.txt", "w");
        double eps1 = 0;
        double eps2 = 0.005; // strToDouble(m0.vals[1]) - предел упругости
        int N = 1000;//1000000;
        for(int i = 0; i < N; i++)
        {
            double eps = eps1+(eps2-eps1)*i/N;
            double sigma = mechMat.sigma(mechMat.elasticParameters0, 0, eps, 0);
            double deps = mechMat.difSigma(mechMat.elasticParameters0, 0, eps, 0);
            fprintf(f, "%le %le\n", eps, sigma);
            fprintf(df, "%le %le\n", eps, deps);
        }
        fclose(df);
        fclose(f);
    }
}
void setPlasticCurve_Yeld(Solid::MechMaterialSource &mechMat, const double k_sigma_eps, const double k_eps_sigma)
{
    using namespace Solid;
    // скорректированный предел текучести. А вообще, лучше задавать кривые и обработку результата сразу подогнанные под коэффициенты Csigma, Ceps
    double elasticSigmaLimit = mechMat.elasticSigmaLimit*sqrt(mechMat.Csigma);//###
    mechMat.elasticSigmaLimit = elasticSigmaLimit;
    mechMat.elasticEpsLimit = elasticSigmaLimit/mechMat.tan_el(mechMat.elasticParameters0);
    // построение сплайнов
    /*
    // eps(sigma)
    {
        POINT3 p1, p2, p3;
        //double k = 10;//0.5;//100000;
        //double k = 10;//100000;
        double E1 = 1/mechMat.tan_el(mechMat.elasticParameters0);//1/(3*mechMat.elasticParameters0.G);
        double E2 = E1*k_eps_sigma;   // неидеальная пластичность
        double percent1 = 0.01;//0.01;
        double percent2 = percent1/k_eps_sigma*2;
        // чуть до текучести
        p1[0] = elasticSigmaLimit*(1 - percent1);
        p1[1] = elasticSigmaLimit*(1 - percent1)*E1;
        // начало текучести
        p2[0] = elasticSigmaLimit;
        p2[1] = elasticSigmaLimit*E1;
        // чуть после текучести
        p3[0] = elasticSigmaLimit + elasticSigmaLimit*percent2;
        p3[1] = elasticSigmaLimit*E1 + elasticSigmaLimit*percent2*E2;

        // Безье
        //mechMat.epsFun.setBezie(E1, p1, p2, p3, E2);
        //mechMat.difEpsFun.setDifBezie(E1, p1, p2, p3, E2);

        // Кусочно-линейная
        {
            POINT3 p1;
            // начало текучести
            p1[0] = elasticSigmaLimit;
            p1[1] = elasticSigmaLimit*E1;
            mechMat.epsFun.setPlasticEps(p1, E1, E2);
            mechMat.difEpsFun.setDifPlasticEps(p1, E1, E2);
        }

        // вывод кривой в файл
        FILE *f = fopen("___plasticEpsSigma.txt", "w");
        FILE *df = fopen("___dif_plasticEpsSigma.txt", "w");
        double sigma1 = 0;
        double sigma2 = elasticSigmaLimit + 20000; // strToDouble(m0.vals[1]) - предел упругости
        int N = 10;//1000000;
        for(int i = 0; i < N; i++)
        {
            double sigma = sigma1+(sigma2-sigma1)*i/N;
            double eps = mechMat.eps(mechMat.elasticParameters0, 0, sigma, 0);
            double deps = mechMat.difEps(mechMat.elasticParameters0, 0, sigma, 0);
            if(eps > 0.003) break;
            //if(sigma >= p1[0] && sigma <= p3[0])
            {
                fprintf(f, "%le %le\n", eps, sigma);
                fprintf(df, "%le %le\n", eps, deps);
            }
            int i_new = i + N/1000;
            double sigma_new = sigma1+(sigma2-sigma1)*i_new/N;
            if(sigma_new <= p3[0])
                i = i_new;
        }
        fclose(df);
        fclose(f);
    }
    */
    // sigma(eps)
    {
        double E1 = mechMat.tan_el(mechMat.elasticParameters0);//3*mechMat.elasticParameters0.G;
        double E2 = E1/k_sigma_eps;
        mechMat.sigmaFun.setPlasticSigma_Yeld(elasticSigmaLimit, E1, E2);
        mechMat.difSigmaFun.setDiffPlasticSigma_Yeld(elasticSigmaLimit, E1, E2);
    }
}
// пластичность с изотропным степенным упрочнением с коэффициентом n
void setPlasticCurve_Yeld_hardening(Solid::MechMaterialSource &mechMat, const double n)
{
    using namespace Solid;
    // скорректированный предел текучести. А вообще, лучше задавать кривые и обработку результата сразу подогнанные под коэффициенты Csigma, Ceps
    double elasticSigmaLimit = mechMat.elasticSigmaLimit*sqrt(mechMat.Csigma);//###
    mechMat.elasticSigmaLimit = elasticSigmaLimit;
    mechMat.elasticEpsLimit = elasticSigmaLimit/mechMat.tan_el(mechMat.elasticParameters0);
    // построение кусочно-линейного сплайна
    // sigma(eps_pl)
    {
        double E = mechMat.elasticParameters0.E;
        mechMat.sigmaFun.setPlasticSigma_Yeld_hardening(mechMat.elasticSigmaLimit, E, n);
        mechMat.difSigmaFun.setDiffPlasticSigma_Yeld_hardening(mechMat.elasticSigmaLimit, E, n);
    }
}

// инициализация кривых ползучего материала
void setCreepMaterialCurve(const double A, const double n, const double m, Solid::MechMaterialSource &mechMat)
{
    mechMat.epsFun.setCreepEps(A, n, m, mechMat.tan_el(mechMat.elasticParameters0));
    mechMat.difEpsFun.setDifCreepEps(A, n, m, mechMat.tan_el(mechMat.elasticParameters0));
    mechMat.sigmaFun.setCreepSigma(A, n, m, mechMat.tan_el(mechMat.elasticParameters0));
    mechMat.difSigmaFun.setDifCreepSigma(A, n, m, mechMat.tan_el(mechMat.elasticParameters0));
    // вывод кривых в файл
    // eps(sigma)
    {
        FILE *f = fopen("___creepEpsCurve.txt", "w");
        int N = 1000;//1000000;
        double t1 = 0;
        double t2 = 1000;
        for(int i = 0; i < N; i++)
        {
            double t = t1+(t2-t1)*i/N;
            double sigmaEqv = 200;
            //fun.setArgumentValue(0, sigmaEqv);
            //fun.setArgumentValue(1, t);
            //double eps = fun.solve();
            double eps = mechMat.eps(mechMat.elasticParameters0, 0, sigmaEqv, t);
            fprintf(f, "%le %le %le\n", sigmaEqv, t, eps);
        }
        fclose(f);
    }
    // sigma(eps)
    {
        FILE *f = fopen("___creepSigmaCurve.txt", "w");
        int N = 1000;//1000000;
        double t1 = 0;
        double t2 = 1000;
        for(int i = 0; i < N; i++)
        {
            double t = t1+(t2-t1)*i/N;
            double epsEqv = 0.002;
            //fun.setArgumentValue(0, sigmaEqv);
            //fun.setArgumentValue(1, t);
            //double eps = fun.solve();
            double sigma = mechMat.sigma(mechMat.elasticParameters0, 0, epsEqv, t);
            fprintf(f, "%le %le %le\n", epsEqv, t, sigma);
        }
        fclose(f);
    }
}
// инициализировать упругий или упруго-пластичный материал
void setPlasticMaterialCurve(const int plasticPlot, Solid::MechMaterialSource &mechMat, const double k_sigma_eps, const double k_eps_sigma)
{
    using namespace Solid;
    // упругость
    if(plasticPlot == 0)
    {
        mechMat.plasticityMethodType = MechPlasticityMethodType::Elasticity;
        return;
    }
    // кусочно-линейная sigma(eps) (и дополнительно Безье eps(sigma))
    if(plasticPlot == 1)
    {
        setPlasticCurve(mechMat, k_sigma_eps, k_eps_sigma);
        mechMat.plasticityMethodType = MechPlasticityMethodType::D_pl;
        mechMat.PCDependenceType = MechPlasticityCurveDependenceType::None;
        mechMat.PCUnloadingType = MechPlasticityCurveUnloadingType::Elastic;
        mechMat.PCType = MechPlasticityCurveType::Sigma_eps;
    }
    // Безье eps(sigma) (и дополнительно кусочно-линейная sigma(eps))
    if(plasticPlot == 2)
    {
        setPlasticCurve(mechMat, k_sigma_eps, k_eps_sigma);
        mechMat.plasticityMethodType = MechPlasticityMethodType::D_pl;
        mechMat.PCDependenceType = MechPlasticityCurveDependenceType::None;
        mechMat.PCUnloadingType = MechPlasticityCurveUnloadingType::Elastic;
        mechMat.PCType = MechPlasticityCurveType::Eps_sigma;
    }
}
void setPlasticMaterialCurve_Yeld(const int plasticPlot, MechMaterialSource &mechMat, const double k_sigma_eps, const double k_eps_sigma)
{
    using namespace Solid;
    // упругость
    if(plasticPlot == 0)
    {
        mechMat.plasticityMethodType = MechPlasticityMethodType::Elasticity;
        return;
    }
    // кусочно-линейная sigma(eps) (и дополнительно Безье eps(sigma))
    if(plasticPlot == 1)
    {
        setPlasticCurve_Yeld(mechMat, k_sigma_eps, k_eps_sigma);
        mechMat.plasticityMethodType = MechPlasticityMethodType::D_pl;
        mechMat.PCDependenceType = MechPlasticityCurveDependenceType::None;
        mechMat.PCUnloadingType = MechPlasticityCurveUnloadingType::Elastic;
        mechMat.PCType = MechPlasticityCurveType::Sigma_eps;
    }
    // Безье eps(sigma) (и дополнительно кусочно-линейная sigma(eps))
    if(plasticPlot == 2)
    {
        setPlasticCurve_Yeld(mechMat, k_sigma_eps, k_eps_sigma);
        mechMat.plasticityMethodType = MechPlasticityMethodType::D_pl;
        mechMat.PCDependenceType = MechPlasticityCurveDependenceType::None;
        mechMat.PCUnloadingType = MechPlasticityCurveUnloadingType::Elastic;
        mechMat.PCType = MechPlasticityCurveType::Eps_sigma;
    }
}

void setPlasticMaterialCurve_Yeld_hardening(const int plasticPlot, MechMaterialSource &mechMat, const double n)
{
    using namespace Solid;
    // упругость
    if(plasticPlot == 0)
    {
        mechMat.plasticityMethodType = MechPlasticityMethodType::Elasticity;
        return;
    }
    // кусочно-линейная sigma(eps) (и дополнительно Безье eps(sigma))
    if(plasticPlot == 1)
    {
        setPlasticCurve_Yeld_hardening(mechMat, n);
        mechMat.plasticityMethodType = MechPlasticityMethodType::D_pl;
        mechMat.PCDependenceType = MechPlasticityCurveDependenceType::None;
        mechMat.PCUnloadingType = MechPlasticityCurveUnloadingType::Elastic;
        mechMat.PCType = MechPlasticityCurveType::Sigma_eps;
    }
}

// отсутствие 2-х краевых
void setBc2_none(std::vector<Solid::MechBoundaryCondition2Source_base *> *&bc2Source)
{
    using namespace Solid;
    bc2Source = new std::vector<MechBoundaryCondition2Source_base *>;
    bc2Source->clear();
}
// внутреннее давление на полый шар линейной функцией
void setBc2Sphere(std::vector<Solid::MechBoundaryCondition2Source_base *> *&bc2Source, double P1, double P2, double t1, double t2)
{
    bc2Source = new std::vector<MechBoundaryCondition2Source_base *>(2);
    // линейная функция
    setLinearBc2_el(P1, P2, t1, t2,
                    (*bc2Source)[0]);
    (*bc2Source)[1] = new MechBoundaryCondition2Source_None;
}

void setBcSphereT_1(std::vector<Heat::ThermBoundaryCondition1Source> *&thermBc1, std::vector<Heat::ThermBoundaryCondition2Source> *&thermBc2,
                    double T1, double T2)
{
    using namespace Heat;
    // первые краевые условия
    thermBc1 = new std::vector<ThermBoundaryCondition1Source>(5);
    (*thermBc1)[0] = {-1, -1};
    (*thermBc1)[1] = {-1, -1};
    (*thermBc1)[2] = {-1, -1};
    (*thermBc1)[3] = {0, T1};
    (*thermBc1)[4] = {0, T2};
    // вторые краевые условия
    thermBc2 = new std::vector<ThermBoundaryCondition2Source>(1);
    ThermBoundaryCondition2Source &bc2_el = (*thermBc2)[0];
    bc2_el.mode = ThermBoundaryCondition2Type::scalar;// к/у присутствует
    bc2_el.hi = 0;  // коэффициент конвективного теплообмена
    bc2_el.Ta = 0;  // температура окружающей среды
    bc2_el.q = 0*10;   // текущая плотность набегающего теплового потока
}
void setBcSphereT_2(std::vector<Heat::ThermBoundaryCondition1Source> *&thermBc1, std::vector<Heat::ThermBoundaryCondition2Source> *&thermBc2,
                    double T1, double Ta, double hi)
{
    using namespace Heat;
    // первые краевые условия
    thermBc1 = new std::vector<ThermBoundaryCondition1Source>(5);
    (*thermBc1)[0] = {-1, -1};
    (*thermBc1)[1] = {-1, -1};
    (*thermBc1)[2] = {-1, -1};
    (*thermBc1)[3] = {0, T1};
    (*thermBc1)[4] = {-1, -1};
    // вторые краевые условия
    thermBc2 = new std::vector<ThermBoundaryCondition2Source>(2);
    (*thermBc2)[0].mode = ThermBoundaryCondition2Type::none;
    ThermBoundaryCondition2Source &bc2_el = (*thermBc2)[1];
    bc2_el.mode = ThermBoundaryCondition2Type::scalar;    // к/у присутствует
    bc2_el.hi = hi;     // коэффициент конвективного теплообмена
    bc2_el.Ta = Ta;     // температура окружающей среды
    bc2_el.q = 0*10;    // текущая плотность набегающего теплового потока
}
void setBcSphereT_3(std::vector<Heat::ThermBoundaryCondition1Source> *&thermBc1, std::vector<Heat::ThermBoundaryCondition2Source> *&thermBc2,
                    double Ta1, double hi1, double Ta2, double hi2)
{
    using namespace Heat;
    // первые краевые условия
    thermBc1 = new std::vector<ThermBoundaryCondition1Source>(5);
    (*thermBc1)[0] = {-1, -1};
    (*thermBc1)[1] = {-1, -1};
    (*thermBc1)[2] = {-1, -1};
    (*thermBc1)[3] = {-1, -1};
    (*thermBc1)[4] = {-1, -1};
    // вторые краевые условия
    thermBc2 = new std::vector<ThermBoundaryCondition2Source>(2);
    {
        ThermBoundaryCondition2Source &bc2_el = (*thermBc2)[0];
        bc2_el.mode = ThermBoundaryCondition2Type::scalar;    // к/у присутствует
        bc2_el.hi = hi1;    // коэффициент конвективного теплообмена
        bc2_el.Ta = Ta1;    // температура окружающей среды
        bc2_el.q = 0;       // текущая плотность набегающего теплового потока
    }
    {
        ThermBoundaryCondition2Source &bc2_el = (*thermBc2)[1];
        bc2_el.mode = ThermBoundaryCondition2Type::scalar;    // к/у присутствует
        bc2_el.hi = hi2;    // коэффициент конвективного теплообмена
        bc2_el.Ta = Ta2;    // температура окружающей среды
        bc2_el.q = 0;       // текущая плотность набегающего теплового потока
    }
}
void setBcSphereT_4(std::vector<Heat::ThermBoundaryCondition1Source> *&thermBc1, std::vector<Heat::ThermBoundaryCondition2Source> *&thermBc2, Heat::ThermMaterialSource &thermMat,
                    double q, double L, double T2)
{
    using namespace Heat;
    // тензор теплопроводности
    thermMat.L.initDiag(L);
    // первые краевые условия
    thermBc1 = new std::vector<ThermBoundaryCondition1Source>(5);
    (*thermBc1)[0] = {-1, -1};
    (*thermBc1)[1] = {-1, -1};
    (*thermBc1)[2] = {-1, -1};
    (*thermBc1)[3] = {-1, -1};
    (*thermBc1)[4] = {0, T2};
    // вторые краевые условия
    thermBc2 = new std::vector<ThermBoundaryCondition2Source>(2);
    {
        ThermBoundaryCondition2Source &bc2_el = (*thermBc2)[0];
        bc2_el.mode = ThermBoundaryCondition2Type::scalar;    // к/у присутствует
        bc2_el.hi = 0;    // коэффициент конвективного теплообмена
        bc2_el.Ta = 0;    // температура окружающей среды
        bc2_el.q = q;       // текущая плотность набегающего теплового потока
    }
    {
        ThermBoundaryCondition2Source &bc2_el = (*thermBc2)[1];
        bc2_el.mode = ThermBoundaryCondition2Type::none;
    }
}
void setBcSphereT_5(std::vector<Heat::ThermBoundaryCondition1Source> *&thermBc1, std::vector<Heat::ThermBoundaryCondition2Source> *&thermBc2, Heat::ThermMaterialSource &thermMat,
                    double q, double L)
{
    using namespace Heat;
    // тензор теплопроводности
    thermMat.L.initDiag(L);
    // первые краевые условия
    thermBc1 = new std::vector<ThermBoundaryCondition1Source>(5);
    (*thermBc1)[0] = {-1, -1};
    (*thermBc1)[1] = {-1, -1};
    (*thermBc1)[2] = {-1, -1};
    (*thermBc1)[3] = {-1, -1};
    (*thermBc1)[4] = {-1, -1};
    // вторые краевые условия
    thermBc2 = new std::vector<ThermBoundaryCondition2Source>(2);
    {
        ThermBoundaryCondition2Source &bc2_el = (*thermBc2)[0];
        bc2_el.mode = ThermBoundaryCondition2Type::scalar;    // к/у присутствует
        bc2_el.hi = 0;    // коэффициент конвективного теплообмена
        bc2_el.Ta = 0;    // температура окружающей среды
        bc2_el.q = q;       // текущая плотность набегающего теплового потока
    }
    {
        ThermBoundaryCondition2Source &bc2_el = (*thermBc2)[1];
        bc2_el.mode = ThermBoundaryCondition2Type::none;
    }
}



void sortms(const Elementary::VECTOR3 &ms, int *ind)
{
    ind[0] = 0;
    ind[1] = 1;
    ind[2] = 2;
    double mindelta = 1.e100;
    int ind0, ind1;
    for(int i0 = 0; i0 < 3; i0++)
        for(int i1 = i0+1; i1 < 3; i1++)
        {
            double delta = fabs((ms[i0] - ms[i1])/(MIN(fabs(ms[i0]), fabs(ms[i1]))));
            if(delta < mindelta)
            {
                mindelta = delta;
                ind0 = i0;
                ind1 = i1;
            }
        }
    if(ind0 == 0 && ind1 == 1)
    {
        ind[0] = 2;
        ind[1] = 0;
        ind[2] = 1;
    }
    if(ind0 == 0 && ind1 == 2)
    {
        ind[0] = 1;
        ind[1] = 0;
        ind[2] = 2;
    }

}

void sortms_maxToMin(VECTOR3 &ms)
{
    for(int i0 = 0; i0 < 3; i0++)
        for(int i1 = i0 + 1; i1 < 3; i1++)
        {
            if(ms[i0] < ms[i1])
            {
                std::swap(ms[i0], ms[i1]);
            }
        }
}

void calcAnalit(double a,
                 double b,
                 double R,
                 double P0,
                 double elasticSigmaLimit,
                 double c,
                 double &sigma_r,
                 double &sigma_fi,
                 double &sigma_r_0,
                 double &sigma_fi_0
                 )
{
    double sigma_r_upr;
    double sigma_fi_upr;
    double q = 0;
    // упругое решение (вспомогательное)
    {
        double PP = P0 * pow(a,3) / (pow(b,3) - pow(a,3));
        sigma_r_upr = PP * (1 - pow(b/R, 3));
        sigma_fi_upr = PP * (1 + 0.5 * pow(b/R, 3));
    }
    if(c == -1) // упругое либо пластичное решение
    {
        double PnnUprMax = elasticSigmaLimit * 2./3.*(pow(b,3) - pow(a,3))/pow(b,3);
        if(P0 <= PnnUprMax)
        {
            // упругий закон
            sigma_r = sigma_r_upr;
            sigma_fi = sigma_fi_upr;
        }
        else
        {
            // пластичность
            sigma_r = 2.*elasticSigmaLimit*log(R/a) - P0;
            sigma_fi = sigma_r+elasticSigmaLimit;
        }
    }
    else    // упруго-пластичное решение
    {
        q = fabs(2*elasticSigmaLimit*log(c/a) - P0);
        if(R > c)
        {
            // упругий закон
            double PP1 = q * pow(c,3.) / (pow(b,3) - pow(c,3));
            sigma_r = PP1 * (1 - pow(b/R, 3));
            sigma_fi = PP1 * (1 + 0.5 * pow(b/R, 3.));
        }
        else
        {
            // пластичность
            sigma_r = 2.*elasticSigmaLimit*log(R/a) - P0;
            sigma_fi = sigma_r + elasticSigmaLimit;
        }
    }
    sigma_r_0 = sigma_r - sigma_r_upr;
    sigma_fi_0 = sigma_fi - sigma_fi_upr;
}
double calc_c(double a, double b, double P0, double elasticSigmaLimit)
{
    double c;
    double c1 = a, c2 = b;
    double LOW = 1.e-15;
    for(;;)
    {
        c = (c1+c2)/2;
        if(c2-c1 < LOW) break;
        double d = log(c/a)-pow(c/b,3.)/3. - P0/2./elasticSigmaLimit+1./3.;
        if(fabs(d) < LOW) break;
        if(d < 0)
            c1 = c;
        else
            c2 = c;
    }
    if(c2-c1 < LOW)
        // упругое либо пластичное решение
        c = -1;
    return c;
}


// Деформация пластины с отверстием
/*
void setBc2Kirch(std::vector<Solid::MechBoundaryCondition2Source_base *> *&bc2Source, double P1, double P2, double t1, double t2, int mode)
{
    using namespace Solid;
     bc2Source = new std::vector<Solid::MechBoundaryCondition2Source_base *>(5);
     //(*bc2Source)[0] = new MechBoundaryCondition2Source_None;
     //(*bc2Source)[1] = new MechBoundaryCondition2Source_None;
     //(*bc2Source)[2] = new MechBoundaryCondition2Source_None;
     //(*bc2Source)[3] = new MechBoundaryCondition2Source_None;
     //(*bc2Source)[4] = new MechBoundaryCondition2Source_None;
     // линейная функция
     FunParser::Function *fun = new FunParser::Function;
     char Px[1000];
     double PA, PB, PC;
     solveCoef(t1, t2, P1, P2, PA, PB, PC);
     if(PA > 0)
     {
         if(PB > 0)
             sprintf(Px, "(t-%.16le)*%.16le+%.16le\0", PA, PB, PC);
         else
             sprintf(Px, "(t-%.16le)*(0-%.16le)+%.16le\0", PA, fabs(PB), PC);
     }
     else
     {
         if(PB > 0)
             sprintf(Px, "(t+%.16le)*%.16le+%.16le\0", fabs(PA), PB, PC);
         else
             sprintf(Px, "(t+%.16le)*(0-%.16le)+%.16le\0", fabs(PA), fabs(PB), PC);
     }
     (*fun).setExpression(std::string(Px));
     (*fun).args.clear();
     (*fun).addArgument("t\0");
     (*fun).parse();
     (*fun).setName("P(t)");

     if(mode == 1)
     {
         //(*bc2Source)[0].dPdt = {dPdt, 0, 0};
         //(*bc2Source)[1].dPdt = {dPdt, 0, 0};
         //(*bc2Source)[2].dPdt = {0, 0, 0};
         //(*bc2Source)[3].dPdt = {0, 0, 0};
         //(*bc2Source)[4].dPdt = {0, 0, 0};
     }
     if(mode == 2)
     {
         for(int i = 0; i < 4; i++)
         {
             MechBoundaryCondition2Source_ScalarFunction *mbc2Source_sf = new MechBoundaryCondition2Source_ScalarFunction;
             mbc2Source_sf->init(fun);
             (*bc2Source)[i] = mbc2Source_sf;
         }
         (*bc2Source)[4] = new MechBoundaryCondition2Source_None;
     }
     if(mode == 3)
     {
         for(int i = 0; i < 4; i++)
         {
             (*bc2Source)[i] =  new MechBoundaryCondition2Source_None;
         }
         MechBoundaryCondition2Source_ScalarFunction *mbc2Source_sf = new MechBoundaryCondition2Source_ScalarFunction;
         mbc2Source_sf->init(fun);
         (*bc2Source)[4] = mbc2Source_sf;
     }
}
*/
/*
void genTestKirch(TaskTest &task)
{
    // индекс теста
    task.testIndex = 2;
    // подрубка
    task.mechTask.enabled = true;
    task.thermTask.enabled = false;
    // общие шаги и сетка
    std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
    Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
    // сетка
    KirschParameters &kp = task.kirch;
    kp.a = 100;
    kp.b = 100;
    kp.c = 1;
    kp.r = 1;
    kp.N_r = 32;
    kp.N_fi = kp.N_r;
    kp.N_z = 1;
    kp.q = 1. + 1. / 4.;
    kp.Pmode = 2;               // 1 - одноосное растяжение, 2 - давление снаружи
                                // 3 - давление изнутри
    // 32-1/4
    // 64-1/10
    // 96-1/16
    grid->genKirch(task.kirch);



    // параметры
    int fixGrid = 0;        // 0 - подвижная сетка, 1 - зафиксировать сетку
    int plasticPlot = 2;    // 0 - упругость, 1 - сплайн, 2 - Безье
    int plasticityCurveMode = 1;    // 0 - eps(sigma), 1 - sigma(eps)
    bool usingSecantMethod = true;
    IncForsesMode incForsesMode = IncForsesMode::IncrementP;         // 1 - dSigma = integral(...), 2 - в приращениях (P->dP)
    double epsResidualLimit = 1.e-14;
    int StepsNumber1 = 160/2;    // шаги
    int StepsNumber2 = 160/2;
    int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny


    double Time = 1;
    // шаги
    GlobalStep s;       // шаг
    // первый шаг (1 временной слой, заведомо упругое нагружение P1)
    s.t1 = 0;                          // (имеет значение только для первого глобального шага)
    s.t2 = Time*1;
    s.dt0 = Time/StepsNumber1;
    step->push_back(s);
    // второй шаг (пластичное нагружение P2)
    s.t1 = Time*1;
    s.t2 = Time*2;
    s.dt0 = Time/StepsNumber2;
    step->push_back(s);
    // сетка и шаги общие
    task.grid = grid;
    task.step = step;


    // МДТТ
    if(task.mechTask.enabled)
    {
        using namespace Solid;
        // сетка
        task.mechTask.grid = grid;
        // общие параметры шагов
        task.mechTask.step = step;
        // материал
        task.mechTask.material = new std::vector<MechMaterialSource>(1);
        MechMaterialSource &mechMat = (*task.mechTask.material)[0];
        mechMat.elasticSigmaLimit = 2.e7;
        mechMat.Ceps = 4./9.;
        mechMat.Csigma = 1;
        mechMat.set_E_NU(1.e10,0.3);
        mechMat.set_D_isotropic_XY();
        mechMat.set_M_sigma();
        mechMat.F = {0,0,0};
        mechMat.ro = 0;//1e10;
        setPlasticMaterialCurve(plasticPlot, plasticityCurveMode, Solid::MechStressStrainMethod::variableD, usingSecantMethod, mechMat);
        // шаги (на каждом шаге задаются вторые краевые условия)
        std::vector<MechGlobalStep> *mechStep = new std::vector<MechGlobalStep>;
        // нагружение
        if(kp.Pmode == 3)
            task.beamTask.P = 2.0e7;//1*mechMat.elasticSigmaLimit;  // мах = 1.154700538sT = 2.309401077
        if(kp.Pmode == 2)
            task.beamTask.P = 0.995*mechMat.elasticSigmaLimit;
        //task.beamTask.P = 1000;
        double P = task.beamTask.P;
        double P1;
        if(kp.Pmode == 3)
            P1 = P/2.;//+(mat.elasticSigmaLimit/2.); // /sqrt(3)
        if(kp.Pmode == 2)
            P1 = P/2.;//+(mat.elasticSigmaLimit/2.); // /sqrt(3)
        double P2 = P;
        MechGlobalStep ms;  // шаг для МДТТ
        GlobalStep s;       // шаг
        std::vector<MechBoundaryCondition2Source> *bc2Source;
        ms.slausolverParameters = task.slausolver_parameters;
        ms.timeMode = 0;    // 0 - квазистатическая задача
        ms.incForsesMode = incForsesMode;  // 1 - dSigma = integral(...), 2 - приращения
        ms.fixGrid = fixGrid;     // 1 - зафиксировать сетку
        ms.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
        ms.controlMode = 0; // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
        ms.iterLimit = 128;
        ms.slauResidualLimit = 1;
        ms.epsResidualLimit = epsResidualLimit;
        ms.sigmaResidualLimit = 1;
        // первый шаг (1 временной слой, заведомо упругое нагружение P1)
        s = (*step)[0];
        setBc2Kirch(bc2Source, 0, P1, s.t1, s.t2, kp.Pmode);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);
        // второй шаг (пластичное нагружение P2)
        s = (*step)[1];
        setBc2Kirch(bc2Source, P1, P2, s.t1, s.t2, kp.Pmode);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);
        task.mechTask.mechStep = mechStep;
        // первые краевые условия
        std::vector<MechBoundaryCondition1Source> *mechBc1 = new std::vector<MechBoundaryCondition1Source>(4);
        (*mechBc1)[0].mode = {{0, -1, -1}};
        (*mechBc1)[0].u0 = {{0, -1, -1}};
        (*mechBc1)[1].mode = {{-1, 0, -1}};
        (*mechBc1)[1].u0 = {{-1, 0, -1}};
        (*mechBc1)[2].mode = {{-1, -1, 0}};
        (*mechBc1)[2].u0 = {{-1, -1, 0}};
        (*mechBc1)[3].mode = {{-1, -1, -1}};
        (*mechBc1)[3].u0 = {{-1, -1, -1}};
        task.mechTask.bc1Source = mechBc1;
        // инициализация начальных деформаций и напряжений (=0)
        task.mechTask.fe = new std::vector<MechFeData>;
        task.mechTask.fe->resize(grid->fe.size());
        for (size_t i = 0; i < grid->fe.size(); i++)
        {
            (*task.mechTask.fe)[i].init(BasisType_1L,
                     HomogenyMode, Integration::IntegrationType::Gauss3);
        }
        // инициализация начальных перемещений (=0)
        task.mechTask.vertex = new std::vector<MechVertexData>;
        for (size_t i = 0; i < grid->vertex.size(); i++)
        {
            MechVertexData v_el;
            for (int j = 0; j < 3; j++)
            {
                v_el.sum_du[j] = 0;
                v_el.du[j] = 0; // не имеет значения
            }
            task.mechTask.vertex->push_back(v_el);
        }
        // инициализация начальных скоростей (=0) и ускорений (=0)
        task.mechTask.V0 = new Vector(grid->vertex.size()*3);
        task.mechTask.dV0 = new Vector(grid->vertex.size()*3);
        for (size_t i = 0; i < task.mechTask.V0->size(); i++)
        {
            (*task.mechTask.V0)[i] = 0;
            (*task.mechTask.dV0)[i] = 0;
        }
    }
}*/

/*
void BuldGridRegularParallelepiped(const GridRegularParallelepiped_parameters &rpp, Grid3D *&g)
{
    using namespace Operations;

    //int reindex1[6] =
    //{
    //    0, 1, 2,
    //    4, 5, 6
    //};
    //int reindex2[6] =
    //{
    //    2, 1, 3,
    //    6, 5, 7
    //};
    g = new Grid3D;
    int i0, i1, i2;
    int vind0;
    int offset0 = 1;
    int offset1 = rpp.N[0] + 1;
    int offset2 = (rpp.N[0] + 1)*(rpp.N[1] + 1);
    // построение точек
    //g->vertex.resize((rpp.N[0] + 1) * (rpp.N[1] + 1) * (rpp.N[2] + 1));		// количество точек
    for (i2 = 0; i2 <= rpp.N[2]; i2++)
        for (i1 = 0; i1 <= rpp.N[1]; i1++)
            for (i0 = 0; i0 <= rpp.N[0]; i0++)
            {
                POINT3 p;
                findPointOnTheLine_1d(rpp.p1.x[0], rpp.p2.x[0], rpp.N[0], 1, i0, p.x[0]);
                findPointOnTheLine_1d(rpp.p1.x[1], rpp.p2.x[1], rpp.N[1], 1, i1, p.x[1]);
                findPointOnTheLine_1d(rpp.p1.x[2], rpp.p2.x[2], rpp.N[2], 1, i2, p.x[2]);
                // ind = i0*(i1+1)*(i2+1)+i1*(i2+1)+i0 - индекс точки i0,i1,i2
                g->vertex.push_back(p);
            }
    for (i2 = 0; i2 < rpp.N[2]; i2++)
        for (i1 = 0; i1 < rpp.N[1]; i1++)
            for (i0 = 0; i0 < rpp.N[0]; i0++)
            {
                vind0 = i0 + i1*offset1 + i2*offset2;		// индекс вершины, определяющей шестигранник
                int vi[8] =
                {
                    MOVE(vind0, 0, 0, 0), MOVE(vind0, 1, 0, 0), MOVE(vind0, 0, 1, 0), MOVE(vind0, 1, 1, 0),
                    MOVE(vind0, 0, 0, 1), MOVE(vind0, 1, 0, 1), MOVE(vind0, 0, 1, 1), MOVE(vind0, 1, 1, 1)
                };
                FiniteElement fe_el;
                // шестигранник
                fe_el.mi = 0;
                for (int t = 0; t < 8; t++)
                    fe_el.vi[t] = vi[t];
                g->fe.push_back(fe_el);
                //if (rpp.fe_type == GridFEType::Hexagon)
                //{
                //}
                //else
                //if (rpp.fe_type == GridFEType::Prism)
                //{
                //    // 2 призмы
                //    int count = g->fe.push();
                //    g->fe[count].init(GridFEType::Prism, 0);
                //    for (int t = 0; t < 6; t++)
                //        g->fe[count].vi[t] = vi[reindex1[t]];
                //    count = g->fe.push();
                //    g->fe[count].init(GridFEType::Prism, 0);
                //    for (int t = 0; t < 6; t++)
                //        g->fe[count].vi[t] = vi[reindex2[t]];
                //}

                // второе краевое условие (продолжение)
                if (i0 == 0)
                {// сила и грань для 2 краевого слева (bc2_x1)
                    g->bc2.push_back(BoundaryCondition2_4gonal{{
                                        MOVE(vind0, 0, 0, 0),
                                        MOVE(vind0, 0, 1, 0),
                                        MOVE(vind0, 0, 0, 1),
                                        MOVE(vind0, 0, 1, 1),},
                                        0,});
                }
                if (i0 == rpp.N[0] - 1)
                {// сила и грань для 2 краевого справа (bc2_x2)
                    g->bc2.push_back(BoundaryCondition2_4gonal{{
                                        MOVE(vind0, 1, 1, 0),
                                        MOVE(vind0, 1, 0, 0),
                                        MOVE(vind0, 1, 1, 1),
                                        MOVE(vind0, 1, 0, 1),},
                                        1,});
                }
                if (i1 == 0)
                {// сила и грань для 2 краевого снизу (bc2_y1)
                    g->bc2.push_back(BoundaryCondition2_4gonal{{
                                        MOVE(vind0, 0, 0, 0),
                                        MOVE(vind0, 1, 0, 0),
                                        MOVE(vind0, 0, 0, 1),
                                        MOVE(vind0, 1, 0, 1),},
                                        2,});
                }
                if (i1 == rpp.N[1] - 1)
                {// сила и грань для 2 краевого сверху (bc2_y2)
                    g->bc2.push_back(BoundaryCondition2_4gonal{{
                                        MOVE(vind0, 0, 1, 0),
                                        MOVE(vind0, 1, 1, 0),
                                        MOVE(vind0, 0, 1, 1),
                                        MOVE(vind0, 1, 1, 1),},
                                        3,});
                }
                if (i2 == 0)
                {// сила и грань для 2 краевого вблизи (bc2_z1)

                    g->bc2.push_back(BoundaryCondition2_4gonal{{
                                        MOVE(vind0, 0, 0, 0),
                                        MOVE(vind0, 1, 0, 0),
                                        MOVE(vind0, 0, 1, 0),
                                        MOVE(vind0, 1, 1, 0),},
                                        4,});

                    // ударяющая сила (bc2_bang)
                    if(i0 >= rpp.bang_x1 && i0 < rpp.bang_x2 &&
                       i1 >= rpp.bang_y1 && i1 < rpp.bang_y2)
                    {
                        g->bc2.push_back(BoundaryCondition2_4gonal{{
                                            MOVE(vind0, 0, 0, 0),
                                            MOVE(vind0, 1, 0, 0),
                                            MOVE(vind0, 0, 1, 0),
                                            MOVE(vind0, 1, 1, 0),},
                                            6,});
                    }
                }
                if (i2 == rpp.N[2] - 1)
                {// сила и грань для 2 краевого вдали (bc2_z2)
                    g->bc2.push_back(BoundaryCondition2_4gonal{{
                                        MOVE(vind0, 1, 0, 1),
                                        MOVE(vind0, 0, 0, 1),
                                        MOVE(vind0, 1, 1, 1),
                                        MOVE(vind0, 0, 1, 1),},
                                        5,});
                }
            }
    // первые краевые
    //g->bc1.resize((rpp.N[0] + 1)*(rpp.N[1] + 1) + (rpp.N[0] + 1)*(rpp.N[2] + 1) + (rpp.N[1] + 1)*(rpp.N[2] + 1) + 1);
    // фиксируем x на плоскости x = 0
    i0 = rpp.N[0] / 2;
    for (i2 = 0; i2 <= rpp.N[2]; i2++)
        for (i1 = 0; i1 <= rpp.N[1]; i1++)
        {
            BoundaryCondition1 bc1_el;
            bc1_el.vi = i0 + i1*offset1 + i2*offset2;
            bc1_el.si = 0;
            g->bc1.push_back(bc1_el);
        }
    // фиксируем y на плоскости y = 0
    i1 = rpp.N[1] / 2;
    for (i0 = 0; i0 <= rpp.N[0]; i0++)
        for (i2 = 0; i2 <= rpp.N[2]; i2++)
        {
            BoundaryCondition1 bc1_el;
            bc1_el.vi = i0 + i1*offset1 + i2*offset2;
            bc1_el.si = 1;
            g->bc1.push_back(bc1_el);
        }
    // фиксируем z на плоскости z = 0
    i2 = rpp.N[2] / 2;
    for (i0 = 0; i0 <= rpp.N[0]; i0++)
        for (i1 = 0; i1 <= rpp.N[1]; i1++)
        {
            BoundaryCondition1 bc1_el;
            bc1_el.vi = i0 + i1*offset1 + i2*offset2;
            bc1_el.si = 2;
            g->bc1.push_back(bc1_el);
        }
    // фиксируем x,y,z в центре
    BoundaryCondition1 bc1_el;
    i0 = rpp.N[0] / 2;
    i1 = rpp.N[1] / 2;
    i2 = rpp.N[2] / 2;
    bc1_el.vi = i0 + i1*offset1 + i2*offset2;
    bc1_el.si = 3;
    g->bc1.push_back(bc1_el);
    // фиксируем x на плоскости x = min
    i0 = 0;
    for (i2 = 0; i2 <= rpp.N[2]; i2++)
        for (i1 = 0; i1 <= rpp.N[1]; i1++)
        {
            BoundaryCondition1 bc1_el;
            bc1_el.vi = i0 + i1*offset1 + i2*offset2;
            bc1_el.si = 4;
            g->bc1.push_back(bc1_el);
        }
    // фиксируем x на плоскости x = max
    i0 = rpp.N[0];
    for (i2 = 0; i2 <= rpp.N[2]; i2++)
        for (i1 = 0; i1 <= rpp.N[1]; i1++)
        {
            BoundaryCondition1 bc1_el;
            bc1_el.vi = i0 + i1*offset1 + i2*offset2;
            bc1_el.si = 5;
            g->bc1.push_back(bc1_el);
        }
    // фиксируем y на плоскости y = min
    i1 = 0;
    for (i0 = 0; i0 <= rpp.N[0]; i0++)
        for (i2 = 0; i2 <= rpp.N[2]; i2++)
        {
            BoundaryCondition1 bc1_el;
            bc1_el.vi = i0 + i1*offset1 + i2*offset2;
            bc1_el.si = 6;
            g->bc1.push_back(bc1_el);
        }
    // фиксируем y на плоскости y = max
    i1 = rpp.N[1];
    for (i0 = 0; i0 <= rpp.N[0]; i0++)
        for (i2 = 0; i2 <= rpp.N[2]; i2++)
        {
            BoundaryCondition1 bc1_el;
            bc1_el.vi = i0 + i1*offset1 + i2*offset2;
            bc1_el.si = 7;
            g->bc1.push_back(bc1_el);
        }
    // фиксируем z на плоскости z = min
    i2 = 0;
    for (i0 = 0; i0 <= rpp.N[0]; i0++)
        for (i1 = 0; i1 <= rpp.N[1]; i1++)
        {
            BoundaryCondition1 bc1_el;
            bc1_el.vi = i0 + i1*offset1 + i2*offset2;
            bc1_el.si = 8;
            g->bc1.push_back(bc1_el);
        }
    // фиксируем z на плоскости z = max
    i2 = rpp.N[2];
    for (i0 = 0; i0 <= rpp.N[0]; i0++)
        for (i1 = 0; i1 <= rpp.N[1]; i1++)
        {
            BoundaryCondition1 bc1_el;
            bc1_el.vi = i0 + i1*offset1 + i2*offset2;
            bc1_el.si = 9;
            g->bc1.push_back(bc1_el);
        }

}
*/


}   // namespace Tests
