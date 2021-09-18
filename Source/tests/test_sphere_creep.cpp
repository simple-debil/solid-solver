#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"
#include <fstream>

#include "tests.h"
#include "console.h"

using namespace Elementary;
using namespace Solid;
using namespace Grid;

namespace Tests
{
Test_sphere_creep::Test_sphere_creep()
{

    // дирректория данных теста
    dir = "./tests/hollowSphere_creep/";
    // имена gnuplot файлов(+пути к файлам в опциях) и соответствующих картинок
    fn_filesList = dir + "files.txt";
    resGraph.read_gnuplot_from_file(fn_filesList);
    // инициализация графиков с невязками
    initStepsGraphs();

}
Test_base::Type Test_sphere_creep::get_type() const
{
    return Type::Sphere_creep;
}
void Test_sphere_creep::initTask(Task &task)
{
    const int creepCurve = 0;
    const int plasticCurve = 1;

    // кривая (ползучесть/пластичность)
    int creepCurveMode = creepCurve;
    //int creepCurveMode = plasticCurve;

    // параметры кривой ползучести
    //gp.n = 3;
    gp.n = 2;
    //gp.n = 1.5;   // 1 - упругость, ->бесконечности - пластичность

    gp.B = 1.e-14;

    // метод (начальных напряжений, начальных деформаций, пластичность)
    //MechPlasticityMethodType creepMode = MechPlasticityMethodType::D_pl;
    //MechPlasticityMethodType creepMode = MechPlasticityMethodType::InitialEps;
    //MechPlasticityMethodType creepMode = MechPlasticityMethodType::Combo_D_pl_InitialEps;
    //MechPlasticityMethodType creepMode = MechPlasticityMethodType::InitialSigma;
    MechPlasticityMethodType creepMode = MechPlasticityMethodType::Combo_D_pl_InitialSigma;
    //### кривые не для параметра упрочнения Одквиста

    //IncForsesMode incForsesMode = IncForsesMode::MinusIntegral;
    IncForsesMode incForsesMode = IncForsesMode::bPlusR;
    //IncForsesMode incForsesMode = IncForsesMode::IncrementP;

    bool terminateIfAccuracyIsNotAchieving = false;
    double w = 1;
    int controlMode = 0; // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
    // подрубка
    task.thermTask.enabled = false;
    task.mechTask.enabled = true;
    // общие шаги и сетка
    std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
    Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
    // сетка
    gp.r1 = 1;
    gp.r2 = 4;
    gp.N = 8;   // 4 промежутка на четвертьокружностях    //8, 32
    gp.Nparts = gp.N*4;
    gp.q = 1. + 1./(gp.Nparts/2);//1.0625;//1.0625;//1.03125;
    gp.curvilinear = 1;    // отображение (0 - линейное, 1 - квадратичное)
    gp.buildingMethod = 0; // способ построения (0 - делим дуги, 1 - отображение на куб)
    grid->genSphere(gp);
    // параметры
    sigmaSolvingType = 0;  // способ расчёта решения (0 - главные напряжения, 1 - проекции площадки)
    int fixGrid = 0;            // 0 - подвижная сетка, 1 - зафиксировать сетку
    int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny
    int StepsNumber1 = 1;      // шаги
    int StepsNumber2 = 100;//100;
    double Time0 = 1;
    double Time1 = 1.0000001;
    double Time2 = 100;//400
    double P = 1.0e5;
    //double P = 2.0e7;
    double P1;
    double P2;

    double NU = 0.3; //0.4
    double E = 1.e10;//1.e10;
    int nonlinearIterLimit = 1000;//2        // ограничение на количество итераций
    MaxSigmaEpsResidual plasticResidualLimit;   // желаемая невязка по эквиволентным деформациям

    StepsNumber1 = 1;// 1 // для пластичности
    StepsNumber2 = 10;//10 // для пластичности

    // метод начальных напряжений
    if(creepMode == MechPlasticityMethodType::InitialSigma ||
       creepMode == MechPlasticityMethodType::Combo_D_pl_InitialSigma)
    {
        plasticResidualLimit.sigma = 1.e-4;//1000;//1.e-14;//1.e-13;
        plasticResidualLimit.eps = 1.e-14;
        nonlinearIterLimit = 1000;
        // пластичность
        if(creepCurveMode == plasticCurve)
        {
            gp.n = 1.5;
            Time0 = 0;
            Time1 = 1;
            Time2 = 2;
            P1 = 0;//1.30e7;
            P2 = 2.0e7;
        }
        // ползучесть, n = 2
        if(creepCurveMode == creepCurve && gp.n == 2)
        {
            //plasticResidualLimit.eps = 1.e-5;
            Time0 = 0;
            Time1 = 1.e-7;
            Time2 = 10;//1000;//100;//10;//1;
            StepsNumber1 = 1;
            StepsNumber2 = 100;//10000;//1000;//100;//10;
            P1 = P;
            P2 = P;
        }
        // ползучесть, n = 3
        if(creepCurveMode == creepCurve && gp.n == 3)
        {
            Time0 = 0;
            Time1 = 1.e-7;//1.e-10;
            Time2 = .0001;//.0001;
            StepsNumber1 = 1;
            StepsNumber2 = 100;//100;//2000
            P1 = P;
            P2 = P;
        }
        // ползучесть, n = 1.5
        if(creepCurveMode == creepCurve && gp.n == 1.5)
        {
            Time0 = 0;
            Time1 = 0.0000001;
            Time2 = 1000;//500//400
            StepsNumber1 = 1;
            StepsNumber2 = 100;//2000
            P1 = P;
            P2 = P;
        }
    }
    // метод начальных деформаций
    /*
    if(creepMode == MechPlasticityMethodType::Combo_D_pl_InitialEps)
    {
        plasticResidualLimit.eps = 1.e-14;
        plasticResidualLimit.sigma = 1000;
        if(creepCurveMode == creepCurve && gp.n == 1.5)
        {
            terminateIfAccuracyIsNotAchieving = true;
            //incForsesMode = IncForsesMode::MinusIntegral;
            incForsesMode = IncForsesMode::IncrementP; // в скоростях
            nonlinearIterLimit = 200;
            gp.n = 1.5;
            Time0 = 0;
            Time1 = 0.0000001;
            Time2 = 100;//500//100//400
            StepsNumber1 = 1;
            StepsNumber2 = 100;//500;//100
            P1 = P;
            P2 = P;
        }
        if(creepCurveMode == plasticCurve)
        {
            terminateIfAccuracyIsNotAchieving = true;
            w = 0.1;// 0.1/3;// сходится если брать коэффициент наклона k = 10 и StepsNumber2 = 1000
            //incForsesMode = IncForsesMode::IncrementP; // в скоростях
            nonlinearIterLimit = 1000;
            gp.n = 1.5;
            Time0 = 0;
            Time1 = 1;
            Time2 = 10;//500//100//400
            P1 = 0;//1.30e7;
            P2 = 2.0e7;
        }
    }
    */
    // пластичность
    if(creepMode == MechPlasticityMethodType::D_pl)
    {
        nonlinearIterLimit = 100;//2        // ограничение на количество итераций
        plasticResidualLimit.eps = 1.e-14;       // желаемая невязка по эквиволентным деформациям
        plasticResidualLimit.sigma = 1000;//1.e-6;
        controlMode = 0;//1;

        /*
        // приращения
        if(incForsesMode == 1)
        {

        }
        gp.n = 1.5;
        Time1 = 1;
        Time2 = Time1 + 1000;//100//60
        StepsNumber1 = 1;
        StepsNumber2 = 1000;//2000//600
        P1 = 0.001*P;//0.3*P;
        P2 = P;
        */

        // пластичность
        if(creepCurveMode == plasticCurve)
        {
            gp.n = 1.5;
            Time0 = 0;
            Time1 = 1;
            Time2 = 2;
            P1 = 0;//1.30e7;
            //P2 = 2.5e7;
            P2 = 2.0e7;
        }
        // ползучесть, n = 2
        if(creepCurveMode == creepCurve && gp.n == 2)
        {
            Time0 = 0;
            Time1 = 0.0001;
            Time2 = Time1 + 100;
            StepsNumber1 = 10;
            StepsNumber2 = 2000;    //сходится
            //StepsNumber2 = 1000;    //расходится
            P1 = 0.0001*P;
            P2 = P;
        }
        // ползучесть, n = 1.5
        if(creepCurveMode == creepCurve && gp.n == 1.5)
        {
            Time0 = 0;
            Time1 = 0.0001;
            Time2 = Time1 + 1000;
            StepsNumber1 = 1;
            StepsNumber2 = 100;
            P1 = 0.0001*P;
            P2 = P;
        }

        /*
        //сходится
        Time1 = 1;
        Time2 = Time1 + 1000;//100//60
        StepsNumber1 = 1;
        StepsNumber2 = 100;//2000//600
        P1 = 0.001*P;//0.3*P;
        P2 = P;
        */
        /*
        // расходится
        Time1 = 30;
        Time2 = Time1 + 50;//100//60
        StepsNumber1 = 100;
        StepsNumber2 = 10000;//2000//600
        P1 = 0.8*P;//0.3*P;
        P2 = P;
        */
        /*
        Time1 = 10;
        Time2 = Time1 + 50;//100//60
        StepsNumber1 = 1;
        StepsNumber2 = 1000;//2000//600
        P1 = 0.5*P;//0.3*P;
        P2 = P;
        */
        //^ 1000 - не сходится, 100 - сходится
        //P1 = P;
        //P2 = P;

        /*
           incForsesMode = IncForsesMode::MinusIntegral;
        // скорости
        if(incForsesMode == IncForsesMode::IncrementP)
        {
            gp.n = 1.5;
            Time1 = 0.0001;
            Time2 = 1000;//400
            StepsNumber1 = 1;
            StepsNumber2 = 1000;
            P1 = 0.0000001*P;
            P2 = P;
        }
        */
    }

    // шаги
    GlobalStep s;       // шаг
    // первый шаг - 1 временной слой, заведомо упругое нагружение P
    s.t_start = Time0;                          // (имеет значение только для первого глобального шага)
    s.t_finish = Time1;
    s.dt0 = (Time1 - Time0)/StepsNumber1;
    step->push_back(s);
    // второй шаг - ползучесть)
    s.t_start = Time1;
    s.t_finish = Time2;
    s.dt0 = (Time2 - Time1)/StepsNumber2;
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
        std::vector<MechMaterialSource> *material;                      // материалы
        std::vector<MechBoundaryCondition1Source> *bc1Source;           // 1-e краевые условия
        std::vector<MechBoundaryCondition2Source_base *> *bc2Source;    // 2-е краевые
        std::vector<Grid::MechBoundaryCondition2> *bc2;                 // 2-е краевые: индекс КЭ-поверхности + индекс MechBoundaryCondition2Source_base
        // материал
        material = new std::vector<MechMaterialSource>(1);
        MechMaterialSource &mechMat = (*material)[0];
        mechMat.Ceps = 4./9.;//4./3.;//4./9.;
        mechMat.Csigma = 1.;//1./3.;//1;
        mechMat.elasticSigmaLimit = 2.e7;//2.e70;//2.e7;
        mechMat.set_E_NU(E, NU);
        mechMat.elasticParameters0.Talpha = 0;
        mechMat.set_M_sigma();
        mechMat.F = {0,0,0};
        mechMat.elasticParameters0.ro = 0;//1e10;
        // кривая

        double k_sigma_eps = 1000;
        double k_eps_sigma = 10;
        if(creepCurveMode == plasticCurve)
        {
            setPlasticCurve(mechMat, k_sigma_eps, k_eps_sigma);
            mechMat.PCDependenceType = MechPlasticityCurveDependenceType::None;
        }
        if(creepCurveMode == creepCurve)
        {
            setCreepMaterialCurve(gp.B, gp.n, 1, mechMat);
            mechMat.PCDependenceType = MechPlasticityCurveDependenceType::Time;
        }
        // метод
        if(creepMode == MechPlasticityMethodType::InitialSigma ||
           creepMode == MechPlasticityMethodType::Combo_D_pl_InitialSigma)
        {
            // метод начальных напряжений
            mechMat.plasticityMethodType = creepMode;
            mechMat.PCUnloadingType = MechPlasticityCurveUnloadingType::Curve;
            mechMat.PCType = MechPlasticityCurveType::Sigma_eps;
            mechMat.w_midPoint = 0;
            mechMat.w_project = w;
        }
        /*
        if(creepMode == MechPlasticityMethodType::Combo_D_pl_InitialEps)
        {
            // метод начальных деформаций
            mechMat.plasticityMethodType = creepMode;
            mechMat.PCUnloadingType = MechPlasticityCurveUnloadingType::Curve;
            mechMat.PCType = MechPlasticityCurveType::Eps_sigma;
            mechMat.w = w;  // 0.1
        }*/
        if(creepMode == MechPlasticityMethodType::D_pl)
        {
            // пластичность
            mechMat.plasticityMethodType = creepMode;
            mechMat.PCUnloadingType = MechPlasticityCurveUnloadingType::Curve;
            mechMat.PCType = MechPlasticityCurveType::Sigma_eps;
            mechMat.w_midPoint = 0;
            mechMat.w_project = w;
        }
        // первые краевые условия
        bc1Source = new std::vector<MechBoundaryCondition1Source>(5);
        (*bc1Source)[0].mode = {{ 0, -1, -1}};
        (*bc1Source)[0].u0 =   {{ 0, -1, -1}};
        (*bc1Source)[1].mode = {{-1,  0, -1}};
        (*bc1Source)[1].u0 =   {{-1,  0, -1}};
        (*bc1Source)[2].mode = {{-1, -1,  0}};
        (*bc1Source)[2].u0 =   {{-1, -1,  0}};
        (*bc1Source)[3].mode = {{-1, -1, -1}};
        (*bc1Source)[3].u0 =   {{-1, -1, -1}};
        (*bc1Source)[4].mode = {{-1, -1, -1}};
        (*bc1Source)[4].u0 =   {{-1, -1, -1}};
        //(*bc1Source)[4].mode = {{0, 0, 0}};
        //(*bc1Source)[4].u0 =   {{0, 0, 0}};
        // вторые краевые
        bc2 = new std::vector<MechBoundaryCondition2>;
        MechBoundaryCondition2 bc2_el;
        bc2_el.FEsurfaceInd = 0;
        bc2_el.bc2SourceIndex = 0;
        (*bc2).push_back(bc2_el);
        bc2_el.FEsurfaceInd = 1;
        bc2_el.bc2SourceIndex = 1;
        (*bc2).push_back(bc2_el);


        // шаги (на каждом шаге задаются вторые краевые условия)
        std::vector<MechGlobalStepParameters> *mechStep = new std::vector<MechGlobalStepParameters>;
        // нагружение - разгрузка
        MechGlobalStepParameters ms;  // шаг для МДТТ
        GlobalStep s;       // шаг
        ms.slausolverParameters = slausolver_parameters;
        ms.timeMode = 0;    // 0 - квазистатическая задача
        ms.incForsesMode = incForsesMode;  // 1 - dSigma = integral(...), 2 - приращения
        ms.fixGrid = fixGrid;     // 1 - зафиксировать сетку
        ms.switchIterationsMode = SwitchIterationsMode::Serial;
        ms.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
        ms.controlMode = controlMode; // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
        ms.terminateIfAccuracyIsNotAchieving = terminateIfAccuracyIsNotAchieving;
        ms.iterLimit = nonlinearIterLimit;
        ms.slauResidualLimit = 1000;
        ms.plasticResidualLimit = plasticResidualLimit;
        ms.material = material;
        ms.bc1Source = bc1Source;
        //ms.bc2Source = bc2Source; - будет меняться
        ms.bc2 = bc2;
        // первый шаг - 1 временной слой, заведомо упругое нагружение P1
        s = (*step)[0];
        setBc2Sphere(bc2Source, 0, P1, s.t_start, s.t_finish);
        //setBc2Sphere(bc2Source, 0, P, s.t1, s.t2);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);
        // второй шаг - ползучесть
        s = (*step)[1];
        setBc2Sphere(bc2Source, P1, P2, s.t_start, s.t_finish);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);

        task.mechTask.mechStep = mechStep;

        // Контакт (отсутствует)
        // поверхность: жёсткий неподвижный цилиндр
        task.mechTask.rigidSurface = new std::vector<Surface_base *>;
        // контакты КЭ-поверхность - аналитическая поверхность
        task.mechTask.CC_FE_Rigid = new std::vector<ContactCondition_FE_Rigid>;
        // инициализация нулями начальных данных
        task.mechTask.initNull(Integration::IntegrationType::Gauss3);
    }

    // Температура
    if(task.thermTask.enabled)
    {
    }
}
void Test_sphere_creep::writeResults(const Task &task, const OutData &out)
{
    using namespace Interpolation;
    int globalStepNumber = (int)(*out.mechOut).size() - 2;

    // пути к файлам
    {
    // подпапка для данного шага
    std::string subdir = dir + std::to_string(globalStepNumber);
    {
        OS::Console c;
        c.exec(("mkdir " + subdir).c_str());
    }
    subdir += "/";
    // пути к текстовым файлам с данными (в корневой директории, общие для всех глобальных шагов)
    //fn_Fn = dir + "Fn.txt";
    // пути к текстовым файлам с данными (для каждого глобального шага в отдельной директории)
    fn_sigma_r_analit = subdir + "sigma_r_analit.txt";
    fn_sigma_r = subdir + "sigma_r.txt";
    fn_sigma_fi_analit = subdir + "sigma_fi_analit.txt";
    fn_sigma_fi = subdir + "sigma_fi.txt";
    fn_sigma_eqv = subdir + "sigma_eqv.txt";
    fn_eps_el_eqv = subdir + "eps_el_eqv.txt";
    fn_eps_pl_eqv = subdir + "eps_pl_eqv.txt";
    // пути к файлам для отображения графиков
{
    std::ofstream f;
    if(globalStepNumber == 0)
    {
        f.open(fn_filesList, std::ofstream::out);
        resGraph.generalGraphsNomber = 0;
        resGraph.graphsPerGlobalStepNomber = 5;
        f << resGraph.generalGraphsNomber << "\n";
        f << resGraph.graphsPerGlobalStepNomber << "\n";
        /*
        // Fn
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "Fn.png") +
             genGnuplotCommandLineParameter("fn_Fn", fn_Fn) +
             genGnuplotCommandLineParameter("fn_Fn_analit", fn_Fn_analit) +
             "\" " +
             dir + "Fn.gnu" << "\n";
        f << dir + "Fn.png" << "\n";
        f << "Fn.png" << "\n";
        */
    }
    else
    {
        f.open(fn_filesList, std::ofstream::app);
    }
    f << std::to_string(globalStepNumber) + "\n";
    // sigma_fi
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "sigma_fi.png") +
         genGnuplotCommandLineParameter("fn_sigma_fi", subdir + "sigma_fi.txt") +
         genGnuplotCommandLineParameter("fn_sigma_fi_analit", subdir + "sigma_fi_analit.txt") +
         "\" " +
         dir + "sigma_fi.gnu" << "\n";
    f << subdir + "sigma_fi.png" << "\n";
    f << "sigma_fi.png" << "\n";
    // fn_sigma_r
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "sigma_r.png") +
         genGnuplotCommandLineParameter("fn_sigma_r", subdir + "sigma_r.txt") +
         genGnuplotCommandLineParameter("fn_sigma_r_analit", subdir + "sigma_r_analit.txt") +
         "\" " +
         dir + "sigma_r.gnu" << "\n";
    f << subdir + "sigma_r.png" << "\n";
    f << "sigma_r.png" << "\n";
    // sigma_eqv
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "sigma_eqv.png") +
         genGnuplotCommandLineParameter("fn_sigma_eqv", subdir + "sigma_eqv.txt") +
         "\" " +
         dir + "sigma_eqv.gnu" << "\n";
    f << subdir + "sigma_eqv.png" << "\n";
    f << "sigma_eqv.png" << "\n";
    // eps_el_eqv
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "eps_el_eqv.png") +
         genGnuplotCommandLineParameter("fn_eps_el_eqv", subdir + "eps_el_eqv.txt") +
         "\" " +
         dir + "eps_el_eqv.gnu" << "\n";
    f << subdir + "eps_el_eqv.png" << "\n";
    f << "eps_el_eqv.png" << "\n";
    // eps_pl_eqv
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "eps_pl_eqv.png") +
         genGnuplotCommandLineParameter("fn_eps_pl_eqv", subdir + "eps_pl_eqv.txt") +
         "\" " +
         dir + "eps_pl_eqv.gnu" << "\n";
    f << subdir + "eps_pl_eqv.png" << "\n";
    f << "eps_pl_eqv.png" << "\n";
    f.close();
}
    resGraph.read_gnuplot_from_file(fn_filesList);
    }

    FILE *fsigma_r_analit = fopen(fn_sigma_r_analit.c_str(), "w");;
    FILE *fsigma_r = fopen(fn_sigma_r.c_str(), "w");
    FILE *fsigma_fi_analit = fopen(fn_sigma_fi_analit.c_str(), "w");
    FILE *fsigma_fi = fopen(fn_sigma_fi.c_str(), "w");
    FILE *fsigma_eqv = fopen(fn_sigma_eqv.c_str(), "w");
    FILE *feps_el_eqv = fopen(fn_eps_el_eqv.c_str(), "w");
    FILE *feps_pl_eqv = fopen(fn_eps_pl_eqv.c_str(), "w");

    Grid::Grid3D *grid = task.grid;

    {
        // вывод данных МДТТ
        if(task.mechTask.enabled)
        {
            // механическое решение
            int numLocalStep = (int)(*out.mechOut)[globalStepNumber].step.size() - 1;
            MechOutStepData &mechOutStep_el = (*out.mechOut)[globalStepNumber].step[numLocalStep];
            MechBoundaryCondition2Source_ScalarFunction *mbc2Source_sf = (MechBoundaryCondition2Source_ScalarFunction *)(*(*task.mechTask.mechStep)[globalStepNumber].bc2Source)[0];
            double P_i = mbc2Source_sf->value(mechOutStep_el.t0);    // внутреннее давление
            double P_o = 0;                                         // внешнее давление
            POINT3 v1 = grid->vertex[0];
            int N = gp.N/2+1;
            POINT3 v2 = grid->vertex[(3*N*N-3*N)*(gp.Nparts+1)];
            double a = v1.abs();
            double b = v2.abs();
            MechIterInf nlInf = mechOutStep_el.nlInf.iterInf.back();
            /*
            int globalStep = -1;
            if(mechOutStep_el.t == 0.0000001)
                globalStep = 0;
            if(globalStepNumber == (int)task.step->size() - 1)
                globalStep = 1;
            if(globalStep != -1)
            */
            {
                MechMaterialSource &m0 = (*(*task.mechTask.mechStep)[globalStepNumber].material)[0];
                // расчет численного решения в конечных элементах
                //if(globalStep == 1)
                for(int j0 = 0; j0 < gp.Nparts; j0++)//int j0 = gp.Nparts - 1;
                {
                    double R = 0;
                    double sigma_r = 0;
                    double sigma_fi = 0;
                    int resSize = (int)(*out.mechOut)[globalStepNumber].fe.size();
                    std::vector<ShperePoint> sp(resSize / gp.Nparts * 27);
                    int spInd = 0;
                    for(int j = j0; j < resSize; j += gp.Nparts)// выбераем КЭ на одном уровне по R
                    {
                        for (size_t pInd = 0; pInd < (*out.mechOut)[globalStepNumber].fe[j].pd.size(); pInd++)	// индекс точки Гаусса внутри конечного элемента
                        {
                            MechOutFePointData &r = (*out.mechOut)[globalStepNumber].fe[j].pd[pInd];
                            POINT3 C = r.p;
                            R = C.abs();
                            VECTOR3 ms;
                            if(r.sumSigma*r.sumSigma == 0)
                            {
                                ms[0] = 0;
                                ms[1] = 0;
                                ms[2] = 0;
                            }
                            else
                            {
                                if(sigmaSolvingType == 0)
                                {
                                    // главные напряжения
                                    r.mainStresses(ms);
                                }
                                if(sigmaSolvingType == 1)
                                {
                                    // напряжения на площадке ###
                                    VECTOR3 norm = -r.p;
                                    norm = norm / norm.abs();  // направление к центру
                                    MATR3x3 sigma3x3;
                                    VECTOR3 vectSigma;
                                    r.sumSigma.ToMATR3x3_sigma(sigma3x3); // тензор напряжений
                                    {
                                        vectSigma = sigma3x3*norm;// ms - вектор напряжений на площадке
                                        double absSigma = vectSigma*norm;  // проекция на нормаль к площадке
                                        ms[0] = absSigma;
                                    }
                                    VECTOR3 norm1 = {0, norm[2], -norm[1]};
                                    norm1 = norm1 / norm1.abs();
                                    {
                                        vectSigma = sigma3x3*norm1;// ms - вектор напряжений на площадке
                                        double absSigma = vectSigma*norm1;  // проекция на нормаль к площадке
                                        ms[1] = absSigma;
                                        ms[2] = absSigma;
                                    }
                                }
                            }

                            int ind[3];
                            sortms(ms, ind);

                            sigma_r = ms[ind[0]];
                            sigma_fi = ms[ind[1]];//(ms[ind[1]] + ms[ind[2]]) / 2.;

                            sp[spInd].R = R;
                            sp[spInd].sigma_r = sigma_r;
                            sp[spInd].sigma_fi = sigma_fi;
                            double c = calc_c(a, b, P_i, m0.elasticSigmaLimit);
                            double sigma_r_load;
                            double sigma_fi_load;
                            double sigma_r_unload;
                            double sigma_fi_unload;
                            calcAnalit(a,b,R,P_i,m0.elasticSigmaLimit,c,sigma_r_load,sigma_fi_load,sigma_r_unload,sigma_fi_unload);

                            //if(mechOutStep_el.t == t2)
                            {
                                sp[spInd].sigma_r_analit = sigma_r_load;
                                sp[spInd].sigma_fi_analit = sigma_fi_load;
                                sp[spInd].sigma_r_pogr_abs = -(sp[spInd].sigma_r_analit - sp[spInd].sigma_r);
                                sp[spInd].sigma_fi_pogr_abs = -(sp[spInd].sigma_fi_analit - sp[spInd].sigma_fi);
                                sp[spInd].sigma_r_pogr = -((sp[spInd].sigma_r_analit - sp[spInd].sigma_r)/fabs(sp[spInd].sigma_r_analit));
                                sp[spInd].sigma_fi_pogr = -((sp[spInd].sigma_fi_analit - sp[spInd].sigma_fi)/fabs(sp[spInd].sigma_fi_analit));

                                /*
                                if(fabs(sp[spInd].sigma_fi_pogr_abs) > fi_pogr_max_abs)
                                    fi_pogr_max_abs = fabs(sp[spInd].sigma_fi_pogr_abs);
                                if(fabs(sp[spInd].sigma_r_pogr_abs) > r_pogr_max_abs)
                                    r_pogr_max_abs = fabs(sp[spInd].sigma_r_pogr_abs);
                                if(fabs(sp[spInd].sigma_fi_pogr) > fi_pogr_max)
                                    fi_pogr_max = fabs(sp[spInd].sigma_fi_pogr);
                                if(fabs(sp[spInd].sigma_r_pogr) > r_pogr_max)
                                    r_pogr_max = fabs(sp[spInd].sigma_r_pogr);
                                */
                                fprintf(fsigma_r, "%le %le\n", R, sp[spInd].sigma_r);
                                fprintf(fsigma_fi, "%le %le\n", R, sp[spInd].sigma_fi);
                                //fprintf(fsigma_r_pogr, "%le %le\n", R, fabs(sp[spInd].sigma_r_pogr)*100);
                                //fprintf(fsigma_fi_pogr, "%le %le\n", R, fabs(sp[spInd].sigma_fi_pogr)*100);
                                fprintf(fsigma_eqv, "%le %le\n", R, r.sumSigma.eqv_SIGMA());
                                //fprintf(feps_eqv, "%le %le\n", R, m0.epsEqv(r.sumEpsElastoplastic));

                                // упруго-пластические деформации
                                fprintf(feps_el_eqv, "%le %le\n", R, r.sumEpsElastic.eqv_EPS());
                                // ползучие/пластичные деформации
                                fprintf(feps_pl_eqv, "%le %le\n", R, r.sumEpsPlastic.eqv_EPS());
                            }

                            spInd++;
                        }
                    }
                }

                // расчет аналитическог решения
                //if(globalStep == 1)
                {
                    int NumPoints = 1000;
                    double c = calc_c(a, b, P_i, m0.elasticSigmaLimit);
                    for(int j = 0; j <= NumPoints; j++)
                    {
                        // упругое решение
                        double R = a + (b-a)*j/NumPoints;
                        double sigma_r_load;
                        double sigma_fi_load;
                        double sigma_r_unload;
                        double sigma_fi_unload;
                        calcAnalit(a,b,R,P_i,m0.elasticSigmaLimit,c,sigma_r_load,sigma_fi_load,sigma_r_unload,sigma_fi_unload);
                        //if(mechOutStep_el.t == t2)
                        {
                            //fprintf(fsigma_r, "%le %le\n", P0, m0.elasticSigmaLimit);
                            //fprintf(fsigma_r, "%le %le\n", R, sigma_r);
                            //fprintf(fsigma_fi, "%le %le\n", R, sigma_fi);
                        }
                        // ползучее решение
                        double n = gp.n;
                        //double k = gp.r2/gp.r1;
                        //double _R = R/gp.r1;
                        double k = b/a;
                        double _R = R/a;
                        double kk = pow(k, -3./n);
                        double RR = pow(_R, -3./n);
                        double sigma_r_creep = P_i/(kk - 1)*(RR - kk);
                        double sigma_fi_creep = P_i/(kk - 1)*((1 - 3./(2.*n))*RR - kk);
                        fprintf(fsigma_r_analit, "%le %le %le\n", R, sigma_r_load, sigma_r_creep);
                        fprintf(fsigma_fi_analit, "%le %le %le\n", R, sigma_fi_load, sigma_fi_creep);
                    }
                }
            }

        }
    }

    fclose(fsigma_r_analit);
    fclose(fsigma_r);
    fclose(fsigma_fi_analit);
    fclose(fsigma_fi);
    fclose(fsigma_eqv);
    fclose(feps_el_eqv);
    fclose(feps_pl_eqv);
    // построение графиков
    makeGraphsAfterGlobalStep(globalStepNumber, task);
}
}
