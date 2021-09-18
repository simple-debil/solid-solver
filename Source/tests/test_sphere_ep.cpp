#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"
#include <fstream>

#include "tests.h"
#include "console.h"

#include "test_sphere_ep.h"

using namespace Elementary;
using namespace Solid;
using namespace Grid;

namespace Tests
{
Test_sphere_ep::Test_sphere_ep()
{
    // дирректория данных теста
    dir = "./tests/hollowSphere_ep/";
    // имена gnuplot файлов(+пути к файлам в опциях) и соответствующих картинок
    fn_filesList = dir + "files.txt";
    resGraph.read_gnuplot_from_file(fn_filesList);
    // инициализация графиков с невязками
    initStepsGraphs();
}
Test_base::Type Test_sphere_ep::get_type() const
{
    return Type::Sphere_ep;
}
void Test_sphere_ep::initTask(Task &task)
{
    fn_inf = dir + "_inf.txt";
    FILE *f_inf = fopen(fn_inf.c_str(), "w");
    //bool fidesys = true;
    bool fidesys = false;
    // подрубка
    task.thermTask.enabled = false;
    task.mechTask.enabled = true;
    // общие шаги и сетка
    std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
    Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
    // сетка
    if(fidesys)
    {
        gp.r1 = 2.5;
        gp.r2 = 5;
    }
    else
    {
        gp.r1 = 1;
        gp.r2 = 4;
    }
    gp.N = 8;   // 4 промежутка на четвертьокружностях    //8, 32
    gp.Nparts = 32;
    if(fidesys)
        gp.q = 1. + 0*1./(gp.Nparts/2);//1.0625;//1.0625;//1.03125;
    else
        gp.q = 1. + 1*1./(gp.Nparts/2);//1.0625;//1.0625;//1.03125;
    gp.curvilinear = 0;    // отображение (0 - линейное, 1 - квадратичное)
    gp.buildingMethod = 0; // способ построения (0 - делим дуги, 1 - отображение на куб)
    grid->genSphere(gp);
    double NU = 0.3; //0.4

    // параметры
    //task.slausolver_parameters.preconditioning = SlauSolving::Preconditioning::SlauPreconditioning_LLT;

    sigmaSolvingType = 0;  // способ расчёта решения (0 - главные напряжения, 1 - проекции площадки)
    int fixGrid;                // 0 - подвижная сетка, 1 - зафиксировать сетку

    if(fidesys)
    {
        fixGrid = 0;
    }
    else
    {
        fixGrid = 0;
    }
    int plasticPlot = 1;            // 0 - упругость, 1 - Безье/кусочно-линейная
    //MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::Elasticity;
    MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::D_pl;
    //MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::D_pl_Solov;
    //MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::InitialSigma;
    //MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::Combo_D_pl_InitialSigma;
    double w_midPoint = 1;
    double w_project = 1;

    MechPlasticityCurveUnloadingType PCUnloadingType = MechPlasticityCurveUnloadingType::Elastic;
    //MechPlasticityCurveUnloadingType PCUnloadingType = MechPlasticityCurveUnloadingType::Curve;
    //double k_sigma_eps = -10;
    double k_sigma_eps = 1000.*10000*10.e100;
    double k_eps_sigma = 10;
    double n = 0.0; // 0;

    IncForsesMode incForsesMode = IncForsesMode::IncrementP;
    //IncForsesMode incForsesMode = IncForsesMode::bPlusR;
    //IncForsesMode incForsesMode = IncForsesMode::MinusIntegral;

    MaxSigmaEpsResidual plasticResidualLimit;
    plasticResidualLimit.eps = 1.e-10;   // желаемая невязка по эквиволентным деформациям
    plasticResidualLimit.sigma = 1;//1.e-10;//1.e-10;//1000;

    int nonlinearIterLimit = 50;        // ограничение на количество итераций
    int StepsNumber1 = 1;//8;//8;//128/2;      // шаги
    int StepsNumber2 = 4;//300;//512;//40;//128/2;
    int StepsNumber3 = 1;
    int StepsNumber4 = 1;
    double P = 30;//2.0e7;//2.5e7;//1000;//2.5e7;
    if(fidesys)
        P = 30;
    else
        P = 2.0e7;
    double P1 = 0.1*P;//P/2;//0.4*P;
    double P2 = P;
    double P3 = P*0.1;//P/2.;//P-((P2-P1)/StepsNumber2)/100;//-P*1/2;//-(P2/StepsNumber);
    double P4 = 0;//-P*1/2;//-(P - P2/StepsNumber);


    // нагружение за 1 шаг
    StepsNumber1 = 1;
    StepsNumber2 = 64;
    StepsNumber3 = 1;
    StepsNumber4 = 16;
    P1 = 0;//0.1*P;
    P2 = P;
    P3 = P;//P*0.99; //P*0.9
    P4 = 0;




    //P = 2.e7;
    //epsResidualLimit = 1.e-7;
    //P1 = 0.001*P;
    //P2 = P;
    //P3 = P;
    //P4 = P;
    /*
    double P1 = P*1/2;
    double P2 = P;
    double P3 = P*1/2;
    double P4 = 0;
    */
    int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny

    double Time = 1;
    // шаги
    GlobalStep s;       // шаг
    // первый шаг (1 временной слой, заведомо упругое нагружение P1)
    s.t_start = 0;                          // (имеет значение только для первого глобального шага)
    s.t_finish = Time*1;
    s.dt0 = Time/StepsNumber1;
    step->push_back(s);
    // второй шаг (пластичное нагружение P2)
    s.t_start = Time*1;
    s.t_finish = Time*2;
    s.dt0 = Time/StepsNumber2;
    step->push_back(s);
    // третий шаг (1 временной слой, начало разгрузки)
    s.t_start = Time*2;
    s.t_finish = Time*3;
    s.dt0 = Time/StepsNumber3;
    step->push_back(s);
    // четвертый шаг (1 временной слой, заведомо упругая разгрузка до 0)
    s.t_start = Time*3;
    s.t_finish = Time*4;
    s.dt0 = Time/StepsNumber4;
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
        if(fidesys)
            mechMat.elasticSigmaLimit = 24;//2.e7;
        else
            mechMat.elasticSigmaLimit = 2.e7;//2.e20;//2.e7;
        mechMat.Ceps = 4./9.;
        mechMat.Csigma = 1;

        if(fidesys)
        {
            mechMat.set_E_NU(21000., NU);        //1.e10;
            mechMat.elasticParameters0.Talpha = 0;
        }
        else
        {
            mechMat.set_E_NU(1.e10, NU);
            mechMat.elasticParameters0.Talpha = 0;
        }
        //mechMat.set_K_G(70.0e9, 26.0e9);
        mechMat.set_M_sigma();
        mechMat.F = {0,0,0};
        mechMat.elasticParameters0.ro = 0;//1e10;
        //setPlasticMaterialCurve(plasticPlot, mechMat, k_sigma_eps, k_eps_sigma);
        if(n == 0)
            setPlasticMaterialCurve_Yeld(plasticPlot, mechMat, k_sigma_eps, k_eps_sigma);
        else
            setPlasticMaterialCurve_Yeld_hardening(plasticPlot, mechMat, n);
        mechMat.w_midPoint = w_midPoint;
        mechMat.w_project = w_project;
        mechMat.plasticityMethodType = plasticityMethodType;
        mechMat.PCUnloadingType = PCUnloadingType;
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
        ms.controlMode = 0; // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
        ms.terminateIfAccuracyIsNotAchieving = false;
        ms.iterLimit = nonlinearIterLimit;
        ms.slauResidualLimit = 1000;
        ms.plasticResidualLimit = plasticResidualLimit;
        ms.material = material;
        ms.bc1Source = bc1Source;
        //ms.bc2Source = bc2Source; - будет меняться
        ms.bc2 = bc2;
        // первый шаг (1 временной слой, заведомо упругое нагружение P1)
        s = (*step)[0];
        setBc2Sphere(bc2Source, 0, P1, s.t_start, s.t_finish);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);
        // второй шаг (пластичное нагружение P2)
        s = (*step)[1];
        setBc2Sphere(bc2Source, P1, P2, s.t_start, s.t_finish);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);
        // третий шаг (1 временной слой, начало разгрузки)
        s = (*step)[2];
        setBc2Sphere(bc2Source, P2, P3, s.t_start, s.t_finish);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);
        // четвертый шаг (1 временной слой, заведомо упругая разгрузка до 0)
        s = (*step)[3];
        setBc2Sphere(bc2Source, P3, P4, s.t_start, s.t_finish);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);

        task.mechTask.mechStep = mechStep;

        // Контакт (отсутствует)
        // жёсткие поверхности
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

    // сохренение информации о сетке в файл
    fprintf(f_inf, "r1 = %lf, r2 = %lf\n", gp.r1, gp.r2);
    fprintf(f_inf, "Grid: %dx%dx%d, q = %le\n", gp.N, gp.N, gp.Nparts, gp.q);
    if(gp.curvilinear == 0)
        fprintf(f_inf, "Linear\n");
    if(gp.curvilinear == 1)
        fprintf(f_inf, "Quadratic\n");
    if(gp.buildingMethod == 0)
        fprintf(f_inf, "buildingMethod = Divide arcs\n");
    if(gp.buildingMethod == 1)
        fprintf(f_inf, "Quadratic = Cube to sphere\n");
    fprintf(f_inf, "FENumber = %d\n", (int)grid->fe.size());
    fprintf(f_inf, "VertexesNumber = %d\n", (int)grid->vertex.size());
    fprintf(f_inf, "VertexesForCurvatureNumber = %d\n", (int)grid->vertexForCurvature.size());
    fprintf(f_inf, "bc1Number = %d\n", (int)grid->bc1.size());
    fprintf(f_inf, "FEsurface = %d\n", (int)grid->FESurface.size());
    for (size_t FEsurfaceInd = 0; FEsurfaceInd < grid->FESurface.size(); FEsurfaceInd++)
    {
        fprintf(f_inf, "FEsurface[%d]Size = %d\n", (int)FEsurfaceInd, (int)grid->FESurface[FEsurfaceInd].face.size());
    }
    //grid->buldOpenSCADModel("setka.scad");
    // сохренение информации о параметрах в файл
    fprintf(f_inf, "fixGrid = %d\n", fixGrid);
    fprintf(f_inf, "plasticPlot = %d\n", plasticPlot);
    fprintf(f_inf, "epsResidualLimit = %le\n", plasticResidualLimit.eps);
    fprintf(f_inf, "nonlinearIterLimit = %d\n", nonlinearIterLimit);
    fprintf(f_inf, "StepsNumber1 = %d\n", StepsNumber1);
    fprintf(f_inf, "StepsNumber2 = %d\n", StepsNumber2);
    fprintf(f_inf, "StepsNumber3 = %d\n", StepsNumber3);
    fprintf(f_inf, "StepsNumber4 = %d\n", StepsNumber4);
    fprintf(f_inf, "P = %lf\n", P);
    fprintf(f_inf, "P1 = %lf\n", P1);
    fprintf(f_inf, "P2 = %lf\n", P2);
    fprintf(f_inf, "P3 = %lf\n", P3);
    fprintf(f_inf, "P4 = %lf\n", P4);
    fprintf(f_inf, "incForsesMode = %d\n", (int)incForsesMode);
    fprintf(f_inf, "HomogenyMode = %d\n", HomogenyMode);

    fclose(f_inf);

}
void Test_sphere_ep::writeResults(const Task &task, const OutData &out)
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
    fn_inf = subdir + "_inf.txt";

    fn_curve = dir + "curve.txt";
    fn_dcurve = dir + "dcurve.txt";

    fn_sigma_r_analit = subdir + "sigma_r_analit.txt";
    fn_sigma_r = subdir + "sigma_r.txt";
    fn_sigma_fi_analit = subdir + "sigma_fi_analit.txt";
    fn_sigma_fi = subdir + "sigma_fi.txt";

    fn_sigma_r_pogr_abs = subdir + "sigma_r_pogr_abs.txt";
    fn_sigma_r_pogr = subdir + "sigma_r_pogr.txt";
    fn_sigma_fi_pogr_abs = subdir + "sigma_fi_pogr_abs.txt";
    fn_sigma_fi_pogr = subdir + "sigma_fi_pogr.txt";

    fn_res_maxpogr = subdir + "_maxpogr.txt";
    fn_s_nonlinearStateFENumber = subdir + "_nonlinearStateFENumber.txt";
    // пути к файлам для отображения графиков
{
    std::ofstream f;
    if(globalStepNumber == 0)
    {
        f.open(fn_filesList, std::ofstream::out);
        resGraph.generalGraphsNomber = 4;
        resGraph.graphsPerGlobalStepNomber = 2;
        f << resGraph.generalGraphsNomber << "\n";
        f << resGraph.graphsPerGlobalStepNomber << "\n";
        // curve
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "curve.png") +
             genGnuplotCommandLineParameter("fn_curve", dir + "curve.txt") +
             "\" " +
             dir + "curve.gnu" << "\n";
        f << dir + "curve.png" << "\n";
        f << "curve.png" << "\n";
        // dcurve
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "dcurve.png") +
             genGnuplotCommandLineParameter("fn_dcurve", dir + "dcurve.txt") +
             "\" " +
             dir + "dcurve.gnu" << "\n";
        f << dir + "dcurve.png" << "\n";
        f << "dcurve.png" << "\n";
        // sigma_fi_13
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "sigma_fi_13.png") +
             genGnuplotCommandLineParameter("fn_load_sigma_fi", dir + "1/sigma_fi.txt") +
             genGnuplotCommandLineParameter("fn_unload_sigma_fi", dir + "3/sigma_fi.txt") +
             genGnuplotCommandLineParameter("fn_load_sigma_fi_analit", dir + "1/sigma_fi_analit.txt") +
             genGnuplotCommandLineParameter("fn_unload_sigma_fi_analit", dir + "3/sigma_fi_analit.txt") +
             "\" " +
             dir + "sigma_fi_13.gnu" << "\n";
        f << dir + "sigma_fi_13.png" << "\n";
        f << "sigma_fi_13.png" << "\n";
        // sigma_r_13
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "sigma_r_13.png") +
             genGnuplotCommandLineParameter("fn_load_sigma_r", dir + "1/sigma_r.txt") +
             genGnuplotCommandLineParameter("fn_unload_sigma_r", dir + "3/sigma_r.txt") +
             genGnuplotCommandLineParameter("fn_load_sigma_r_analit", dir + "1/sigma_r_analit.txt") +
             genGnuplotCommandLineParameter("fn_unload_sigma_r_analit", dir + "3/sigma_r_analit.txt") +
             "\" " +
             dir + "sigma_r_13.gnu" << "\n";
        f << dir + "sigma_r_13.png" << "\n";
        f << "sigma_r_13.png" << "\n";
    }
    else
    {
        f.open(fn_filesList, std::ofstream::app);
    }
    f << std::to_string(globalStepNumber) + "\n";
    // sigma_fi_pogr
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "sigma_fi_pogr.png") +
         genGnuplotCommandLineParameter("fn_sigma_fi_pogr", subdir + "sigma_fi_pogr.txt") +
         "\" " +
         dir + "sigma_fi_pogr.gnu" << "\n";
    f << subdir + "sigma_fi_pogr.png" << "\n";
    f << "sigma_fi_pogr.png" << "\n";
    // sigma_r_pogr
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "sigma_r_pogr.png") +
         genGnuplotCommandLineParameter("fn_sigma_r_pogr", subdir + "sigma_r_pogr.txt") +
         "\" " +
         dir + "sigma_r_pogr.gnu" << "\n";
    f << subdir + "sigma_r_pogr.png" << "\n";
    f << "sigma_r_pogr.png" << "\n";
    f.close();
}
    resGraph.read_gnuplot_from_file(fn_filesList);
    }

    FILE *fsigma_r_analit = fopen(fn_sigma_r_analit.c_str(), "w");
    FILE *fsigma_r = fopen(fn_sigma_r.c_str(), "w");
    FILE *fsigma_fi_analit = fopen(fn_sigma_fi_analit.c_str(), "w");
    FILE *fsigma_fi = fopen(fn_sigma_fi.c_str(), "w");

    FILE *fsigma_r_pogr_abs = fopen(fn_sigma_r_pogr_abs.c_str(), "w");
    FILE *fsigma_fi_pogr_abs = fopen(fn_sigma_fi_pogr_abs.c_str(), "w");
    FILE *fsigma_r_pogr = fopen(fn_sigma_r_pogr.c_str(), "w");
    FILE *fsigma_fi_pogr = fopen(fn_sigma_fi_pogr.c_str(), "w");

    FILE *fres_maxpogr = fopen(fn_res_maxpogr.c_str(), "w");
    FILE *fs_nonlinearStateFENumber = fopen(fn_s_nonlinearStateFENumber.c_str(), "w");

    Grid::Grid3D *grid = task.grid;

    // кривая пластичности
    if(globalStepNumber == 0)
    {
        FILE *fcurve = fopen(fn_curve.c_str(), "w");
        FILE *fdcurve = fopen(fn_dcurve.c_str(), "w");
        MechMaterialSource &m0 = (*(*task.mechTask.mechStep)[globalStepNumber].material)[0];
        double eps1 = 0;
        double eps2 = 0.004; // strToDouble(m0.vals[1]) - предел упругости
        int N = 1000;//1000000;
        for(int i = 0; i <= N; i++)
        {
            double eps = eps1+(eps2-eps1)*i/N;
            double sigma = m0.sigma(m0.elasticParameters0, 0, eps, 0);
            double deps = m0.difSigma(m0.elasticParameters0, 0, eps, 0);
            fprintf(fcurve, "%le %le\n", eps, sigma);
            fprintf(fdcurve, "%le %le\n", eps, deps);
        }
        fclose(fcurve);
        fclose(fdcurve);
    }

    {
        double fi_pogr_max_abs = 0;
        double r_pogr_max_abs = 0;
        double sigma_r_max = 0;
        double sigma_fi_max = 0;
        double fi_pogr_max = 0;
        double r_pogr_max = 0;

        // вывод данных МДТТ
        if(task.mechTask.enabled)
        {
            int numLocalStep = (int)(*out.mechOut)[globalStepNumber].step.size() - 1;
            MechOutStepData &mechOutStep_el = (*out.mechOut)[globalStepNumber].step[numLocalStep];
            // количество КЭ участвующих в итерационном процессе
            MechIterInf nlInf = mechOutStep_el.nlInf.iterInf.back();
            fprintf(fs_nonlinearStateFENumber, "%le %d\n", mechOutStep_el.t0, nlInf.plastic.nonlinearStateFENumber);

            //MechBoundaryCondition2Source_ScalarFunction *mbc2Source_sf = (MechBoundaryCondition2Source_ScalarFunction *)(*(*task.mechTask.mechStep)[globalStepNumber].bc2Source)[0];
            //double P0 = mbc2Source_sf->value(mechOutStep_el.t);
            bool load = globalStepNumber == 1; // нагружение
            bool unload = globalStepNumber == 3; // разгрузка
            if(load || unload)
            {
                // давление после нагружения
                double P0 = 0;
                {
                    int gInd = 1;// индекс глобального шага, в результате которого тело нагружается
                    int lInd = (int)(*out.mechOut)[gInd].step.size() - 1;
                    double t = (*out.mechOut)[gInd].step[lInd].t0;
                    MechBoundaryCondition2Source_ScalarFunction *mbc2Source_sf = (MechBoundaryCondition2Source_ScalarFunction *)(*(*task.mechTask.mechStep)[gInd].bc2Source)[0];//####
                    P0 = mbc2Source_sf->value(t);
                }
                // радиусы
                POINT3 v1 = grid->vertex[0];
                int N = gp.N/2+1;
                POINT3 v2 = grid->vertex[(3*N*N-3*N)*(gp.Nparts+1)];
                double a = v1.abs();
                double b = v2.abs();

                MechMaterialSource &m0 = (*(*task.mechTask.mechStep)[globalStepNumber].material)[0];
                // расчет численного решения в конечных элементах
                for(int j0 = 0; j0 < gp.Nparts; j0++)//int j0 = gp.Nparts - 1;
                {
                    int resSize = (int)(*out.mechOut)[globalStepNumber].fe.size();
                    std::vector<ShperePoint> sp(resSize / gp.Nparts * 27);
                    int spInd = 0;
                    for(int j = j0; j < resSize; j += gp.Nparts)// выбераем КЭ на одном уровне по R
                    {
                        for (size_t pInd = 0; pInd < (*out.mechOut)[globalStepNumber].fe[j].pd.size(); pInd++)	// индекс точки Гаусса внутри конечного элемента
                        {
                            MechOutFePointData &r = (*out.mechOut)[globalStepNumber].fe[j].pd[pInd];
                            POINT3 C = r.p;
                            double R = C.abs();
                            double sigma_r = 0;
                            double sigma_fi = 0;

                            VECTOR3 ms;
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
                                /*VECTOR3 norm2 = norm*norm1;
                                {
                                    vectSigma = sigma3x3*norm2;// ms - вектор напряжений на площадке
                                    double absSigma = vectSigma*norm2;  // проекция на нормаль к площадке
                                    ms[2] = absSigma;
                                }*/
                            }
                            int ind[3];
                            sortms(ms, ind);
                            sigma_r = ms[ind[0]];
                            sigma_fi = ms[ind[1]];//(ms[ind[1]] + ms[ind[2]]) / 2.;
                            sp[spInd].R = R;
                            sp[spInd].sigma_r = sigma_r;
                            sp[spInd].sigma_fi = sigma_fi;
                            // аналитическое решение
                            {
                                double c = calc_c(a, b, P0, m0.elasticSigmaLimit);
                                double sigma_r_load;
                                double sigma_fi_load;
                                double sigma_r_unload;
                                double sigma_fi_unload;
                                calcAnalit(a,b,R,P0,m0.elasticSigmaLimit,c,sigma_r_load,sigma_fi_load,sigma_r_unload,sigma_fi_unload);
                                if(load)
                                {
                                    sp[spInd].sigma_r_analit = sigma_r_load;
                                    sp[spInd].sigma_fi_analit = sigma_fi_load;
                                }
                                if(unload)
                                {
                                    sp[spInd].sigma_r_analit = sigma_r_unload;
                                    sp[spInd].sigma_fi_analit = sigma_fi_unload;
                                }
                            }
                            sp[spInd].sigma_r_pogr_abs = -(sp[spInd].sigma_r_analit - sp[spInd].sigma_r);
                            sp[spInd].sigma_fi_pogr_abs = -(sp[spInd].sigma_fi_analit - sp[spInd].sigma_fi);
                            sp[spInd].sigma_r_pogr = -((sp[spInd].sigma_r_analit - sp[spInd].sigma_r)/fabs(sp[spInd].sigma_r_analit));
                            sp[spInd].sigma_fi_pogr = -((sp[spInd].sigma_fi_analit - sp[spInd].sigma_fi)/fabs(sp[spInd].sigma_fi_analit));
                            if(fabs(sp[spInd].sigma_fi_pogr_abs) > fi_pogr_max_abs)
                                fi_pogr_max_abs = fabs(sp[spInd].sigma_fi_pogr_abs);
                            if(fabs(sp[spInd].sigma_r_pogr_abs) > r_pogr_max_abs)
                                r_pogr_max_abs = fabs(sp[spInd].sigma_r_pogr_abs);
                            if(fabs(sp[spInd].sigma_fi_pogr) > fi_pogr_max)
                                fi_pogr_max = fabs(sp[spInd].sigma_fi_pogr);
                            if(fabs(sp[spInd].sigma_r_pogr) > r_pogr_max)
                                r_pogr_max = fabs(sp[spInd].sigma_r_pogr);
                            fprintf(fsigma_r, "%le %le\n", R, sp[spInd].sigma_r);
                            fprintf(fsigma_fi, "%le %le\n", R, sp[spInd].sigma_fi);
                            fprintf(fsigma_r_pogr_abs, "%le %le\n", R, sp[spInd].sigma_r_pogr_abs*100);
                            fprintf(fsigma_fi_pogr_abs, "%le %le\n", R, sp[spInd].sigma_fi_pogr_abs*100);
                            fprintf(fsigma_r_pogr, "%le %le\n", R, fabs(sp[spInd].sigma_r_pogr)*100);
                            fprintf(fsigma_fi_pogr, "%le %le\n", R, fabs(sp[spInd].sigma_fi_pogr)*100);
                            spInd++;
                        }
                    }
                }

                // расчет аналитического решения
                int NumPoints = 1000;
                double c = calc_c(a, b, P0, m0.elasticSigmaLimit);
                for(int j = 0; j <= NumPoints; j++)
                {
                    double R = a + (b-a)*j/NumPoints;
                    double sigma_r;
                    double sigma_fi;
                    {
                        double sigma_r_load;
                        double sigma_fi_load;
                        double sigma_r_unload;
                        double sigma_fi_unload;
                        calcAnalit(a,b,R,P0,m0.elasticSigmaLimit,c,sigma_r_load,sigma_fi_load,sigma_r_unload,sigma_fi_unload);
                        if(load)
                        {
                            sigma_r = sigma_r_load;
                            sigma_fi = sigma_fi_load;
                        }
                        if(unload)
                        {
                            sigma_r = sigma_r_unload;
                            sigma_fi = sigma_fi_unload;
                        }
                    }
                    fprintf(fsigma_r_analit, "%le %le\n", R, sigma_r);
                    fprintf(fsigma_fi_analit, "%le %le\n", R, sigma_fi);
                    if(fabs(sigma_r) > sigma_r_max)
                        sigma_r_max = fabs(sigma_r);
                    if(fabs(sigma_fi) > sigma_fi_max)
                        sigma_fi_max = fabs(sigma_fi);
                }
            }
        }
        // максимальные погрешности
        fprintf(fres_maxpogr, "r_pogr_max_abs = %le\n", r_pogr_max_abs/sigma_r_max*100);
        fprintf(fres_maxpogr, "fi_pogr_max_abs = %le\n", fi_pogr_max_abs/sigma_fi_max*100);
        fprintf(fres_maxpogr, "r_pogr_max = %le\n", r_pogr_max*100);
        fprintf(fres_maxpogr, "fi_pogr_max = %le\n", fi_pogr_max*100);
    }

    fclose(fsigma_r_analit);
    fclose(fsigma_r);
    fclose(fsigma_fi_analit);
    fclose(fsigma_fi);

    fclose(fsigma_r_pogr_abs);
    fclose(fsigma_fi_pogr_abs);
    fclose(fsigma_r_pogr);
    fclose(fsigma_fi_pogr);

    fclose(fs_nonlinearStateFENumber);
    fclose(fres_maxpogr);
    // построение графиков
    makeGraphsAfterGlobalStep(globalStepNumber, task);
}
}
