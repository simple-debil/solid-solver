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
void Test_creep::setBc2Creep(double Px1, double Px2,
                 double Py1, double Py2,
                 double t1, double t2, std::vector<Solid::MechBoundaryCondition2Source_base *> *&bc2Source)
{
    bc2Source = new std::vector<MechBoundaryCondition2Source_base *>(2);
    // линейная функция
    setLinearBc2_el(Px1, Px2, t1, t2,
                    (*bc2Source)[0]);
    setLinearBc2_el(Py1, Py2, t1, t2,
                    (*bc2Source)[1]);
}

Test_creep::Test_creep()
{
    // дирректория данных теста
    dir = "./tests/creep/";
    // имена gnuplot файлов(+пути к файлам в опциях) и соответствующих картинок
    fn_filesList = dir + "files.txt";
    resGraph.read_gnuplot_from_file(fn_filesList);
    // инициализация графиков с невязками
    initStepsGraphs();
}
Test_base::Type Test_creep::get_type() const
{
    return Type::Creep;
}
void Test_creep::initTask(Solid::Task &task)
{
    //nonlinearIterLimit = 2
    //incForsesMode = IncForsesMode::IncrementP;
    //mechMat.plasticityMethodType = MechPlasticityMethodType::CreepInitialEps;
    //stepsNumber[2] = 4;
    //mechMat.w = 1;
    // подрубка
    task.mechTask.enabled = true;
    task.thermTask.enabled = false;
    // общие шаги и сетка
    std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
    Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
    // сетка

    // X
    cp.cube.i[0][0] = -2;
    cp.cube.i[0][1] = +2;
    // Y
    cp.cube.i[1][0] = -2;
    cp.cube.i[1][1] = +2;
    // Z
    cp.cube.i[2][0] = -2;
    cp.cube.i[2][1] = +2;

    cp.N[0] = 8;
    cp.N[1] = 8;
    cp.N[2] = 8;

    cp.E = 200000;//1.e8;//1.e5;//200*1.e3*100*100;
    cp.Nu = 0.3;

    cp.sigma0 = 200;//1.e4;    //double P = 100;//200*100*100 * 1.e0;     // давление
    cp.alpha = -1;//-0.1;      //double alpha = -1;
    cp.A = 3.125*1.e-14;
    cp.n = 5;
    cp.m = 0.5;

    double creepResidualLimit = 1.e-10;           // желаемая относительная погрешность по деформациям ползучести
    int nonlinearIterLimit;//2                // ограничение на количество итераций

    grid->genCreep(cp);
    //grid->buldOpenSCADModel("setka.scad");
    // параметры
    int fixGrid = 0;                // 0 - подвижная сетка, 1 - зафиксировать сетку
    int plasticityCurveMode = 1;    // 0 - eps(sigma), 1 - sigma(eps)
    bool usingSecantMethod = false;
    MaxSigmaEpsResidual plasticResidualLimit;
    plasticResidualLimit.eps = 1.e-10;            // желаемая относительная погрешность по эквиволентным деформациям
    plasticResidualLimit.sigma = 10000;

    std::vector<int> stepsNumber_data(5 + 1);
    int *stepsNumber = &(stepsNumber_data[1]);   // - массив из 5 элементов, нумерация начинается с -1

    std::vector<double> Time_data(5 + 1);
    double *Time = &(Time_data[1]);   // - массив из 5 элементов, нумерация начинается с -1
    // ползучесть и релаксация
/*
    cp.timeCreep = 1000;   // 1000/1;
    cp.timeRelaxation1 = 10;
    cp.timeRelaxation2 = 1000 - cp.timeRelaxation1;
    stepsNumber[0] = 1;                        // пустой шаг
    stepsNumber[1] = 1;                        // количество шагов нагружения
    stepsNumber[2] = 2;// количество шагов ползучести
    stepsNumber[3] = 100;// количество шагов разгрузки1
    stepsNumber[4] = 100;// количество шагов разгрузки2

    Time[-1] = -0.000000002;          // начало процесса
    Time[0] = -0.000000001;           // время завершения пустого шага
    Time[1] = 0;                      // время завершения нагружения
    Time[2] = Time[1] + cp.timeCreep;      // время завершения ползучести
    Time[3] = Time[2] + cp.timeRelaxation1;// время разгрузки1
    Time[4] = Time[3] + cp.timeRelaxation2;// время разгрузки2
*/

    // только ползучесть
// /*
    // для этого тесто нужно поставить w = 0##
    nonlinearIterLimit = 2;
    cp.timeCreep = 1000;   // 1000/1;
    cp.timeRelaxation1 = 0;
    cp.timeRelaxation2 = 0;
    stepsNumber[0] = 1;                        // пустой шаг
    stepsNumber[1] = 1;                        // количество шагов нагружения
    stepsNumber[2] = 4;// количество шагов ползучести
    stepsNumber[3] = 0;// количество шагов разгрузки1
    stepsNumber[4] = 0;// количество шагов разгрузки2

    Time[-1] = -0.000000002;          // начало процесса
    Time[0] = -0.000000001;           // время завершения пустого шага
    Time[1] = 0;                      // время завершения нагружения
    Time[2] = Time[1] + cp.timeCreep;      // время завершения ползучести
    Time[3] = Time[2] + cp.timeRelaxation1;// время разгрузки1
    Time[4] = Time[3] + cp.timeRelaxation2;// время разгрузки2
// */
    // только релаксация
/*
    nonlinearIterLimit = 2;
    cp.timeCreep = 0;

    //cp.timeRelaxation1 = 0.01;
    //cp.timeRelaxation2 = 1000 - cp.timeRelaxation1;
    //stepsNumber[0] = 1;                        // пустой шаг
    //stepsNumber[1] = 1;                        // количество шагов нагружения
    //stepsNumber[2] = 0; // количество шагов ползучести
    //stepsNumber[3] = 100;// количество шагов разгрузки1
    //stepsNumber[4] = 100;// количество шагов разгрузки2

    cp.timeRelaxation1 = 0.01;
    cp.timeRelaxation2 = 10 - cp.timeRelaxation1;
    stepsNumber[0] = 1;                        // пустой шаг
    stepsNumber[1] = 1;                        // количество шагов нагружения
    stepsNumber[2] = 0; // количество шагов ползучести
    stepsNumber[3] = 50;// количество шагов разгрузки1
    stepsNumber[4] = 50;// количество шагов разгрузки2



    Time[-1] = -0.000000002;          // начало процесса
    Time[0] = -0.000000001;           // время завершения пустого шага
    Time[1] = 0;                      // время завершения нагружения
    Time[2] = Time[1] + cp.timeCreep;    // время завершения ползучести
    Time[3] = Time[2] + cp.timeRelaxation1; // время разгрузки1
    Time[4] = Time[3] + cp.timeRelaxation2; // время разгрузки1
*/


    double Px1 = 0;
    double Px2 = -cp.sigma0;
    double Py1 = 0;
    double Py2 = -cp.alpha*cp.sigma0;
    double elasticSigmaLimit = 3.e100;
    IncForsesMode incForsesMode = IncForsesMode::IncrementP;                         // 1 - dSigma = integral(...), 2 - в приращениях (P->dP)

    // шаги
    GlobalStep s;       // шаг

    for(int globalStep = 0; globalStep < 5; globalStep++)
    {
        if(globalStep == 3)
        {
            // не равномерный шаг по времени
            double t1 = Time[globalStep - 1];
            double t2 = Time[globalStep];
            double N = stepsNumber[globalStep];
            double q = 1.1;
            for(int i = 0; i < N; i++)
            {
                Elementary::Operations::findPointOnTheLine_1d(t1, t2, N, q, i, s.t_start);
                Elementary::Operations::findPointOnTheLine_1d(t1, t2, N, q, i + 1, s.t_finish);
                s.dt0 = s.t_finish - s.t_start;
                step->push_back(s);
            }
        }else
        if(globalStep == 4)
        {
            // не равномерный шаг по времени
            double t1 = Time[globalStep - 1];
            double t2 = Time[globalStep];
            double N = stepsNumber[globalStep];
            double q = 1.1;
            for(int i = 0; i < N; i++)
            {
                Elementary::Operations::findPointOnTheLine_1d(t1, t2, N, q, i, s.t_start);
                Elementary::Operations::findPointOnTheLine_1d(t1, t2, N, q, i + 1, s.t_finish);
                s.dt0 = s.t_finish - s.t_start;
                step->push_back(s);
            }
        }else
        {
            // равномерный шаг по времени
            for(int i = 0; i < stepsNumber[globalStep]; i++)
            {
                s.t_start = Time[globalStep - 1] + (Time[globalStep] - Time[globalStep - 1])*i/stepsNumber[globalStep];
                s.t_finish = Time[globalStep - 1] + (Time[globalStep] - Time[globalStep - 1])*(i + 1)/stepsNumber[globalStep];
                s.dt0 = (Time[globalStep] - Time[globalStep - 1])/stepsNumber[globalStep];
                step->push_back(s);
            }
        }
    }

    /*
    // пустой шаг
    for(int i = 0; i < stepsNumber0; i++)
    {
        s.t1 = Time_start + Time0*i/stepsNumber0;
        s.t2 = Time_start + Time0*(i + 1)/stepsNumber0;
        s.dt0 = Time0/stepsNumber0;
        step->push_back(s);
    }
    // нагружение
    for(int i = 0; i < stepsNumber1; i++)
    {
        s.t1 = Time_start + Time0 + Time1*i/stepsNumber1;
        s.t2 = Time_start + Time0 + Time1*(i + 1)/stepsNumber1;
        s.dt0 = Time1/stepsNumber1;
        step->push_back(s);
    }
    // ползучесть
    for(int i = 0; i < stepsNumber2; i++)
    {
        s.t1 = Time_start + Time0 + Time1 + Time2*i/stepsNumber2;
        s.t2 = Time_start + Time0 + Time1 + Time2*(i + 1)/stepsNumber2;
        s.dt0 = Time2/stepsNumber2;
        step->push_back(s);
    }
    // разгрузка
    for(int i = 0; i < stepsNumber3; i++)
    {
        s.t1 = Time_start + Time0 + Time1 + Time2 + Time3*i/stepsNumber2;
        s.t2 = Time_start + Time0 + Time1 + Time2*(i + 1)/stepsNumber2;
        s.dt0 = Time2/stepsNumber2;
        step->push_back(s);
    }*/
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
        mechMat.elasticSigmaLimit = elasticSigmaLimit;
        mechMat.Ceps = 4./9.;
        mechMat.Csigma = 1;
        mechMat.set_E_NU(cp.E, cp.Nu);
        mechMat.set_M_sigma();
        mechMat.F = {0,0,0};// объёмные силы
        mechMat.elasticParameters0.ro = 0;     // плотность
            mechMat.elasticParameters0.Talpha = 0*1.e-5; // Коэффициент линейного расширения
        setCreepMaterialCurve(cp.A, cp.n, cp.m, mechMat);
        //mechMat.plasticityMethodType = MechPlasticityMethodType::InitialEps;
        //mechMat.plasticityMethodType = MechPlasticityMethodType::Plasticity;
        //mechMat.plasticityMethodType = MechPlasticityMethodType::Combo_D_pl_InitialSigma;
        mechMat.plasticityMethodType = MechPlasticityMethodType::InitialSigma;
        //### убрал InitialEps
        //### кривые не для параметра упрочнения Одквиста
        mechMat.PCDependenceType = MechPlasticityCurveDependenceType::Time;
        mechMat.PCUnloadingType = MechPlasticityCurveUnloadingType::Curve;
        mechMat.PCType = MechPlasticityCurveType::Eps_sigma;
        mechMat.w_midPoint = 0;
        mechMat.w_project = 1;
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
        MechGlobalStepParameters ms;  // шаг для МДТТ
        GlobalStep s;       // шаг
        ms.slausolverParameters = slausolver_parameters;
        ms.timeMode = 0;    // 0 - квазистатическая задача
        ms.incForsesMode = incForsesMode; // 1 - dSigma = integral(...), 2 - приращения
        ms.fixGrid = fixGrid;     // 1 - зафиксировать сетку
        ms.switchIterationsMode = SwitchIterationsMode::Serial;
        ms.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
        ms.controlMode = 0;     // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
        ms.terminateIfAccuracyIsNotAchieving = false;
        ms.iterLimit = nonlinearIterLimit;
        ms.slauResidualLimit = 10000;
        ms.plasticResidualLimit = plasticResidualLimit;
        ms.contactEndPointResidualLimit = 1000;//1.e-10;//1.e-13;
        ms.material = material;
        ms.bc1Source = bc1Source;
        //ms.bc2Source = bc2Source; - будет меняться
        ms.bc2 = bc2;
        int stepCount = 0;
        // пустой шаг
        for(int i = 0; i < stepsNumber[0]; i++)
        {
            s = (*step)[stepCount];
            setBc2Creep(0, 0,
                        0, 0,
                        s.t_start, s.t_finish, bc2Source);
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
            stepCount++;
        }
        // нагружение
        for(int i = 0; i < stepsNumber[1]; i++)
        {
            s = (*step)[stepCount];
            //s.t1 = Time*0 + Time*i/stepsNumber;
            //  s.t2 = Time*0 + Time*(i + 1)/stepsNumber;
            setBc2Creep(Px1 + (Px2 - Px1)*i/stepsNumber[1], Px1 + (Px2 - Px1)*(i + 1)/stepsNumber[1],
                        Py1 + (Py2 - Py1)*i/stepsNumber[1], Py1 + (Py2 - Py1)*(i + 1)/stepsNumber[1],
                        s.t_start, s.t_finish, bc2Source);
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
            stepCount++;
        }
        // ползучесть
        for(int i = 0; i < stepsNumber[2]; i++)
        {
            s = (*step)[stepCount];
            setBc2Creep(Px2, Px2,
                        Py2, Py2,
                        s.t_start, s.t_finish, bc2Source);
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
            stepCount++;
        }
        // первые краевые условия для разгрузки
        bc1Source = new std::vector<MechBoundaryCondition1Source>(5);
        (*bc1Source)[0].mode = {{ 0, -1, -1}};
        (*bc1Source)[0].u0 =   {{ 0, -1, -1}};
        (*bc1Source)[1].mode = {{-1,  0, -1}};
        (*bc1Source)[1].u0 =   {{-1,  0, -1}};
        (*bc1Source)[2].mode = {{-1, -1,  0}};
        (*bc1Source)[2].u0 =   {{-1, -1,  0}};
        (*bc1Source)[3].mode = {{-1, -1, -1}};
        (*bc1Source)[3].u0 =   {{-1, -1, -1}};
        (*bc1Source)[4].mode = {{0, 0, 0}};
        (*bc1Source)[4].u0 =   {{0, 0, 0}};
        // разгрузка1
        for(int i = 0; i < stepsNumber[3]; i++)
        {
            s = (*step)[stepCount];
            setBc2Creep(Px2, Px2,
                        Py2, Py2,
                        s.t_start, s.t_finish, bc2Source);
            ms.bc1Source = bc1Source;
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
            stepCount++;
        }
        // разгрузка2
        for(int i = 0; i < stepsNumber[4]; i++)
        {
            s = (*step)[stepCount];
            setBc2Creep(Px2, Px2,
                        Py2, Py2,
                        s.t_start, s.t_finish, bc2Source);
            ms.bc1Source = bc1Source;
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
            stepCount++;
        }

        task.mechTask.mechStep = mechStep;

        // Контакт (отсутствует)

        // поверхности
        task.mechTask.rigidSurface = new std::vector<Surface_base *>;
        // контакты КЭ-поверхность - жёсткая поверхность
        task.mechTask.CC_FE_Rigid = new std::vector<ContactCondition_FE_Rigid>;
        // инициализация нулями начальных данных
        task.mechTask.initNull(Integration::IntegrationType::Gauss3);
    }

    // Температура
    if(task.thermTask.enabled)
    {
    }
}
void Test_creep::writeResults(const Task &task, const Solid::OutData &out)
{
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
    fn_eps11 = subdir + "eps11.txt";
    fn_eps22 = subdir + "eps22.txt";
    fn_eps33 = subdir + "eps33.txt";
    fn_eps11_analit = subdir + "eps11_analit.txt";
    fn_eps22_analit = subdir + "eps22_analit.txt";
    fn_eps33_analit = subdir + "eps33_analit.txt";
    fn_sigmaEqv = subdir + "sigmaEqv.txt";
    // пути к файлам для отображения графиков
{
    std::ofstream f;
    if(globalStepNumber == 0)
    {
        f.open(fn_filesList, std::ofstream::out);
        resGraph.generalGraphsNomber = 0;
        resGraph.graphsPerGlobalStepNomber = 2;
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
    // eps1122
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "eps1122.png") +
         genGnuplotCommandLineParameter("fn_eps11", subdir + "eps11.txt") +
         genGnuplotCommandLineParameter("fn_eps11_analit", subdir + "eps11_analit.txt") +
         genGnuplotCommandLineParameter("fn_eps22", subdir + "eps22.txt") +
         genGnuplotCommandLineParameter("fn_eps22_analit", subdir + "eps22_analit.txt") +
         genGnuplotCommandLineParameter("fn_eps33", subdir + "eps33.txt") +
         genGnuplotCommandLineParameter("fn_eps33_analit", subdir + "eps33_analit.txt") +
         "\" " +
         dir + "eps1122.gnu" << "\n";
    f << subdir + "eps1122.png" << "\n";
    f << "eps1122.png" << "\n";
    // sigmaEqv
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "sigmaEqv.png") +
         genGnuplotCommandLineParameter("fn_sigmaEqv", subdir + "sigmaEqv.txt") +
         "\" " +
         dir + "sigmaEqv.gnu" << "\n";
    f << subdir + "sigmaEqv.png" << "\n";
    f << "sigmaEqv.png" << "\n";
    f.close();
}
    resGraph.read_gnuplot_from_file(fn_filesList);
    }

    FILE *f_eps11 = fopen(fn_eps11.c_str(), "w");
    FILE *f_eps22 = fopen(fn_eps22.c_str(), "w");
    FILE *f_eps33 = fopen(fn_eps33.c_str(), "w");
    FILE *f_eps11_analit = fopen(fn_eps11_analit.c_str(), "w");
    FILE *f_eps22_analit = fopen(fn_eps22_analit.c_str(), "w");
    FILE *f_eps33_analit = fopen(fn_eps33_analit.c_str(), "w");
    FILE *f_sigmaEqv = fopen(fn_sigmaEqv.c_str(), "w");
    {
        double time = (*out.mechOut)[globalStepNumber].step[0].t0;
        // напряжения (для проверки релаксации)
        //if(time >= cp.timeCreep && time <= cp.timeCreep + cp.timeRelaxation1 + cp.timeRelaxation2)
        double sigmaEqv_min = 1.e100, sigmaEqv_max = -1.e100;
        for(int feInd = 0; feInd  < (int)(*out.mechOut)[globalStepNumber].fe.size(); feInd++)
        {
            MechOutFePointData &r = (*out.mechOut)[globalStepNumber].fe[feInd].pd[0];
            MechMaterialSource &m0 = (*(*task.mechTask.mechStep)[globalStepNumber].material)[0];
            double sigmaEqv = r.sumSigma.eqv_SIGMA();
            if(sigmaEqv < sigmaEqv_min)
                sigmaEqv_min = sigmaEqv;
            if(sigmaEqv > sigmaEqv_max)
                sigmaEqv_max = sigmaEqv;
        }
        fprintf(f_sigmaEqv, "%le %le %le\n", time, sigmaEqv_min, sigmaEqv_max);
        // деформации (для проверки ползучести)
        if(time > 0 && time <= cp.timeCreep)
        {
            double creepEps11_min = 1.e100, creepEps11_max = -1.e100;
            double creepEps22_min = 1.e100, creepEps22_max = -1.e100;
            double creepEps33_min = 1.e100, creepEps33_max = -1.e100;
            for(int feInd = 0; feInd  < (*out.mechOut)[globalStepNumber].fe.size(); feInd++)
            {
                MechOutFePointData &r = (*out.mechOut)[globalStepNumber].fe[feInd].pd[0];
                //MechMaterialSource &m0 = (*(*task.mechTask.mechStep)[globalStepNumber].material)[0];
                if(r.sumEpsPlastic.m[0][0] < creepEps11_min)
                    creepEps11_min = r.sumEpsPlastic.m[0][0];
                if(r.sumEpsPlastic.m[0][0] > creepEps11_max)
                    creepEps11_max = r.sumEpsPlastic.m[0][0];
                if(r.sumEpsPlastic.m[1][1] < creepEps22_min)
                    creepEps22_min = r.sumEpsPlastic.m[1][1];
                if(r.sumEpsPlastic.m[1][1] > creepEps22_max)
                    creepEps22_max = r.sumEpsPlastic.m[1][1];
                if(r.sumEpsPlastic.m[2][2] < creepEps33_min)
                    creepEps33_min = r.sumEpsPlastic.m[2][2];
                if(r.sumEpsPlastic.m[2][2] > creepEps33_max)
                    creepEps33_max = r.sumEpsPlastic.m[2][2];
            }
            fprintf(f_eps11, "%le %le %le\n", time, creepEps11_min, creepEps11_max);
            fprintf(f_eps22, "%le %le %le\n", time, creepEps22_min, creepEps22_max);
            fprintf(f_eps33, "%le %le %le\n", time, creepEps33_min, creepEps33_max);
        }
    }
    // аналитическое решение
    // деформации (для проверки ползучести)
    for(double time = 0; time < cp.timeCreep; time += cp.timeCreep/1000)
    {
        double creepEps11 = cp.A*pow(cp.sigma0, cp.n)*pow(time, cp.m) *
                (0.5*(2 - cp.alpha)*pow(cp.alpha*cp.alpha - cp.alpha + 1, 0.5*(cp.n - 1)));
        double creepEps22 = cp.A*pow(cp.sigma0, cp.n)*pow(time, cp.m) *
                (0.5*(2*cp.alpha - 1)*pow(cp.alpha*cp.alpha - cp.alpha + 1, 0.5*(cp.n - 1)));
        double creepEps33 = cp.A*pow(cp.sigma0, cp.n)*pow(time, cp.m) *
                (0.5*(cp.alpha + 1)*pow(cp.alpha*cp.alpha - cp.alpha + 1, 0.5*(cp.n - 1)));
        fprintf(f_eps11_analit, "%le %le\n", time, creepEps11);
        fprintf(f_eps22_analit, "%le %le\n", time, creepEps22);
        fprintf(f_eps33_analit, "%le %le\n", time, creepEps33);
    }


    /*
    // деформации
    double creepEps11_min = 1.e100, creepEps11_max = -1.e100;
    double creepEps22_min = 1.e100, creepEps22_max = -1.e100;
    double creepEps33_min = 1.e100, creepEps33_max = -1.e100;
    int globalStepNumber = (int)(*out.mechOut).size() - 2;
    {
        for(int feInd = 0; feInd  < (*out.mechOut)[globalStepNumber].fe.size(); feInd++)
        {
            MechOutFePointData &r = (*out.mechOut)[globalStepNumber].fe[feInd].pd[0];
            MechMaterialSource &m0 = (*(*task.mechTask.mechStep)[globalStepNumber].material)[0];
            if(r.sumEpsCreep[0] < creepEps11_min)
                creepEps11_min = r.sumEpsCreep[0];
            if(r.sumEpsCreep[0] > creepEps11_max)
                creepEps11_max = r.sumEpsCreep[0];
            if(r.sumEpsCreep[1] < creepEps22_min)
                creepEps22_min = r.sumEpsCreep[1];
            if(r.sumEpsCreep[1] > creepEps22_max)
                creepEps22_max = r.sumEpsCreep[1];
            if(r.sumEpsCreep[2] < creepEps33_min)
                creepEps33_min = r.sumEpsCreep[2];
            if(r.sumEpsCreep[2] > creepEps33_max)
                creepEps33_max = r.sumEpsCreep[2];
        }
    }
    double creepEps11 = cp.A*pow(cp.sigma0, cp.n)*pow(cp.timeCreep, cp.m) *
            (0.5*(2 - cp.alpha)*pow(cp.alpha*cp.alpha - cp.alpha + 1, 0.5*(cp.n - 1)));
    double creepEps22 = cp.A*pow(cp.sigma0, cp.n)*pow(cp.timeCreep, cp.m) *
            (0.5*(2*cp.alpha - 1)*pow(cp.alpha*cp.alpha - cp.alpha + 1, 0.5*(cp.n - 1)));
    double creepEps33 = cp.A*pow(cp.sigma0, cp.n)*pow(cp.timeCreep, cp.m) *
            (0.5*(cp.alpha + 1)*pow(cp.alpha*cp.alpha - cp.alpha + 1, 0.5*(cp.n - 1)));
    fprintf(f_eps11, "numb = %le .. %le analit = %le\n", creepEps11_min, creepEps11_max, creepEps11);
    fprintf(f_eps22, "numb = %le .. %le analit = %le\n", creepEps22_min, creepEps22_max, creepEps22);
    fprintf(f_eps33, "numb = %le .. %le analit = %le\n", creepEps33_min, creepEps33_max, creepEps33);
    */
/*
    int globalStepNumber = (int)(*out.mechOut).size() - 2;
    MechOutFePointData &r = (*out.mechOut)[globalStepNumber].fe[0].pd[0];
    MechMaterialSource &m0 = (*(*task.mechTask.mechStep)[globalStepNumber].material)[0];
    double creepEps11 = cp.A*pow(cp.sigma0, cp.n)*pow(cp.timeCreep, cp.m) *
            (0.5*(2 - cp.alpha)*pow(cp.alpha*cp.alpha - cp.alpha + 1, 0.5*(cp.n - 1)));
    double creepEps22 = cp.A*pow(cp.sigma0, cp.n)*pow(cp.timeCreep, cp.m) *
            (0.5*(2*cp.alpha - 1)*pow(cp.alpha*cp.alpha - cp.alpha + 1, 0.5*(cp.n - 1)));
    double creepEps33 = cp.A*pow(cp.sigma0, cp.n)*pow(cp.timeCreep, cp.m) *
            (0.5*(cp.alpha + 1)*pow(cp.alpha*cp.alpha - cp.alpha + 1, 0.5*(cp.n - 1)));
    fprintf(f_creepEps11, "numb = %le analit = %le\n", r.sumEpsCreep[0], creepEps11);
    fprintf(f_creepEps22, "numb = %le analit = %le\n", r.sumEpsCreep[1], creepEps22);
    fprintf(f_creepEps33, "numb = %le analit = %le\n", r.sumEpsCreep[2], creepEps33);
*/
    fclose(f_eps11);
    fclose(f_eps22);
    fclose(f_eps33);
    fclose(f_eps11_analit);
    fclose(f_eps22_analit);
    fclose(f_eps33_analit);
    fclose(f_sigmaEqv);
    // построение графиков
    makeGraphsAfterGlobalStep(globalStepNumber, task);
}
bool Test_creep::possibleToShow2d() const
{
    return true;
}
bool Test_creep::getContactSurfaceCircle(const Task &task, const int globalStepIndex, POINT2 &c, double &R) const
{
    return false;
}
void Test_creep::needToDrawFe(const Solid::Task &, const OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const
{
    faceState = {1, 0, 0, 0, 0, 0};
}
bool Test_creep::needToPaintFiniteElementSurface() const
{
    return false;
}
}


namespace Grid
{
void Grid3D::genCreep(const CreepParameters &cp)
{
    #define MOVE(v_ind, d0, d1, d2)	((v_ind + offset0*d0 + offset1*d1 + offset2*d2))
    using namespace Operations;
    fe.clear();
    vertex.clear();
    vertexForCurvature.clear();
    bc1.clear();
    FESurface.clear();
    analiticalSurface.clear();
    ISurface.clear();
    {
        int N0 = cp.N[0];
        int N1 = cp.N[1];
        int N2 = cp.N[2];
        int offset2 = 1;
        int offset1 = N2 + 1;
        int offset0 = (N2 + 1)*(N1 + 1);
        // построение вершин и первых краевых условий
        for (int i0 = 0; i0 <= N0; i0++)            // x
            for (int i1 = 0; i1 <= N1; i1++)        // y
                for (int i2 = 0; i2 <= N2; i2++)    // z
                {
                    POINT3 p;
                    findPointOnTheLine_1d(cp.cube.i[0][0], cp.cube.i[0][1], N0, 1, i0, p[0]);
                    findPointOnTheLine_1d(cp.cube.i[1][0], cp.cube.i[1][1], N1, 1, i1, p[1]);
                    findPointOnTheLine_1d(cp.cube.i[2][0], cp.cube.i[2][1], N2, 1, i2, p[2]);
                    vertex.push_back(p);
                    // фиксируем по x слева
                    if(i0 == 0)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 0;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    // фиксируем по x справа (для разгрузки) или не фиксируем
                    //if(i0 == N0)
                    // фиксируем все границы (для разгрузки) или не фиксируем
                    if(i0 == 0 || i0 == N0
                    || i1 == 0 || i1 == N1
                    || i2 == 0 || i2 == N2)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 4;// подрубается при разгрузке
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    // фиксируем по y снизу
                    if(i1 == 0)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 1;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    // фиксируем по z посередине
                    //if(i2 == N2/2)
                    //if(true)
                    if(i2 == 0)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 2;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    /*
                    // фиксируем по z везде
                    if(true)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.si = 2;
                        bc1_el.vi = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }*/
                    /*
                    // температура снизу
                    if(i1 == 0)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.si = 3;
                        bc1_el.vi = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    // температура сверху
                    if(i1 == N1)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.si = 4;
                        bc1_el.vi = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    */
                }
        // построение конечных элементов и поверхностей
        // будет 2 поверхности (длля 2-х краевых)
        FESurface.resize(2);
        for (size_t FEsurfaceInd = 0; FEsurfaceInd < FESurface.size(); FEsurfaceInd++)
            FESurface[FEsurfaceInd].face.clear();
        for (int i0 = 0; i0 < N0; i0++)            // x
            for (int i1 = 0; i1 < N1; i1++)        // y
                for (int i2 = 0; i2 < N2; i2++)    // z
                {
                    // конечный элемент
                    int vind0 = i0*offset0 + i1*offset1 + i2*offset2;       // индекс вершины, определяющей шестигранник
                    int vi[8] =
                    {
                        MOVE(vind0, 0, 0, 0), MOVE(vind0, 1, 0, 0), MOVE(vind0, 0, 1, 0), MOVE(vind0, 1, 1, 0),
                        MOVE(vind0, 0, 0, 1), MOVE(vind0, 1, 0, 1), MOVE(vind0, 0, 1, 1), MOVE(vind0, 1, 1, 1)
                    };
                    FE_LinearHexagon *fe_el = new FE_LinearHexagon;   // шестигранник, линейное отображение
                    fe_el->mi = 0;                       // материал с индексом 0
                    for (int t = 0; t < 8; t++)
                        fe_el->vi[t] = vi[t];
                    fe.push_back(fe_el);
                    // поверхности для 2-х краевых условий
                    FEFace face;
                    // справа
                    if (i0 == N0 - 1)
                    {
                        face.feIndex = (int)fe.size() - 1;
                        face.faceIndex = 3; // X = +1
                        FESurface[0].face.push_back(face);
                    }
                    // сверху
                    if (i1 == N1 - 1)
                    {
                        face.feIndex = (int)fe.size() - 1;
                        face.faceIndex = 5; // Y = +1
                        FESurface[1].face.push_back(face);
                    }
                }
    }

    // параметры индексации пространства
    {
        Sphere s0;
        s0.O = {0, 0, 0};
        s0.R = 13;
        regionIndexationParameters.q0.initBySphere(s0);
    }
    // нумерация базисных функций
    DOFs = new Grid::GlobalDOFs;
    DOFs->init(vertex.size());
    // ф-и первого порядка
    for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
    {
        Grid::FuncID funcID;
        funcID.type = Grid::FuncType::linear;
        funcID.crackIndex = -1;
        DOFs->addDOFs(vertexIndex, funcID);
    }
}
}
