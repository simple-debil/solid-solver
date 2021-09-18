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
Test_sphere_Te::Test_sphere_Te()
{
    // дирректория данных теста
    dir = "./tests/hollowSphere_Te/";
    // имена gnuplot файлов(+пути к файлам в опциях) и соответствующих картинок
    fn_filesList = dir + "files.txt";
    resGraph.read_gnuplot_from_file(fn_filesList);
    // инициализация графиков с невязками
    initStepsGraphs();
}
Test_base::Type Test_sphere_Te::get_type() const
{
    return Type::Sphere_Te;
}
void Test_sphere_Te::initTask(Task &task)
{
    // индекс теста решателя теплопроводности
    thermTestIndex = 1;             // 1 - t(a)=t1, t(b)=t2
                                    // 2 - t(a)=t1, конвекция на b
                                    // 3 - конвекция на r=a, конвекция на r=b
                                    // 4 - подогрев на r=a, t(b)=t2
    task.thermTask.enabled = true;
    task.mechTask.enabled = true;
    // общие шаги и сетка
    std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
    Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
    gp.r1 = 1;
    gp.r2 = 4;
    posledovatelnost = 0;   //0 +температура -> +сила          -> -сила -> -температура     -эталонное решение
                            //1 *+температура+сила*            -> *-температура-сила*
                            //2 +нагружение  -> *+температура* -> -нагружение    -> *?-температура*?
                            //3 +нагружение  -> *+температура* -> *-температура* -> -нагружение
    IncForsesMode incForsesMode = IncForsesMode::IncrementP;
    //IncForsesMode incForsesMode = IncForsesMode::MinusIntegral;
    //IncForsesMode incForsesMode = IncForsesMode::bPlusR;
    int Steps_number = 8;//1 //8 //64
    gp.N = 8;   // 4 промежутка на четвертьокружностях    //8, 32
    gp.Nparts = 32;
    gp.q = 1 + 0./(gp.Nparts/2);
    gp.curvilinear = 0;    // отображение (0 - линейное, 1 - квадратичное)
    gp.buildingMethod = 0; // способ построения (0 - делим дуги, 1 - отображение на куб)
    grid->genSphere(gp);

    // параметры
    sigmaSolvingType = 0;  // способ расчёта решения (0 - главные напряжения, 1 - проекции площадки)
     double P00 = 1*1.0e2;//100  // нагрузка
     double Talpha = 0*1.e-5;    // коэффициент температурного расширения
    double E0 = 1.e10;
    double NU0 = 0.3; //0.4
    int fixGrid = 0;            // 0 - подвижная сетка, 1 - зафиксировать сетку
    int StepsNumber[4];   // шаги
    std::vector<double> P_data(4 + 1);
    double *P = &(P_data[1]);   // - массив из 5 элементов, нумерация начинается с -1

    // тест1
    double t1_Ta00 = 100;
    double t1_Tb00 = 1;

    std::vector<double> t1_Ta_data(4 + 1);
    std::vector<double> t1_Tb_data(4 + 1);
    double *t1_Ta = &(t1_Ta_data[1]);   // - массив из 5 элементов, нумерация начинается с -1
    double *t1_Tb = &(t1_Tb_data[1]);   // - массив из 5 элементов, нумерация начинается с -1

    // правильность аналитических графиков(сплашной рисовать или пунктиром)
    // T                       да              да                да                да
    // sigma                   да              нет               да                да
    if(posledovatelnost == 0)
    {
        StepsNumber[0] = Steps_number;//1;
        StepsNumber[1] = Steps_number;//1;
        StepsNumber[2] = Steps_number;//1;
        StepsNumber[3] = Steps_number;//1;
        P[-1] = 0;
        P[0] = 0;
        P[1] = P00;
        P[2] = 0;
        P[3] = 0;
        t1_Ta[-1] = 0;
        t1_Tb[-1] = 0;
        t1_Ta[0] = t1_Ta00;
        t1_Tb[0] = t1_Tb00;
        t1_Ta[1] = t1_Ta00;
        t1_Tb[1] = t1_Tb00;
        t1_Ta[2] = t1_Ta00;
        t1_Tb[2] = t1_Tb00;
        t1_Ta[3] = 0;
        t1_Tb[3] = 0;
    }
    if(posledovatelnost == 1)
    {
        StepsNumber[0] = Steps_number;
        StepsNumber[1] = Steps_number;//1;
        StepsNumber[2] = Steps_number;
        StepsNumber[3] = Steps_number;//1;
        P[-1] = 0;
        P[0] = P00;
        P[1] = P00;
        P[2] = 0;
        P[3] = 0;
        t1_Ta[-1] = 0;
        t1_Tb[-1] = 0;
        t1_Ta[0] = t1_Ta00;
        t1_Tb[0] = t1_Tb00;
        t1_Ta[1] = t1_Ta00;
        t1_Tb[1] = t1_Tb00;
        t1_Ta[2] = 0;
        t1_Tb[2] = 0;
        t1_Ta[3] = 0;
        t1_Tb[3] = 0;
    }
    if(posledovatelnost == 2)
    {
        StepsNumber[0] = Steps_number;//1;
        StepsNumber[1] = Steps_number;
        StepsNumber[2] = Steps_number;//1;
        StepsNumber[3] = Steps_number;
        P[-1] = 0;
        P[0] = P00;
        P[1] = P00;
        P[2] = 0;
        P[3] = 0;
        t1_Ta[-1] = 0;
        t1_Tb[-1] = 0;
        t1_Ta[0] = 0;
        t1_Tb[0] = 0;
        t1_Ta[1] = t1_Ta00;
        t1_Tb[1] = t1_Tb00;
        t1_Ta[2] = t1_Ta00;
        t1_Tb[2] = t1_Tb00;
        t1_Ta[3] = 0;
        t1_Tb[3] = 0;
    }
    if(posledovatelnost == 3)
    {
        StepsNumber[0] = Steps_number;//1;
        StepsNumber[1] = Steps_number;
        StepsNumber[2] = Steps_number;
        StepsNumber[3] = Steps_number;//1;
        P[-1] = 0;
        P[0] = P00;
        P[1] = P00;
        P[2] = P00;
        P[3] = 0;
        t1_Ta[-1] = 0;
        t1_Tb[-1] = 0;
        t1_Ta[0] = 0;
        t1_Tb[0] = 0;
        t1_Ta[1] = t1_Ta00;
        t1_Tb[1] = t1_Tb00;
        t1_Ta[2] = 0;
        t1_Tb[2] = 0;
        t1_Ta[3] = 0;
        t1_Tb[3] = 0;
    }

    int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny
    double Time = 1;

    // шаги
    GlobalStep s;       // шаг
    for(int globalStep = 0; globalStep < 4; globalStep++)
    {
        for(int i = 0; i < StepsNumber[globalStep]; i++)
        {
            s.t_start = Time*globalStep + Time*i/StepsNumber[globalStep];
            s.t_finish = Time*globalStep + Time*(i + 1)/StepsNumber[globalStep];
            s.dt0 = Time/StepsNumber[globalStep];
            step->push_back(s);
        }
    }

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
        mechMat.elasticSigmaLimit = 1.e10;//2.e7;
        mechMat.Ceps = 4./9.;
        mechMat.Csigma = 1;
        mechMat.set_E_NU(E0, NU0);
        mechMat.elasticParameters0.Talpha = Talpha;
        mechMat.set_M_sigma();
        mechMat.F = {0,0,0};
        mechMat.elasticParameters0.ro = 0;//1e10;
        mechMat.plasticityMethodType = MechPlasticityMethodType::Elasticity;
        setLinearFunction(E0, E0*0.1, 0, 100, mechMat.E_Fun);
        setLinearFunction(NU0, NU0, 0, 100, mechMat.NU_Fun);
        mechMat.temperatureDependence = true;
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


        // шаги
        std::vector<MechGlobalStepParameters> *mechStep = new std::vector<MechGlobalStepParameters>;
        // нагружение - разгрузка
        MechGlobalStepParameters ms;  // шаг для МДТТ
        GlobalStep s;       // шаг
        ms.slausolverParameters = slausolver_parameters;
        ms.timeMode = 0;    // 0 - квазистатическая задача
        ms.incForsesMode = incForsesMode;  // 1 - dSigma = integral(...), 2 - приращения
        ms.fixGrid = fixGrid;     // 1 - зафиксировать сетку
        ms.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
        ms.controlMode = 0; // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
        ms.terminateIfAccuracyIsNotAchieving = false;
        ms.iterLimit = 1;
        ms.slauResidualLimit = 1000;
        ms.plasticResidualLimit.eps = 1000;
        ms.plasticResidualLimit.sigma = 1000;
        ms.material = material;
        ms.bc1Source = bc1Source;
        //ms.bc2Source = bc2Source; - будет меняться
        ms.bc2 = bc2;


        int stepCount = 0;
        for(int globalStep = 0; globalStep < 4; globalStep++)
        {
            for(int i = 0; i < StepsNumber[globalStep]; i++)
            {
                s = (*step)[stepCount];
                setBc2Sphere(bc2Source,
                        P[globalStep - 1] + (P[globalStep] - P[globalStep - 1])*i/StepsNumber[globalStep],
                        P[globalStep - 1] + (P[globalStep] - P[globalStep - 1])*(i + 1)/StepsNumber[globalStep],
                        s.t_start, s.t_finish);
                ms.bc2Source = bc2Source;
                mechStep->push_back(ms);
                stepCount++;
            }
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
        using namespace Heat;
        // температура
        double T0 = 0;
        // сетка
        task.thermTask.grid = grid;
        // общие параметры шагов
        task.thermTask.step = step;
        std::vector<ThermMaterialSource> *material;             // параметры материалов
        std::vector<ThermBoundaryCondition1Source> *bc1SourceT; // 1-е краевые
        std::vector<ThermBoundaryCondition2Source> *bc2SourceT; // 2-е краевые
        // материал
        material = new std::vector<ThermMaterialSource>(1);
        ThermMaterialSource &thermMat = (*material)[0];
        thermMat.ro = 1;     // плотность
        thermMat.c = 5;      // удельная теплоёмкость
        thermMat.f = 0;      // мощность внутренних объёмных источников (стоков) тепла
        double L0 = 10;//10;
        thermMat.L.initDiag(L0); // тензор теплопроводности
        int timeMode = 0;


        // шаги
        std::vector<ThermGlobalStep> *thermStep = new std::vector<ThermGlobalStep>;
        ThermGlobalStep ts;  // шаг для теплопроводности
        ts.slausolverParameters = slausolver_parameters;
        ts.timeMode = timeMode;    // 0 - статика, 1 - динамика
        ts.material = material;
        //ts.bc1SourceT меняется;
        //ts.bc2SourceT меняется;
        int stepCount = 0;
        for(int globalStep = 0; globalStep < 4; globalStep++)
        {
            for(int i = 0; i < StepsNumber[globalStep]; i++)
            {
                s = (*step)[stepCount];
                // Тест1
                if(thermTestIndex == 1)
                    setBcSphereT_1(bc1SourceT, bc2SourceT,
                            t1_Ta[globalStep - 1] + (t1_Ta[globalStep] - t1_Ta[globalStep - 1])*(i + 1)/StepsNumber[globalStep],    // температура на r = a
                            t1_Tb[globalStep - 1] + (t1_Tb[globalStep] - t1_Tb[globalStep - 1])*(i + 1)/StepsNumber[globalStep]     // температура на r = b
                                   );
                // Тест2
                if(thermTestIndex == 2)
                    setBcSphereT_2(bc1SourceT, bc2SourceT,
                                   1,   // температура на r = a
                                   2,   // температура окружающей среды r > b
                                   1    // коэффициент конвективного теплообмена
                                   );
                // Тест3
                if(thermTestIndex == 3)
                    setBcSphereT_3(bc1SourceT, bc2SourceT,
                                   1,   // температура окружающей среды r < a
                                   10,  // коэффициент конвективного теплообмена
                                   2,   // температура окружающей среды r > b
                                   20   // коэффициент конвективного теплообмена
                                   );
                // Тест4
                if(thermTestIndex == 4)
                {
                    setBcSphereT_4(bc1SourceT, bc2SourceT, thermMat,
                                   100, // плотность набегающего теплового потока на r = a
                                   15,  // коэффициент теплопроводности
                                   1    // температура на r = b
                                   );
                }
                // Тест5
                if(thermTestIndex == 5)
                {
                    setBcSphereT_5(bc1SourceT, bc2SourceT, thermMat,
                                   10, // плотность набегающего теплового потока на r = a
                                   1   // коэффициент теплопроводности
                                   );
                    timeMode = 1;
                }
                ts.bc1SourceT = bc1SourceT;
                ts.bc2SourceT = bc2SourceT;
                thermStep->push_back(ts);
                stepCount++;
            }
        }

        task.thermTask.thermStep = thermStep;
        // инициализация начальных температур внутри кэ
        task.thermTask.fe = new std::vector<ThermFeData>;
        task.thermTask.fe->resize(grid->fe.size());
        for (size_t i = 0; i < grid->fe.size(); i++)
        {
            (*task.thermTask.fe)[i].init(BasisType_1L,
                                         HomogenyMode,
                                         Integration::IntegrationType::Gauss3,
                                         T0);
        }
        // инициализация начальных температур в вершинах
        task.thermTask.vertex = new std::vector<ThermVertexData>;
        for (size_t i = 0; i < grid->vertex.size(); i++)
        {
            ThermVertexData v_el;
            if(thermTestIndex == 5)
            {
                double a = task.grid->vertex[0].abs();
                double R = grid->vertex[i].abs();
                double T_a_ch = 100;
                double Q0 = (*bc2SourceT)[0].q * 4*PI*SQR(a);
                T0 = T_a_ch + Q0/(4*PI*(*material)[0].L.m[0][0])*(1./R - 1./a);
            }
            v_el.T = T0;
            v_el.newT = T0;
            task.thermTask.vertex->push_back(v_el);
        }
    }
}
void Test_sphere_Te::writeResults(const Task &task, const OutData &out)
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
    fn_T_analit = subdir + "T_analit.txt";
    fn_T = subdir + "T.txt";
    fn_T_pogr = subdir + "T_pogr.txt";
    fn_u_analit = subdir + "u_analit.txt";
    fn_u = subdir + "u.txt";
    fn_u_pogr = subdir + "u_pogr.txt";
    fn_sigma_r_analit = subdir + "sigma_r_analit.txt";
    fn_sigma_r = subdir + "sigma_r.txt";
    fn_sigma_r_pogr = subdir + "sigma_r_pogr.txt";
    fn_sigma_fi_analit = subdir + "sigma_fi_analit.txt";
    fn_sigma_fi = subdir + "sigma_fi.txt";
    fn_sigma_fi_pogr = subdir + "sigma_fi_pogr.txt";
    fn_max_pogr = subdir + "_max_pogr.txt";
    // пути к файлам для отображения графиков
{
    std::ofstream f;
    if(globalStepNumber == 0)
    {
        f.open(fn_filesList, std::ofstream::out);
        resGraph.generalGraphsNomber = 0;
        resGraph.graphsPerGlobalStepNomber = 3;
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
    // sigma_r
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "sigma_r.png") +
         genGnuplotCommandLineParameter("fn_sigma_r", subdir + "sigma_r.txt") +
         genGnuplotCommandLineParameter("fn_sigma_r_analit", subdir + "sigma_r_analit.txt") +
         "\" " +
         dir + "sigma_r.gnu" << "\n";
    f << subdir + "sigma_r.png" << "\n";
    f << "sigma_r.png" << "\n";
    // T
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "T.png") +
         genGnuplotCommandLineParameter("fn_T", subdir + "T.txt") +
         genGnuplotCommandLineParameter("fn_T_analit", subdir + "T_analit.txt") +
         "\" " +
         dir + "T.gnu" << "\n";
    f << subdir + "T.png" << "\n";
    f << "T.png" << "\n";
    f.close();
}
    resGraph.read_gnuplot_from_file(fn_filesList);
    }

    FILE *fT_analit = fopen(fn_T_analit.c_str(), "w");
    FILE *fT = fopen(fn_T.c_str(), "w");
    FILE *fT_pogr = fopen(fn_T_pogr.c_str(), "w");
    FILE *fsigma_r_analit = fopen(fn_sigma_r_analit.c_str(), "w");
    FILE *fsigma_r = fopen(fn_sigma_r.c_str(), "w");
    FILE *fsigma_r_pogr = fopen(fn_sigma_r_pogr.c_str(), "w");
    FILE *fsigma_fi_analit = fopen(fn_sigma_fi_analit.c_str(), "w");
    FILE *fsigma_fi = fopen(fn_sigma_fi.c_str(), "w");
    FILE *fsigma_fi_pogr = fopen(fn_sigma_fi_pogr.c_str(), "w");

    FILE *fmax_pogr = fopen(fn_max_pogr.c_str(), "w");

    std::vector<GlobalStep> *step = task.step;
    Grid::Grid3D *grid = task.grid;

    {
        double fi_pogr_max_abs = 0;
        double r_pogr_max_abs = 0;
        double sigma_r_max = 0;
        double sigma_fi_max = 0;
        double fi_pogr_max = 0;
        double r_pogr_max = 0;

        double max_pogr_sigma_fi = 0;
        double max_pogr_sigma_r = 0;

        // вывод данных МДТТ
        if(task.mechTask.enabled)
        {
            int numLocalStep = (int)(*out.mechOut)[globalStepNumber].step.size() - 1;
            // механическое решение
            MechMaterialSource &m0 = (*(*task.mechTask.mechStep)[globalStepNumber].material)[0];
            MechOutStepData &mechOutStep_el = (*out.mechOut)[globalStepNumber].step[numLocalStep];
            MechBoundaryCondition2Source_ScalarFunction *mbc2Source_sf = (MechBoundaryCondition2Source_ScalarFunction *)(*(*task.mechTask.mechStep)[globalStepNumber].bc2Source)[0];
            double P0 = mbc2Source_sf->value(mechOutStep_el.t0);
            POINT3 v1 = grid->vertex[0];
            int N = gp.N/2+1;
            POINT3 v2 = grid->vertex[(3*N*N-3*N)*(gp.Nparts+1)];
            double a = v1.abs();
            double b = v2.abs();
            double c = calc_c(a, b, P0, m0.elasticSigmaLimit);
            MechIterInf nlInf = mechOutStep_el.nlInf.iterInf.back();
            //fprintf(fsigma_r, "t1 = %le t2 = %le t = %le\n", t1, t2, mechOutStep_el.t);

            fi_pogr_max_abs = 0;
            r_pogr_max_abs = 0;

            // расчет численного решения в конечных элементах
            for(int j0 = 0; j0 < gp.Nparts; j0++)//int j0 = gp.Nparts - 1;
            {
                double R = 0;
                double sigma_r_numb = 0;
                double sigma_fi_numb= 0;
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

                        sigma_r_numb = ms[ind[0]];
                        sigma_fi_numb = ms[ind[1]];//(ms[ind[1]] + ms[ind[2]]) / 2.;

                        sp[spInd].R = R;
                        sp[spInd].sigma_r = sigma_r_numb;
                        sp[spInd].sigma_fi = sigma_fi_numb;
                        double sigma_r_load;
                        double sigma_fi_load;
                        double sigma_r_unload;
                        double sigma_fi_unload;
                        calcAnalit(a,b,R,P0,m0.elasticSigmaLimit,c,sigma_r_load,sigma_fi_load,sigma_r_unload,sigma_fi_unload);

                        {
                            sp[spInd].sigma_r_analit = sigma_r_load;
                            sp[spInd].sigma_fi_analit = sigma_fi_load;
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
                            fprintf(fsigma_r_pogr, "%le %le\n", R, fabs(sp[spInd].sigma_r_pogr)*100);
                            fprintf(fsigma_fi_pogr, "%le %le\n", R, fabs(sp[spInd].sigma_fi_pogr)*100);
                        }

                        spInd++;
                    }
                }
            }

            //  максимальная погрешность шага 3
            /*if(globalStepNumber == 3)
            {
                fprintf(fmax_pogr, "%le\n%le\n", fi_pogr_max_abs, r_pogr_max_abs);
            }*/
            fprintf(fmax_pogr, "a = %le\n", a);
            fprintf(fmax_pogr, "b = %le\n", b);
            fprintf(fmax_pogr, "c = %le\n", c);
            fprintf(fmax_pogr, "P0 = %le\n", P0);

            // расчет аналитического решения
            int NumPoints = 1000;
            for(int j = 0; j <= NumPoints; j++)
            {
                double R = a + (b-a)*j/NumPoints;
                double sigma_r_load;
                double sigma_fi_load;
                double sigma_r_unload;
                double sigma_fi_unload;
                calcAnalit(a,b,R,P0,m0.elasticSigmaLimit,c,sigma_r_load,sigma_fi_load,sigma_r_unload,sigma_fi_unload);
                //if(mechOutStep_el.t == t2)
                {
                    fprintf(fsigma_r_analit, "%le %le\n", R, sigma_r_load);
                    fprintf(fsigma_fi_analit, "%le %le\n", R, sigma_fi_load);
                    if(fabs(sigma_r_load) > sigma_r_max)
                        sigma_r_max = fabs(sigma_r_load);
                    if(fabs(sigma_fi_load) > sigma_fi_max)
                        sigma_fi_max = fabs(sigma_fi_load);
                }
            }
        }


        // вывод данных задачи теплопроводности
        double T_pogr_max = 0;
        double T_pogr_max_abs = 0;
        Heat::ThermGlobalStep &ThermGlobalStep_el = (*task.thermTask.thermStep)[globalStepNumber];
        if(task.thermTask.enabled)
        {
            using namespace Heat;
            int numLocalStep = (int)(*out.thermOut)[globalStepNumber].step.size() - 1;
            ThermOutStepData &thermOutStep_el = (*out.thermOut)[globalStepNumber].step[numLocalStep];
            POINT3 v1 = grid->vertex[0];
            int N = gp.N/2+1;
            int N_b = (3*N*N-3*N)*(gp.Nparts+1);
            POINT3 v2 = grid->vertex[N_b];
            double a = v1.abs();
            double b = v2.abs();
            ThermMaterialSource &m0 = (*ThermGlobalStep_el.material)[0];
            // температура
            double T1 = (*ThermGlobalStep_el.bc1SourceT)[3].T0;
            double T2 = (*ThermGlobalStep_el.bc1SourceT)[4].T0;
            double Ta1 = (*ThermGlobalStep_el.bc2SourceT)[0].Ta;
            double hi1 = (*ThermGlobalStep_el.bc2SourceT)[0].hi;
            double h1 = hi1/m0.L.m[0][0];
            double Ta2 = (*ThermGlobalStep_el.bc2SourceT)[1].Ta;
            double hi2 = (*ThermGlobalStep_el.bc2SourceT)[1].hi;
            double h2 = hi2/m0.L.m[0][0];
            double T_b_ch = (*out.thermOut)[globalStepNumber].vertex[N_b].T;
            double T_a_ch = (*out.thermOut)[globalStepNumber].vertex[0].T;
            double Q0 = (*ThermGlobalStep_el.bc2SourceT)[0].q * 4*PI*SQR(a);
            // расчет численного решения и его погрешности
            for(size_t j = 0; j < (*out.thermOut)[globalStepNumber].vertex.size(); j++)
            {
                ThermOutVertexData &v = (*out.thermOut)[globalStepNumber].vertex[j];
                double R = grid->vertex[j].abs();
                //if(R >= 1.05)
                {
                double T_analit;
                if(thermTestIndex == 1)
                {
                    T_analit = (a*T1*(b-R)+b*T2*(R-a)) / (R*(b-a));
                }
                if(thermTestIndex == 2)
                {
                    T_analit = (a*T1*(h2*SQR(b) + R*(1-h2*b))+h2*SQR(b)*Ta2*(R-a)) /
                            (R*(h2*SQR(b) + a*(1-h2*b)));
                }
                if(thermTestIndex == 3)
                {
                    T_analit = (Ta1*SQR(a)*h1*(SQR(b)*h2 - R*(b*h2-1)) + Ta2*SQR(b)*h2*(R*(a*h1+1)-SQR(a)*h1)) /
                            (R*(SQR(b)*h2*(a*h1+1)-SQR(a)*h1*(b*h2-1)));
                }
                if(thermTestIndex == 4)
                {
                    T_analit = T_b_ch + Q0/(4*PI*m0.L.m[0][0])*(1./R - 1./b);
                }
                if(thermTestIndex == 5)
                {
                    T_analit = T_a_ch + Q0/(4*PI*m0.L.m[0][0])*(1./R - 1./a);
                }
                double T_pogr_abs = v.T-T_analit;
                double T_pogr = (v.T-T_analit)/fabs(T_analit);
                if(fabs(T_pogr) > T_pogr_max)
                    T_pogr_max = fabs(T_pogr);
                if(fabs(T_pogr_abs) > T_pogr_max_abs)
                    T_pogr_max_abs = fabs(T_pogr_abs);
                fprintf(fT, "%le %le\n", R, v.T);
                fprintf(fT_pogr, "%le %le\n", R, fabs(T_pogr)*100);
                //fprintf(fT_pogr_abs, "%le %le\n", R, T_pogr_abs*100);
                }
            }
            // расчет аналитического решения
            int NumPoints = 1000;
            for(int j = 0; j <= NumPoints; j++)
            {
                double R = a + (b-a)*j/NumPoints;
                double T_analit;
                if(thermTestIndex == 1)
                {
                    T_analit = (a*T1*(b-R)+b*T2*(R-a)) / (R*(b-a));
                }
                if(thermTestIndex == 2)
                {
                    T_analit = (a*T1*(h2*SQR(b) + R*(1-h2*b))+h2*SQR(b)*Ta2*(R-a)) /
                            (R*(h2*SQR(b) + a*(1-h2*b)));
                }
                if(thermTestIndex == 3)
                {
                    T_analit = (Ta1*SQR(a)*h1*(SQR(b)*h2 - R*(b*h2-1)) + Ta2*SQR(b)*h2*(R*(a*h1+1)-SQR(a)*h1)) /
                            (R*(SQR(b)*h2*(a*h1+1)-SQR(a)*h1*(b*h2-1)));
                }
                if(thermTestIndex == 4)
                {
                    T_analit = T_b_ch + Q0/(4*PI*m0.L.m[0][0])*(1./R - 1./b);
                }
                if(thermTestIndex == 5)
                {
                    T_analit = T_a_ch + Q0/(4*PI*m0.L.m[0][0])*(1./R - 1./a);
                }
                fprintf(fT_analit, "%le %le\n", R, T_analit);
            }

        }

    }

    fclose(fT_analit);
    fclose(fT);
    fclose(fT_pogr);
    fclose(fsigma_r_analit);
    fclose(fsigma_r);
    fclose(fsigma_r_pogr);
    fclose(fsigma_fi_analit);
    fclose(fsigma_fi);
    fclose(fsigma_fi_pogr);
    fclose(fmax_pogr);
    // построение графиков
    makeGraphsAfterGlobalStep(globalStepNumber, task);
}
}
