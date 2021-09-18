#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"
#include <fstream>

#include "tests.h"
#include "console.h"

#include "test_Kirsch.h"

using namespace Elementary;
using namespace Solid;
using namespace Grid;

namespace Tests
{
Test_Kirsch::Test_Kirsch()
{
    // дирректория данных теста
    dir = "./tests/Kirch/";
    // имена gnuplot файлов(+пути к файлам в опциях) и соответствующих картинок
    fn_filesList = dir + "files.txt";
    resGraph.read_gnuplot_from_file(fn_filesList);
    // инициализация графиков с невязками
    initStepsGraphs();
}
Test_base::Type Test_Kirsch::get_type() const
{
    return Type::Kirsch;
}
void Test_Kirsch::initTask(Task &task)
{
    // чтение параметров задачи
    fn_parametors = dir + "_input.in";
    FILE *f_parametors = fopen(fn_parametors.c_str(), "r");
    //fscanf(f_parametors, "%lf", &gp.NU);
    gp.XFEM = true;
    gp.NU = 0.4;
    gp.E = 1.e6;
    gp.Lx = 10;
    gp.Ly = 10;
    gp.Lz = 10;
    gp.r = 1;
    gp.Px = -1;
    gp.Py = 0;
    gp.N_fi = 64;
    gp.N_r = 32;
    gp.N_z = 1;//2
    gp.q = 1 + 1./16;
    /*
    gp.NU = 0.4;
    gp.E = 1.e6;
    gp.Lx = 100;
    gp.Ly = 100;
    gp.Lz = 100;
    gp.r = 1;
    gp.Px = -1;
    gp.Py = 0;
    gp.N_fi = 128;
    gp.N_r = 64;
    gp.N_z = 1;//2
    gp.q = 1 + 1./16;
    */
    fclose(f_parametors);
    // подрубка
    task.thermTask.enabled = false;
    task.mechTask.enabled = true;
    // общие шаги и сетка
    std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
    Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
    // сетка
    grid->genKirch(gp);
    MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::Elasticity;
    //MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::D_pl;
    //MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::D_pl_Solov;
    //MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::InitialSigma;
    //MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::Combo_D_pl_InitialSigma;
    double w_midPoint = 1;
    double w_project = 1;
    MechPlasticityCurveUnloadingType PCUnloadingType = MechPlasticityCurveUnloadingType::Elastic;
    //MechPlasticityCurveUnloadingType PCUnloadingType = MechPlasticityCurveUnloadingType::Curve;

    IncForsesMode incForsesMode = IncForsesMode::IncrementP;
    //IncForsesMode incForsesMode = IncForsesMode::bPlusR;
    //IncForsesMode incForsesMode = IncForsesMode::MinusIntegral;

    MaxSigmaEpsResidual plasticResidualLimit;
    plasticResidualLimit.eps = 1.e-10;
    plasticResidualLimit.sigma = 1.e-10;
    int nonlinearIterLimit = 50;    // ограничение на количество итераций
    int StepsNumber = 1;            // шаги
    double Time = 1;

    // шаги
    GlobalStep s;       // шаг
    // первый шаг (1 временной слой, заведомо упругое нагружение P)
    s.t_start = 0;   // (имеет значение только для первого глобального шага)
    s.t_finish = Time;
    s.dt0 = Time/StepsNumber;
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
        mechMat.elasticSigmaLimit = 1.e10;// ни на что не влияет
        mechMat.Ceps = 4./9.;
        mechMat.Csigma = 1;
        mechMat.set_E_NU(gp.E, gp.NU);
        mechMat.elasticParameters0.Talpha = 0;
        mechMat.set_M_sigma();
        mechMat.F = {0,0,0};
        mechMat.elasticParameters0.ro = 0;
        mechMat.w_midPoint = w_midPoint;
        mechMat.w_project = w_project;
        mechMat.plasticityMethodType = plasticityMethodType;
        mechMat.PCUnloadingType = PCUnloadingType;
        // первые краевые условия
        bc1Source = new std::vector<MechBoundaryCondition1Source>(4);
        (*bc1Source)[0].mode = {{ 0, -1, -1}};
        (*bc1Source)[0].u0 =   {{ 0, -1, -1}};
        (*bc1Source)[1].mode = {{-1,  0, -1}};
        (*bc1Source)[1].u0 =   {{-1,  0, -1}};
        (*bc1Source)[2].mode = {{-1, -1,  0}};
        (*bc1Source)[2].u0 =   {{-1, -1,  0}};
        (*bc1Source)[3].mode = {{-1, -1, -1}};
        (*bc1Source)[3].u0 =   {{-1, -1, -1}};
        // вторые краевые
        bc2 = new std::vector<MechBoundaryCondition2>;
        MechBoundaryCondition2 bc2_el;
        bc2_el.FEsurfaceInd = 0;
        bc2_el.bc2SourceIndex = 0;
        (*bc2).push_back(bc2_el);
        bc2_el.FEsurfaceInd = 1;
        bc2_el.bc2SourceIndex = 1;
        (*bc2).push_back(bc2_el);
        bc2_el.FEsurfaceInd = 2;
        bc2_el.bc2SourceIndex = 2;
        (*bc2).push_back(bc2_el);
        bc2_el.FEsurfaceInd = 3;
        bc2_el.bc2SourceIndex = 3;
        (*bc2).push_back(bc2_el);
        bc2_el.FEsurfaceInd = 4;
        bc2_el.bc2SourceIndex = 4;
        (*bc2).push_back(bc2_el);



        // шаги (на каждом шаге задаются вторые краевые условия)
        std::vector<MechGlobalStepParameters> *mechStep = new std::vector<MechGlobalStepParameters>;
        // нагружение - разгрузка
        MechGlobalStepParameters ms;  // шаг для МДТТ
        GlobalStep s;       // шаг
        ms.slausolverParameters = slausolver_parameters;
        ms.timeMode = 0;    // 0 - квазистатическая задача
        ms.incForsesMode = incForsesMode;  // 1 - dSigma = integral(...), 2 - приращения
        ms.fixGrid = 0;     // 1 - зафиксировать сетку
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
        // первый шаг (1 временной слой, заведомо упругое нагружение Px, Py)
        s = (*step)[0];
        bc2Source = new std::vector<MechBoundaryCondition2Source_base *>(5);
        setLinearBc2_el(0, gp.Px, s.t_start, s.t_finish,
                        (*bc2Source)[0]);
        setLinearBc2_el(0, gp.Px, s.t_start, s.t_finish,
                        (*bc2Source)[1]);
        setLinearBc2_el(0, gp.Py, s.t_start, s.t_finish,
                        (*bc2Source)[2]);
        setLinearBc2_el(0, gp.Py, s.t_start, s.t_finish,
                        (*bc2Source)[3]);
        setLinearBc2_el(0, 0, s.t_start, s.t_finish,
                        (*bc2Source)[4]);
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
}
void Test_Kirsch::writeResults(const Task &task, const OutData &out)
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
    fn_general_inf = dir + "_inf.txt";
    fn_sigma_fi = dir + "sigma_fi.txt";
    // пути к текстовым файлам с данными (для каждого глобального шага в отдельной директории)
    // пути к файлам для отображения графиков
    std::ofstream f;
    if(globalStepNumber == 0)
    {
        f.open(fn_filesList, std::ofstream::out);
        resGraph.generalGraphsNomber = 1;
        resGraph.graphsPerGlobalStepNomber = 0;
        f << resGraph.generalGraphsNomber << "\n";
        f << resGraph.graphsPerGlobalStepNomber << "\n";
        // sigma_fi
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "sigma_fi.png") +
             genGnuplotCommandLineParameter("fn_sigma_fi", dir + "sigma_fi.txt") +
             "\" " +
             dir + "sigma_fi.gnu" << "\n";
        f << dir + "sigma_fi.png" << "\n";
        f << "sigma_fi.png" << "\n";
    }
    else
    {
        f.open(fn_filesList, std::ofstream::app);
    }
    f.close();
    resGraph.read_gnuplot_from_file(fn_filesList);
    }

    FILE *fgeneral_inf = fopen(fn_general_inf.c_str(), "w");
    FILE *fsigma_fi = fopen(fn_sigma_fi.c_str(), "w");

    Grid::Grid3D *grid = task.grid;

    {
        int N_r = gp.N_r + 1; // количество узлов на 1 больше чем количество отрезков по координате r
        int N_fi = gp.N_fi;
        int N_z = gp.N_z;
        int feIndex = 0;
        for (int i = 0; i < N_fi; i++)              // по углу
        {
            for (int j = 0; j < N_r - 1; j++)       // по радиусу
            {
                for (int k = 0; k < N_z; k++)   // по толщине пластинки
                {
                    // КЭ около отверстия
                    if(j == 0)
                    {
                        MechOutFePointData &fePointData = (*out.mechOut)[globalStepNumber].fe[feIndex].pd[0];
                        MechMaterialSource &m0 = (*(*task.mechTask.mechStep)[globalStepNumber].material)[0];
                        POINT3 p = fePointData.p; // центр КЭ
                        p[2] = 0;
                        POINT3 p0(1, 0, 0);
                        double fi;
                        if(p[1] >= 0)
                        {
                            fi = acos((p0 * p) / (p0.abs() * p.abs()));
                        }
                        else
                        {
                            fi = 2*PI - acos((p0 * p) / (p0.abs() * p.abs()));
                        }
                        double sigma_fi =
                                  fePointData.sumSigma.m[0][0] * SQR(sin(fi))
                                + fePointData.sumSigma.m[1][1] * SQR(cos(fi))
                                - 2 * fePointData.sumSigma.m[0][1] * sin(fi) * cos(fi);
                        double sigma_fi_analit = 1 - 2*cos(2*fi);
                        double sigma_fi_residual = sigma_fi - sigma_fi_analit;
                        fprintf(fsigma_fi, "%le\t%le\t%le\t%le\n", fi, sigma_fi, sigma_fi_analit, sigma_fi_residual*10);
                        //fprintf(fsigma_fi, "%le\t%le\t%le\t%le\t%le\t%le\n", sigma_fi, sigma_fi_analit, sigma_fi_residual, p[0], p[1], p[2]);
                    }
                    feIndex++;
                }
            }
        }
    }


    /*
    gp.NU = 0.3;
    gp.E = 1.e10;
    gp.Lx = 10;
    gp.Ly = 10;
    gp.Lz = 10;
    gp.r = 1;
    gp.Px = 1;
    gp.Py = 0;
    gp.N_fi = 8;
    gp.N_r = 8;
    gp.N_z = 2;
    gp.q = 1;
    // сохренение информации о сетке в файл
    fprintf(f_inf, "r1 = %lf, r2 = %lf\n", gp.r1, gp.r2);

    fprintf(f_inf, "r1 = %lf, r2 = %lf\n", gp.r1, gp.r2);
    fprintf(f_inf, "Grid: %dx%dx%d, q = %le\n", gp.N, gp.N, gp.Nparts, gp.q);
    if(gp.curvilinear == 0)
        fprintf(f_parametors, "Linear\n");
    if(gp.curvilinear == 1)
        fprintf(f_parametors, "Quadratic\n");
    if(gp.buildingMethod == 0)
        fprintf(f_parametors, "buildingMethod = Divide arcs\n");
    if(gp.buildingMethod == 1)
        fprintf(f_parametors, "Quadratic = Cube to sphere\n");
    fprintf(f_parametors, "FENumber = %d\n", (int)grid->fe.size());
    fprintf(f_parametors, "VertexesNumber = %d\n", (int)grid->vertex.size());
    fprintf(f_parametors, "VertexesForCurvatureNumber = %d\n", (int)grid->vertexForCurvature.size());
    fprintf(f_parametors, "bc1Number = %d\n", (int)grid->bc1.size());
    fprintf(f_parametors, "FEsurface = %d\n", (int)grid->FESurface.size());
    for (size_t FEsurfaceInd = 0; FEsurfaceInd < grid->FESurface.size(); FEsurfaceInd++)
    {
        fprintf(f_parametors, "FEsurface[%d]Size = %d\n", (int)FEsurfaceInd, (int)grid->FESurface[FEsurfaceInd].face.size());
    }
    //grid->buldOpenSCADModel("setka.scad");
    // сохренение информации о параметрах в файл
    fprintf(f_inf, "fixGrid = %d\n", fixGrid);
    fprintf(f_inf, "plasticPlot = %d\n", plasticPlot);
    fprintf(f_parametors, "epsResidualLimit = %le\n", plasticResidualLimit.eps);
    fprintf(f_parametors, "nonlinearIterLimit = %d\n", nonlinearIterLimit);
    fprintf(f_inf, "StepsNumber1 = %d\n", StepsNumber1);
    fprintf(f_inf, "StepsNumber2 = %d\n", StepsNumber2);
    fprintf(f_inf, "StepsNumber3 = %d\n", StepsNumber3);
    fprintf(f_inf, "StepsNumber4 = %d\n", StepsNumber4);
    fprintf(f_inf, "P = %lf\n", P);
    fprintf(f_inf, "P1 = %lf\n", P1);
    fprintf(f_inf, "P2 = %lf\n", P2);
    fprintf(f_inf, "P3 = %lf\n", P3);
    fprintf(f_inf, "P4 = %lf\n", P4);
    fprintf(f_parametors, "incForsesMode = %d\n", (int)incForsesMode);
    fprintf(f_inf, "HomogenyMode = %d\n", HomogenyMode);
*/

    /*
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
    }*/

    fclose(fgeneral_inf);
    fclose(fsigma_fi);

    // построение графиков
    makeGraphsAfterGlobalStep(globalStepNumber, task);
}
bool Test_Kirsch::possibleToShow2d() const
{
    return true;
}
bool Test_Kirsch::getContactSurfaceCircle(const Task &task, const int globalStepIndex, POINT2 &c, double &R) const
{
    return false;
}
void Test_Kirsch::needToDrawFe(const Solid::Task &task, const OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const
{
    if(feInd%gp.N_z == gp.N_z / 2)
    {
        faceState = {0, 0, 1, 0, 0, 0};
    }
    else
    {
        faceState = {0, 0, 0, 0, 0, 0};
    }
    /*
    if(!(feInd%gp.N[2] == (size_t)z_index ||
            feInd >= gp.N[0] * gp.N[1] * gp.N[2]))
    {
        faceState = {0, 0, 0, 0, 0, 0};
    }
    else
    {
        faceState = {1, 0, 0, 0, 0, 0};
    }*/
}
bool Test_Kirsch::needToPaintFiniteElementSurface() const
{
    return false;
}

}

void calc_POINT3(double fi, double c, POINT3 &p, double a, double b, double r)
{
    POINT3 p1;
    POINT3 p2;
    p1[0] = r*cos(fi);
    p1[1] = r*sin(fi);
    // круглая пластинка
    //p2[0] = a*cos(fi);
    //p2[1] = b*sin(fi);
    // прямоугольная пластинка
    // верхняя сторона
    if(fi >= PI/4 && fi <= 3*PI/4)
    {
        p2[1] = b/2;
        p2[0] = p2[1]*cos(fi)/sin(fi);
    }else
    // левая сторона
    if(fi > 3*PI/4 && fi < 5*PI/4)
    {
        p2[0] = -a/2;
        p2[1] = p2[0]*sin(fi)/cos(fi);
    }else
    // нижняя сторона
    if(fi >= 5*PI/4 && fi <= 7*PI/4)
    {
        p2[1] = -b/2;
        p2[0] = p2[1]*cos(fi)/sin(fi);
    }else
    // правая сторона
    if(fi > 7*PI/4 || fi < PI/4)
    {
        p2[0] = a/2;
        p2[1] = p2[0]*sin(fi)/cos(fi);
    }
    else
    {
        p2[0] = 0;
        p2[1] = 0;

    }
    // задан отрезок p1-p2
    // нужно взять от него часть c (0..1)
    p[0] = p1[0] + (p2[0] - p1[0])*c;
    p[1] = p1[1] + (p2[1] - p1[1])*c;
}

namespace Grid
{
void Grid3D::genKirch(const KirschParameters &kp)
{
    int N_r = kp.N_r + 1; // количество узлов на 1 больше чем количество отрезков по координате r
    int N_fi = kp.N_fi;
    int N_z = kp.N_z;
    double fi;		// полярные координаты
    double t = 0;	// доля отрезка которую отсечет очередная точка
    double b;
    POINT3 p;
    fe.clear();
    vertex.clear();
    vertexForCurvature.clear();
    bc1.clear();
    FESurface.clear();

    // вершины
    b = 1 / (pow(kp.q, N_r - 1) - 1);
    for (int i = 0; i < N_fi; i++)  // индекс по углу
    {
        fi = 2 * PI*i / N_fi;	// угол
        t = b;
        for (int j = 0; j < N_r; j++)   // индекс по радиусу
        {
            if (kp.q == 1)
            {
                t = (double)j / (N_r - 1);	// равномерное разбиение
                calc_POINT3(fi, t, p, kp.Lx, kp.Ly, kp.r);		// найдена точка на отрезке
            }
            else
            {								// геометрическая прогрессия
                calc_POINT3(fi, t - b, p, kp.Lx, kp.Ly, kp.r);	// найдена точка на отрезке
                t = b*pow(kp.q, j+1);
                //t *= kp.q;					// увеличиваем шаг
            }

            for (int k = 0; k <= N_z; k++)   // индекс по толщине пластинки
            {
                POINT3 vertex_el;
                vertex_el[0] = p[0];
                vertex_el[1] = p[1];
                vertex_el[2] = -kp.Lz / 2 + kp.Lz * (((double)k)/N_z);
                if(k == N_z/2 && N_z%2 == 0)
                {
                    vertex_el[2] = 0;
                }
                vertex.push_back(vertex_el);
            }
        }
    }
    // 6-гранники и краевые условия 2 рода
    // будет 5 поверхностей
    FESurface.resize(5);
    for (size_t FEsurfaceInd = 0; FEsurfaceInd < FESurface.size(); FEsurfaceInd++)
        FESurface[FEsurfaceInd].face.clear();
    int offset_r = N_z + 1;     // прибавление смещает на 1 по r
    int offset_fi = N_r * offset_r;	// прибавление offset переводит на индекс той же точки на следующем отрезке (для чуть большего fi)

    for (int i = 0; i < N_fi; i++)              // по углу
    {
        for (int j = 0; j < N_r - 1; j++)       // по радиусу
        {
            // ind1,2 Для z=min, r=min
            int ind1 = j * offset_r + i*offset_fi;         // индекс текущей вершины
            int ind2 = j * offset_r + (i + 1)*offset_fi;   // индекс вершины, смещенной от текущей на 1 шаг по fi
            if (i == N_fi - 1)
                ind2 = j * offset_r;    // (замыкание по углу)
            for (int k = 0; k < N_z; k++)   // по толщине пластинки
            {
                int vi[8] =
                {
                    ind1 + k,            ind1 + (k+1),          ind2 + k,          ind2 + (k+1),
                    ind1 + k + offset_r, ind1 + (k+1)+offset_r, ind2 + k+offset_r, ind2 + (k+1)+offset_r
                };
                // шестигранник
                FE_LinearHexagon *fe_el = new FE_LinearHexagon;
                fe_el->mi = 0;
                for (int t = 0; t < 8; t++)
                    fe_el->vi[t] = vi[t];
                fe.push_back(fe_el);
                // второе краевое условие
                FEFace face;
                if (j == 0)
                {// сила и грань для 2 краевого изнутри

                    face.feIndex = (int)fe.size() - 1;
                    face.faceIndex = 1;   // Z = 1
                    FESurface[4].face.push_back(face);
                }
                if (j == N_r - 2 && (i >= N_fi / 8 * 1 && i < N_fi / 8 * 3))
                {// сила и грань для 2 краевого сверху

                    face.feIndex = (int)fe.size() - 1;
                    face.faceIndex = 0;   // Z = -1
                    FESurface[2].face.push_back(face);
                }
                if (j == N_r - 2 && (i >= N_fi / 8 * 3 && i < N_fi / 8 * 5))
                {// сила и грань для 2 краевого слева
                    face.feIndex = (int)fe.size() - 1;
                    face.faceIndex = 0;   // Z = -1
                    FESurface[0].face.push_back(face);
                }
                if (j == N_r - 2 && (i >= N_fi / 8 * 5 && i < N_fi / 8 * 7))
                {// сила и грань для 2 краевого снизу
                    face.feIndex = (int)fe.size() - 1;
                    face.faceIndex = 0;   // Z = -1
                    FESurface[3].face.push_back(face);
                }
                if (j == N_r - 2 && (i >= N_fi / 8 * 7 || i < N_fi / 8 * 1))
                {// сила и грань для 2 краевого справа

                    face.feIndex = (int)fe.size() - 1;
                    face.faceIndex = 0;   // Z = -1
                    FESurface[1].face.push_back(face);
                }
            }
        }
    }
    // первые краевые
    BoundaryCondition1 bc1_el;
    // фиксируем x на центральной вертикали
    for (int i = 0; i < N_r * offset_r; i++)
    {
        bc1_el.bc1SourceIndex = 0;
        bc1_el.vertexIndex = (N_fi * 1 / 4)*offset_fi + i;
        bc1.push_back(bc1_el);
    }
    for (int i = 0; i < N_r * offset_r; i++)
    {
        bc1_el.bc1SourceIndex = 0;
        bc1_el.vertexIndex = (N_fi * 3 / 4)*offset_fi + i;
        bc1.push_back(bc1_el);
    }
    // фиксируем y на центральной горизонтале
    for (int i = 0; i < N_r * offset_r; i++)
    {
        bc1_el.bc1SourceIndex = 1;
        bc1_el.vertexIndex = (N_fi * 0 / 4)*offset_fi + i;
        bc1.push_back(bc1_el);
    }
    for (int i = 0; i < N_r * offset_r; i++)
    {
        bc1_el.bc1SourceIndex = 1;
        bc1_el.vertexIndex = (N_fi * 2 / 4)*offset_fi + i;
        bc1.push_back(bc1_el);
    }
    // фиксируем z во всех вершинах
    for(size_t i = 0; i < vertex.size(); i++)
    {
        bc1_el.bc1SourceIndex = 2;
        bc1_el.vertexIndex = (int)i;
        bc1.push_back(bc1_el);
    }
    // фиксируем z по центру
    /*
    {
        int vertexIndex = 0;
        for (int i = 0; i < N_fi; i++)  // индекс по углу
        {
            for (int j = 0; j < N_r; j++)   // индекс по радиусу
            {
                for (int k = 0; k <= N_z; k++)   // индекс по толщине пластинки
                {
                    if(k == N_z / 2)
                    {
                        bc1_el.bc1SourceIndex = 2;
                        bc1_el.vertexIndex = vertexIndex;
                        bc1.push_back(bc1_el);
                    }
                    vertexIndex++;
                }
            }
        }
    }
    */
    // tt
    /*for(size_t i = 0; i < vertex.size(); i++)
    {
        if(vertex[i][2] == 0)
        {
            bc1_el.bc1SourceIndex = 2;
            bc1_el.vertexIndex = (int)i;
            bc1.push_back(bc1_el);
        }
    }*/
    /*
    BoundaryCondition1 bc1_el;
    // фиксируем x на центральной вертикали
    for (i = 0; i < N_r * 2; i++)
    {
        bc1_el.si = 0;
        bc1_el.vi = (N_fi * 1 / 4)*offset + i;
        bc1.push_back(bc1_el);
    }
    for (i = 0; i < N_r * 2; i++)
    {
        bc1_el.si = 0;
        bc1_el.vi = (N_fi * 3 / 4)*offset + i;
        bc1.push_back(bc1_el);
    }
    // фиксируем y на центральной горизонтале
    for (i = 0; i < N_r * 2; i++)
    {
        bc1_el.si = 1;
        bc1_el.vi = (N_fi * 0 / 4)*offset + i;
        bc1.push_back(bc1_el);
    }
    for (i = 0; i < N_r * 2; i++)
    {
        bc1_el.si = 1;
        bc1_el.vi = (N_fi * 2 / 4)*offset + i;
        bc1.push_back(bc1_el);
    }
    // фиксируем z во всех вершинах
    for(i = 0; i < N_r*N_fi*2; i += 1)
    {
        bc1_el.si = 2;
        bc1_el.vi = i;
        bc1.push_back(bc1_el);
    }*/
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
