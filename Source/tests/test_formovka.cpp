#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"
#include <fstream>

#include "tests.h"
#include "console.h"

#include "interpolation.h"

using namespace Elementary;
using namespace Solid;
using namespace Grid;

void setBcFormovkaT(std::vector<Heat::ThermBoundaryCondition1Source> *&thermBc1, std::vector<Heat::ThermBoundaryCondition2Source> *&thermBc2,
                    double q)
{
    using namespace Heat;
    // первые краевые условия
    thermBc1 = new std::vector<ThermBoundaryCondition1Source>(5);
    (*thermBc1)[0] = {-1, -1};
    (*thermBc1)[1] = {-1, -1};
    (*thermBc1)[2] = {-1, -1};
    (*thermBc1)[3] = {0, 0};        // снизу
    (*thermBc1)[4] = {0, 0};        // сверху
    // вторые краевые условия
    thermBc2 = new std::vector<ThermBoundaryCondition2Source>(2);
    {
        ThermBoundaryCondition2Source &bc2_el = (*thermBc2)[1];
        bc2_el.mode = ThermBoundaryCondition2Type::scalar;    // к/у присутствует
        bc2_el.hi = 0;      // коэффициент конвективного теплообмена
        bc2_el.Ta = 0;      // температура окружающей среды
        bc2_el.q = q;       // текущая плотность набегающего теплового потока
    }
    {
        ThermBoundaryCondition2Source &bc2_el = (*thermBc2)[0];
        bc2_el.mode = ThermBoundaryCondition2Type::none;
    }
}

namespace Tests
{
Test_formovka::Test_formovka()
{
    // дирректория данных теста
    dir = "./tests/formovka/";
    // имена gnuplot файлов(+пути к файлам в опциях) и соответствующих картинок
    fn_filesList = dir + "files.txt";
    resGraph.read_gnuplot_from_file(fn_filesList);
    // инициализация графиков с невязками
    initStepsGraphs();

}
Test_base::Type Test_formovka::get_type() const
{
    return Type::Formovka;
}
void Test_formovka::initTask(Solid::Task &task)
{
    // подрубка
    task.mechTask.enabled = true;
    task.thermTask.enabled = false;
    // общие шаги и сетка
    std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
    Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
    // способ представления поверхности
    //gp.surfaceType = SurfaceType::none;
    gp.surfaceType = SurfaceType::AnaliticalSurface_Cylinder;
    //gp.surfaceType = SurfaceType::InterpolantSurface_Hermite3;
    //gp.surfaceType = SurfaceType::InterpolantSurface_Lagrange3;
        //gp.contact_mode = SurfaceType::FiniteElementSurface;//нет поиска пересечения
    // способ декомпозиции поверхности
    gp.decompositionType = RegionDecompositionType::none;
    //gp.decompositionType = DecompositionType::cubeVertexIndexation;
    //gp.decompositionType = DecompositionType::surfaceAsSpheres;

    // подрубка универсальной оптимизации
    bool noContactRadiusOptimization = false;
    //bool noContactRadiusOptimization = true;


    int NN = 13;//7; //9//5
    gp.p1 = VECTOR3(-10, -1, -1.0);
    gp.p2 = VECTOR3(10, -0.5, +1.0);
    //gp.p1 = VECTOR3(-10, -1, -1);
    //gp.p2 = VECTOR3(10, -0.5, -0.5);
    gp.N = VECTOR3_uint(NN*40, NN, 2);//8
    gp.N_react = 0;//(NN+1)/2;
    gp.N_p = (NN+1)/2;
    /*
    gp.p1 = VECTOR3(-10, -1, -1);
    gp.p2 = VECTOR3(10, 1, 2);
    gp.N = VECTOR3_int(90, 9, 2);
    gp.N_react = 5;
    gp.N_p = 5;
    */

    // параметры подвижного цилиндра
    gp.R0 = 5;
    gp.y0 = -1 - gp.R0;
    gp.R_fun = new Interpolation::Interpolant1D_Lagrange1Unregular;
    gp.x_fun = new Interpolation::Interpolant1D_Lagrange1Unregular;
    gp.y_fun = new Interpolation::Interpolant1D_Lagrange1Unregular;

    Grid::GridRectangleRegular1D gr;
    gp.R_fun->init(gr, 0);// аргументы игнорируются для нерегулярной сетки
    gp.x_fun->init(gr, 0);
    gp.y_fun->init(gr, 0);
    gp.R_fun->addPoint(0, gp.R0);
    gp.R_fun->addPoint(1000, gp.R0);
    gp.R_fun->buildInterpolant();
    gp.x_fun->addPoint(0, 0);
    gp.x_fun->addPoint(1000, 0);
    gp.x_fun->buildInterpolant();
    gp.y_fun->addPoint(0, gp.y0);
    gp.y_fun->addPoint(1000, gp.y0);
    gp.y_fun->buildInterpolant();

    grid->genFormovka(gp);

    //for(;;);
    //grid->buldOpenSCADModel("setka.scad");
    //for(;;);
    // параметры
    int fixGrid = 0;                // 0 - подвижная сетка, 1 - зафиксировать сетку
    MaxSigmaEpsResidual plasticResidualLimit;   // желаемая невязка по эквиволентным деформациям
    double contactDeltaFResidualLimit = 1.e-7;   // желаемая относительная погрешность по силе реакции опоры
    double contactEndPointResidualLimit = 1.e-13;
    double w_stiffness;
    Solid::ContactType constantMetod;
    bool contactNormal;
    //double stiffness_min;            // минимальная(начальная) контактная жёсткость
    //double stiffness_max;            // минимальная(начальная) контактная жёсткость
    //double h_max = 1.e-7;                // ограничение сверху расстояния от узла до поверхности в случае контакта
    int nonlinearIterLimit = 500;        // ограничение на количество итераций
    int StepsNumber1 = 1;        //
    int StepsNumber2 = 40;    // механическое нагружение
    int StepsNumber3 = 1;        // механическая разгрузка
  //int StepsNumber4 = 1;        // температурное нагружение
    // пластичность: 500 шагов, stiffness = 1.e4, P2 = 8.e5*1.5 - вроде в конце падает а мб и нет
    double P1 = 1.e-10;
    double P2 = 1.e7;//8.e5*1.5;//0.10e8/2;//0.20e8/2;//-0.20e8/2;//-0.20e8/3.34;//-0.20e8/3.7(если 400 шагов);//-0.20e8/3.36;(если 40 шагов)  //-0.05e8;//-0.08e7;//-0.77e7;//-0.075e8;//-0.075e8;//-1.5e8;//-1.1e8;//-2.e8;//-1.e7*1;//-5.e7*1; //-5.e7/20*1; - для пластичности
    double Time0;
    double Time1;
    double Time2;
    // 9, 100, -0.20e8/3.34, 2.e7*10
    double elasticSigmaLimit = 2.e7;
    //double P2 = -0.20e8;
    //double elasticSigmaLimit = 2.e7*16;
  //double P3 = P2;
    int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny




    const int creepCurve = 0;
    const int plasticCurve = 1;

    // кривая (ползучесть/пластичность)
    //int creepCurveMode = creepCurve;
    int creepCurveMode = plasticCurve;

    // параметры кривой ползучести
    //gp.n = 3;
    //gp.n = 2;
    gp.n = 1.5;   // 1 - упругость, ->бесконечности - пластичность
    gp.B = 1.e-14;

    // метод (начальных напряжений, начальных деформаций, пластичность)
    //MechPlasticityMethodType creepMode = MechPlasticityMethodType::D_pl;
    //MechPlasticityMethodType creepMode = MechPlasticityMethodType::InitialEps;
    //MechPlasticityMethodType creepMode = MechPlasticityMethodType::InitialSigma;
    MechPlasticityMethodType creepMode = MechPlasticityMethodType::Combo_D_pl_InitialSigma;

    //IncForsesMode incForsesMode = IncForsesMode::bPlusR;
    //IncForsesMode incForsesMode = IncForsesMode::IncrementP;
    IncForsesMode incForsesMode = IncForsesMode::MinusIntegral;

    bool terminateIfAccuracyIsNotAchieving = false;
    double w = 1;
    int controlMode = 0; // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг


    // пластичность (на пределе)
    /*
incForsesMode = IncForsesMode::MinusIntegral;
    creepCurveMode = plasticCurve;
    MechPlasticityMethodType creepMode = MechPlasticityMethodType::Combo_D_pl_InitialSigma;
    //elasticSigmaLimit = 2.e100;
    elasticSigmaLimit = 2.e7;
    gp.n = 1.5;
    //plasticResidualLimit.sigma = 1.e-14;//1.e-13;
    //plasticResidualLimit.eps = 100000;
    plasticResidualLimit.sigma = 1.e-10;//1.e-13;
    plasticResidualLimit.eps = 1000;
    nonlinearIterLimit = 1000;
    Time0 = 0;
    Time1 = 1.e-7;
    Time2 = 10;//10
    StepsNumber1 = 1;
    StepsNumber2 = 80;//'100
    StepsNumber3 = 1;
    P1 = 0;
    P2 = 5.e5;
    */

    // ползучесть
/*
incForsesMode = IncForsesMode::MinusIntegral;
    creepCurveMode = creepCurve;
    creepMode = MechPlasticityMethodType::Combo_D_pl_InitialSigma;
    //stiffness0 = 1.e5;
    stiffness_min = 1.e10;
    stiffness_max = 1.e10;
    h_max = 1.e-7;
    contactDeltaFResidualLimit = 1.e-7;   // желаемая относительная погрешность по силе реакции опоры
    contactEndPointResidualLimit = 1.e-12;
    elasticSigmaLimit = 2.e100;
    gp.n = 1.5;
    //plasticResidualLimit.sigma = 1.e-14;//1.e-13;
    //plasticResidualLimit.eps = 100000;
    plasticResidualLimit.sigma = 1.e-10;//1.e-13;
    plasticResidualLimit.eps = 1000;
    nonlinearIterLimit = 1000;
    Time0 = 0;
    Time1 = 1.e-7;
    Time2 = 500;//500
    StepsNumber1 = 1;
    StepsNumber2 = 200;//200
    StepsNumber3 = 1;
    P1 = 1.e5;
    P2 = 1.e5;
*/


    // упругость
    creepCurveMode = plasticCurve;
    constantMetod = Solid::ContactType::AugmentedLagrange;
    //contactType = Solid::ContactType::Fastest;
    contactNormal = true;//false;
    creepMode = MechPlasticityMethodType::Elasticity;
    //creepMode = MechPlasticityMethodType::D_pl;
    creepMode = MechPlasticityMethodType::Combo_D_pl_InitialSigma;
    w_stiffness = 10;
    contactDeltaFResidualLimit = 1.e-8;//1.e-6;//1;//1.e-7;   // желаемая относительная погрешность по силе реакции опоры
    contactEndPointResidualLimit = 1.e-11;
    //elasticSigmaLimit = 2.e100;
    //elasticSigmaLimit = 2.e7*16;
    gp.n = 1.5;
    plasticResidualLimit.sigma = 1.e-10;//1.e-13;
    plasticResidualLimit.eps = 1000;
    nonlinearIterLimit = 1000;
    Time0 = 0;
    Time1 = 1.e-7;
    Time2 = 500;//10
    StepsNumber1 = 1;
    StepsNumber2 = 200;//50;//'100
    StepsNumber3 = 1;
    P1 = 0;
    P2 = 2.e7;//1.5e7;

/*
    double epsResidualLimit = 1.e-9;            // желаемая относительная погрешность по эквиволентным деформациям
    double contactDeltaFResidualLimit = 1.e-7;   // желаемая относительная погрешность по силе реакции опоры
    double stiffness = 1.e11;            // контактная жёсткость
    int nonlinearIterLimit = 200;        // ограничение на количество итераций
    int StepsNumber1 = 1;        //
    int StepsNumber2 = 40;    // механическое нагружение
    double P1 = 0;
    double P2 = 1.e7;//8.e5*1.5;//0.10e8/2;//0.20e8/2;//-0.20e8/2;//-0.20e8/3.34;//-0.20e8/3.7(если 400 шагов);//-0.20e8/3.36;(если 40 шагов)
*/




    // шаги
    GlobalStep s;       // шаг
    // 1-й шаг (нагревание)
    s.t_start = Time0;                          // (имеет значение только для первого глобального шага)
    s.t_finish = Time1;
    s.dt0 = (Time1 - Time0)/StepsNumber1;
    step->push_back(s);
    // 2-й шаг (нагружение P)
    for(int i = 0; i < StepsNumber2; i++)
    {
        s.t_start = Time1 + (Time2 - Time1)*i/StepsNumber2;
        s.t_finish = Time1 + (Time2 - Time1)*(i + 1)/StepsNumber2;
        s.dt0 = (Time2 - Time1)/StepsNumber2;
        step->push_back(s);
    }
    /*
    // 3-й шаг (разгрузка)
    for(int i = 0; i < StepsNumber3; i++)
    {
        s.t1 = Time*2 + Time*i/StepsNumber3;
        s.t2 = Time*2 + Time*(i + 1)/StepsNumber3;
        s.dt0 = Time/StepsNumber3;
        step->push_back(s);
    }
    // 4-й шаг (остужение)
    s.t1 = Time*3;
    s.t2 = Time*4;
    s.dt0 = Time/StepsNumber4;
    step->push_back(s);
    */
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
        mechMat.set_E_NU(1.e10,0.3);
        //mechMat.set_E_NU(1.e10,0.4);
        //mechMat.set_D_isotropic_XY();
        //mechMat.set_M_sigma_XY();
        //mechMat.set_D_isotropic();
        //mechMat.set_M_sigma();
        mechMat.set_M_sigma();       //mechMat.set_M_sigma_XY();       // эквивалентные напряжения плоские
        mechMat.F = {0,0,0};
        mechMat.elasticParameters0.ro = 0;     // Плотность
        mechMat.elasticParameters0.Talpha = 0;//1.e-5; // Коэффициент линейного расширения



        double k_sigma_eps = 1000;
        double k_eps_sigma = 10;
        // кривая
        if(creepCurveMode == plasticCurve)
        {
            mechMat.elasticSigmaLimit = 1.e7*10;//2.e7*16;
            Tests::setPlasticMaterialCurve_Yeld(1, mechMat, k_sigma_eps, k_eps_sigma);
            //setPlasticCurve(mechMat, k_sigma_eps, k_eps_sigma);
            mechMat.PCDependenceType = MechPlasticityCurveDependenceType::None;
        }
        if(creepCurveMode == creepCurve)
        {
            setCreepMaterialCurve(gp.B, gp.n, 1, mechMat);
            mechMat.PCDependenceType = MechPlasticityCurveDependenceType::Time;
        }
        // метод
        if(creepMode == MechPlasticityMethodType::Elasticity)
        {
            // упругость
            mechMat.plasticityMethodType = creepMode;
            mechMat.PCUnloadingType = MechPlasticityCurveUnloadingType::Curve;
        }
        if(creepMode == MechPlasticityMethodType::InitialSigma ||
           creepMode == MechPlasticityMethodType::Combo_D_pl_InitialSigma)
        {
            // метод начальных напряжений
            mechMat.plasticityMethodType = creepMode;
            mechMat.PCUnloadingType = MechPlasticityCurveUnloadingType::Curve;
            mechMat.PCType = MechPlasticityCurveType::Sigma_eps;
            mechMat.w_midPoint = 0.5;
            mechMat.w_project = w;
        }
        if(creepMode == MechPlasticityMethodType::D_pl)
        {
            // пластичность
            mechMat.plasticityMethodType = creepMode;
            mechMat.PCUnloadingType = MechPlasticityCurveUnloadingType::Curve;
            mechMat.PCType = MechPlasticityCurveType::Sigma_eps;
            mechMat.w_midPoint = 0.5;
            mechMat.w_project = w;
        }


        // первые краевые условия
        bc1Source = new std::vector<MechBoundaryCondition1Source>(5);
        (*bc1Source)[0].mode = {{ 0, -1, -1}};
        (*bc1Source)[0].u0 =   {{ 0, -1, -1}};
        //(*bc1Source)[1].mode = {{0,  0, 0}};
        //(*bc1Source)[1].u0 =   {{0,  0, 0}};
        (*bc1Source)[1].mode = {{-1,  0, -1}};
        (*bc1Source)[1].u0 =   {{-1,  0, -1}};
        //(*bc1Source)[2].mode = {{-1, -1,  -1}};
        //(*bc1Source)[2].u0 =   {{-1, -1,  -1}};
        (*bc1Source)[2].mode = {{-1, -1,  0}};
        (*bc1Source)[2].u0 =   {{-1, -1,  0}};

        (*bc1Source)[3].mode = {{-1, -1,  -1}};
        (*bc1Source)[3].u0 =   {{-1, -1,  -1}};
        (*bc1Source)[4].mode = {{-1, -1,  -1}};
        (*bc1Source)[4].u0 =   {{-1, -1,  -1}};
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
        // нагревание - нагружение - разгрузка
        MechGlobalStepParameters ms;  // шаг для МДТТ
        GlobalStep s;       // шаг
        ms.slausolverParameters = slausolver_parameters;
        ms.timeMode = 0;    // 0 - квазистатическая задача
        ms.incForsesMode = incForsesMode; // 1 - dSigma = integral(...), 2 - приращения
        ms.fixGrid = fixGrid;     // 1 - зафиксировать сетку
        //ms.switchIterationsMode = SwitchIterationsMode::Serial;
        ms.switchIterationsMode = SwitchIterationsMode::Parallel;
        ms.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
        ms.controlMode = controlMode;     // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
        ms.terminateIfAccuracyIsNotAchieving = terminateIfAccuracyIsNotAchieving;
        ms.iterLimit = nonlinearIterLimit;
        ms.slauResidualLimit = 10000;
        ms.plasticResidualLimit = plasticResidualLimit;
        ms.contactEndPointResidualLimit = contactEndPointResidualLimit;//1.e-10;//1.e-13;
        ms.contactDeltaFResidualLimit = contactDeltaFResidualLimit;//1.e-13;
        ms.material = material;
        ms.bc1Source = bc1Source;
        //ms.bc2Source = bc2Source; - будет меняться
        ms.bc2 = bc2;
        // 1-й шаг (нагревание)
        s = (*step)[0];
        setBc2Sphere(bc2Source, 0, P1, s.t_start, s.t_finish);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);
        // 2-й шаг (нагружение)
        for(int i = 0; i < StepsNumber2; i++)
        {
            s.t_start = Time1 + (Time2 - Time1)*i/StepsNumber2;
            s.t_finish = Time1 + (Time2 - Time1)*(i + 1)/StepsNumber2;
            setBc2Sphere(bc2Source, P1 + (P2 - P1)*i/StepsNumber2, P1 + (P2 - P1)*(i + 1)/StepsNumber2, s.t_start, s.t_finish);
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
        }
        /*
        // 3-й шаг (разгрузка)
        for(int i = 0; i < StepsNumber3; i++)
        {
            s.t1 = Time*2 + Time*i/StepsNumber3;
            s.t2 = Time*2 + Time*(i + 1)/StepsNumber3;
            //setBc2Sphere(bc2Source, P2, P2, s.t1, s.t2);
            setBc2Sphere(bc2Source, P2 + (P3 - P2)*i/StepsNumber3, P2 + (P3 - P2)*(i + 1)/StepsNumber3, s.t1, s.t2);
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
        }
        // 4-й шаг (остужение)
        s = (*step)[1 + StepsNumber2 + StepsNumber3];
        setBc2Sphere(bc2Source, 0, 0, s.t1, s.t2);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);
        */
        task.mechTask.mechStep = mechStep;


        /*
        // инициализация начальных деформаций и напряжений (=0)
        task.mechTask.fe = new std::vector<MechFeData_base *>;
        task.mechTask.fe->resize(grid->fe.size());
        for (size_t i = 0; i < grid->fe.size(); i++)
        {
            MechFeData_LinearHexagon_homogeny_E_NLE_EP *feDataEl = new MechFeData_LinearHexagon_homogeny_E_NLE_EP;
            feDataEl->init(Integration::IntegrationType::Gauss3);
            (*task.mechTask.fe)[i] = feDataEl;
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
        // для криволинейных границ ###
        task.mechTask.vertexForCurvature = new std::vector<MechVertexData>;
        for (size_t i = 0; i < grid->vertexForCurvature.size(); i++)
        {
            MechVertexData v_el;
            for (int j = 0; j < 3; j++)
            {
                v_el.sum_du[j] = 0;
                v_el.du[j] = 0; // не имеет значения
            }
            task.mechTask.vertexForCurvature->push_back(v_el);
        }
        // инициализация начальных скоростей (=0) и ускорений (=0)
        task.mechTask.V0 = new Vector(grid->vertex.size()*3);
        task.mechTask.dV0 = new Vector(grid->vertex.size()*3);
        for (size_t i = 0; i < task.mechTask.V0->size(); i++)
        {
            (*task.mechTask.V0)[i] = 0;
            (*task.mechTask.dV0)[i] = 0;
        }*/


        // Контакт

        // поверхность: жёсткий цилиндр
        task.mechTask.rigidSurface = new std::vector<Surface_base *>;
        (*task.mechTask.rigidSurface).resize(3);
        (*task.mechTask.rigidSurface)[0] = grid->analiticalSurface[0];   // поверхность заданная аналитически
        (*task.mechTask.rigidSurface)[1] = grid->ISurface[0];            // поверхность заданная интерполянтом
        (*task.mechTask.rigidSurface)[2] = &(grid->FESurface[3]);        // поверхность заданная сеткой (если contact_mode != 3 то эта дополнительная сетка не строится)



        // 1)цилиндр задан аналитически
        //AnaliticalSurface_Cylinder *cylinder = new AnaliticalSurface_Cylinder;
        //cylinder->initCylinderZ({0, -1 - 5, 0}, 5);
        //(*task.mechTask.rigidSurface).push_back(cylinder);

        // 2) цилиндр задан сеткой
        //grid->FEsurface[3].init(true);  // неподвижна - обноврять декомпозицию пространства не нужно
        //(*task.mechTask.rigidSurface).push_back(&grid->FESurface[3]);


        // контакты КЭ-поверхность - жёсткая поверхность
        task.mechTask.CC_FE_Rigid = new std::vector<ContactCondition_FE_Rigid>;
        if(gp.surfaceType != SurfaceType::none)
        {
            ContactCondition_FE_Rigid CC_FE_Rigid_el;
            CC_FE_Rigid_el.w_stiffness = w_stiffness;
            CC_FE_Rigid_el.method = constantMetod;
            CC_FE_Rigid_el.constantNormal = contactNormal;
            CC_FE_Rigid_el.FESurfaceInd = 2;
            if(gp.surfaceType == SurfaceType::AnaliticalSurface_Cylinder)
                CC_FE_Rigid_el.RigidSurfaceInd = 0;
            if(gp.surfaceType == SurfaceType::InterpolantSurface_Hermite3 ||
               gp.surfaceType == SurfaceType::InterpolantSurface_Lagrange3)
                CC_FE_Rigid_el.RigidSurfaceInd = 1;
            if(gp.surfaceType == SurfaceType::FiniteElementSurface)
                CC_FE_Rigid_el.RigidSurfaceInd = 2;
            CC_FE_Rigid_el.noContactRadiusOptimization = noContactRadiusOptimization;
            (*task.mechTask.CC_FE_Rigid).push_back(CC_FE_Rigid_el);
        }

        // инициализация нулями начальных данных
        task.mechTask.initNull(Integration::IntegrationType::Gauss3);

        // задание 1-ной контактной точки по центру, чтобы не требовались 1-е краевые по веритикали
        for(size_t CV_FE_RigidInd = 0; CV_FE_RigidInd < task.mechTask.contactVertex_FE_Rigid.size(); CV_FE_RigidInd++)
        {
            ContactSurfaceVertexData_FE_Rigid &csv = task.mechTask.contactVertex_FE_Rigid[CV_FE_RigidInd];
            POINT3 &x = grid->vertex[csv.vertexIndex];
            if(x[0] == 0 && fabs(x[2] - 0) < 1.e-6)
            {
                POINT3 x2_nearestPoint;
                POINT3 x2_normal;
                int x2_side;
                bool x2_onBorder;
                Grid::Surface_base *sb = (*task.mechTask.rigidSurface)[(*task.mechTask.CC_FE_Rigid)[csv.si].RigidSurfaceInd];
                sb->findNearestPoint(x, x, csv.ssd2,
                                     csv.ssd2, x2_nearestPoint, x2_normal, x2_side, x2_onBorder);
                csv.contact = true;
                csv.it.contactNew = true;
                csv.it.updateEndPoint(x, x);
                //csv.it.setNormal(x2_normal);
                csv.it.normal = x2_normal;
                csv.it_new = csv.it;
            }
        }

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
        thermMat.ro = 2200;//1;//2200;     // плотность
        thermMat.c = 1;//5;//1000;      // удельная теплоёмкость
        thermMat.f = 0;        // мощность внутренних объёмных источников (стоков) тепла
        double L0 = 10;//0.7;//10;
        thermMat.L.initDiag(L0); // тензор теплопроводности
        int timeMode = 1;       // 0 - статика, 1 - динамика
        setBcFormovkaT(bc1SourceT, bc2SourceT,
                       10000000*98 // плотность набегающего теплового потока
                       );
        // шаги
        std::vector<ThermGlobalStep> *thermStep = new std::vector<ThermGlobalStep>;
        ThermGlobalStep ts;         // шаг для теплопроводности
        ts.slausolverParameters = slausolver_parameters;
        ts.timeMode = timeMode;     // статика/динамика
        ts.material = material;
        ts.bc1SourceT = bc1SourceT;
        ts.bc2SourceT = bc2SourceT;
        thermStep->push_back(ts);
        for(int i = 0; i < StepsNumber2; i++)
            thermStep->push_back(ts);
        for(int i = 0; i < StepsNumber3; i++)
            thermStep->push_back(ts);
        thermStep->push_back(ts);
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
            v_el.T = T0;
            v_el.newT = T0;
            task.thermTask.vertex->push_back(v_el);
        }

        // ## 2-е краевые нужно переделать
    }
}
void Test_formovka::writeResults(const Task &task, const Solid::OutData &out)
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
    fn_F = subdir + "F.txt";
    fn_h = subdir + "h.txt";
    fn_sigmaTop1 = subdir + "sigmaTop1.txt";
    fn_sigmaTop2 = subdir + "sigmaTop2.txt";
    fn_sigmaTop3 = subdir + "sigmaTop3.txt";
    fn_sigmaMiddle1 = subdir + "sigmaMiddle1.txt";
    fn_sigmaMiddle2 = subdir + "sigmaMiddle2.txt";
    fn_sigmaMiddle3 = subdir + "sigmaMiddle3.txt";
    fn_sigmaBottom1 = subdir + "sigmaBottom1.txt";
    fn_sigmaBottom2 = subdir + "sigmaBottom2.txt";
    fn_sigmaBottom3 = subdir + "sigmaBottom3.txt";
    fn_T = subdir + "T.txt";
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
    // F
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "F.png") +
         genGnuplotCommandLineParameter("fn_F", subdir + "F.txt") +
         "\" " +
         dir + "F.gnu" << "\n";
    f << subdir + "F.png" << "\n";
    f << "F.png" << "\n";
    // h
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "h.png") +
         genGnuplotCommandLineParameter("fn_h", subdir + "h.txt") +
         "\" " +
         dir + "h.gnu" << "\n";
    f << subdir + "h.png" << "\n";
    f << "h.png" << "\n";
    f.close();
}
    resGraph.read_gnuplot_from_file(fn_filesList);
    }

    //task.grid->buldOpenSCADModel("setka_.scad");
    FILE *f_F = fopen(fn_F.c_str(), "w");
    FILE *f_h = fopen(fn_h.c_str(), "w");
    FILE *f_sigmaTop1 = fopen(fn_sigmaTop1.c_str(), "w");
    FILE *f_sigmaTop2 = fopen(fn_sigmaTop2.c_str(), "w");
    FILE *f_sigmaTop3 = fopen(fn_sigmaTop3.c_str(), "w");
    FILE *f_sigmaMiddle1 = fopen(fn_sigmaMiddle1.c_str(), "w");
    FILE *f_sigmaMiddle2 = fopen(fn_sigmaMiddle2.c_str(), "w");
    FILE *f_sigmaMiddle3 = fopen(fn_sigmaMiddle3.c_str(), "w");
    FILE *f_sigmaBottom1 = fopen(fn_sigmaBottom1.c_str(), "w");
    FILE *f_sigmaBottom2 = fopen(fn_sigmaBottom2.c_str(), "w");
    FILE *f_sigmaBottom3 = fopen(fn_sigmaBottom3.c_str(), "w");
    FILE *f_T = fopen(fn_T.c_str(), "w");
    size_t i0, i1, i2;

    //for (int globalStepNumber = 0; globalStepNumber < (int)task.step->size(); globalStepNumber++)
    //int globalStepNumber = 1;  // нагрузили
    //int globalStepNumber = 2;  // +охладили
    int FeInd = 0;  // индекс конечного элемента
    if(task.mechTask.enabled)
    for (i0 = 0; i0 < gp.N[0]; i0++)            // x
        for (i1 = 0; i1 < gp.N[1]; i1++)        // y
            for (i2 = 0; i2 < gp.N[2]; i2++)    // z
            {
                MechOutFePointData &r = (*out.mechOut)[globalStepNumber].fe[FeInd].pd[0];
                VECTOR3 ms;
                r.mainStressesXY(ms);
                double sigma1 = ms[0];
                double sigma2 = ms[1];
                double sigma3 = ms[2];
                // главные напряжения
                /*
                VECTOR3 ms;
                r.sumSigma.solveMain(ms);
                int ind[3];
                sortms(ms, ind);
                double sigma1 = ms[ind[0]];
                double sigma2 = ms[ind[1]];
                double sigma3 = ms[ind[2]];
                */
                /*
                //VECTOR3 ms;
                //r.sumSigma.solveMain(ms);
                //double sigma1 = ms[0];
                //double sigma2 = ms[1];
                //double sigma3 = ms[2];
                */
                if(i1 == 0)
                {
                    fprintf(f_sigmaBottom1, "%lu\t%le\n", i0, sigma1);
                    fprintf(f_sigmaBottom2, "%lu\t%le\n", i0, sigma2);
                    fprintf(f_sigmaBottom3, "%lu\t%le\n", i0, sigma3);
                }
                if(i1 == gp.N[1] - 1)
                {
                    fprintf(f_sigmaTop1, "%lu\t%le\n", i0, sigma1);
                    fprintf(f_sigmaTop2, "%lu\t%le\n", i0, sigma2);
                    fprintf(f_sigmaTop3, "%lu\t%le\n", i0, sigma3);
                }
                if(i1 == gp.N[1]/2)
                {
                    fprintf(f_sigmaMiddle1, "%lu\t%le\n", i0, sigma1);
                    fprintf(f_sigmaMiddle2, "%lu\t%le\n", i0, sigma2);
                    fprintf(f_sigmaMiddle3, "%lu\t%le\n", i0, sigma3);
                }
                FeInd++;
            }
    // температура и реакции опоры
    int vertexInd = 0;  // индекс вершины
    for (i0 = 0; i0 <= gp.N[0]; i0++)            // x
        for (i1 = 0; i1 <= gp.N[1]; i1++)        // y
            for (i2 = 0; i2 <= gp.N[2]; i2++)    // z
            {
                if(task.thermTask.enabled)
                {
                    Heat::ThermOutVertexData &r = (*out.thermOut)[globalStepNumber].vertex[vertexInd];
                    //POINT3 v = task.grid->vertex[vertexInd];
                    //fprintf(f_T, "%le %le\n", v[1], r.T);
                    fprintf(f_T, "%le\t%le\n", (double)i1, r.T);
                }
                bool contacted = (*out.mechOut)[globalStepNumber].vertex[vertexInd].contact;
                if(contacted)
                //if(i1 == 0)
                {
                    POINT3 p0 = (*out.mechOut)[globalStepNumber].vertex[vertexInd].p;
                    VECTOR3 F = (*out.mechOut)[globalStepNumber].vertex[vertexInd].F_sum;
                    double h = (*out.mechOut)[globalStepNumber].vertex[vertexInd].h;
                    //if(sumF != 0)
                    fprintf(f_F, "%lf\t%le\n", p0[0], F.abs());
                    //if(F.abs() != 0)
                        fprintf(f_h, "%lf\t%le\n", p0[0], h);
                    //VECTOR3 sumF = (*out.mechOut)[globalStepNumber].vertex[vertexInd].sumF_vector;
                    //fprintf(f_F, "%le %le\n", (double)i0, sumF.abs());
                }
                vertexInd++;
            }
    fclose(f_F);
    fclose(f_h);
    fclose(f_sigmaTop1);
    fclose(f_sigmaTop2);
    fclose(f_sigmaTop3);
    fclose(f_sigmaMiddle1);
    fclose(f_sigmaMiddle2);
    fclose(f_sigmaMiddle3);
    fclose(f_sigmaBottom1);
    fclose(f_sigmaBottom2);
    fclose(f_sigmaBottom3);
    fclose(f_T);
    // построение графиков
    makeGraphsAfterGlobalStep(globalStepNumber, task);

}
bool Test_formovka::possibleToShow2d() const
{
    return true;
}
bool Test_formovka::getContactSurfaceCircle(const Task &task, const int globalStepIndex, POINT2 &c, double &R) const
{
    if(gp.surfaceType != Grid::SurfaceType::none)
    {
        const auto &s = (*task.mechTask.step)[globalStepIndex];
        (*task.mechTask.rigidSurface)[0]->update(s.t_finish, nullptr);
        R = ((Interpolation::AnaliticalSurface_Cylinder*)(*task.mechTask.rigidSurface)[0])->R;
        c[0] = ((Interpolation::AnaliticalSurface_Cylinder*)(*task.mechTask.rigidSurface)[0])->C[0];
        c[1] = ((Interpolation::AnaliticalSurface_Cylinder*)(*task.mechTask.rigidSurface)[0])->C[1];
        return true;
    }
    else
        return false;
}
void Test_formovka::needToDrawFe(const Solid::Task &, const OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const
{
    if(!(feInd%gp.N[2] == (size_t)z_index ||
            feInd >= gp.N[0] * gp.N[1] * gp.N[2]))
    {
        faceState = {0, 0, 0, 0, 0, 0};
    }
    else
    {
        faceState = {1, 0, 0, 0, 0, 0};
    }
}
bool Test_formovka::needToPaintFiniteElementSurface() const
{
    return gp.surfaceType == Grid::SurfaceType::FiniteElementSurface;
}
}


void solveDiskPoint(const double x0, const double y0, const double r, const double fi, POINT3 &p)
{
    p[0] = x0 + r*cos(fi);
    p[1] = y0 + r*sin(fi);
}

namespace Grid
{
void Grid3D::genFormovka(const FormovkaParameters &gp)
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
        int N0 = gp.N[0];
        int N1 = gp.N[1];
        int N2 = gp.N[2];
        int offset2 = 1;
        int offset1 = N2 + 1;
        int offset0 = (N2 + 1)*(N1 + 1);
        //g->vertex.resize((rpp.N[0] + 1) * (rpp.N[1] + 1) * (rpp.N[2] + 1));		// количество точек
        // построение вершин и первых краевых условий
        for (int i0 = 0; i0 <= N0; i0++)            // x
            for (int i1 = 0; i1 <= N1; i1++)        // y
                for (int i2 = 0; i2 <= N2; i2++)    // z
                {
                    POINT3 p;
                    findPointOnTheLine_1d_conc(gp.p1[0], gp.p2[0], N0, 1.- 1./64*0, i0, p[0]);
                    findPointOnTheLine_1d(gp.p1[1], gp.p2[1], N1, 1, i1, p[1]);
                    findPointOnTheLine_1d(gp.p1[2], gp.p2[2], N2, 1, i2, p[2]);
                    // ind = i0*(i1+1)*(i2+1)+i1*(i2+1)+i0 - индекс точки i0,i1,i2
                    vertex.push_back(p);
                    // фиксируем по x посередине
                    if(i0 == N0/2)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 0;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    // фиксируем по y на центральном участке снизу
                    /*
                    int i0_seredina = N0/2;
                    if(i1 == 0 && i0 >= i0_seredina - gp.N_react && i0 <= i0_seredina + gp.N_react)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 1;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }*/
                    // фиксируем по z посередине
                    if(i2 == N2/2)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 2;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    // фиксируем по z везде
                    if(false)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 2;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    // температура снизу
                    if(i1 == 0)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 3;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    // температура сверху
                    if(i1 == N1)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 4;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                }

        // построение конечных элементов и поверхностей
        // будет 4 поверхности (2 длля 2-х краевых и 2 для контакта)
        FESurface.resize(4);
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
                    if (i1 == N1 - 1 && (i0 < gp.N_p || i0 >= N0 - gp.N_p))
                    {
                        face.feIndex = (int)fe.size() - 1;
                        face.faceIndex = 5; // 2    Y = +1
                        FESurface[0].face.push_back(face);
                    }
                    if (i1 == N1 - 1)
                    {
                        face.feIndex = (int)fe.size() - 1;
                        face.faceIndex = 5; // 2    Y = +1
                        FESurface[1].face.push_back(face);
                    }
                    // поверхность для контакта (активная)
                    if(i1 == 0)
                    {
                        Grid::FEFace face;
                        face.feIndex = (int)fe.size() - 1;
                        face.faceIndex = 4; //-2    Y = -1
                        FESurface[2].face.push_back(face);
                    }
                }
    }



    // параметры индексации пространства
    {
        regionIndexationParameters.h_min = 2;
        regionIndexationParameters.h_max = 6;
        Sphere s0;
        s0.O = {0, -1, 0};
        s0.R = 1;
        //s0.O = {0, -6, 0};
        //s0.R = 10;
        regionIndexationParameters.q0.initBySphere(s0);
        regionIndexationParameters.div = 1;
    }

    // поверхность заданная сеткой
    // дополнительная неподвижная сетка для задания жёсткой контактной поверхности
    // (если contact_mode != 3 то эта дополнительная сетка не строится)
    if(gp.surfaceType == SurfaceType::FiniteElementSurface)
    {
        int vind_start0 =  (int)vertex.size();
        int N0 = gp.N[0];
        int N1 = 1;
        int N2 = 1;
        int offset2 = 1;
        int offset1 = N2 + 1;
        int offset0 = (N2 + 1)*(N1 + 1);
        // построение точек
        for (int i0 = 0; i0 <= N0; i0++)               // x (fi)
            for (int i1 = 0; i1 <= N1; i1++)           // y (r)
                for (int i2 = 0; i2 <= N2; i2++)       // z
                {
                    double r1 = 4.8;
                    double r2 = 5;
                    double fi1 = PI + PI/4;
                    double fi2 = -PI/4;
                    double r = r1 + (r2 - r1)*((double)i1)/N1;         // r1..r2
                    double fi = fi1 + (fi2 - fi1)*((double)i0)/N0;     // fi1..fi2
                    double x0 = 0, y0 = -1 - r2;        // центр круга

                    POINT3 p;
                    solveDiskPoint(x0, y0, r, fi, p);                              // x, y
                    double z0 = (gp.p1[2] + gp.p2[2]) / 2;
                    double z1 = z0 + (gp.p1[2] - z0)*2 / 1;
                    double z2 = z0 + (gp.p2[2] - z0)*2 / 1;
                    findPointOnTheLine_1d(z1, z2, N2, 1, i2, p[2]);    // z ось направлена на нас!
                    //findPointOnTheLine_1d(gp.p1[2], gp.p2[2], N2, 1, i2, p[2]);
                    //findPointOnTheLine_1d(gp.p1[2], gp.p2[2], N2, 1, i2, p[2]);    // z
                    vertex.push_back(p);
                    BoundaryCondition1 bc1_el;
                    // фиксируем по x
                    bc1_el.bc1SourceIndex = 0;
                    bc1_el.vertexIndex = (int)vertex.size() - 1;
                    bc1.push_back(bc1_el);
                    // фиксируем по y
                    bc1_el.bc1SourceIndex = 1;
                    bc1_el.vertexIndex = (int)vertex.size() - 1;
                    bc1.push_back(bc1_el);
                    // фиксируем по z
                    bc1_el.bc1SourceIndex = 2;
                    bc1_el.vertexIndex = (int)vertex.size() - 1;
                    bc1.push_back(bc1_el);
                }

        // построение конечных элементов и поверхностей
        for (int i0 = 0; i0 < N0; i0++)                // x (fi)
            for (int i1 = 0; i1 < N1; i1++)            // y (r)
                for (int i2 = 0; i2 < N2; i2++)        // z
                {
                    // конечный элемент
                    int vind0 = vind_start0 + i0*offset0 + i1*offset1 + i2*offset2;       // индекс вершины, определяющей шестигранник
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
                    // поверхность для контакта (пассивная, жёсткая)
                    Grid::FEFace face;
                    face.feIndex = (int)fe.size() - 1;
                    face.faceIndex = 5; //2 // Y = 1
                    FESurface[3].face.push_back(face);
                }
    }

    GridRectangleRegular2D it_grid;
    int N[2] = {80, 40};
    if(gp.surfaceType == SurfaceType::InterpolantSurface_Lagrange3)    // интерполянт лагранжа
    {
        //N[1] = 1024;
        if(gp.decompositionType == RegionDecompositionType::cubeVertexIndexation)
            //it_grid.init(0 - PI/2, PI + PI/2, gp.p1[2] - 20, gp.p2[2] + 20, N[0], N[1]);
            it_grid.init(0 - PI/2, PI + PI/2, gp.p1[2] - 20, gp.p2[2] + 20, N[0], N[1]);
        else
            //it_grid.init(0 - PI/2, PI + PI/2, gp.p1[2] - 2, gp.p2[2] + 2, N[0], N[1]);
            it_grid.init(0 - PI/2, PI + PI/2, gp.p1[2] - 20, gp.p2[2] + 20, N[0], N[1]);

    }
    if(gp.surfaceType == SurfaceType::InterpolantSurface_Hermite3)    // интерполянт эрмита
    {
        if(gp.decompositionType == RegionDecompositionType::cubeVertexIndexation)
            it_grid.init(0 - PI/2, PI + PI/2, gp.p1[2] - 20, gp.p2[2] + 20, N[0], N[1]);
        else
            it_grid.init(0 - PI/2, PI + PI/2, gp.p1[2] - 2, gp.p2[2] + 2, N[0], N[1]);
    }
    if(gp.surfaceType == SurfaceType::AnaliticalSurface_Cylinder)    // аналитическая поверхность
    {
        it_grid.init(0 - PI/2, PI + PI/2, gp.p1[2] - 2, gp.p2[2] + 2, N[0], N[1]);
    }

    // поверхность заданная аналитически
    analiticalSurface.resize(1);
    VECTOR3 V(0, -0.001*0, 0);
    if(gp.surfaceType == SurfaceType::AnaliticalSurface_Cylinder)
    {
        Interpolation::AnaliticalSurface_Cylinder *cylinder =
                new Interpolation::AnaliticalSurface_Cylinder(
                    gp.decompositionType,
                    it_grid,
                    gp.R_fun,
                    gp.x_fun,
                    gp.y_fun,
                    0, //0
                    0, //0
                    0, //0
                    0, //0
                    1.e-8);
        analiticalSurface[0] = cylinder;
    }
    else
    {
        // аналитическую поверхность инициализируем в любом случае (без декомпозиции)
        Interpolation::AnaliticalSurface_Cylinder *cylinder =
                new Interpolation::AnaliticalSurface_Cylinder(
                    RegionDecompositionType::none,
                    it_grid,
                    gp.R_fun,
                    gp.x_fun,
                    gp.y_fun,
                    0, //0
                    0, //0
                    0, //0
                    0, //0
                    1.e-8);
        analiticalSurface[0] = cylinder;
    }
    // поверхность заданная интерполянтом
    ISurface.resize(1);
    if(gp.surfaceType == SurfaceType::InterpolantSurface_Hermite3 ||
       gp.surfaceType == SurfaceType::InterpolantSurface_Lagrange3)    // интерполянт
    {
        //int NN[2] = {N[0]*3, N[1]*3};   // без регуляризации
        //int NN[2] = {N[0], N[1]};     // с регуляризацией
        //double CoefAlpha = 0.001;
        //it_grid.init(0, PI, 0, 5., N[0], N[1]);
        //it_grid.init(0, PI, gp.p1[2] - 11, gp.p2[2] + 11, N[0], N[1]);
        //it_grid.init(0, PI, -20, 20, N[0], N[1]);
        //it_grid.init(-PI/4, PI/4, gp.p1[2] - 1, gp.p2[2] + 1, N[0], N[1]);
        //it_grid.init(0 - PI/2, PI + PI/2, -20, 20, N[0], N[1]);
        //it_grid.init(0 - PI/2, PI + PI/2, gp.p1[2] - 1, gp.p2[2] + 1, N[0], N[1]);
        //it_grid.init(0 - PI/2, PI + PI/2, gp.p1[2] - 20, gp.p2[2] + 20, N[0], N[1]);
        double CoefAlpha;
        Interpolation::Interpolant2D_base *itx;
        Interpolation::Interpolant2D_base *ity;
        Interpolation::Interpolant2D_base *itz;
        if(gp.surfaceType == SurfaceType::InterpolantSurface_Hermite3)    // интерполянт эрмита
        {
            CoefAlpha = 0.0001;
            itx = new Interpolation::Interpolant2D_Hermite3;
            ity = new Interpolation::Interpolant2D_Hermite3;
            itz = new Interpolation::Interpolant2D_Hermite3;
        }
        if(gp.surfaceType == SurfaceType::InterpolantSurface_Lagrange3)    // интерполянт лагранжа
        {
            CoefAlpha = 0.0000;
            itx = new Interpolation::Interpolant2D_Lagrange3;
            ity = new Interpolation::Interpolant2D_Lagrange3;
            itz = new Interpolation::Interpolant2D_Lagrange3;
        }

        itx->init(it_grid, CoefAlpha);
        ity->init(it_grid, CoefAlpha);
        itz->init(it_grid, CoefAlpha);

        std::vector<POINT2> coordinates;
        itx->getNodesCoordinates(coordinates);
        Vector nodeValuex, nodeValuey, nodeValuez;

        // определение ф-и в наборе точек
        nodeValuex.resize(coordinates.size());
        nodeValuey.resize(coordinates.size());
        nodeValuez.resize(coordinates.size());
        for(size_t i = 0; i < coordinates.size(); i++)
        {
            double alpha = coordinates[i][0];
            double zz = coordinates[i][1];
            nodeValuex[i] = gp.R0*cos(alpha);
            nodeValuey[i] = gp.R0*sin(alpha) + gp.y0;
            nodeValuez[i] = zz;
            //itx.addPoint(coordinates[i], nodeValuex[i]);
            //ity.addPoint(coordinates[i], nodeValuey[i]);
            //itz.addPoint(coordinates[i], nodeValuez[i]);
        }

        itx->buildInterpolantByAllNodes(nodeValuex);
        ity->buildInterpolantByAllNodes(nodeValuey);
        itz->buildInterpolantByAllNodes(nodeValuez);

        // интерполянт эрмита
        if(gp.surfaceType == SurfaceType::InterpolantSurface_Hermite3)
        {
            Interpolation::InterpolantSurface_Hermite3 *it =
                    new Interpolation::InterpolantSurface_Hermite3(
                        gp.decompositionType,
                        (Interpolation::Interpolant2D_Hermite3 *)itx,
                        (Interpolation::Interpolant2D_Hermite3 *)ity,
                        (Interpolation::Interpolant2D_Hermite3 *)itz,
                        20,
                        1.e-7,
                        1.e-12,
                        1.e-20,
                        1.e-8);
            ISurface[0] = it;
        }
        /*200,
        1.e-10,
        1.e-13,
        1.e-20,
        1.e-10);*/
        // интерполянт лагранжа
        if(gp.surfaceType == SurfaceType::InterpolantSurface_Lagrange3)
        {
            Interpolation::InterpolantSurface_Lagrange3 *it =
                    new Interpolation::InterpolantSurface_Lagrange3(
                        gp.decompositionType,
                        (Interpolation::Interpolant2D_Lagrange3 *)itx,
                        (Interpolation::Interpolant2D_Lagrange3 *)ity,
                        (Interpolation::Interpolant2D_Lagrange3 *)itz,
                        20,//100,
                        1.e-7,//1.e-12,
                        1.e-12,
                        1.e-20,
                        1.e-8);//1.e-12
            ISurface[0] = it;
        }






        FILE *fx = fopen("fx", "w");
        FILE *fy = fopen("fy", "w");
        FILE *fz = fopen("fz", "w");
        FILE *fpr = fopen("fpr", "w");

        int color = 3;

        itx->save("x.bmp", 800, 800, color, 1);
        ity->save("y.bmp", 800, 800, color, 1);
        itz->save("z.bmp", 800, 800, color, 1);

        for(size_t i = 0; i < coordinates.size(); i++)
        {
            double alpha = coordinates[i][0];
            double zz = coordinates[i][1];
            fprintf(fx, "alpha = %lf, zz = %lf, Fx = %lf, ip = %lf\n", alpha, zz, nodeValuex[i], itx->fun({alpha, zz}));
            fprintf(fy, "alpha = %lf, zz = %lf, Fy = %lf, ip = %lf\n", alpha, zz, nodeValuey[i], ity->fun({alpha, zz}));
            fprintf(fz, "alpha = %lf, zz = %lf, Fz = %lf, ip = %lf\n", alpha, zz, nodeValuez[i], itz->fun({alpha, zz}));
        }


        // проекции на поверхность
        for(int i = -1; i <= 101; i++)
        {
            double alpha = PI*i/100;
            double zz = 0.1;  // -2 .. 0.5
            POINT3 p0;
            //p0 = {0, R + 1, 0};
            p0[0] = (gp.R0 + 1)*cos(alpha);
            p0[1] = (gp.R0 + 1)*sin(alpha) - gp.R0 - 1;
            p0[2] = zz;
            POINT3 nearestPointAn;
            nearestPointAn[0] = (gp.R0)*cos(alpha);
            nearestPointAn[1] = (gp.R0)*sin(alpha) + gp.y0;
            nearestPointAn[2] = zz;
            POINT3 normalAn = (p0 - nearestPointAn) / (p0 - nearestPointAn).abs();

            POINT3 nearestPoint;
            POINT3 normal;
            int side;
            bool onBorder;
            SurfacePositionData ssd1, ssd2;
            ssd1.index = -1;    // данные о начальном приближении отсутствуют
            ISurface[0]->findNearestPoint(p0, {0,0,0}, ssd1,
                                          ssd2, nearestPoint, normal, side, onBorder);
            fprintf(fpr, "p0 = (%lf, %lf, %lf):\n\tnearestPointAn = (%lf, %lf, %lf)\n\tnearestPoint   = (%lf, %lf, %lf)\n\tnormalAn = (%lf, %lf, %lf)\n\tnormal   = (%lf, %lf, %lf)\n\tside = %d, onBorder = %d\n",
                    p0[0], p0[1], p0[2],
                    nearestPointAn[0], nearestPointAn[1], nearestPointAn[2],
                    nearestPoint[0], nearestPoint[1], nearestPoint[2],
                    normalAn[0], normalAn[1], normalAn[2],
                    normal[0], normal[1], normal[2],
                    side, (int)onBorder);
        }

        fclose(fx);
        fclose(fy);
        fclose(fz);
        fclose(fpr);
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
