#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"
#include <fstream>

#include "tests.h"
#include "post.h"
#include "console.h"

#include "interpolation.h"

using namespace Elementary;
using namespace Solid;
using namespace Grid;

namespace Tests
{
Test_contactSphere::Test_contactSphere()
{
    // дирректория данных теста
    dir = "./tests/contactSphere/";
    // имена gnuplot файлов(+пути к файлам в опциях) и соответствующих картинок
    fn_filesList = dir + "files.txt";
    resGraph.read_gnuplot_from_file(fn_filesList);
    // инициализация графиков с невязками
    initStepsGraphs();
}
Test_base::Type Test_contactSphere::get_type() const
{
    return Type::ContactSphere;
}
// расстояние от центра поверхности контакта(круга) до точки с этой поверхности
// считаем просто расстояние а не длину дуги
double Test_contactSphere::calc_r(const POINT3 &p)
{
    return sqrt(SQR(p[0]) + SQR(p[1]));
}
double Test_contactSphere::calc_P(const double d, const Grid::ContactSphereParameters &csp)
{
    double E = csp.E / (1 - SQR(csp.Nu));
    double R = csp.R0;
    double P;
    double P1 = 0, P2 = 1.e11;
    // функция d(P) растёт
    for(int i = 0; i < 100; i++)
    {
        P = (P1 + P2) / 2;
        double d_trial = pow(3./(4.*E)*P, 2./3.) / pow(R, 1./3.);
        if(d_trial < d)
        {
            P1 = P;
        }
        else
        {
            P2 = P;
        }
    }
    return P;
}
void Test_contactSphere::initTask(Solid::Task &task)
{
    // подрубка
    task.mechTask.enabled = true;
    task.thermTask.enabled = false;
    // общие шаги и сетка
    std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
    Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
    // способ представления поверхности
    ccp.surfaceType = SurfaceType::AnaliticalSurface_Sphere;

    // способ декомпозиции поверхности
    ccp.decompositionType = RegionDecompositionType::none;

    // подрубка универсальной оптимизации
    bool noContactRadiusOptimization = false;
    //bool noContactRadiusOptimization = true;

    ccp.E = 1.e10;
    ccp.Nu = 0.3;

// при N = 12, n = 0
    // было: 1922.7
// подрубил -march=native и стало: 1768.2
// подрубил -O3 -march=native -mtune=native и стало: 1512.9
// подрубил -O3 -march=native -mtune=native -funroll-loops и стало: 1501.4
// память 3200 - стало 1470.3
// обновление линукса - стало 1510.0

// при N = 16, n = 0
// 7230.5 с = 2 часа ровно, жрёт 2гб, то есть за 2 часа можно решить 6 задач, с тратой 12гб памяти?
    ccp.NN = 8;      // подробность сетки
    ccp.bYeldCoeff = 550; //I) 110; 220; *550*; 1100;
    ccp.bd_max = 110;     //II) 1; 3; 6; 10; 15; 20; 30; 50; 70; 90; *110*;
    ccp.n = 0.0;      // степенное упрочнение III) *0*; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8;
    ccp.R0 = 1;
    ccp.grid_mnojitel1 = 1.5;
    ccp.grid_mnojitel2 = 15;
    ccp.w_project = 1;

    ccp.w_stiffness = 1;
    ccp.sigmaResidualLimit = 1.e-10;          // желаемая относительная погрешность по эквиволентным напряжениям
    ccp.epsResidualLimit = 1.e-10;            // желаемая относительная погрешность по эквиволентным деформациям
    ccp.contactDeltaFResidualLimit = 1.e-10;;//1.e-12;//1.e-7;  // желаемая относительная погрешность по силе реакции опоры
    ccp.contactEndPointResidualLimit = 1.e-14;//1.e-15

    // ввод параметров задачи из файла
    fn_general_inf = dir + "_input.in";
    FILE *f_parametors = fopen(fn_general_inf.c_str(), "r");
    fscanf(f_parametors, "%lf", &ccp.n);
    fscanf(f_parametors, "%lf", &ccp.bd_max);
    fscanf(f_parametors, "%lf", &ccp.bYeldCoeff);

    fscanf(f_parametors, "%d", &ccp.NN);
    fscanf(f_parametors, "%d", &ccp.STEPS_PUSH);
    fscanf(f_parametors, "%lf", &ccp.grid_mnojitel1);
    fscanf(f_parametors, "%lf", &ccp.grid_mnojitel2);


    fscanf(f_parametors, "%lf", &ccp.R0);
    fscanf(f_parametors, "%lf", &ccp.E);
    fscanf(f_parametors, "%lf", &ccp.Nu);


    fscanf(f_parametors, "%d", &ccp.plasticityMethodType); // 4 - sigma0, 5 - C_pl + sigma0, 1 - C_pl, 2 - C_pl алгоритм Соловейчика
    fscanf(f_parametors, "%d", &ccp.incForsesMode);        // 1 - безумство, 2 - без безумства
    fscanf(f_parametors, "%d", &ccp.constantNormal);


    fscanf(f_parametors, "%lf", &ccp.w_project);
    fscanf(f_parametors, "%lf", &ccp.w_midPoint);
    fscanf(f_parametors, "%lf", &ccp.w_stiffness);
    fscanf(f_parametors, "%lf", &ccp.cosTettaMin);

    fscanf(f_parametors, "%lf", &ccp.sigmaResidualLimit);
    fscanf(f_parametors, "%lf", &ccp.epsResidualLimit);
    fscanf(f_parametors, "%lf", &ccp.contactDeltaFResidualLimit);
    fscanf(f_parametors, "%lf", &ccp.contactEndPointResidualLimit);
    fclose(f_parametors);


    int &NN = ccp.NN;
    double &c1 = ccp.c1;
    double &c2 = ccp.c2;
    double &c3 = ccp.c3;
    double &q = ccp.q;

    c1 = 0.8;// 1.5; // 2;
    c2 = 10; // 20
    c3 = 20; // 20
    q = pow(1./40., 1./NN); //1 - 1./4.; // сгущение
    // рассчёт предела упругости, глубины, количества шагов, сгущения
    {
        double E = ccp.E / (1 - SQR(ccp.Nu));   // приведённый модуль упругости
        ccp.elasticSigmaLimit = E/ccp.bYeldCoeff;
        double c = 1.08;

        //double d_max_const = 0.05;
        //ccp.R0 = d_max_const/ccp.bd_max*16./9./SQR(PI)/SQR(c)*SQR(E/ccp.elasticSigmaLimit);

        ccp.P_y = 9./16.*pow(PI, 3)*pow(c, 3)*SQR(ccp.R0)*ccp.elasticSigmaLimit/SQR(E/ccp.elasticSigmaLimit);
        ccp.A_y = 9./16.*pow(PI, 3)*pow(c, 2)*SQR(ccp.R0)/SQR(E/ccp.elasticSigmaLimit);
        ccp.d_y = 9./16.*pow(PI, 2)*pow(c, 2)*ccp.R0/SQR(E/ccp.elasticSigmaLimit);
        ccp.d0 = ccp.bd_max*ccp.d_y;
        //ccp.STEPS_PUSH = MAX((int)round(ccp.bd_max*1.0), 110); // количество шагов вдавливания и вынимания
        // сгущение сетки в зависимости от глубины
        double L_contact = sqrt(2*ccp.R0*ccp.d0 - SQR(ccp.d0));
        c1 = ccp.grid_mnojitel1*L_contact;
        //c1 = 3*sqrt(2*ccp.R0*ccp.d0 - SQR(ccp.d0));
        c2 = ccp.grid_mnojitel2*L_contact;//c1*10;     //c1*20;     //c1*10
        c3 = c2*2;//c1*10*2;   //c1*20*2;   //c1*10*2

            //c2 = c1*20;     //c1*20;     //c1*10
            //c3 = c1*20*2;   //c1*20*2;   //c1*10*2
        // коэффициент разрядки
        {
            double L0 = c1/NN;  // желаемый размер самого разряжённого КЭ
            double l = c2 - c1;
            double y = (l - L0) / l;
            // решение уравнения f(q) = y
            // функция f(q) убывает
            double q1 = 0;
            double q2 = 0.9999;
            for(int i = 0; i < 100; i++)
            {
                q = (q1 + q2) / 2;
                double f_trial = (1. - pow(q, NN - 1.)) / (1. - pow(q, NN));
                if(f_trial > y)
                {
                    q1 = q;
                }
                else
                {
                    q2 = q;
                }
            }
            //q = pow(1./80., 1./NN); // сгущение // c1 = 0.8 - q = pow(1./40., 1./NN);
        }
    }

    ccp.sortNodes = true;   // переупорядочение узлов сетки
    int unload = 1;   // 1 - есть разгрузка, 0 - нет разгрузки
    double k_unload_step = 0.5;//1.e-10;//0.00001;//0.5;
    ccp.CORR = 0;        // количество шагов коррекции после каждого шага (1 или 0)



    /*
    ccp.elasticSigmaLimit = 1.e100;
    ccp.STEPS_PUSH = 1;
    unload = 0;
    */




    /*
    int NN = 6;      // подробность сетки
    double c1 = 0.8;// 1.5; // 2;
    double q = pow(1./40., 1./NN); //1 - 1./4.; // сгущение
    ccp.sortNodes = true;   // переупорядочение узлов сетки

    ccp.d = 0.01;//+1.0;//0.5//0.1
    ccp.STEPS_PUSH = 30; // количество шагов вдавливания и вынимания
    int unload = 1;   // 1 - есть разгрузка, 0 - нет разгрузки
    double k_unload_step = 0.5;//1.e-10;//0.00001;//0.5;
    ccp.CORR = 0;        // количество шагов коррекции после каждого шага (1 или 0)

    ccp.E = 1.e10;
    ccp.Nu = 0.3;
    ccp.elasticSigmaLimit = 1.e7;//1.e8;//3.e8;

    ccp.d = 0.05;//+1.0;//0.5//0.1
    ccp.STEPS_PUSH = 64; // количество шагов вдавливания и вынимания
    int unload = 1;   // 1 - есть разгрузка, 0 - нет разгрузки
    double k_unload_step = 0.5;//1.e-10;//0.00001;//0.5;
    ccp.elasticSigmaLimit = 5.e7;//1.e8;//3.e8;
    */

    //упругость
    int plasticityMethodType = ccp.plasticityMethodType;
    //plasticityMethodType = (int)MechPlasticityMethodType::Elasticity;
    // пластичность
    //plasticityMethodType = (int)MechPlasticityMethodType::D_pl;
    //plasticityMethodType = (int)MechPlasticityMethodType::InitialSigma;
     //plasticityMethodType = (int)MechPlasticityMethodType::Combo_D_pl_InitialSigma;
    double w_project = ccp.w_project;//0.1;//0.5;//1;
    double w_midPoint = ccp.w_midPoint;

    MechPlasticityCurveUnloadingType PCUnloadingType = MechPlasticityCurveUnloadingType::Elastic;
    //MechPlasticityCurveUnloadingType PCUnloadingType = MechPlasticityCurveUnloadingType::Curve;
    double k_sigma_eps = 1.e50;//1000.*10000;//1000
    double k_eps_sigma = 10;
    ccp.plasticOut = (plasticityMethodType != (int)MechPlasticityMethodType::Elasticity);

    Solid::ContactType contactType = Solid::ContactType::AugmentedLagrange;
    //Solid::ContactType contactType = Solid::ContactTy::Fastest;
    //ccp.constantNormal = 1
//    ccp.stiffnessMultiplier = 1;
    int plasticPlot = 1;
//    ccp.sigmaResidualLimit = 1.e-10;            // желаемая относительная погрешность по эквиволентным напряжениям
//    ccp.epsResidualLimit = 1.e-10;            // желаемая относительная погрешность по эквиволентным деформациям


    IncForsesMode incForsesMode = (IncForsesMode)ccp.incForsesMode;


    int stepsNumber = ccp.STEPS_PUSH*(1 + ccp.CORR)*(1);//128;//64*5;//40/10;                        // количество шагов
    int nonlinearIterLimit = 100;                // ограничение на количество итераций
//    ccp.contactDeltaFResidualLimit = 1.e-10;;//1.e-12;//1.e-7;  // желаемая относительная погрешность по силе реакции опоры
//    ccp.contactEndPointResidualLimit = 1.e-14;//1.e-15



    //plasticityMethodType = MechPlasticityMethodType::InitialSigma;

    nonlinearIterLimit = 100;








    // для упругого и упругопластического тестов на одном графике
    //ccp.elasticSigmaLimit = 1.e20;
    //ccp.plasticOut = false;










    double Time = stepsNumber;  // 50 - только нагружение, 100 - нагружение и разгрузка

    // параметры сетки
{
    std::vector<POINT3> &av = ccp.av;
    std::vector<HexagonArea> &ah = ccp.ah;
    // вершины
    av.resize(19);
    // z = -20 (-c3)
    av[0] = POINT3{-c2, -c2, -c3};
    av[1] = POINT3{  0, -c2, -c3};
    av[2] = POINT3{-c2,   0, -c3};
    av[3] = POINT3{  0,   0, -c3};
    // z = -10
    av[4] = POINT3{-c2, -c2, -c2};
    av[5] = POINT3{  0, -c2, -c2};
    av[6] = POINT3{-c2,   0, -c2};
    av[7] = POINT3{  0,   0, -c2};
    // z = -1 (-c1)
    av[8] = POINT3{-c1, -c1, -c1};
    av[9] = POINT3{ 0, -c1,  -c1};
    av[10] = POINT3{-c1,  0, -c1};
    av[11] = POINT3{ 0,  0,  -c1};
    // z = 0
    av[12] = POINT3{-c2, -c2,  0};
    av[13] = POINT3{  0, -c2,  0};
     av[14] = POINT3{-c1, -c1, 0};
     av[15] = POINT3{ 0, -c1,  0};
    av[16] = POINT3{-c2,   0,  0};
     av[17] = POINT3{-c1,  0,  0};
    av[18] = POINT3{  0,   0,  0};
    // 6-гранники
    ah.resize(5);
    ah[0].vi = {0, 1, 2, 3, 4, 5, 6, 7};
    ah[1].vi = {4, 5, 6, 7, 8, 9, 10, 11};
    ah[2].vi = {4, 5, 8, 9, 12, 13, 14, 15};
    ah[3].vi = {4, 8, 6, 10, 12, 14, 16, 17};
    ah[4].vi = {8, 9, 10, 11, 14, 15, 17, 18};
    for(int i = 0; i < 5; i++)
    {
        // разбиения
        ah[i].N = {NN, NN, NN};
        // сгущение
        ah[i].condensation_coord = -1; // нет сгущения
        //материал
        ah[i].mi = 0;
        // поверхности(пока нету)
        ah[i].surfaceIndex = {-1, -1, -1, -1, -1, -1};
        // 1-е краевые(пока нету)
        ah[i].bc1Index = {-1, -1, -1, -1, -1, -1};
    }
    // разбиения
    ah[0].N = {NN, NN, NN/2};
    // сгущения
    //k = 1;
    ah[1].condensation_q = q;
    ah[1].condensation_coord = 2;//2
    ah[2].condensation_q = q;
    ah[2].condensation_coord = 1;//1
    ah[3].condensation_q = q;
    ah[3].condensation_coord = 0;//0
    // поверхность контакта
    ah[4].surfaceIndex = {-1, 0, -1, -1, -1, -1};   // Z = +1
    // 1-е краевые
    // фиксация по X (X = +1)
    ah[0].bc1Index = {-1, -1, -1, 0, -1, -1};
    ah[1].bc1Index = {-1, -1, -1, 0, -1, -1};
    ah[2].bc1Index = {-1, -1, -1, 0, -1, -1};
    ah[4].bc1Index = {-1, -1, -1, 0, -1, -1};
    // фиксация по Y (Y = +1)
    ah[0].bc1Index = {-1, -1, -1, 0, -1, 1};// дополнение предыдущего
    ah[1].bc1Index = {-1, -1, -1, 0, -1, 1};// дополнение предыдущего
    ah[3].bc1Index = {-1, -1, -1, -1, -1, 1};
    ah[4].bc1Index = {-1, -1, -1, 0, -1, 1};// дополнение предыдущего
    // фиксация по Z (Z = -1)
    ah[0].bc1Index = {2, -1, -1, 0, -1, 1};// дополнение предыдущего
}
    // параметры подвижного шара
{
    ccp.z0 = ccp.R0;

    ccp.R_fun = new Interpolation::Interpolant1D_Lagrange1Unregular;
    ccp.x_fun = new Interpolation::Interpolant1D_Lagrange1Unregular;
    ccp.y_fun = new Interpolation::Interpolant1D_Lagrange1Unregular;
    ccp.z_fun = new Interpolation::Interpolant1D_Lagrange1Unregular;

    Grid::GridRectangleRegular1D gr;
    ccp.R_fun->init(gr, 0);// аргументы игнорируются для нерегулярной сетки
    ccp.x_fun->init(gr, 0);
    ccp.y_fun->init(gr, 0);
    ccp.z_fun->init(gr, 0);
    ccp.R_fun->addPoint(0, ccp.R0);
    ccp.R_fun->addPoint(1000, ccp.R0);
    ccp.R_fun->buildInterpolant();
    ccp.x_fun->addPoint(0, 0);
    ccp.x_fun->addPoint(1000, 0);
    ccp.x_fun->buildInterpolant();
    ccp.y_fun->addPoint(0, 0);
    ccp.y_fun->addPoint(1000, 0);
    ccp.y_fun->buildInterpolant();
}

    grid->genContactSphere(ccp);

    //grid->buldOpenSCADModel("setka.scad");
    // параметры
    int fixGrid = 0;                // 0 - подвижная сетка, 1 - зафиксировать сетку
    //int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny

    // шаги
    GlobalStep s;       // шаг
    // вдавливание

    //int firstTimes = 0;
    /*
    int firstTimes = 1;
    double firstStepTime = Time;
    double firstStep_d = ccp.d / stepsNumber / 10;
    s.t_start = 0;
    s.t_finish = firstStepTime;
    s.dt0 = firstStepTime;
    step->push_back(s);
    ccp.z_fun->addPoint(s.t_start, ccp.z0 - firstStep_d);
    */
    ccp.z_fun->addPoint(0, ccp.z0);
    for(int i = 0; i < stepsNumber; i++)
    {
        s.t_start = Time*i/stepsNumber;
        s.t_finish = Time*(i + 1)/stepsNumber;
        s.dt0 = Time/stepsNumber;
        ccp.z_fun->addPoint(s.t_finish, ccp.z0 - s.t_finish/Time*ccp.d0);
        /*
        ccp.z_fun->addPoint(s.t_start, ccp.z0 - firstStep_d - (s.t_start - firstStepTime)/Time*(ccp.d - firstStep_d));
        if(i == stepsNumber - 1)
            ccp.z_fun->addPoint(s.t_finish, ccp.z0 - firstStep_d - (s.t_finish - firstStepTime)/Time*(ccp.d - firstStep_d));
        */
        step->push_back(s);
    }
    if(unload == 1)
    {
        double firstStepTime = Time;
        double firstStep_d = ccp.d0;
        for(int i = 0; i <= stepsNumber; i++)
        {
            if(i == 0)
            {
                s.t_start = firstStepTime;
                s.t_finish = firstStepTime + Time/stepsNumber;
                s.dt0 = Time/stepsNumber;
                ccp.z_fun->addPoint(s.t_finish, ccp.z0 - firstStep_d + ccp.d0/stepsNumber*k_unload_step);
                step->push_back(s);
            }
            else
            {
                s.t_start = firstStepTime + Time*i/stepsNumber;
                s.t_finish = firstStepTime + Time*(i + 1)/stepsNumber;
                s.dt0 = Time/stepsNumber;
                if(i == stepsNumber)
                    ccp.z_fun->addPoint(s.t_finish, ccp.z0 - firstStep_d + (s.t_start - firstStepTime)/Time*(ccp.d0) + ccp.d0);
                else
                    ccp.z_fun->addPoint(s.t_finish, ccp.z0 - firstStep_d + (s.t_start - firstStepTime)/Time*(ccp.d0));
                /*
                ccp.z_fun->addPoint(s.t_start, ccp.z0 - firstStep_d + (s.t_start - firstStepTime)/Time*(ccp.d));
                if(i == stepsNumber)
                    ccp.z_fun->addPoint(s.t_finish, ccp.z0 - firstStep_d + (s.t_finish - firstStepTime)/Time*(ccp.d) + ccp.d);
                */
                step->push_back(s);
            }
        }

    }


    ccp.z_fun->buildInterpolant();

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
        mechMat.elasticSigmaLimit = ccp.elasticSigmaLimit;
        mechMat.Ceps = 4./9.;
        mechMat.Csigma = 1;
        mechMat.set_E_NU(ccp.E, ccp.Nu);
        mechMat.set_M_sigma();
        mechMat.F = {0,0,0};// объёмные силы
        mechMat.elasticParameters0.ro = 0;     // плотность
            mechMat.elasticParameters0.Talpha = 0*1.e-5; // Коэффициент линейного расширения
        if(ccp.n == 0)
            setPlasticMaterialCurve_Yeld(plasticPlot, mechMat, k_sigma_eps, k_eps_sigma);
        else
            setPlasticMaterialCurve_Yeld_hardening(plasticPlot, mechMat, ccp.n);
        mechMat.plasticityMethodType = (MechPlasticityMethodType)plasticityMethodType;
        mechMat.PCUnloadingType = PCUnloadingType;
        mechMat.w_midPoint = w_midPoint;//1.0;//0.5//0
        mechMat.w_project = w_project;
        mechMat.cosTettaMin = ccp.cosTettaMin;

        // первые краевые условия
        bc1Source = new std::vector<MechBoundaryCondition1Source>(4);
        (*bc1Source)[0].mode = {{ 0, -1, -1}};
        (*bc1Source)[0].u0 =   {{ 0, -1, -1}};
        (*bc1Source)[1].mode = {{-1,  0, -1}};
        (*bc1Source)[1].u0 =   {{-1,  0, -1}};
        (*bc1Source)[2].mode = {{-1, -1,  0}};
        (*bc1Source)[2].u0 =   {{-1, -1,  0}};
        (*bc1Source)[3].mode = {{-1, -1,  -1}};
        (*bc1Source)[3].u0 =   {{-1, -1,  -1}};
        // вторые краевые
        bc2 = new std::vector<MechBoundaryCondition2>;
        bc2->clear();
        // шаги (на каждом шаге задаются вторые краевые условия)
        std::vector<MechGlobalStepParameters> *mechStep = new std::vector<MechGlobalStepParameters>;
        MechGlobalStepParameters ms;  // шаг для МДТТ
        ms.slausolverParameters = slausolver_parameters;
        ms.timeMode = 0;    // 0 - квазистатическая задача
        ms.incForsesMode = incForsesMode; // 1 - dSigma = integral(...), 2 - приращения
        ms.fixGrid = fixGrid;     // 1 - зафиксировать сетку
        //ms.switchIterationsMode = SwitchIterationsMode::Serial;
        ms.switchIterationsMode = SwitchIterationsMode::Parallel;
        ms.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
        ms.controlMode = 0;     // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
        ms.terminateIfAccuracyIsNotAchieving = false;
        ms.iterLimit = nonlinearIterLimit;
        ms.slauResidualLimit = 10000;
        ms.plasticResidualLimit.eps = ccp.epsResidualLimit;
        ms.plasticResidualLimit.sigma = ccp.sigmaResidualLimit;
        ms.contactEndPointResidualLimit = ccp.contactEndPointResidualLimit;//1.e-10;//1.e-13;
        ms.contactDeltaFResidualLimit = ccp.contactDeltaFResidualLimit;//1.e-13;
        ms.material = material;
        ms.bc1Source = bc1Source;
        //ms.bc2Source = bc2Source; - будет меняться
        ms.bc2 = bc2;
        // вдавливание шара
        for(size_t i = 0; i < step->size(); i++)
        {
            setBc2_none(bc2Source);
            ms.bc2Source = bc2Source;
            /*
            if(unload && (int)i == stepsNumber)
                ms.forsedElasticD = false;
                //ms.forsedElasticD = true;
            else
                ms.forsedElasticD = false;
            */
            mechStep->push_back(ms);
        }
        task.mechTask.mechStep = mechStep;

        // Контакт

        // поверхность: жёсткая сфера
        task.mechTask.rigidSurface = new std::vector<Surface_base *>;
        (*task.mechTask.rigidSurface).resize(2);
        (*task.mechTask.rigidSurface)[0] = grid->analiticalSurface[0];   // поверхность заданная аналитически
        //(*task.mechTask.rigidSurface)[1] = grid->ISurface[0];            // поверхность заданная интерполянтом
        //(*task.mechTask.rigidSurface)[2] = &(grid->FESurface[3]);        // поверхность заданная сеткой (если contact_mode != 3 то эта дополнительная сетка не строится)

        // контакты КЭ-поверхность - жёсткая поверхность
        task.mechTask.CC_FE_Rigid = new std::vector<ContactCondition_FE_Rigid>;
        if(ccp.surfaceType != SurfaceType::none)
        {
            ContactCondition_FE_Rigid CC_FE_Rigid_el;
            CC_FE_Rigid_el.method = contactType;
            CC_FE_Rigid_el.constantNormal = ccp.constantNormal;
            CC_FE_Rigid_el.w_stiffness = ccp.w_stiffness;
            CC_FE_Rigid_el.FESurfaceInd = 0;
            if(ccp.surfaceType == SurfaceType::AnaliticalSurface_Sphere)
                CC_FE_Rigid_el.RigidSurfaceInd = 0;
            //if(ccp.surfaceType == SurfaceType::InterpolantSurface_Hermite3 ||
            //   ccp.surfaceType == SurfaceType::InterpolantSurface_Lagrange3)
            //    CC_FE_Rigid_el.RigidSurfaceInd = 1;
            //if(ccp.surfaceType == SurfaceType::FiniteElementSurface)
            //    CC_FE_Rigid_el.RigidSurfaceInd = 2;
            CC_FE_Rigid_el.noContactRadiusOptimization = noContactRadiusOptimization;
            (*task.mechTask.CC_FE_Rigid).push_back(CC_FE_Rigid_el);
        }
        // инициализация нулями начальных данных
        task.mechTask.initNull(Integration::IntegrationType::Gauss3);
    }

    // Температура
    if(task.thermTask.enabled)
    {
    }
}

// аналитическое решение (напряжения по отношению к среднему давлению)
void Test_contactSphere::calc_elastic_ms(const Grid::ContactSphereParameters &csp, const double a, const double z, const double r,
                                          VECTOR3 &r_fi_z, VECTOR3 &ms)
{
    double NU = csp.Nu;

    double w = sqrt(0.5*(SQR(r) + SQR(z) - SQR(a) + sqrt(SQR(SQR(r) + SQR(z) - SQR(a)) + 4*SQR(a*z))));
    double sigma_fi = -((1.-2.*NU)/2.*SQR(a/r)*(1 - pow(z/w, 3))
                      + 3./2.*z/w*(2*NU + (1 - NU)*SQR(w)/(SQR(a) + SQR(w)) - (1 + NU)*w/a*atan(a/w)));
    double sigma_z = -3./2.*pow(z/w, 3)*SQR(a*w)/(pow(w, 4) + SQR(a*z));
    double sigma_r = -(sigma_fi + sigma_z) - 3.*z/w*(1 + NU)*(1 - w/a*atan(a/w));
    double sigma_rz = -3./2.*r*w*SQR(z)/(pow(w, 4) + SQR(a*z))*SQR(a)/(SQR(a) + SQR(w));

    ms[0] = (sigma_r + sigma_z)/2. + sqrt(SQR((sigma_r - sigma_z)/2.) + SQR(sigma_rz));
    ms[2] = (sigma_r + sigma_z)/2. - sqrt(SQR((sigma_r - sigma_z)/2.) + SQR(sigma_rz));
    ms[1] = sigma_fi;
    r_fi_z[0] = sigma_r;
    r_fi_z[1] = sigma_fi;
    r_fi_z[2] = sigma_z;

    /*
    double NU = csp.Nu;
    double _rza = SQR(r) + SQR(z) - SQR(a);
    double u = 0.5*(_rza + sqrt(SQR(_rza) + 4*SQR(a*z)));
    double _uaz = SQR(u) + SQR(a*z);
    double _zu12 = z/sqrt(u);
    double _uttt = u*(1 - NU)/(SQR(a) + u);
    double _1ttt = (1+NU)*sqrt(u)/a*atan(a/sqrt(u));
    double sigma_r = 3./2.*((1.-2.*NU)/3.*SQR(a/r)*(1-pow(_zu12, 3))
                            + pow(_zu12, 3)*SQR(a)*u/_uaz
                            + _zu12*(_uttt + _1ttt - 2));   // точно -2 а не -2*NU??
    double sigma_fi = -3./2.*((1.-2.*NU)/3.*SQR(a/r)*(1-pow(_zu12, 3))
                              + _zu12*(2*NU + _uttt - _1ttt));
    double sigma_z = -3./2.*pow(z/sqrt(u), 3)*SQR(a)*u/_uaz;
    double sigma_rz = -3./2.*r*SQR(z)/_uaz*SQR(a)*sqrt(u)/(SQR(a) + u);
    ms[0] = (sigma_r + sigma_z)/2. + sqrt(SQR((sigma_r - sigma_z)/2.) + SQR(sigma_rz));
    ms[2] = (sigma_r + sigma_z)/2. - sqrt(SQR((sigma_r - sigma_z)/2.) + SQR(sigma_rz));
    ms[1] = sigma_fi;
    r_fi_z[0] = sigma_r;
    r_fi_z[1] = sigma_fi;
    r_fi_z[2] = sigma_z;
    */
}

// численное решение (напряжения по отношению к (аналитическому) среднему давлению)
void Test_contactSphere::calc_number_ms(const double a_analit, const double P_m_analit, const MechOutFePointData &fePointData,
                                         double &z, double &r, VECTOR3 &r_fi_z, VECTOR3 &ms)
{
    r = calc_r(fePointData.p);
    z = fePointData.p[2];
    // 1) главные напряжения
    {
        fePointData.mainStresses(ms);
        sortms_maxToMin(ms);
    }
    // 2) напряжения на площадках
    {
        MATR3x3 sigma3x3;
        fePointData.sumSigma.ToMATR3x3_sigma(sigma3x3); // тензор напряжений
        // sigma_r
        {
            VECTOR3 norm_r = VECTOR3(fePointData.p[0], fePointData.p[1], 0);    // направление r
            norm_r = norm_r / norm_r.abs();
            {
                VECTOR3 vectSigma = sigma3x3*norm_r;    // вектор напряжений на площадке
                r_fi_z[0] = vectSigma*norm_r;           // проекция на нормаль к площадке
            }
        }
        // sigma_fi
        {
            VECTOR3 norm_fi = VECTOR3(-fePointData.p[1], fePointData.p[0], 0);    // направление fi
            norm_fi = norm_fi / norm_fi.abs();
            {
                VECTOR3 vectSigma = sigma3x3*norm_fi;   // вектор напряжений на площадке
                r_fi_z[1] = vectSigma*norm_fi;          // проекция на нормаль к площадке
            }
        }
        // sigma_z
        {
            VECTOR3 norm_z = VECTOR3(0, 0, 1);    // направление z
            norm_z = norm_z / norm_z.abs();
            {
                VECTOR3 vectSigma = sigma3x3*norm_z;   // вектор напряжений на площадке
                r_fi_z[2] = vectSigma*norm_z;          // проекция на нормаль к площадке
            }
        }
    }
    ms /= P_m_analit;
    r_fi_z /= P_m_analit;
}

void Test_contactSphere::writeResults(const Task &task, const Solid::OutData &out)
{
    int globalStepNumber = (int)(*out.mechOut).size() - 2;
    //bool isLoadLastGlobalStep = (globalStepNumber == (int)task.mechTask.step->size() - 1);
    //bool isLoading = (globalStepNumber <= ccp.STEPS_PUSH*(1 + ccp.CORR) - 1);
    bool isUnloading = (globalStepNumber >= ccp.STEPS_PUSH*(1 + ccp.CORR) - 1);
    bool isLoadLastGlobalStep = (globalStepNumber == ccp.STEPS_PUSH*(1 + ccp.CORR) - 1);
     bool isUnLoadLastGlobalStep = (globalStepNumber == ccp.STEPS_PUSH*(1 + ccp.CORR)*2 + 1 - 1);
    const std::vector<MechOutVertexData> &vertexData = (*out.mechOut)[globalStepNumber].vertex;
    FiniteElementSurface &contactSurface = task.mechTask.grid->FESurface[0];

    // пути к файлам
    std::string subdir;
    {
    // подпапка для данного шага
    subdir = dir + std::to_string(globalStepNumber);
    {
        OS::Console c;
        c.exec(("mkdir " + subdir).c_str());
    }
    subdir += "/";
    // пути к текстовым файлам с данными (в корневой директории, общие для всех глобальных шагов)
    {
        fn_curve_diagramma = dir +"curve_diagramma.txt";
        fn_curve = dir + "curve.txt";
        fn_dcurve = dir + "dcurve.txt";

        fn_Fn = dir + "Fn.txt";
        fn_Fn_analit = dir + "Fn_analit.txt";
        fn_a = dir + "a.txt";
        fn_a_analit = dir + "a_analit.txt";
        fn_Pmax = dir + "Pmax.txt";
        fn_Pmax_analit = dir + "Pmax_analit.txt";

         fn_bP_bd = dir + "bP_bd.txt";
          fn_bP_bd_analit = dir + "bP_bd_analit.txt";
         fn_unload_bP_bd = dir + "unload_bP_bd.txt";
         fn_unload_bP_bd_analit = dir + "unload_bP_bd_analit.txt";

         fn_bA_bd_true = dir + "bA_bd_true.txt";
         fn_bA_bd_geom = dir + "bA_bd_geom.txt";
         fn_bA_bd_shamanstvo = dir + "bA_bd_shamanstvo.txt";
         fn_unload_bA_bd = dir + "unload_bA_bd.txt";
         fn_unload_bA_bd_analit = dir + "unload_bA_bd_analit.txt";

        // пути к текстовым файлам с данными (для каждого глобального шага в отдельной директории)
        fn_inf = subdir + "_inf.txt";
        fn_general_inf = dir + "_general_inf.txt";
        fn_P = subdir + "P.txt";
        fn_P_analit = subdir + "P_analit.txt";
        fn_h = subdir + "h.txt";
        fn_F = subdir + "F.txt";

        fn_z0_sigma_r = subdir + "z0_sigma_r.txt";
        fn_z0_sigma_fi = subdir + "z0_sigma_fi.txt";
        fn_z0_sigma_z = subdir + "z0_sigma_z.txt";
        fn_z0_sigma1 = subdir + "z0_sigma1.txt";
        fn_z0_sigma2 = subdir + "z0_sigma2.txt";
        fn_z0_sigma3 = subdir + "z0_sigma3.txt";

        /*
        fn_2d_sigma_r_analit = dir + "2d_sigma_r_analit.txt";
        fn_2d_sigma_fi_analit = dir + "2d_sigma_fi_analit.txt";
        fn_2d_sigma_z_analit = dir + "2d_sigma_z_analit.txt";
        fn_2d_sigma1_analit = dir + "2d_sigma1_analit.txt";
        fn_2d_sigma2_analit = dir + "2d_sigma2_analit.txt";
        fn_2d_sigma3_analit = dir + "2d_sigma3_analit.txt";
        */
        fn_2d_sigma_r = dir + "2d_sigma_r.txt";
        fn_2d_sigma_fi = dir + "2d_sigma_fi.txt";
        fn_2d_sigma_z = dir + "2d_sigma_z.txt";
        fn_2d_sigmaEqv = dir + "2d_sigmaEqv.txt";
        fn_2d_sigma1 = dir + "2d_sigma1.txt";
        fn_2d_sigma2 = dir + "2d_sigma2.txt";
        fn_2d_sigma3 = dir + "2d_sigma3.txt";
        fn_2d_sigmaEqv_unload = dir + "2d_sigmaEqv_unload.txt";
        fn_2d_sigma1_unload = dir + "2d_sigma1_unload.txt";
        fn_2d_sigma2_unload = dir + "2d_sigma2_unload.txt";
        fn_2d_sigma3_unload = dir + "2d_sigma3_unload.txt";
    }

    // пути к файлам для отображения графиков
{
    std::ofstream f;
    if(globalStepNumber == 0)
    {
        f.open(fn_filesList, std::ofstream::out);
        resGraph.generalGraphsNomber = 12 + 4*2;
        resGraph.graphsPerGlobalStepNomber = 3;
        f << resGraph.generalGraphsNomber << "\n";
        f << resGraph.graphsPerGlobalStepNomber << "\n";

        // curve_diagramma
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "curve_diagramma.png") +
             genGnuplotCommandLineParameter("fn_curve_diagramma", dir + "curve_diagramma.txt") +
             "\" " +
             dir + "curve_diagramma.gnu" << "\n";
        f << dir + "curve_diagramma.png" << "\n";
        f << "curve_diagramma.png" << "\n";
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

        // fn_unload_bA_bd
        if(ccp.plasticOut)
        {
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "unload_bA_bd.png") +
                 genGnuplotCommandLineParameter("fn_unload_bA_bd", fn_unload_bA_bd) +
                 genGnuplotCommandLineParameter("fn_unload_bA_bd_analit", fn_unload_bA_bd_analit) +
                 "\" " +
                 dir + "unload_bA_bd_pl.gnu" << "\n";
            f << dir + "unload_bA_bd.png" << "\n";
            f << "unload_bA_bd.png" << "\n";
        }
        else
        {
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "unload_bA_bd.png") +
                 genGnuplotCommandLineParameter("fn_unload_bA_bd", fn_unload_bA_bd) +
                 genGnuplotCommandLineParameter("fn_unload_bA_bd_analit", fn_unload_bA_bd_analit) +
                 "\" " +
                 dir + "unload_bA_bd_el.gnu" << "\n";
            f << dir + "unload_bA_bd.png" << "\n";
            f << "unload_bA_bd.png" << "\n";
        }

        // fn_bA_bd_true
        if(ccp.plasticOut)
        {
            std::string fn_out = "bA_bd_true.png";
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + fn_out) +
                 genGnuplotCommandLineParameter("fn_bA_bd", fn_bA_bd_true) +
                 "\" " +
                 dir + "bA_bd_pl.gnu" << "\n";
            f << dir + fn_out << "\n";
            f << fn_out << "\n";
        }
        else
        {
            std::string fn_out = "bA_bd_true.png";
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + fn_out) +
                 genGnuplotCommandLineParameter("fn_bA_bd", fn_bA_bd_true) +
                 "\" " +
                 dir + "bA_bd_el.gnu" << "\n";
            f << dir + fn_out << "\n";
            f << fn_out << "\n";
        }

        // fn_bA_bd_geom
        if(ccp.plasticOut)
        {
            std::string fn_out = "bA_bd_geom.png";
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + fn_out) +
                 genGnuplotCommandLineParameter("fn_bA_bd", fn_bA_bd_geom) +
                 "\" " +
                 dir + "bA_bd_pl.gnu" << "\n";
            f << dir + fn_out << "\n";
            f << fn_out << "\n";
        }
        else
        {
            std::string fn_out = "bA_bd_geom.png";
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + fn_out) +
                 genGnuplotCommandLineParameter("fn_bA_bd", fn_bA_bd_geom) +
                 "\" " +
                 dir + "bA_bd_el.gnu" << "\n";
            f << dir + fn_out << "\n";
            f << fn_out << "\n";
        }

        // fn_bA_bd_shamanstvo
        if(ccp.plasticOut)
        {
            std::string fn_out = "bA_bd_shamanstvo.png";
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + fn_out) +
                 genGnuplotCommandLineParameter("fn_bA_bd", fn_bA_bd_shamanstvo) +
                 "\" " +
                 dir + "bA_bd_pl.gnu" << "\n";
            f << dir + fn_out << "\n";
            f << fn_out << "\n";
        }
        else
        {
            std::string fn_out = "bA_bd_shamanstvo.png";
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + fn_out) +
                 genGnuplotCommandLineParameter("fn_bA_bd", fn_bA_bd_shamanstvo) +
                 "\" " +
                 dir + "bA_bd_el.gnu" << "\n";
            f << dir + fn_out << "\n";
            f << fn_out << "\n";
        }

        // fn_unload_bP_bd
        if(ccp.plasticOut)
        {
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "unload_bP_bd.png") +
                 genGnuplotCommandLineParameter("fn_unload_bP_bd", fn_unload_bP_bd) +
                 genGnuplotCommandLineParameter("fn_unload_bP_bd_analit", fn_unload_bP_bd_analit) +
                 "\" " +
                 dir + "unload_bP_bd_pl.gnu" << "\n";
            f << dir + "unload_bP_bd.png" << "\n";
            f << "unload_bP_bd.png" << "\n";
        }
        else
        {
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "unload_bP_bd.png") +
                 genGnuplotCommandLineParameter("fn_unload_bP_bd", fn_unload_bP_bd) +
                 genGnuplotCommandLineParameter("fn_unload_bP_bd_analit", fn_unload_bP_bd_analit) +
                 "\" " +
                 dir + "unload_bP_bd_el.gnu" << "\n";
            f << dir + "unload_bP_bd.png" << "\n";
            f << "unload_bP_bd.png" << "\n";
        }

        // fn_bP_bd
        if(ccp.plasticOut)
        {
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "bP_bd.png") +
                 genGnuplotCommandLineParameter("fn_bP_bd", fn_bP_bd) +
                 genGnuplotCommandLineParameter("fn_bP_bd_analit", fn_bP_bd_analit) +
                 "\" " +
                 dir + "bP_bd_pl.gnu" << "\n";
            f << dir + "bP_bd.png" << "\n";
            f << "bP_bd.png" << "\n";
        }
        else
        {
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "bP_bd.png") +
                 genGnuplotCommandLineParameter("fn_bP_bd", fn_bP_bd) +
                 genGnuplotCommandLineParameter("fn_bP_bd_analit", fn_bP_bd_analit) +
                 "\" " +
                 dir + "bP_bd_el.gnu" << "\n";
            f << dir + "bP_bd.png" << "\n";
            f << "bP_bd.png" << "\n";
        }

        // Fn
        if(ccp.plasticOut)
        {
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "Fn.png") +
                 genGnuplotCommandLineParameter("fn_Fn", fn_Fn) +
                 genGnuplotCommandLineParameter("fn_Fn_analit", fn_Fn_analit) +
                 "\" " +
                 dir + "Fn_pl.gnu" << "\n";
            f << dir + "Fn.png" << "\n";
            f << "Fn.png" << "\n";
        }
        else
        {
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "Fn.png") +
                 genGnuplotCommandLineParameter("fn_Fn", fn_Fn) +
                 genGnuplotCommandLineParameter("fn_Fn_analit", fn_Fn_analit) +
                 "\" " +
                 dir + "Fn.gnu" << "\n";
            f << dir + "Fn.png" << "\n";
            f << "Fn.png" << "\n";
        }
        // a
        if(ccp.plasticOut)
        {
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "a.png") +
                 genGnuplotCommandLineParameter("fn_a", fn_a) +
                 genGnuplotCommandLineParameter("fn_a_analit", fn_a_analit) +
                 "\" " +
                 dir + "a_pl.gnu" << "\n";
            f << dir + "a.png" << "\n";
            f << "a.png" << "\n";
        }
        else
        {
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "a.png") +
                 genGnuplotCommandLineParameter("fn_a", fn_a) +
                 genGnuplotCommandLineParameter("fn_a_analit", fn_a_analit) +
                 "\" " +
                 dir + "a.gnu" << "\n";
            f << dir + "a.png" << "\n";
            f << "a.png" << "\n";
        }
        // Pmax
        if(ccp.plasticOut)
        {
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "Pmax.png") +
                 genGnuplotCommandLineParameter("fn_Pmax", fn_Pmax) +
                 genGnuplotCommandLineParameter("fn_Pmax_analit", fn_Pmax_analit) +
                 "\" " +
                 dir + "Pmax_pl.gnu" << "\n";
            f << dir + "Pmax.png" << "\n";
            f << "Pmax.png" << "\n";
        }
        else
        {
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "Pmax.png") +
                 genGnuplotCommandLineParameter("fn_Pmax", fn_Pmax) +
                 genGnuplotCommandLineParameter("fn_Pmax_analit", fn_Pmax_analit) +
                 "\" " +
                 dir + "Pmax.gnu" << "\n";
            f << dir + "Pmax.png" << "\n";
            f << "Pmax.png" << "\n";
        }

        char NU[100];
        char a[100];
        double a_analit;
        {
            double time = (*out.mechOut)[globalStepNumber].step[0].t0;
            double d = fabs(ccp.z0 - ccp.z_fun->fun(time));  // сближение шара и полупространства
            double P_analit = calc_P(d, ccp);
            double R = ccp.R0;
            double E = ccp.E / (1 - SQR(ccp.Nu));   // приведённый модуль упругости
            a_analit = pow((3*R)/(4*E)*P_analit, 1./3.);
        }
        sprintf(NU, "%le", ccp.Nu);
        sprintf(a, "%le", a_analit);
        // нагружение
        if(ccp.plasticOut)
        {
            // 2d_sigmaEqv
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "2d_sigmaEqv.png") +
                 genGnuplotCommandLineParameter("fn_2d_sigmaEqv", dir + "2d_sigmaEqv.txt") +
                 genGnuplotCommandLineParameter("a", a) +
                 genGnuplotCommandLineParameter("NU", NU) +
                 "\" " +
                 dir + "2d_sigmaEqv_pl.gnu" << "\n";
            f << dir + "2d_sigmaEqv.png" << "\n";
            f << "2d_sigmaEqv.png" << "\n";
            // 2d_sigma1
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma1.png") +
                 genGnuplotCommandLineParameter("fn_2d_sigma1", dir + "2d_sigma1.txt") +
                 genGnuplotCommandLineParameter("a", a) +
                 genGnuplotCommandLineParameter("NU", NU) +
                 "\" " +
                 dir + "2d_sigma1_pl.gnu" << "\n";
            f << dir + "2d_sigma1.png" << "\n";
            f << "2d_sigma1.png" << "\n";
            // 2d_sigma2
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma2.png") +
                 genGnuplotCommandLineParameter("fn_2d_sigma2", dir + "2d_sigma2.txt") +
                 genGnuplotCommandLineParameter("a", a) +
                 genGnuplotCommandLineParameter("NU", NU) +
                 "\" " +
                 dir + "2d_sigma2_pl.gnu" << "\n";
            f << dir + "2d_sigma2.png" << "\n";
            f << "2d_sigma2.png" << "\n";
            // 2d_sigma3
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma3.png") +
                 genGnuplotCommandLineParameter("fn_2d_sigma3", dir + "2d_sigma3.txt") +
                 genGnuplotCommandLineParameter("a", a) +
                 genGnuplotCommandLineParameter("NU", NU) +
                 "\" " +
                 dir + "2d_sigma3_pl.gnu" << "\n";
            f << dir + "2d_sigma3.png" << "\n";
            f << "2d_sigma3.png" << "\n";
        }
        else
        {
            // 2d_sigmaEqv
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "2d_sigmaEqv.png") +
                 genGnuplotCommandLineParameter("fn_2d_sigmaEqv", dir + "2d_sigmaEqv.txt") +
                 genGnuplotCommandLineParameter("a", a) +
                 genGnuplotCommandLineParameter("NU", NU) +
                 "\" " +
                 dir + "2d_sigmaEqv.gnu" << "\n";
            f << dir + "2d_sigmaEqv.png" << "\n";
            f << "2d_sigmaEqv.png" << "\n";
            // 2d_sigma1, 2d_sigma1_analit
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma1.png") +
                 genGnuplotCommandLineParameter("fn_2d_sigma1", dir + "2d_sigma1.txt") +
                 genGnuplotCommandLineParameter("a", a) +
                 genGnuplotCommandLineParameter("NU", NU) +
                 "\" " +
                 dir + "2d_sigma1.gnu" << "\n";
            f << dir + "2d_sigma1.png" << "\n";
            f << "2d_sigma1.png" << "\n";
            // 2d_sigma2, 2d_sigma2_analit
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma2.png") +
                 genGnuplotCommandLineParameter("fn_2d_sigma2", dir + "2d_sigma2.txt") +
                 genGnuplotCommandLineParameter("a", a) +
                 genGnuplotCommandLineParameter("NU", NU) +
                 "\" " +
                 dir + "2d_sigma2.gnu" << "\n";
            f << dir + "2d_sigma2.png" << "\n";
            f << "2d_sigma2.png" << "\n";
            // 2d_sigma3, 2d_sigma3_analit
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma3.png") +
                 genGnuplotCommandLineParameter("fn_2d_sigma3", dir + "2d_sigma3.txt") +
                 genGnuplotCommandLineParameter("a", a) +
                 genGnuplotCommandLineParameter("NU", NU) +
                 "\" " +
                 dir + "2d_sigma3.gnu" << "\n";
            f << dir + "2d_sigma3.png" << "\n";
            f << "2d_sigma3.png" << "\n";
        }
        // разгрузка
        {
            // 2d_sigmaEqv
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "2d_sigmaEqv_unload.png") +
                 genGnuplotCommandLineParameter("fn_2d_sigmaEqv", dir + "2d_sigmaEqv_unload.txt") +
                 genGnuplotCommandLineParameter("a", a) +
                 genGnuplotCommandLineParameter("NU", NU) +
                 "\" " +
                 dir + "2d_sigmaEqv_pl.gnu" << "\n";
            f << dir + "2d_sigmaEqv_unload.png" << "\n";
            f << "2d_sigmaEqv_unload.png" << "\n";
            // 2d_sigma1
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma1_unload.png") +
                 genGnuplotCommandLineParameter("fn_2d_sigma1", dir + "2d_sigma1_unload.txt") +
                 genGnuplotCommandLineParameter("a", a) +
                 genGnuplotCommandLineParameter("NU", NU) +
                 "\" " +
                 dir + "2d_sigma1_pl.gnu" << "\n";
            f << dir + "2d_sigma1_unload.png" << "\n";
            f << "2d_sigma1_unload.png" << "\n";
            // 2d_sigma2
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma2_unload.png") +
                 genGnuplotCommandLineParameter("fn_2d_sigma2", dir + "2d_sigma2_unload.txt") +
                 genGnuplotCommandLineParameter("a", a) +
                 genGnuplotCommandLineParameter("NU", NU) +
                 "\" " +
                 dir + "2d_sigma2_pl.gnu" << "\n";
            f << dir + "2d_sigma2_unload.png" << "\n";
            f << "2d_sigma2_unload.png" << "\n";
            // 2d_sigma3
            f << std::string("-e \"") +
                 genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma3_unload.png") +
                 genGnuplotCommandLineParameter("fn_2d_sigma3", dir + "2d_sigma3_unload.txt") +
                 genGnuplotCommandLineParameter("a", a) +
                 genGnuplotCommandLineParameter("NU", NU) +
                 "\" " +
                 dir + "2d_sigma3_pl.gnu" << "\n";
            f << dir + "2d_sigma3_unload.png" << "\n";
            f << "2d_sigma3_unload.png" << "\n";
        }


/*
// r_fi_z аналитическое решение
        // 2d_sigma_r_analit
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma_r_analit.png") +
             genGnuplotCommandLineParameter("fn_in", dir + "2d_sigma_r_analit.txt") +
             "\" " +
             dir + "2d.gnu" << "\n";
        f << dir + "2d_sigma_r_analit.png" << "\n";
        f << "2d_sigma_r_analit.png" << "\n";
        // 2d_sigma_fi_analit
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma_fi_analit.png") +
             genGnuplotCommandLineParameter("fn_in", dir + "2d_sigma_fi_analit.txt") +
             "\" " +
             dir + "2d.gnu" << "\n";
        f << dir + "2d_sigma_fi_analit.png" << "\n";
        f << "2d_sigma_fi_analit.png" << "\n";
        // 2d_sigma_z_analit
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma_z_analit.png") +
             genGnuplotCommandLineParameter("fn_in", dir + "2d_sigma_z_analit.txt") +
             "\" " +
             dir + "2d.gnu" << "\n";
        f << dir + "2d_sigma_z_analit.png" << "\n";
        f << "2d_sigma_z_analit.png" << "\n";
// r_fi_z численное решение
        // 2d_sigma_r
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma_r.png") +
             genGnuplotCommandLineParameter("fn_in", dir + "2d_sigma_r.txt") +
             "\" " +
             dir + "2d.gnu" << "\n";
        f << dir + "2d_sigma_r.png" << "\n";
        f << "2d_sigma_r.png" << "\n";
        // 2d_sigma_fi
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma_fi.png") +
             genGnuplotCommandLineParameter("fn_in", dir + "2d_sigma_fi.txt") +
             "\" " +
             dir + "2d.gnu" << "\n";
        f << dir + "2d_sigma_fi.png" << "\n";
        f << "2d_sigma_fi.png" << "\n";
        // 2d_sigma_z
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma_z.png") +
             genGnuplotCommandLineParameter("fn_in", dir + "2d_sigma_z.txt") +
             "\" " +
             dir + "2d.gnu" << "\n";
        f << dir + "2d_sigma_z.png" << "\n";
        f << "2d_sigma_z.png" << "\n";
// sigma123 аналитическое решение
        // 2d_sigma1_analit
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma1_analit.png") +
             genGnuplotCommandLineParameter("fn_in", dir + "2d_sigma1_analit.txt") +
             "\" " +
             dir + "2d.gnu" << "\n";
        f << dir + "2d_sigma1_analit.png" << "\n";
        f << "2d_sigma1_analit.png" << "\n";
        // 2d_sigma2_analit
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma2_analit.png") +
             genGnuplotCommandLineParameter("fn_in", dir + "2d_sigma2_analit.txt") +
             "\" " +
             dir + "2d.gnu" << "\n";
        f << dir + "2d_sigma2_analit.png" << "\n";
        f << "2d_sigma2_analit.png" << "\n";
        // 2d_sigma3_analit
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma3_analit.png") +
             genGnuplotCommandLineParameter("fn_in", dir + "2d_sigma3_analit.txt") +
             "\" " +
             dir + "2d.gnu" << "\n";
        f << dir + "2d_sigma3_analit.png" << "\n";
        f << "2d_sigma3_analit.png" << "\n";
// sigma123 численное решение
        // 2d_sigma1
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma1.png") +
             genGnuplotCommandLineParameter("fn_in", dir + "2d_sigma1.txt") +
             "\" " +
             dir + "2d.gnu" << "\n";
        f << dir + "2d_sigma1.png" << "\n";
        f << "2d_sigma1.png" << "\n";
        // 2d_sigma2
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma2.png") +
             genGnuplotCommandLineParameter("fn_in", dir + "2d_sigma2.txt") +
             "\" " +
             dir + "2d.gnu" << "\n";
        f << dir + "2d_sigma2.png" << "\n";
        f << "2d_sigma2.png" << "\n";
        // 2d_sigma3
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "2d_sigma3.png") +
             genGnuplotCommandLineParameter("fn_in", dir + "2d_sigma3.txt") +
             "\" " +
             dir + "2d.gnu" << "\n";
        f << dir + "2d_sigma3.png" << "\n";
        f << "2d_sigma3.png" << "\n";
*/
    }
    else
    {
        f.open(fn_filesList, std::ofstream::app);
    }
    f << std::to_string(globalStepNumber) + "\n";
    // P
    if(ccp.plasticOut)
    {
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", subdir + "P.png") +
             genGnuplotCommandLineParameter("fn_P", subdir + "P.txt") +
             genGnuplotCommandLineParameter("fn_P_analit", subdir + "P_analit.txt") +
             genGnuplotCommandLineParameter("fn_z0_sigma1", subdir + "z0_sigma1.txt") +
             genGnuplotCommandLineParameter("fn_z0_sigma2", subdir + "z0_sigma2.txt") +
             genGnuplotCommandLineParameter("fn_z0_sigma3", subdir + "z0_sigma3.txt") +
             genGnuplotCommandLineParameter("fn_z0_sigma_r", subdir + "z0_sigma_r.txt") +
             genGnuplotCommandLineParameter("fn_z0_sigma_fi", subdir + "z0_sigma_fi.txt") +
             genGnuplotCommandLineParameter("fn_z0_sigma_z", subdir + "z0_sigma_z.txt") +
             "\" " +
             dir + "P_pl.gnu" << "\n";
        f << subdir + "P.png" << "\n";
        f << "P.png" << "\n";
    }
    else
    {
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", subdir + "P.png") +
             genGnuplotCommandLineParameter("fn_P", subdir + "P.txt") +
             genGnuplotCommandLineParameter("fn_P_analit", subdir + "P_analit.txt") +
             genGnuplotCommandLineParameter("fn_z0_sigma1", subdir + "z0_sigma1.txt") +
             genGnuplotCommandLineParameter("fn_z0_sigma2", subdir + "z0_sigma2.txt") +
             genGnuplotCommandLineParameter("fn_z0_sigma3", subdir + "z0_sigma3.txt") +
             genGnuplotCommandLineParameter("fn_z0_sigma_r", subdir + "z0_sigma_r.txt") +
             genGnuplotCommandLineParameter("fn_z0_sigma_fi", subdir + "z0_sigma_fi.txt") +
             genGnuplotCommandLineParameter("fn_z0_sigma_z", subdir + "z0_sigma_z.txt") +
             "\" " +
             dir + "P.gnu" << "\n";
        f << subdir + "P.png" << "\n";
        f << "P.png" << "\n";
    }
    // h
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "h.png") +
         genGnuplotCommandLineParameter("fn_h", subdir + "h.txt") +
         "\" " +
         dir + "h.gnu" << "\n";
    f << subdir + "h.png" << "\n";
    f << "h.png" << "\n";
    // F
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "F.png") +
         genGnuplotCommandLineParameter("fn_F", subdir + "F.txt") +
         "\" " +
         dir + "F.gnu" << "\n";
    f << subdir + "F.png" << "\n";
    f << "F.png" << "\n";

    f.close();
}
    resGraph.read_gnuplot_from_file(fn_filesList);
    }
    // сохранение сетки
    //task.save(subdir);


     FILE *f_inf = fopen(fn_inf.c_str(), "w");
    FILE *f_P = fopen(fn_P.c_str(), "w");
    FILE *f_P_analit = fopen(fn_P_analit.c_str(), "w");
    FILE *f_h = fopen(fn_h.c_str(), "w");
    FILE *f_F = fopen(fn_F.c_str(), "w");

    FILE *f_z0_sigma_r = fopen(fn_z0_sigma_r.c_str(), "w");
    FILE *f_z0_sigma_fi = fopen(fn_z0_sigma_fi.c_str(), "w");
    FILE *f_z0_sigma_z = fopen(fn_z0_sigma_z.c_str(), "w");
    FILE *f_z0_sigma1 = fopen(fn_z0_sigma1.c_str(), "w");
    FILE *f_z0_sigma2 = fopen(fn_z0_sigma2.c_str(), "w");
    FILE *f_z0_sigma3 = fopen(fn_z0_sigma3.c_str(), "w");

    /*
    FILE *f_2d_sigma_r_analit = fopen(fn_2d_sigma_r_analit.c_str(), "w");
    FILE *f_2d_sigma_fi_analit = fopen(fn_2d_sigma_fi_analit.c_str(), "w");
    FILE *f_2d_sigma_z_analit = fopen(fn_2d_sigma_z_analit.c_str(), "w");
    FILE *f_2d_sigma1_analit = fopen(fn_2d_sigma1_analit.c_str(), "w");
    FILE *f_2d_sigma2_analit = fopen(fn_2d_sigma2_analit.c_str(), "w");
    FILE *f_2d_sigma3_analit = fopen(fn_2d_sigma3_analit.c_str(), "w");
    */

    FILE *f_2d_sigma_r;
    FILE *f_2d_sigma_fi;
    FILE *f_2d_sigma_z;
    FILE *f_2d_sigmaEqv;
    FILE *f_2d_sigma1;
    FILE *f_2d_sigma2;
    FILE *f_2d_sigma3;
    FILE *f_2d_sigmaEqv_unload;
    FILE *f_2d_sigma1_unload;
    FILE *f_2d_sigma2_unload;
    FILE *f_2d_sigma3_unload;
    if(isLoadLastGlobalStep)
    {
        f_2d_sigma_r = fopen(fn_2d_sigma_r.c_str(), "w");
        f_2d_sigma_fi = fopen(fn_2d_sigma_fi.c_str(), "w");
        f_2d_sigma_z = fopen(fn_2d_sigma_z.c_str(), "w");
        f_2d_sigmaEqv = fopen(fn_2d_sigmaEqv.c_str(), "w");
        f_2d_sigma1 = fopen(fn_2d_sigma1.c_str(), "w");
        f_2d_sigma2 = fopen(fn_2d_sigma2.c_str(), "w");
        f_2d_sigma3 = fopen(fn_2d_sigma3.c_str(), "w");
    }
    if(isUnLoadLastGlobalStep)
    {
        f_2d_sigmaEqv_unload = fopen(fn_2d_sigmaEqv_unload.c_str(), "w");
        f_2d_sigma1_unload = fopen(fn_2d_sigma1_unload.c_str(), "w");
        f_2d_sigma2_unload = fopen(fn_2d_sigma2_unload.c_str(), "w");
        f_2d_sigma3_unload = fopen(fn_2d_sigma3_unload.c_str(), "w");
    }


    FILE *f_Fn;
    FILE *f_a;
    FILE *f_Pmax;
    FILE *f_bP_bd;
    FILE *f_unload_bP_bd;

    FILE *f_bA_bd_true;
    FILE *f_bA_bd_geom;
    FILE *f_bA_bd_shamanstvo;
    FILE *f_unload_bA_bd;

    std::string fileOpenMode;
    if(globalStepNumber == 0)
    {
        fileOpenMode = "w";
    }
    else
    {
        fileOpenMode = "a";
    }
    f_Fn = fopen(fn_Fn.c_str(), fileOpenMode.c_str());
    f_a = fopen(fn_a.c_str(), fileOpenMode.c_str());
    f_Pmax = fopen(fn_Pmax.c_str(), fileOpenMode.c_str());
    f_bP_bd = fopen(fn_bP_bd.c_str(), fileOpenMode.c_str());
    f_unload_bP_bd = fopen(fn_unload_bP_bd.c_str(), fileOpenMode.c_str());

    f_bA_bd_true = fopen(fn_bA_bd_true.c_str(), fileOpenMode.c_str());
    f_bA_bd_geom = fopen(fn_bA_bd_geom.c_str(), fileOpenMode.c_str());
    f_bA_bd_shamanstvo = fopen(fn_bA_bd_shamanstvo.c_str(), fileOpenMode.c_str());
    f_unload_bA_bd = fopen(fn_unload_bA_bd.c_str(), fileOpenMode.c_str());


    // кривая пластичности
    if(globalStepNumber == 0)
    {
        MechMaterialSource &m0 = (*(*task.mechTask.mechStep)[globalStepNumber].material)[0];
        // кривая пластичности
        {

            FILE *fcurve = fopen(fn_curve.c_str(), "w");
            FILE *fdcurve = fopen(fn_dcurve.c_str(), "w");
            double q1 = 0;
            double q2 = 0.004; // strToDouble(m0.vals[1]) - предел упругости
            int N = 1000;//1000000;
            for(int i = 0; i <= N; i++)
            {
                double q = q1+(q2-q1)*i/N;
                double F = m0.sigma_Yeld(m0.elasticParameters0, q, 0);
                double dF = m0.difSigma_Yeld(m0.elasticParameters0, q, 0);
                fprintf(fcurve, "%le %le\n", q, F);
                fprintf(fdcurve, "%le %le\n", q, dF);
            }
            fclose(fcurve);
            fclose(fdcurve);
        }
        // диаграмма одноосного деформирования
        {
            FILE *fcurve_diagramma= fopen(fn_curve_diagramma.c_str(), "w");
            double eps1 = 0;
            double eps2 = 0.04; // strToDouble(m0.vals[1]) - предел упругости
            int N = 1000;//1000000;
            for(int i = 0; i <= N; i++)
            {
                double eps = eps1+(eps2-eps1)*i/N;
                double E = m0.elasticParameters0.E;
                double sigma_y = m0.elasticSigmaLimit;
                double eps_y = sigma_y/E;
                double diagr;
                if(eps < eps_y)
                {
                    diagr = E*eps;
                }
                else
                {
                    diagr = sigma_y*pow(eps/eps_y, ccp.n);
                }
                fprintf(fcurve_diagramma, "%le %le\n", eps, diagr);
            }
            fclose(fcurve_diagramma);
        }
    }

    // реакции опоры и расстояния от уонтактных узлов до поверхности
    VECTOR3 Fn_sum = VECTOR3_NULL;
    for(size_t vertexInd = 0; vertexInd < vertexData.size(); vertexInd++)
    {
        bool contacted = vertexData[vertexInd].contact;
        if(contacted)
        {
            double r = calc_r(task.grid->vertex[vertexInd]);
            VECTOR3 F = vertexData[vertexInd].F_sum;
            double h = vertexData[vertexInd].h;
            fprintf(f_F, "%lf\t%le\n", r, F.abs());
            if(F.abs() != 0)
                fprintf(f_h, "%lf\t%le\n", r, h);
            Fn_sum += vertexData[vertexInd].F_sum;
        }
    }
    // площадь зоны контакта
    double S_true;          // учёт количества контактных узлов на граничных ячейках
    double S_geom;          // учёт пересечения жёсткой поверхности с граничными ячейками
    double S_shamanstvo;    // учёт величины давления на граничных ячейках
    std::vector<VECTOR3> vertexP;   // давления в узлах
    contactSurface.grid = task.mechTask.grid;
    Post::calcContactArea(contactSurface, *task.mechTask.grid, vertexData,
                           vertexP, S_true, S_geom, S_shamanstvo);
    S_true *= 4;
    S_geom *= 4;
    S_shamanstvo *= 4;
    // вывод результата
    double time = (*out.mechOut)[globalStepNumber].step[0].t0;
    double d = fabs(ccp.z_fun->fun(time) - ccp.z0);  // сближение шара и полупространства
    //double d = ccp.z0 - ccp.z_fun->fun(time);  // сближение шара и полупространства
    double R = ccp.R0;
    double E = ccp.E / (1 - SQR(ccp.Nu));   // приведённый модуль упругости
    // область которая будет выводиться на 2d графике
    double z_out_max;
    double r_out_max;
    // аналитическое решение
    double P_m_analit;// среднее давление
    double a_analit;
    // множитель для вывода нормализованных значений в случае пластичности
    double mnoj = 0;
    {
        double P_analit;
        if(!ccp.plasticOut && isUnLoadLastGlobalStep)
        {
            d = 0;
            P_analit = 0;
        }
        else
        {
            P_analit = calc_P(d, ccp);
        }
        a_analit = pow((3*R)/(4*E)*P_analit, 1./3.);
        double Pmax_analit = 3*P_analit/(2*PI*SQR(a_analit));
        P_m_analit = P_analit / (PI * SQR(a_analit));
        if(ccp.plasticOut)
        {
            mnoj = P_m_analit / ccp.elasticSigmaLimit;
        }
        else
        {
            if(isUnLoadLastGlobalStep)
                mnoj = P_m_analit;  // после упругой разгрузки интересуют абсолютные значения
            else
                mnoj = 1;
        }
        z_out_max = a_analit*5;
        r_out_max = a_analit*5;

        // давление на z = 0
        int rN = 1000;
        for(int ri = 0; ri <= rN; ri++)
        {
            double r = a_analit * ri / rN;
            fprintf(f_P_analit, "%le\t%le\n", r/a_analit, -Pmax_analit*sqrt(1 - SQR(r/a_analit)) / P_m_analit * mnoj);
        }
        // 2d решение
        /*
        {
            int rN = 100;
            int zN = 100;
            for(int zi = 0; zi <= zN; zi++)
            {
                double z = -(z_out_min + (z_out_max - z_out_min) * zi / zN);
                for(int ri = 0; ri <= rN; ri++)
                {
                    double r = r_out_min + (r_out_max - r_out_min) * ri / rN;

                    //double rza = SQR(r) + SQR(z) - SQR(a_analit);
                    //double u = 0.5*(rza + sqrt(SQR(rza) + 4*SQR(a_analit*z)));
                    //double sigma_z = -3./2.*pow(z/sqrt(u), 3)*SQR(a_analit)*u/(SQR(u)+SQR(a_analit*z));
                    //fprintf(f_2d_sigma_z, "%le\t%le\t%le\n", r/a_analit, z/a_analit, sigma_z);

                    //if(r/a_analit >= 1 && z != 0)
                    if(z != 0 && r != 0)
                    {
                        VECTOR3 ms;
                        VECTOR3 r_fi_z;
                        calc_elastic_ms(ccp, a_analit, -z, r, r_fi_z, ms);
                        fprintf(f_2d_sigma_r_analit, "%le\t%le\t%le\n", r/a_analit, z/a_analit, r_fi_z[0]);
                        fprintf(f_2d_sigma_fi_analit, "%le\t%le\t%le\n", r/a_analit, z/a_analit, r_fi_z[1]);
                        fprintf(f_2d_sigma_z_analit, "%le\t%le\t%le\n", r/a_analit, z/a_analit, r_fi_z[2]);
                        fprintf(f_2d_sigma1_analit, "%le\t%le\t%le\n", r/a_analit, z/a_analit, ms[0]);
                        fprintf(f_2d_sigma2_analit, "%le\t%le\t%le\n", r/a_analit, z/a_analit, ms[1]);
                        fprintf(f_2d_sigma3_analit, "%le\t%le\t%le\n", r/a_analit, z/a_analit, ms[2]);
                    }
                }
            }
        }*/
    }

    double Pmax_numb = -1.e10;

    // численное решение на приповерхностных КЭ
    for(size_t faceInd = 0; faceInd < contactSurface.face.size(); faceInd++)
    {
        const Grid::FEFace &faceEl = contactSurface.face[faceInd];
        // главные напряжения в приповерхностных КЭ
        {
            MechOutFePointData &fePointData = (*out.mechOut)[globalStepNumber].fe[faceEl.feIndex].pd[0];
            double r;
            double z;
            VECTOR3 ms;
            VECTOR3 r_fi_z;
            calc_number_ms(a_analit, P_m_analit, fePointData, z, r, r_fi_z, ms);
            // z ~ 0
            // 1) главные напряжения
            if(r/a_analit < 2)
            {
                fprintf(f_z0_sigma1, "%le\t%le\n", r/a_analit, ms[0]*mnoj);
                fprintf(f_z0_sigma2, "%le\t%le\n", r/a_analit, ms[1]*mnoj);
                fprintf(f_z0_sigma3, "%le\t%le\n", r/a_analit, ms[2]*mnoj);
            }
            // 2) напряжения на площадках
            if(r/a_analit < 2)
            {
                fprintf(f_z0_sigma_r, "%le\t%le\n", r/a_analit, r_fi_z[0]*mnoj);
                fprintf(f_z0_sigma_fi, "%le\t%le\n", r/a_analit, r_fi_z[1]*mnoj);
                fprintf(f_z0_sigma_z, "%le\t%le\n", r/a_analit, r_fi_z[2]*mnoj);
            }

            // расчёт максимального напряжения sigma_z в приповерхностных КЭ
            //if(fabs(r_fi_z[2]*mnoj) > Pmax_numb)
            //    Pmax_numb = fabs(r_fi_z[2]*mnoj);
        }
    }

    // давление по Y в контактных узлах
    for(size_t vertexInd = 0; vertexInd < vertexData.size(); vertexInd++)
    {
        VECTOR3 P_numb;
        P_numb = vertexP[vertexInd];
        if(P_numb.abs() > 0)
        {
            double r = calc_r(task.grid->vertex[vertexInd]);
            fprintf(f_P, "%le\t%le\n", r/a_analit, P_numb[2] / P_m_analit * mnoj);
            //fprintf(f_P, "%le\t%le\n", r/a_analit, -P_numb.abs() / P_m_analit);
            // расчёт максимального давления
            //if(P_numb.abs() > Pmax_numb)
            //    Pmax_numb = P_numb.abs();
            if(r < 1.e-10)// центральный узел
            {
                if(ccp.plasticOut)
                {
                    Pmax_numb = P_numb[2] / ccp.elasticSigmaLimit;
                }
                else
                {
                    Pmax_numb = -P_numb[2];
                }
            }
        }
    }

    // численное решение в центрах всех КЭ
    if(isLoadLastGlobalStep || isUnLoadLastGlobalStep)
    for(size_t feIndex = 0; feIndex < task.mechTask.grid->fe.size(); feIndex++)
    {
        // главные напряжения в центре КЭ
        MechOutFePointData &fePointData = (*out.mechOut)[globalStepNumber].fe[feIndex].pd[0];
        MechMaterialSource &m0 = (*(*task.mechTask.mechStep)[globalStepNumber].material)[0];
        double r;
        double z;
        VECTOR3 ms;
        VECTOR3 r_fi_z;
        calc_number_ms(a_analit, P_m_analit, fePointData, z, r, r_fi_z, ms);
        if(r < r_out_max && -z < z_out_max)
        {
            if(isLoadLastGlobalStep)
            {
                // 1) главные и эквивалентные напряжения
                fprintf(f_2d_sigmaEqv, "%le\t%le\t%le\n", r/a_analit, z/a_analit, fePointData.sumSigma.eqv_SIGMA()/P_m_analit*mnoj);
                fprintf(f_2d_sigma1, "%le\t%le\t%le\n", r/a_analit, z/a_analit, ms[0]*mnoj);
                fprintf(f_2d_sigma2, "%le\t%le\t%le\n", r/a_analit, z/a_analit, ms[1]*mnoj);
                fprintf(f_2d_sigma3, "%le\t%le\t%le\n", r/a_analit, z/a_analit, ms[2]*mnoj);
                // 2) напряжения на площадках
                fprintf(f_2d_sigma_r, "%le\t%le\t%le\n", r/a_analit, z/a_analit, r_fi_z[0]*mnoj);
                fprintf(f_2d_sigma_fi, "%le\t%le\t%le\n", r/a_analit, z/a_analit, r_fi_z[1]*mnoj);
                fprintf(f_2d_sigma_z, "%le\t%le\t%le\n", r/a_analit, z/a_analit, r_fi_z[2]*mnoj);
            }
            if(isUnLoadLastGlobalStep)
            {
                // главные и эквивалентные напряжения
                fprintf(f_2d_sigmaEqv_unload, "%le\t%le\t%le\n", r/a_analit, z/a_analit, fePointData.sumSigma.eqv_SIGMA()/P_m_analit*mnoj);
                fprintf(f_2d_sigma1_unload, "%le\t%le\t%le\n", r/a_analit, z/a_analit, ms[0]*mnoj);
                fprintf(f_2d_sigma2_unload, "%le\t%le\t%le\n", r/a_analit, z/a_analit, ms[1]*mnoj);
                fprintf(f_2d_sigma3_unload, "%le\t%le\t%le\n", r/a_analit, z/a_analit, ms[2]*mnoj);
            }
        }
    }

    // максимальные значения эквивалентных напряжений и деформаций, и параметра Одквиста
    {
        if(globalStepNumber == 0)
        {
            ccp.max_sigma_eqv = 0;
            ccp.max_epsElastoPlastic_eqv = 0;
            ccp.max_q = 0;
        }
        for(size_t feIndex = 0; feIndex < task.mechTask.grid->fe.size(); feIndex++)
        {
            MechOutFePointData &fePointData = (*out.mechOut)[globalStepNumber].fe[feIndex].pd[0];
            MechMaterialSource &m0 = (*(*task.mechTask.mechStep)[globalStepNumber].material)[0];
            ccp.max_sigma_eqv = MAX(ccp.max_sigma_eqv, fePointData.sumSigma.eqv_SIGMA());
            ccp.max_epsElastoPlastic_eqv = MAX(ccp.max_epsElastoPlastic_eqv, (fePointData.sumEpsElastic + fePointData.sumEpsPlastic).eqv_EPS());
            ccp.max_q = MAX(ccp.max_q, fePointData.q);
        }
    }


    // вывод без всяких обезразмериваний
    // численное решение
    {
        double P_numb = fabs(Fn_sum[2])*4;
        FILE *f_P_d;
        FILE *f_A_d_true;
        FILE *f_A_d_geom;
        FILE *f_A_d_shamanstvo;
        std::string f_mode;
        if(globalStepNumber == 0)
        {
             f_mode = "w";
        }
        else
        {
            f_mode = "a";
        }
        f_P_d = fopen((dir + "P_d.txt").c_str(), f_mode.c_str());
        f_A_d_true = fopen((dir + "A_d_true.txt").c_str(), f_mode.c_str());
        f_A_d_geom = fopen((dir + "A_d_geom.txt").c_str(), f_mode.c_str());
        f_A_d_shamanstvo = fopen((dir + "A_d_shamanstvo.txt").c_str(), f_mode.c_str());
        if(P_numb != 0)
        {
            fprintf(f_P_d, "%le\t%le\n", d, P_numb);
            fprintf(f_A_d_true, "%le\t%le\n", d, S_true);
            fprintf(f_A_d_geom, "%le\t%le\n", d, S_geom);
            fprintf(f_A_d_shamanstvo, "%le\t%le\n", d, S_shamanstvo);
        }
        fclose(f_P_d);
        fclose(f_A_d_true);
        fclose(f_A_d_geom);
        fclose(f_A_d_shamanstvo);
    }
    // аналитические решения для разной глубины вдавливания
    if(globalStepNumber == 0)
    {
        FILE *f_P_analit;
        FILE *f_A_analit;
        f_P_analit = fopen((dir + "P_d_analit.txt").c_str(), "w");
        f_A_analit = fopen((dir + "A_d_analit.txt").c_str(), "w");
        int dN = 1000;
        for(int di = 0; di <= dN; di++)
        {
            double d = ccp.d0*di/dN;
            double P_analit = calc_P(d, ccp);  // нормальная сила контактного взаимодействия
            double a_analit = pow((3*R)/(4*E)*P_analit, 1./3.);// радиус области контакта
            double A_analit = PI*SQR(a_analit); // площадь области контакта
            fprintf(f_P_analit, "%le\t%le\n", d, P_analit);
            fprintf(f_A_analit, "%le\t%le\n", d, A_analit);
        }
        fclose(f_P_analit);
        fclose(f_A_analit);
    }


    // общее
    {
        double P_numb = fabs(Fn_sum[2]);    //==Fn_numb
        //double a_numb = r_max;
        double a_numb = sqrt(S_shamanstvo/PI);//r_average;
        P_numb *= 4;    // симметрия: рассматриваем четвертинку
        //Pmax_numb расчитано выше
        fprintf(f_Fn, "%le\t%le\n", d, P_numb);
        fprintf(f_a, "%le\t%le\n", d, a_numb);
        if(Pmax_numb == -1.e10)
        {
            fprintf(f_Pmax, "%le\t%le\n", d, 0.);
        }
        else
        {
            fprintf(f_Pmax, "%le\t%le\n", d, Pmax_numb);
        }
        // безразмерные сила и безразмерный радиус зоны контакта в зависимости от безразмерной глубины
        {
            // нагружение численное решение
            double bd = d / ccp.d_y;
            double bP_numb = P_numb / ccp.P_y;
            double bA_numb_true = S_true / ccp.A_y;
            double bA_numb_geom = S_geom / ccp.A_y;
            double bA_numb_shamanstvo = S_shamanstvo / ccp.A_y;
            //double bA_numb = (PI * SQR(a_numb)) / ccp.A_y;
            if(bP_numb != 0)
            {
                fprintf(f_bP_bd, "%le\t%le\n", bd, bP_numb);
                fprintf(f_bA_bd_true, "%le\t%le\n", bd, bA_numb_true);
                fprintf(f_bA_bd_geom, "%le\t%le\n", bd, bA_numb_geom);
                fprintf(f_bA_bd_shamanstvo, "%le\t%le\n", bd, bA_numb_shamanstvo);
            }
            if(isLoadLastGlobalStep)
            {
                ccp.bP_max_numb = bP_numb;
                ccp.bA_max_numb_shamanstvo = bA_numb_shamanstvo;
            }
            // (не совпадает)нагружение аналитическое решение для разной глубины вдавливания
            if(globalStepNumber == 0)
            {
                FILE *f_bP_bd_analit;
                f_bP_bd_analit = fopen(fn_bP_bd_analit.c_str(), fileOpenMode.c_str());
                int dN = 1000;
                for(int di = 0; di <= dN; di++)
                {
                    double d = ccp.d0*di/dN;
                    double bd = d / ccp.d_y;
                    double bP_analit = (4.0374 + 0.62*pow(ccp.n, 0.45)) * pow(bd, 1.0929 + 0.4*ccp.n);
                    if(ccp.plasticOut)
                    {
                        fprintf(f_bP_bd_analit, "%le\t%le\n", bd, bP_analit);
                    }
                    else
                    {
                        //fprintf(f_bP_bd_analit, "%le\t%le\n", bd, bP_analit);
                    }
                }
                fclose(f_bP_bd_analit);
            }
            // разгрузка численное решение
            double unload_bd_res;
            /*
            if(ccp.n == 0)
            {
                unload_bd_res = ccp.bd_max * (1 - pow(ccp.bd_max, -1./3.)) * (1 - pow(ccp.bd_max, -2./3.));
            }
            else
                */
            {
                double dmax = ccp.bd_max * ccp.d_y;
                double ksi = 3.*PI*1.08/4./sqrt(2.)*pow(0.3, ccp.n);
                double tt = SQR(ksi)*pow(dmax/2./R, ccp.n - 1.)*pow(ccp.bYeldCoeff, 2.*(ccp.n - 1.));
                unload_bd_res = ccp.bd_max * (1. - pow(tt, 1./3.)) * (1. - pow(tt, 2./3.));
            }
            double unload_bd = (bd - unload_bd_res) / (ccp.bd_max - unload_bd_res);
            double unload_bP_numb = bP_numb / ccp.bP_max_numb;
            double unload_bA_numb_shamanstvo = bA_numb_shamanstvo / ccp.bA_max_numb_shamanstvo;
            if(Pmax_numb != -1.e10 && isUnloading)
            {
                fprintf(f_unload_bP_bd, "%le\t%le\n", unload_bd, unload_bP_numb);
                fprintf(f_unload_bA_bd, "%le\t%le\n", unload_bd, unload_bA_numb_shamanstvo);
            }
            if(Pmax_numb == -1.e10 && isUnloading)
            {
                // при первои отлипании сохраняется глубина
                if(ccp.unload_bd_res_numb == -1)
                {
                    // отлипание ещё не происходило
                    ccp.unload_bd_res_numb = bd / ccp.bd_max;
                }
            }
            // разгрузка аналитическое решение для разной глубины вдавливания
            if(globalStepNumber == 0)
            {
                FILE *f_unload_bP_bd_analit = fopen(fn_unload_bP_bd_analit.c_str(), fileOpenMode.c_str());
                FILE *f_unload_bA_bd_analit = fopen(fn_unload_bA_bd_analit.c_str(), fileOpenMode.c_str());;
                int bdN = 1000;

                for(int bdi = 0; bdi <= bdN; bdi++)
                {
                    // unload_bd_res ... ccp.bd_max
                    double bd = unload_bd_res + (ccp.bd_max - unload_bd_res)*bdi/bdN;
                    double unload_bd = (bd - unload_bd_res) / (ccp.bd_max - unload_bd_res);
                    double Alpha = 3./2.*pow(ccp.bd_max, -0.036);
                    double unload_bP_analit = pow(unload_bd, Alpha);
                    double Betta = 0.9*(1. - pow(ccp.bd_max, -0.036));
                    double unload_bA_analit = 1./Betta * (1. - pow(1. - Betta, unload_bd));
                    if(ccp.plasticOut)
                    {
                        fprintf(f_unload_bP_bd_analit, "%le\t%le\n", unload_bd, unload_bP_analit);
                        fprintf(f_unload_bA_bd_analit, "%le\t%le\n", unload_bd, unload_bA_analit);
                    }
                    else
                    {
                        //fprintf(f_unload_bP_bd_analit, "%le\t%le\n", unload_bd, unload_bP_analit);
                    }
                }
                fclose(f_unload_bP_bd_analit);
                fclose(f_unload_bA_bd_analit);
            }
            // в отдельный файл выводятся параметры и глубина отлипания по отношению к максимальной глубине
            if(isUnLoadLastGlobalStep)
            {
                FILE *f_general_inf = fopen(fn_general_inf.c_str(), "w");
                fprintf(f_general_inf, "n = %le\n", ccp.n);
                fprintf(f_general_inf, "bd_max = %le\n", ccp.bd_max);
                fprintf(f_general_inf, "bYeldCoeff = %le\n", ccp.bYeldCoeff);
                fprintf(f_general_inf, "\n");

                fprintf(f_general_inf, "NN = %d\n", ccp.NN);
                fprintf(f_general_inf, "STEPS_PUSH = %d\n", ccp.STEPS_PUSH);
                fprintf(f_general_inf, "grid_mnojitel1 = %le\n", ccp.grid_mnojitel1);
                fprintf(f_general_inf, "grid_mnojitel2 = %le\n", ccp.grid_mnojitel2);
                fprintf(f_general_inf, "\n");


                fprintf(f_general_inf, "R0 = %le\n", ccp.R0);
                fprintf(f_general_inf, "E = %le\n", ccp.E);
                fprintf(f_general_inf, "Nu = %le\n", ccp.Nu);


                fprintf(f_general_inf, "plasticityMethodType = %d\n", ccp.plasticityMethodType);
                fprintf(f_general_inf, "incForsesMode = %d\n", ccp.incForsesMode);
                fprintf(f_general_inf, "constantNormal = %d\n", ccp.constantNormal);

                fprintf(f_general_inf, "w_project = %le\n", ccp.w_project);
                fprintf(f_general_inf, "w_midPoint = %le\n", ccp.w_midPoint);
                fprintf(f_general_inf, "w_stiffness = %le\n", ccp.w_stiffness);
                fprintf(f_general_inf, "cosTettaMin = %le\n", ccp.cosTettaMin);
                fprintf(f_general_inf, "\n");

                fprintf(f_general_inf, "sigmaResidualLimit = %le\n", ccp.sigmaResidualLimit);
                fprintf(f_general_inf, "epsResidualLimit = %le\n", ccp.epsResidualLimit);
                fprintf(f_general_inf, "contactDeltaFResidualLimit = %le\n", ccp.contactDeltaFResidualLimit);
                fprintf(f_general_inf, "contactEndPointResidualLimit = %le\n", ccp.contactEndPointResidualLimit);
                fprintf(f_general_inf, "\n");
                fprintf(f_general_inf, "____________________\n");

                fprintf(f_general_inf, "c1 = %le\n", ccp.c1);
                fprintf(f_general_inf, "c2 = %le\n", ccp.c2);
                fprintf(f_general_inf, "c3 = %le\n", ccp.c3);
                fprintf(f_general_inf, "q = %le\n", ccp.q);
                fprintf(f_general_inf, "\n");

                fprintf(f_general_inf, "elasticSigmaLimit = %le\n", ccp.elasticSigmaLimit);
                fprintf(f_general_inf, "max_sigma_eqv = %le\n", ccp.max_sigma_eqv);
                fprintf(f_general_inf, "max_epsElastoPlastic_eqv = %le\n", ccp.max_epsElastoPlastic_eqv);
                fprintf(f_general_inf, "max_q = %le\n", ccp.max_q);

                fprintf(f_general_inf, "d0 = %le\n", ccp.d0);
                fprintf(f_general_inf, "d_y = %le\n", ccp.d_y);
                fprintf(f_general_inf, "P_y = %le\n", ccp.P_y);
                fprintf(f_general_inf, "a_y = %le\n", ccp.A_y);

                 fprintf(f_general_inf, "\tbP_max_numb = %le\n", ccp.bP_max_numb);
                 fprintf(f_general_inf, "\tbd_residual = %le\n", ccp.unload_bd_res_numb);

                fclose(f_general_inf);
            }
        }


        // аналитические решения для разной глубины вдавливания
        if(globalStepNumber == 0)
        {
            FILE *f_Fn_analit;
            FILE *f_a_analit;
            FILE *f_Pmax_analit;
            f_Fn_analit = fopen(fn_Fn_analit.c_str(), "w");
            f_a_analit = fopen(fn_a_analit.c_str(), "w");
            f_Pmax_analit = fopen(fn_Pmax_analit.c_str(), "w");
            int dN = 1000;
            for(int di = 0; di <= dN; di++)
            {
                double d = ccp.d0*di/dN;
                double P_analit = calc_P(d, ccp);  // нормальная сила контактного взаимодействия
                double a_analit = pow((3*R)/(4*E)*P_analit, 1./3.);// радиус области контакта
                double P_m_analit = P_analit / (PI * SQR(a_analit));
                double Pmax_analit = 3*P_analit/(2*PI*SQR(a_analit));
                fprintf(f_Fn_analit, "%le\t%le\n", d, P_analit);
                fprintf(f_a_analit, "%le\t%le\n", d, a_analit);
                if(ccp.plasticOut)
                {
                    fprintf(f_Pmax_analit, "%le\t%le\n", d, -Pmax_analit / ccp.elasticSigmaLimit);
                }
                else
                {
                    fprintf(f_Pmax_analit, "%le\t%le\n", d, Pmax_analit);
                }
            }
            fclose(f_Fn_analit);
            fclose(f_a_analit);
            fclose(f_Pmax_analit);
        }

        double P_analit = calc_P(d, ccp);
        double a_analit = pow((3*R)/(4*E)*P_analit, 1./3.);
        double Pmax_analit = 3*P_analit/(2*PI*SQR(a_analit));
        /*
        int rN = 1000;
        for(int ri = 0; ri <= rN; ri++)
        {
            double r = a_analit * ri / rN;
            fprintf(f_P_analit, "%le\t%le\n", r, Pmax_analit*sqrt(1 - SQR(r/a_analit)));
        }*/
        fprintf(f_inf, "d = %le\n", d);
        fprintf(f_inf, "a = %le\n", a_analit);
        fprintf(f_inf, "Pmax = %le\n", Pmax_analit);
        fprintf(f_inf, "Fn = %le\n", P_analit);
        fprintf(f_inf, "\nFn_sum = %le, %le, %le\n", Fn_sum[0], Fn_sum[1], Fn_sum[2]);

        fprintf(f_inf, "\n");
    }

    fclose(f_Fn);
    fclose(f_a);
    fclose(f_Pmax);
    fclose(f_bP_bd);
    fclose(f_unload_bP_bd);
    fclose(f_bA_bd_true);
    fclose(f_bA_bd_geom);
    fclose(f_bA_bd_shamanstvo);
    fclose(f_unload_bA_bd);

    fclose(f_inf);
    fclose(f_P);
    fclose(f_P_analit);
    fclose(f_h);
    fclose(f_F);

    fclose(f_z0_sigma_r);
    fclose(f_z0_sigma_fi);
    fclose(f_z0_sigma_z);
    fclose(f_z0_sigma1);
    fclose(f_z0_sigma2);
    fclose(f_z0_sigma3);

    /*
    fclose(f_2d_sigma_r_analit);
    fclose(f_2d_sigma_fi_analit);
    fclose(f_2d_sigma_z_analit);
    fclose(f_2d_sigma1_analit);
    fclose(f_2d_sigma2_analit);
    fclose(f_2d_sigma3_analit);
    */

    if(isLoadLastGlobalStep)
    {
        fclose(f_2d_sigma_r);
        fclose(f_2d_sigma_fi);
        fclose(f_2d_sigma_z);
        fclose(f_2d_sigmaEqv);
        fclose(f_2d_sigma1);
        fclose(f_2d_sigma2);
        fclose(f_2d_sigma3);
    }
    if(isUnLoadLastGlobalStep)
    {
        fclose(f_2d_sigmaEqv_unload);
        fclose(f_2d_sigma1_unload);
        fclose(f_2d_sigma2_unload);
        fclose(f_2d_sigma3_unload);
    }

    // построение графиков
    makeGraphsAfterGlobalStep(globalStepNumber, task);
}
bool Test_contactSphere::possibleToShow2d() const
{
    return true;
}
bool Test_contactSphere::getContactSurfaceCircle(const Task &task, const int globalStepIndex, POINT2 &c, double &R) const
{
    const auto &s = (*task.mechTask.step)[globalStepIndex];
    (*task.mechTask.rigidSurface)[0]->update(s.t_finish, nullptr);
    R = ((Interpolation::AnaliticalSurface_Sphere*)(*task.mechTask.rigidSurface)[0])->R;
    c[0] = ((Interpolation::AnaliticalSurface_Sphere*)(*task.mechTask.rigidSurface)[0])->C[0];
    c[1] = ((Interpolation::AnaliticalSurface_Sphere*)(*task.mechTask.rigidSurface)[0])->C[2];
    return true;
}
void Test_contactSphere::notFixedCoordinates(int &cn1, int &cn2)
{
    cn1 = 0;
    cn2 = 2;
}
void Test_contactSphere::needToDrawFe(const Solid::Task &task, const OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const
{
    //FE_LinearHexagon *fe_el = (*task.mechTask.fe)[feInd];
    //MechFeData_base *fe_el = (*task.mechTask.fe)[feInd];
    FE_base *fe_el = task.mechTask.grid->fe[feInd];
    int vi_8[8];
    fe_el->getVertexIndexes(vi_8);
    for(int i = 0; i < 8; i++)
    {
        //POINT3 p0 = task.mechTask.grid->vertex[vi_8[i]];
        POINT3 p0 = (*out.mechOut)[0].vertex[vi_8[i]].p;
        if(fabs(p0[1] - 0) < 1.e-10 /*&& feInd == task.mechTask.grid->fe.size() - 3*/)
        {
            // найден КЭ который примыкает к плоскости y = 0
            faceState = {0, 0, 0, 0, 0, 1};
            return;
        }
    }
    faceState = {0, 0, 0, 0, 0, 0};
}
bool Test_contactSphere::needToPaintFiniteElementSurface() const
{
    return ccp.surfaceType == Grid::SurfaceType::FiniteElementSurface;
}
}


namespace Grid
{

void Grid3D::genContactSphere(const ContactSphereParameters &ccp)
{
    using namespace Operations;
    build(ccp.av, ccp.ah);
    // параметры индексации пространства
    {
        regionIndexationParameters.h_min = 2;
        regionIndexationParameters.h_max = 6;
        regionIndexationParameters.div = 1;
        regionIndexationParameters.q0.i[0][0] = -11;
        regionIndexationParameters.q0.i[0][1] = +11;
        regionIndexationParameters.q0.i[1][0] = -21;
        regionIndexationParameters.q0.i[1][1] = +10;
        regionIndexationParameters.q0.i[2][0] = -11;
        regionIndexationParameters.q0.i[2][1] = +11;
    }

    // поверхность заданная аналитически
    analiticalSurface.resize(1);
    if(ccp.surfaceType == SurfaceType::AnaliticalSurface_Sphere)
    {
        GridRectangleRegular2D it_grid;// не нужна
        Interpolation::AnaliticalSurface_Sphere *sphere =
                new Interpolation::AnaliticalSurface_Sphere(
                    ccp.decompositionType,
                    it_grid,
                    ccp.R_fun,
                    ccp.x_fun,
                    ccp.y_fun,
                    ccp.z_fun,
                    0, //0
                    0, //0
                    0, //0
                    0, //0
                    1.e-8);
        analiticalSurface[0] = sphere;
    }
    // поверхность заданная интерполянтом
    ISurface.resize(0);

    if(ccp.sortNodes)
    {
        std::vector<int> contactFESurfaceIndex;
        contactFESurfaceIndex.push_back(0); // КЭ поверхность с индексом 0 - контактсная
        //cuthillMcKee(contactFESurfaceIndex);
        nestedDissection(contactFESurfaceIndex);
    }
}

}
// давления в центрах граней
/*
for(size_t faceInd = 0; faceInd < contactSurface.face.size(); faceInd++)
{
    const Grid::FEFace &faceEl = contactSurface.face[faceInd];
    const Grid::FE_base *feEl = task.mechTask.grid->fe[faceEl.feIndex];
    int vi[4];
    feEl->getFaceVertexIndexes(faceEl.faceIndex, vi);
    POINT3 c = VECTOR3_NULL;       // центр поверхности
    VECTOR3 P_numb = VECTOR3_NULL; // давление в центре поверхности
    bool needOut = true;
    for(int i = 0; i < 4; i++)
    {
        c += task.grid->vertex[vi[i]];
        P_numb += vertexP[vi[i]];
        if(vertexForce[vi[i]].abs() == 0)
            needOut = false;
    }
    c /= 4;
    P_numb /= 4;
    needOut = true;// крайние тоже выводить
    if(P_numb.abs() > 0 && needOut)
    {
        if(P_numb.abs() > Pmax_numb)
            Pmax_numb = P_numb.abs();
        double r = calc_r(c);
        fprintf(f_P, "%le\t%le\n", r, P_numb.abs());
    }
}
*/
