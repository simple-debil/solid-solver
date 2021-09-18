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
Test_contactCylinder::Test_contactCylinder()
{
    // дирректория данных теста
    dir = "./tests/contactCylinder/";
    // имена gnuplot файлов(+пути к файлам в опциях) и соответствующих картинок
    fn_filesList = dir + "files.txt";
    resGraph.read_gnuplot_from_file(fn_filesList);
    // инициализация графиков с невязками
    initStepsGraphs();
}
Test_base::Type Test_contactCylinder::get_type() const
{
    return Type::ContactCylinder;
}
// расстояние от центра поверхности контакта(круга) до точки с этой поверхности
// считаем просто координату x а не длину дуги
double Test_contactCylinder::calc_r(const POINT3 &p)
{
    return fabs(p[0]);
}
double Test_contactCylinder::calc_P(const double d, const Grid::ContactCylinderParameters &ccp)
{
    double E = ccp.E / (1 - SQR(ccp.Nu));
    double P;
    double P1 = 0, P2 = 1.e10;
    // функция d(P) растёт
    for(int i = 0; i < 100; i++)
    {
        P = (P1 + P2) / 2;
        //double d_trial = 2*P*(1-ccp.Nu*ccp.Nu)/PI/ccp.E * (-0.5*log(P) + log(2*ccp.Ly*sqrt(PI*E/ccp.R0)) - ccp.Nu/(2.*(1.-ccp.Nu)));
        double d_trial = P*(1-SQR(ccp.Nu))/PI/ccp.E * (log(SQR(2*ccp.Ly)*PI*E/ccp.R0/P) - ccp.Nu/(1.-ccp.Nu));
        //double d_trial = P*(1-SQR(ccp.Nu))/PI/ccp.E * (log(SQR(2*ccp.Ly)*PI*E/ccp.R0/P) - 1);
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
void Test_contactCylinder::initTask(Solid::Task &task)
{
    // подрубка
    task.mechTask.enabled = true;
    task.thermTask.enabled = false;
    // общие шаги и сетка
    std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
    Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
    // способ представления поверхности
    //gp.surfaceType = SurfaceType::none;
    ccp.surfaceType = SurfaceType::AnaliticalSurface_Cylinder;
    //gp.surfaceType = SurfaceType::InterpolantSurface_Hermite3;
    //gp.surfaceType = SurfaceType::InterpolantSurface_Lagrange3;
    //gp.contact_mode = SurfaceType::FiniteElementSurface;//нет поиска пересечения

    // способ декомпозиции поверхности
    ccp.decompositionType = RegionDecompositionType::none;
    //gp.decompositionType = DecompositionType::cubeVertexIndexation;
    //gp.decompositionType = DecompositionType::surfaceAsSpheres;

    // подрубка универсальной оптимизации
    bool noContactRadiusOptimization = false;
    //bool noContactRadiusOptimization = true;


    int unload = 0;   // 1 - есть разгрузка, 0 - нет разгрузки
    int NN = 16;      // подробность сетки
    double c1 = 0.8;
    double q = pow(1./40., 1./NN); //1 - 1./4.; // сгущение
    ccp.halhGrid = true;
    //ccp.halhGrid = false;

    ccp.d = 0.05;//+1.0;//0.5//0.1
    ccp.STEPS_PUSH = 1;  // количество шагов вдавливания и вынимания
    ccp.CORR = 0;        // количество шагов коррекции после каждого шага (1 или 0)

    ccp.E = 1.e10;
    ccp.Nu = 0.45;

    // 3D
    bool is2Dxy = false;
    // 2D
    //bool is2Dxy = true;

    //упругость
    MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::Elasticity;
    // пластичность
    //MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::D_pl;
    //MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::InitialSigma;
     //MechPlasticityMethodType plasticityMethodType = MechPlasticityMethodType::Combo_D_pl_InitialSigma;

    MechPlasticityCurveUnloadingType PCUnloadingType = MechPlasticityCurveUnloadingType::Elastic;
    //MechPlasticityCurveUnloadingType PCUnloadingType = MechPlasticityCurveUnloadingType::Curve;
    double k_sigma_eps = 1000.*10000;
    double k_eps_sigma = 10;


    int plasticPlot = 1;            // 0 - упругость, 1 - Безье/кусочно-линейная
    double epsResidualLimit = 1.e-10;            // желаемая относительная погрешность по эквиволентным деформациям
    double elasticSigmaLimit = 5.e7;//1.5e8;//2.0e8;//3.e8;
    // течение начинается на глубине 0.78a кога p0 > 1.67*elasticSigmaLimit
    //1.67*2.2=3.674

    Solid::ContactType contactType = Solid::ContactType::AugmentedLagrange;
    //Solid::ContactType contactType = Solid::ContactType::Fastest;
    ccp.w_stiffness = 1;//10;

    //IncForsesMode incForsesMode = IncForsesMode::IncrementP;
    //IncForsesMode incForsesMode = IncForsesMode::bPlusR;
    IncForsesMode incForsesMode = IncForsesMode::MinusIntegral;


    // параметры
    int fixGrid = 0;                // 0 - подвижная сетка, 1 - зафиксировать сетку
    //int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny


    int stepsNumber = ccp.STEPS_PUSH*(1 + ccp.CORR)*(1 + unload);//128;//64*5;//40/10;                        // количество шагов
    int nonlinearIterLimit = 300;                // ограничение на количество итераций
    double contactDeltaFResidualLimit = 1.e-10;;//1.e-12;//1.e-7;  // желаемая относительная погрешность по силе реакции опоры
    double contactEndPointResidualLimit = 1.e-14;//1.e-15

    double Time = stepsNumber;  // 50 - только нагружение, 100 - нагружение и разгрузка

    double z0 = -0.1*1;//-0.1;
    double z1 = +0.1*1;//+0.1;
    ccp.L = fabs(z1 - z0);


{
    std::vector<POINT3> &av = ccp.av;
    std::vector<HexagonArea> &ah = ccp.ah;
    const int s = 14;
    ccp.Nz = 1;
    av.resize(s*2);
    for(int i = 0; i <= 1; i++)
    {
        double z;
        int di;
        if(i == 0)
        {
            z = z0;
            di = 0;
        }
        else
        {
            z = z1;
            di = s;
        }
        //пробую: ####
        const double xscale = 1;
        const double yscale = 1;

        const double x1 = 10*xscale;    //10*xscale;
        const double x2 = c1*xscale;   //2*xscale;
        const double x3 = 0*xscale;     //0*xscale;
         const double y1 = -20*yscale;  //-20*yscale;
        const double y2 = -10*yscale;   //-10*yscale;
        const double y3 = -c1*yscale;  //-2*yscale;
        const double y4 = 0*yscale;     //0*yscale;
        av[di + 0] = POINT3{-x1, y1, z};
        av[di + 1] = POINT3{ x3, y1, z};
        av[di + 2] = POINT3{+x1, y1, z};
        av[di + 3] = POINT3{-x1, y2, z};
        av[di + 4] = POINT3{ x3, y2, z};
        av[di + 5] = POINT3{+x1, y2, z};
        av[di + 6] = POINT3{-x2, y3, z};
        av[di + 7] = POINT3{ x3, y3, z};
        av[di + 8] = POINT3{+x2, y3, z};
        av[di +  9] = POINT3{-x1, y4, z};
        av[di + 10] = POINT3{-x2, y4, z};
        av[di + 11] = POINT3{ x3, y4, z};
        av[di + 12] = POINT3{+x2, y4, z};
        av[di + 13] = POINT3{+x1, y4, z};
        ccp.Ly = fabs(y1 - y4);
    }
    ah.resize(8);
    // вершины
    // порядок вершин такой, что
    // z -> x
    // x -> y
    // y -> z
    // (поверхность Y = +1 -> Z = +1)
    ah[0].vi = {0, 0+s, 1, 1+s, 3, 3+s, 4, 4+s};
    ah[1].vi = {3, 3+s, 4, 4+s, 6, 6+s, 7, 7+s};
    ah[2].vi = {3, 3+s, 6, 6+s, 9, 9+s, 10, 10+s};
    ah[3].vi = {6, 6+s, 7, 7+s, 10, 10+s, 11, 11+s};

    ah[4].vi = {7, 7+s, 8, 8+s, 11, 11+s, 12, 12+s};
    ah[5].vi = {8, 8+s, 5, 5+s, 12, 12+s, 13, 13+s};
    ah[6].vi = {4, 4+s, 5, 5+s, 7, 7+s, 8, 8+s};
    ah[7].vi = {1, 1+s, 2, 2+s, 4, 4+s, 5, 5+s};

    for(int i = 0; i < 8; i++)
    {
        // разбиения
        ah[i].N = {ccp.Nz, NN, NN};
        // сгущение
        ah[i].condensation_coord = -1;
        //материал
        ah[i].mi = 0;
        // поверхности(пока нету)
        ah[i].surfaceIndex = {-1, -1, -1, -1, -1, -1};
        // 1-е краевые(фиксация по Z (Z = +1 -> X = +1))
        ah[i].bc1Index = {-1, -1, 2, 2, -1, -1};
    }
    // разбиения
    ah[0].N = {ccp.Nz, NN, NN/2};
    ah[7].N = {ccp.Nz, NN, NN/2};
    // сгущение
    ah[1].condensation_coord = 2;
    ah[1].condensation_q = q;
    ah[2].condensation_coord = 1;
    ah[2].condensation_q = q;

    // поверхность контакта (Y = +1 -> Z = +1)
    ah[3].surfaceIndex = {-1, 0, -1, -1, -1, -1};
    ah[4].surfaceIndex = {-1, 0, -1, -1, -1, -1};
    // 1-е краевые (по Z фиксируем везде)
    // фиксация по X на x=0 (X = +1 -> Y = +1)
    ah[0].bc1Index = {-1, -1, 2, 2, -1, 0};
    ah[1].bc1Index = {-1, -1, 2, 2, -1, 0};
    ah[3].bc1Index = {-1, -1, 2, 2, -1, 0};
    // фиксация по X на x=0 (X = -1 -> Y = -1)
    ah[4].bc1Index = {-1, -1, 2, 2, 0, -1};
    ah[6].bc1Index = {-1, -1, 2, 2, 0, -1};
    ah[7].bc1Index = {-1, -1, 2, 2, 0, -1};
    // фиксация по Y на y=-20 (Y = -1 -> Z = -1)
    ah[0].bc1Index = {1, -1, 2, 2, -1, 0};// дополнение предыдущего
    ah[7].bc1Index = {1, -1, 2, 2, 0, -1};// дополнение предыдущего

    // пробую: ####
    /*
    // фиксация по X и Y на y=-20
    ah[0].bc1Index = {3, -1, 2, 2, -1, 0};// дополнение предыдущего
    ah[7].bc1Index = {3, -1, 2, 2, 0, -1};// дополнение предыдущего
    // фиксация по X на x=-10 (X = -1 -> Y = -1)
    ah[0].bc1Index = {3, -1, 2, 2, 0, 0};// дополнение предыдущего
    ah[2].bc1Index = {-1, -1, 2, 2, 0, -1};
    // фиксация по X на x=+10 (X = -1 -> Y = -1)
    ah[5].bc1Index = {-1, -1, 2, 2, -1, 0};
    ah[7].bc1Index = {3, -1, 2, 2, 0, 0};// дополнение предыдущего
    */

    if(ccp.halhGrid)
    {
        //убираем вторую половинку
        ah.resize(4);
    }
}


    // параметры подвижного цилиндра
    ccp.y0 = 5;
    ccp.R0 = 5;
    ccp.R_fun = new Interpolation::Interpolant1D_Lagrange1Unregular;
    ccp.x_fun = new Interpolation::Interpolant1D_Lagrange1Unregular;
    ccp.y_fun = new Interpolation::Interpolant1D_Lagrange1Unregular;

    Grid::GridRectangleRegular1D gr;
    ccp.R_fun->init(gr, 0);// аргументы игнорируются для нерегулярной сетки
    ccp.x_fun->init(gr, 0);
    ccp.y_fun->init(gr, 0);
    ccp.R_fun->addPoint(0, ccp.R0);
    ccp.R_fun->addPoint(1000, ccp.R0);
    ccp.R_fun->buildInterpolant();
    ccp.x_fun->addPoint(0, 0);
    ccp.x_fun->addPoint(1000, 0);
    ccp.x_fun->buildInterpolant();




    grid->genContactCylinder(ccp);
    //grid->buldOpenSCADModel("setka.scad");

    // шаги
    GlobalStep s;       // шаг
    // вдавливание
    for(int i = 0; i < stepsNumber; i++)
    {
        s.t_start = Time*0 + Time*i/stepsNumber;
        s.t_finish = Time*0 + Time*(i + 1)/stepsNumber;
        s.dt0 = Time/stepsNumber;
        step->push_back(s);
        //ccp.C = {0, 5, 0};
        ccp.y_fun->addPoint(s.t_start, ccp.y0 - s.t_start/Time*ccp.d);
        if(i == stepsNumber - 1)
            ccp.y_fun->addPoint(s.t_finish, ccp.y0 - s.t_finish/Time*ccp.d);
    }
    ccp.y_fun->buildInterpolant();

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
        mechMat.set_E_NU(ccp.E, ccp.Nu);
        mechMat.set_M_sigma();
if(is2Dxy)
mechMat.set_2Dxy();// 2D
        mechMat.F = {0,0,0};// объёмные силы
        mechMat.elasticParameters0.ro = 0;     // плотность
            mechMat.elasticParameters0.Talpha = 0*1.e-5; // Коэффициент линейного расширения
        setPlasticMaterialCurve_Yeld(plasticPlot, mechMat, k_sigma_eps, k_eps_sigma);
        //setPlasticMaterialCurve(plasticPlot, mechMat, k_sigma_eps, k_eps_sigma);
        mechMat.plasticityMethodType = plasticityMethodType;
        mechMat.PCUnloadingType = PCUnloadingType;
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
        (*bc1Source)[3].mode = {{0, 0, 0}};
        (*bc1Source)[3].u0 =   {{0, 0, 0}};
        (*bc1Source)[4].mode = {{-1, -1,  -1}};
        (*bc1Source)[4].u0 =   {{-1, -1,  -1}};
        // вторые краевые
        bc2 = new std::vector<MechBoundaryCondition2>;
        bc2->clear();
        // шаги (на каждом шаге задаются вторые краевые условия)
        std::vector<MechGlobalStepParameters> *mechStep = new std::vector<MechGlobalStepParameters>;
        MechGlobalStepParameters ms;  // шаг для МДТТ
        GlobalStep s;       // шаг
        ms.slausolverParameters = slausolver_parameters;
        ms.timeMode = 0;    // 0 - квазистатическая задача
        ms.incForsesMode = incForsesMode; // 1 - dSigma = integral(...), 2 - приращения
        ms.fixGrid = fixGrid;     // 1 - зафиксировать сетку
        ms.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
        //ms.switchIterationsMode = SwitchIterationsMode::Serial;
        ms.switchIterationsMode = SwitchIterationsMode::Parallel;
        ms.controlMode = 0;     // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
        ms.terminateIfAccuracyIsNotAchieving = false;
        ms.iterLimit = nonlinearIterLimit;
        ms.slauResidualLimit = 10000;
        ms.plasticResidualLimit.eps = epsResidualLimit;
        ms.plasticResidualLimit.sigma = 10000;
        ms.contactEndPointResidualLimit = contactEndPointResidualLimit;//1.e-10;//1.e-13;
        ms.contactDeltaFResidualLimit = contactDeltaFResidualLimit;//1.e-13;
        ms.material = material;
        ms.bc1Source = bc1Source;
        //ms.bc2Source = bc2Source; - будет меняться
        ms.bc2 = bc2;
        // вдавливание цилиндра
        s = (*step)[0];
        for(int i = 0; i < stepsNumber; i++)
        {
            s.t_start = Time*0 + Time*i/stepsNumber;
            s.t_finish = Time*0 + Time*(i + 1)/stepsNumber;
            setBc2_none(bc2Source);
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
        }
        task.mechTask.mechStep = mechStep;

        // Контакт

        // поверхность: жёсткий цилиндр
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
            CC_FE_Rigid_el.constantNormal = false;
            CC_FE_Rigid_el.w_stiffness = ccp.w_stiffness;
            CC_FE_Rigid_el.FESurfaceInd = 0;
            if(ccp.surfaceType == SurfaceType::AnaliticalSurface_Cylinder)
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
void Test_contactCylinder::writeResults(const Task &task, const Solid::OutData &out)
{
    //task.grid->buldOpenSCADModel("setka_.scad");
    int globalStepNumber = (int)(*out.mechOut).size() - 2;
    const std::vector<MechOutVertexData> &vertexData = (*out.mechOut)[globalStepNumber].vertex;
    FiniteElementSurface &contactSurface = task.mechTask.grid->FESurface[0];
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
    fn_Fn = dir + "Fn.txt";
    fn_Fn_analit = dir + "Fn_analit.txt";
    fn_a = dir + "a.txt";
    fn_a_analit = dir + "a_analit.txt";
    fn_Pmax = dir + "Pmax.txt";
    fn_Pmax_analit = dir + "Pmax_analit.txt";
    // пути к текстовым файлам с данными (для каждого глобального шага в отдельной директории)
    fn_inf = subdir + "_inf.txt";
    fn_P = subdir + "P.txt";
    fn_P_analit_d = subdir + "P_analit_d.txt";
    fn_P_analit_F = subdir + "P_analit_F.txt";
    fn_h = subdir + "h.txt";
    fn_Pxyz = subdir + "Pxyz.txt";
    fn_F = subdir + "F.txt";
    fn_sigma1 = subdir + "sigma1.txt";
    fn_sigma2 = subdir + "sigma2.txt";
    fn_sigma3 = subdir + "sigma3.txt";
    fn_stiffness = subdir + "stiffness.txt";
    // пути к файлам для отображения графиков
    {
    std::ofstream f;
    if(globalStepNumber == 0)
    {
        f.open(fn_filesList, std::ofstream::out);
        resGraph.generalGraphsNomber = 3;
        resGraph.graphsPerGlobalStepNomber = 4;
        f << resGraph.generalGraphsNomber << "\n";
        f << resGraph.graphsPerGlobalStepNomber << "\n";
        // Fn
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "Fn.png") +
             genGnuplotCommandLineParameter("fn_Fn", fn_Fn) +
             genGnuplotCommandLineParameter("fn_Fn_analit", fn_Fn_analit) +
             "\" " +
             dir + "Fn.gnu" << "\n";
        f << dir + "Fn.png" << "\n";
        f << "Fn.png" << "\n";
        // a
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "a.png") +
             genGnuplotCommandLineParameter("fn_a", fn_a) +
             genGnuplotCommandLineParameter("fn_a_analit", fn_a_analit) +
             "\" " +
             dir + "a.gnu" << "\n";
        f << dir + "a.png" << "\n";
        f << "a.png" << "\n";
        // Pmax
        f << std::string("-e \"") +
             genGnuplotCommandLineParameter("fn_out", dir + "Pmax.png") +
             genGnuplotCommandLineParameter("fn_Pmax", fn_Pmax) +
             genGnuplotCommandLineParameter("fn_Pmax_analit", fn_Pmax_analit) +
             "\" " +
             dir + "Pmax.gnu" << "\n";
        f << dir + "Pmax.png" << "\n";
        f << "Pmax.png" << "\n";
    }
    else
    {
        f.open(fn_filesList, std::ofstream::app);
    }
    f << std::to_string(globalStepNumber) + "\n";
    // P
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "P.png") +
         genGnuplotCommandLineParameter("fn_P", subdir + "P.txt") +
         genGnuplotCommandLineParameter("fn_sigma1", subdir + "sigma1.txt") +
         genGnuplotCommandLineParameter("fn_sigma2", subdir + "sigma2.txt") +
         genGnuplotCommandLineParameter("fn_sigma3", subdir + "sigma3.txt") +
         genGnuplotCommandLineParameter("fn_P_analit_d", subdir + "P_analit_d.txt") +
         genGnuplotCommandLineParameter("fn_P_analit_F", subdir + "P_analit_F.txt") +
         "\" " +
         dir + "P.gnu" << "\n";
    f << subdir + "P.png" << "\n";
    f << "P.png" << "\n";
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
    // stiffness
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "stiffness.png") +
         genGnuplotCommandLineParameter("fn_stiffness", subdir + "stiffness.txt") +
         "\" " +
         dir + "stiffness.gnu" << "\n";
    f << subdir + "stiffness.png" << "\n";
    f << "stiffness.png" << "\n";
    f.close();
    resGraph.read_gnuplot_from_file(fn_filesList);
    }
    }

    FILE *f_Fn;
    FILE *f_a;
    FILE *f_Pmax;

     FILE *f_inf = fopen(fn_inf.c_str(), "w");
    FILE *f_P = fopen(fn_P.c_str(), "w");
    FILE *f_sigma1 = fopen(fn_sigma1.c_str(), "w");
    FILE *f_sigma2 = fopen(fn_sigma2.c_str(), "w");
    FILE *f_sigma3 = fopen(fn_sigma3.c_str(), "w");
    FILE *f_P_analit_d = fopen(fn_P_analit_d.c_str(), "w");
    FILE *f_P_analit_F = fopen(fn_P_analit_F.c_str(), "w");
    FILE *f_h = fopen(fn_h.c_str(), "w");
     FILE *f_Pxyz = fopen(fn_Pxyz.c_str(), "w");
     FILE *f_F = fopen(fn_F.c_str(), "w");
     FILE *f_stiffness = fopen(fn_stiffness.c_str(), "w");

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
    // реакции опоры
    for(size_t vertexInd = 0; vertexInd < vertexData.size(); vertexInd++)
    {
        if(task.thermTask.enabled)
        {
            //Heat::ThermOutVertexData &r = (*out->thermOut)[globalStepNumber].vertex[vertexInd];
            //POINT3 v = task.grid->vertex[vertexInd];
            //fprintf(f_T, "%le %le\n", v[1], r.T);
            //fprintf(f_T, "%le\t%le\n", (double)i1, r.T);
        }
        bool contacted = vertexData[vertexInd].contact;
        if(contacted)
        {
            //POINT3 p0 = (*out.mechOut)[globalStepNumber].vertex[vertexInd].p;
            POINT3 p0 = task.grid->vertex[vertexInd];
            VECTOR3 F = vertexData[vertexInd].F_sum;
            double h = vertexData[vertexInd].h;
            double stiffness = vertexData[vertexInd].stiffness;
            //if(sumF != 0)
            fprintf(f_F, "%lf\t%le\n", calc_r(p0), F.abs());
            fprintf(f_stiffness, "%lf\t%le\n", calc_r(p0), stiffness);
            if(F.abs() != 0)
                fprintf(f_h, "%lf\t%le\n", calc_r(p0), h);
            //VECTOR3 sumF = (*out->mechOut)[globalStepNumber].vertex[vertexInd].sumF_vector;
            //fprintf(f_F, "%le %le\n", (double)i0, sumF.abs());
        }

    }
    // расчёт поверхностных сил, эквивалентных узловым контактным силам
    std::vector<VECTOR3> vertexForce(task.mechTask.grid->vertex.size());    // силы для каждого узла сетки
    // заполнение вектора контактных сил
    VECTOR3 Fn_sum = VECTOR3_NULL;
    for(size_t vertexInd = 0; vertexInd < vertexData.size(); vertexInd++)
    {
        bool contacted = vertexData[vertexInd].contact;
        if(contacted)
        {
            vertexForce[vertexInd] = vertexData[vertexInd].F_sum;
//vertexForce[vertexInd][0] = 0;// искусственно зануляем силы по x
//vertexForce[vertexInd][1] = 0;// искусственно зануляем силы по y
            Fn_sum += vertexData[vertexInd].F_sum;
        }
        else
        {
            vertexForce[vertexInd] = VECTOR3_NULL;
        }
    }
    // расчёт
    //std::vector<VECTOR3> vertexP; // силы для каждого узла сетки(решение СЛАУ)
    std::vector<VECTOR3> vertexP;   // давления в узлах
    vertexP.clear();
    contactSurface.grid = task.mechTask.grid;
    //fprintf(stderr, "solvePressureByNodeForces..\n");
    {
        // площадь зоны контакта
        double S_true = 0;          // учёт количества контактных узлов на граничных ячейках
        double S_geom = 0;          // учёт пересечения жёсткой поверхности с граничными ячейками
        double S_shamanstvo = 0;    // учёт величины давления на граничных ячейках
        contactSurface.grid = task.mechTask.grid;
        Post::calcContactArea(contactSurface, *task.mechTask.grid, vertexData,
                               vertexP, S_true, S_geom, S_shamanstvo);
        S_true *= 4;
        S_geom *= 4;
        S_shamanstvo *= 4;
        //contactSurface.solvePressureByNodeForces(vertexForce,
        //                                         vertexP);
    }
    //fprintf(stderr, "solvePressureByNodeForces finished\n");
    // вывод результата
    // vertexP[vertexInd] - давление в узле vertexInd с координатами task.grid->vertex[vertexInd]
    double Pmax_numb = -1.e10;

    // подсчёт радиуса области контакта
    double r_summ = 0;
    int r_count = 0;
    for(size_t faceInd = 0; faceInd < contactSurface.face.size(); faceInd++)
    {
        const Grid::FEFace &faceEl = contactSurface.face[faceInd];
        const Grid::FE_base *feEl = task.mechTask.grid->fe[faceEl.feIndex];
        int vi[4];
        feEl->getFaceVertexIndexes(faceEl.faceIndex, vi);
        // подсчёт среднего расстояния до пограничных контактных узлов
        int contactPointsCount = 0;
        for(int i = 0; i < 4; i++)
        {
            if(vertexData[vi[i]].contact)
                contactPointsCount++;
        }
        // граница?
        if(contactPointsCount >= 1 && contactPointsCount <= 3)
        {
            for(int i = 0; i < 4; i++)
            {
                if(vertexData[vi[i]].contact)
                {
                    const POINT3 p = task.grid->vertex[vi[i]];
                    r_summ += calc_r(p);
                    r_count++;
                }
            }
        }
        // главные напряжения в приповерхностных КЭ
        MechOutFePointData &fePointData = (*out.mechOut)[globalStepNumber].fe[faceEl.feIndex].pd[0];
        VECTOR3 ms;
        int ind[3];
        // 1) главные напряжения
        //fePointData.mainStresses(ms);
        //sortms_maxToMin(ms, ind);
        // 2) напряжения на площадках
        //VECTOR3 norm = POINT3(0, ccp.y0, 0) - fePointData.p;    // направление к центру
        VECTOR3 norm = POINT3(0, 1, 0);    // направление y
        norm = norm / norm.abs();
        MATR3x3 sigma3x3;
        VECTOR3 vectSigma;
        fePointData.sumSigma.ToMATR3x3_sigma(sigma3x3); // тензор напряжений
        {
            vectSigma = sigma3x3 * norm;// ms - вектор напряжений на площадке
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
        ind[0] = 0;
        ind[1] = 1;
        ind[2] = 2;
        // вывод в файл
        double r = calc_r(fePointData.p);
        fprintf(f_sigma1, "%le\t%le\n", r, fabs(ms[ind[0]]));
        fprintf(f_sigma2, "%le\t%le\n", r, fabs(ms[ind[1]]));
        fprintf(f_sigma3, "%le\t%le\n", r, fabs(ms[ind[2]]));
        // расчёт максимального главного напряжения в приповерхностных КЭ
        if(fabs(ms[ind[0]]) > Pmax_numb)
            Pmax_numb = fabs(ms[ind[0]]);
    }
    // давление в контактных узлах
    for(size_t vertexInd = 0; vertexInd < vertexData.size(); vertexInd++)
    {
        VECTOR3 P_numb;
        P_numb = vertexP[vertexInd];
        //P_numb[0] = res[vertexInd*3 + 0];
        //P_numb[1] = res[vertexInd*3 + 1];
        //P_numb[2] = res[vertexInd*3 + 2];
        if(P_numb.abs() > 0)
        {
            fprintf(f_P, "%le\t%le\n", calc_r(task.grid->vertex[vertexInd]), -P_numb[1]);
            //fprintf(f_P, "%le\t%le\n", calc_r(task.grid->vertex[vertexInd]), P_numb.abs());
            //fprintf(f_P, "%le\t%le\n", task.grid->vertex[vertexInd][0], P_numb.abs());
            //fprintf(f_P, "%le\t%le\n", task.grid->vertex[vertexInd][0], -P_numb[1]);
            fprintf(f_Pxyz, "%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                    task.grid->vertex[vertexInd][0],
                    task.grid->vertex[vertexInd][1],
                    task.grid->vertex[vertexInd][2],
                    P_numb[0], //task.grid->vertex[vertexInd][0],
                    P_numb[1],
                    P_numb[2],
                    P_numb.abs());
            // расчёт максимального давления
            //if(P_numb.abs() > Pmax_numb)
            //    Pmax_numb = P_numb.abs();
        }
    }

    {
        double time = (*out.mechOut)[globalStepNumber].step[0].t0;
        double d = fabs(ccp.y0 - ccp.y_fun->fun(time));                       // сближение цилиндра и полупространства
        double R = ccp.R0;
        double L = ccp.L;       // длина цилиндра
        double E = ccp.E / (1 - SQR(ccp.Nu));   // приведённый модуль упругости

        double P_numb;
        double Fn_numb;
        double a_numb;
        a_numb = r_summ / r_count;
        if(ccp.halhGrid)
        {
            Fn_sum[1] *= 2;
            P_numb = fabs(Fn_sum[1] / L);
            Fn_numb = fabs(Fn_sum[1]);

        }
        else
        {
            P_numb = fabs(Fn_sum[1] / L);
            Fn_numb = fabs(Fn_sum[1]);
        }
        //Pmax_numb расчитано выше
        fprintf(f_Fn, "%le\t%le\n", d, Fn_numb);
        fprintf(f_a, "%le\t%le\n", d, a_numb);
        fprintf(f_Pmax, "%le\t%le\n", d, Pmax_numb);

        // аналитические решения
        if(globalStepNumber == 0)
        {
            FILE *f_Fn_analit;
            FILE *f_a_analit;
            FILE *f_Pmax_analit;
            f_Fn_analit = fopen(fn_Fn_analit.c_str(), "w");
            f_a_analit = fopen(fn_a_analit.c_str(), "w");
            f_Pmax_analit = fopen(fn_Pmax_analit.c_str(), "w");
            int xN = 1000;
            for(int xi = 0; xi <= xN; xi++)
            {
                double d = ccp.d*xi/xN;
                double P_analit = calc_P(d, ccp);
                double Fn_analit = P_analit*L;                        // нормальная сила контактного взаимодействия
                double a_analit = sqrt(4*R/PI/E*P_analit);            // радиус области контакта
                double Pmax_analit = 2*P_analit/PI/a_analit;//sqrt(E/PI/R*P);
                fprintf(f_Fn_analit, "%le\t%le\n", d, Fn_analit);
                fprintf(f_a_analit, "%le\t%le\n", d, a_analit);
                fprintf(f_Pmax_analit, "%le\t%le\n", d, Pmax_analit);
            }
            fclose(f_Fn_analit);
            fclose(f_a_analit);
            fclose(f_Pmax_analit);
        }

        double P_analit = calc_P(d, ccp);
        double Fn_analit = P_analit*L;                        // нормальная сила контактного взаимодействия
        double a_analit = sqrt(4*R/PI/E*P_analit);            // радиус области контакта
        double Pmax_analit = 2*P_analit/PI/a_analit;//sqrt(E/PI/R*P);
        int xN = 1000;
        for(int xi = 0; xi <= xN; xi++)
        {
            double x1 = 0;
            double x2 = a_analit;
            double xc = 0;//(ccp.p1[0] + ccp.p2[0])/2;
            double x = x1 + (x2 - x1) * xi / xN;
            fprintf(f_P_analit_d, "%le\t%le\n", x, Pmax_analit*sqrt(1 - SQR((x-xc)/a_analit)));
        }
        fprintf(f_inf, "d = %le\n", d);

        fprintf(f_inf, "L = %le\n", L);
        fprintf(f_inf, "a = %le\n", a_analit);
        fprintf(f_inf, "Pmax = %le\n", Pmax_analit);
        fprintf(f_inf, "Fn = %le\n", Fn_analit);
        fprintf(f_inf, "\nFn_sum = %le, %le, %le\n", Fn_sum[0], Fn_sum[1], Fn_sum[2]);
        fprintf(f_inf, "\n");
    }
    fclose(f_Fn);
    fclose(f_a);
    fclose(f_Pmax);

    fclose(f_inf);
    fclose(f_P);
    fclose(f_P_analit_d);
    fclose(f_P_analit_F);
    fclose(f_h);
    fclose(f_Pxyz);
    fclose(f_F);
    fclose(f_sigma1);
    fclose(f_sigma2);
    fclose(f_sigma3);
    fclose(f_stiffness);
    // построение графиков
    makeGraphsAfterGlobalStep(globalStepNumber, task);
}
bool Test_contactCylinder::possibleToShow2d() const
{
    return true;
}
bool Test_contactCylinder::getContactSurfaceCircle(const Task &task, const int globalStepIndex, POINT2 &c, double &R) const
{
    if(ccp.surfaceType != Grid::SurfaceType::none)
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
void Test_contactCylinder::needToDrawFe(const Solid::Task &, const OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const
{
    if(!(feInd%ccp.Nz == (size_t)z_index))
    {
        faceState = {0, 0, 0, 0, 0, 0};
    }
    else
    {
        //faceState = {0, 0, 1, 1, 0, 0};
        faceState = {1, 1, 1, 1, 1, 1};
    }
}
bool Test_contactCylinder::needToPaintFiniteElementSurface() const
{
    return ccp.surfaceType == Grid::SurfaceType::FiniteElementSurface;
}
}


namespace Grid
{
void Grid3D::genContactCylinder(const ContactCylinderParameters &ccp)
{
    //#define MOVE(v_ind, d0, d1, d2)	((v_ind + offset0*d0 + offset1*d1 + offset2*d2))
    using namespace Operations;

    build(ccp.av, ccp.ah);

    // параметры индексации пространства
    {
        regionIndexationParameters.h_min = 2;
        regionIndexationParameters.h_max = 6;
        //Sphere s0;
        //s0.O = {0, 0, 0};
        //s0.R = 20;//MAX(abs(ccp.p2[0] - ccp.p1[0])/2., abs(ccp.p2[1] - ccp.p1[1])/2.);
        //regionIndexationParameters.q0.initBySphere(s0);
        regionIndexationParameters.q0.i[0][0] = -11;
        regionIndexationParameters.q0.i[0][1] = +11;
        regionIndexationParameters.q0.i[1][0] = -21;
        regionIndexationParameters.q0.i[1][1] = +10;
        regionIndexationParameters.q0.i[2][0] = -2;
        regionIndexationParameters.q0.i[2][1] = +2;
        regionIndexationParameters.div = 1;
    }

    GridRectangleRegular2D it_grid;
    int N[2] = {80, 40};
    if(ccp.surfaceType == SurfaceType::AnaliticalSurface_Cylinder)    // аналитическая поверхность
    {
        it_grid.init(0 - PI/2, PI + PI/2, -22, 2, N[0], N[1]);
    }

    // поверхность заданная аналитически
    analiticalSurface.resize(1);
    if(ccp.surfaceType == SurfaceType::AnaliticalSurface_Cylinder)
    {
        Interpolation::AnaliticalSurface_Cylinder *cylinder =
                new Interpolation::AnaliticalSurface_Cylinder(
                    ccp.decompositionType,
                    it_grid,
                    ccp.R_fun,
                    ccp.x_fun,
                    ccp.y_fun,
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
                    ccp.R_fun,
                    ccp.x_fun,
                    ccp.y_fun,
                    0, //0
                    0, //0
                    0, //0
                    0, //0
                    1.e-8);
        analiticalSurface[0] = cylinder;
    }

    // поверхность заданная интерполянтом
    ISurface.resize(0);

    {
        std::vector<int> contactFESurfaceIndex;
        contactFESurfaceIndex.push_back(0); // КЭ поверхность с индексом 0 - контактсная
        //cuthillMcKee(contactFESurfaceIndex);
        nestedDissection(contactFESurfaceIndex);
    }
}

}


/*fe.clear();
vertex.clear();
vertexForCurvature.clear();
bc1.clear();
FESurface.clear();
analiticalSurface.clear();
ISurface.clear();

CC_vertexMap m;
{
    CC_vertexIndex i;
    for(i.areaIndex = 0; i.areaIndex < 5; i.areaIndex++)
    {
        POINT2 a[4];    // область
        a[0] = ccp.a[i.areaIndex][0];
        a[1] = ccp.a[i.areaIndex][1];
        a[2] = ccp.a[i.areaIndex][2];
        a[3] = ccp.a[i.areaIndex][3];
        VECTOR3_uint N = ccp.N[i.areaIndex];

        for (i.i[1] = 0; i.i[1] <= N[1]; i.i[1]++)            // y
            for (i.i[0] = 0; i.i[0] <= N[0]; i.i[0]++)        // x
                for(i.i[2] = 0; i.i[2] <= N[2]; i.i[2]++)     // z
                {
                    int vIndex = find_CC_vertex(ccp, i, m);;
                    POINT3 v;
                    if(vIndex == -1)
                    {
                        // новая вершина
                        // расчёт координат вершины
                        i.getVertex(ccp, v);
                        // добавление вершины в массив вершин
                        vertex.push_back(v);
                        vIndex = (int)vertex.size() - 1;
                        // добавление вершины в хеш таблицу вершин
                        std::pair<const CC_vertexIndex, int> m_value = {i, vIndex};
                        m.insert(m_value);
                        // фиксируем по x посередине
                        if(i.areaIndex == 0 ||
                           i.areaIndex == 1 ||
                           i.areaIndex == 4)
                        {
                            if(i.i[0] == N[0]/2)
                            {
                                BoundaryCondition1 bc1_el;
                                bc1_el.bc1SourceIndex = 0;
                                bc1_el.vertexIndex = vIndex;
                                bc1.push_back(bc1_el);
                            }
                        }
                        // фиксируем по y снизу
                        if(i.areaIndex == 0)
                        {
                            if(i.i[1] == 0)
                            {
                                BoundaryCondition1 bc1_el;
                                bc1_el.bc1SourceIndex = 1;
                                bc1_el.vertexIndex = vIndex;
                                bc1.push_back(bc1_el);
                            }
                        }
                        // фиксируем по z везде
                        {
                            BoundaryCondition1 bc1_el;
                            bc1_el.bc1SourceIndex = 2;
                            bc1_el.vertexIndex = vIndex;
                            bc1.push_back(bc1_el);
                        }
                    }
                    else
                    {
                        // информация о вершине с индексом i найдена
                        v = vertex[vIndex];
                    }
                }
    }

    // построение конечных элементов и поверхностей
    FESurface.resize(1);
    for (size_t FEsurfaceInd = 0; FEsurfaceInd < FESurface.size(); FEsurfaceInd++)
        FESurface[FEsurfaceInd].face.clear();
    for(i.areaIndex = 0; i.areaIndex < 5; i.areaIndex++)
    {
        POINT2 a[4];    // область
        a[0] = ccp.a[i.areaIndex][0];
        a[1] = ccp.a[i.areaIndex][1];
        a[2] = ccp.a[i.areaIndex][2];
        a[3] = ccp.a[i.areaIndex][3];
        VECTOR3_uint N = ccp.N[i.areaIndex];


        for (i.i[1] = 0; i.i[1] < N[1]; i.i[1]++)            // y
            for (i.i[0] = 0; i.i[0] < N[0]; i.i[0]++)        // x
                for(i.i[2] = 0; i.i[2] < N[2]; i.i[2]++)     // z
                {
                    // индексы вершин 6-гранника
                    int vi[8];
                    VECTOR3_uint di;
                    int vi_index = 0;
                    for(di[2] = 0; di[2] <= 1; di[2]++)
                    for(di[1] = 0; di[1] <= 1; di[1]++)
                    for(di[0] = 0; di[0] <= 1; di[0]++)
                    {
                        CC_vertexIndex fe_i = i;
                        fe_i.i = fe_i.i + di;
                        int vIndex = find_CC_vertex(ccp, fe_i, m);;
                        // информация о вершине с индексом fe_i должна существовать в хеш таблице m
                        if(vIndex == -1)
                        {
                            for(;;);
                        }
                        //POINT3 v = vertex[vIndex];
                        vi[vi_index] = vIndex;
                        vi_index++;
                    }
                    // создание шестигранного КЭ
                    FE_LinearHexagon *fe_el = new FE_LinearHexagon;   // шестигранник, линейное отображение
                    fe_el->mi = 0;                                    // материал с индексом 0
                    for (int t = 0; t < 8; t++)
                        fe_el->vi[t] = vi[t];
                    fe.push_back(fe_el);
                    // поверхность для контакта сверху (активная)
                    if(i.areaIndex == 4)
                    {
                        if(i.i[1] == N[1] - 1)
                        {
                            Grid::FEFace face;
                            face.feIndex = (int)fe.size() - 1;
                            face.faceIndex = 5; // Y = +1
                            FESurface[0].face.push_back(face);
                        }
                    }

                }

    }
}*/
//for (int globalStepNumber = 0; globalStepNumber < (int)task.step->size(); globalStepNumber++)
//int globalStepNumber = 1;  // нагрузили
//int globalStepNumber = 2;  // +охладили
/*
int FeInd = 0;  // индекс конечного элемента
if(task.mechTask.enabled)
for (i0 = 0; i0 < ccp.N[0]; i0++)            // x
    for (i1 = 0; i1 < ccp.N[1]; i1++)        // y
        for (i2 = 0; i2 < ccp.N[2]; i2++)    // z
        {
            MechOutFePointData &r = (*out->mechOut)[globalStepNumber].fe[FeInd].pd[0];
            VECTOR3 ms;
            r.mainStressesXY(ms);
            double sigma1 = ms[0];
            double sigma2 = ms[1];
            double sigma3 = ms[2];
            // главные напряжения
            //VECTOR3 ms;
            //r.sumSigma.solveMain(ms);
            //double sigma1 = ms[0];
            //double sigma2 = ms[1];
            //double sigma3 = ms[2];
            if(i1 == 0)
            {
                fprintf(f_sigmaBottom1, "%llu\t%le\n", i0, sigma1);
                fprintf(f_sigmaBottom2, "%llu\t%le\n", i0, sigma2);
                fprintf(f_sigmaBottom3, "%llu\t%le\n", i0, sigma3);
            }
            if(i1 == ccp.N[1] - 1)
            {
                fprintf(f_sigmaTop1, "%llu\t%le\n", i0, sigma1);
                fprintf(f_sigmaTop2, "%llu\t%le\n", i0, sigma2);
                fprintf(f_sigmaTop3, "%llu\t%le\n", i0, sigma3);
            }
            if(i1 == ccp.N[1]/2)
            {
                fprintf(f_sigmaMiddle1, "%llu\t%le\n", i0, sigma1);
                fprintf(f_sigmaMiddle2, "%llu\t%le\n", i0, sigma2);
                fprintf(f_sigmaMiddle3, "%llu\t%le\n", i0, sigma3);
            }
            FeInd++;
        }
*/
// подсчёт среднего расстояния до центров пограничных граней с контактными узлами
/*
POINT3 c = VECTOR3_NULL;       // центр грани
int contactPointsCount = 0;
for(int i = 0; i < 4; i++)
{
    c += task.grid->vertex[vi[i]];
    if(vertexData[vi[i]].contact)
        contactPointsCount++;
}
// граница?
if(contactPointsCount >= 1 && contactPointsCount <= 3)
{
    c /= 4;
    r_summ += fabs(c[0]);// считаем просто координаты x а не длину дуги
    r_count++;
}
*/
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
    for(int i = 0; i < 4; i++)
    {
        c += task.grid->vertex[vi[i]];
        P_numb += vertexP[vi[i]];
    }
    c /= 4;
    P_numb /= 4;
    if(P_numb.abs() > 0)
    {
        if(P_numb.abs() > Pmax_numb)
            Pmax_numb = P_numb.abs();
        fprintf(f_P, "%le\t%le\n", c[0], P_numb.abs());
        fprintf(f_Pxyz, "%le\t%le\t%le\t%le\n",
                c[0],
                c[1],
                c[2],
                P_numb.abs());
    }
}
*/
