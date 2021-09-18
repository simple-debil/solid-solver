#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"
#include <fstream>

#include "tests.h"
#include "console.h"

#include "interpolation.h"

using namespace Elementary;
using namespace Solid;
using namespace Grid;

#define MOVE(v_ind, d0, d1, d2)	((v_ind + offset0*d0 + offset1*d1 + offset2*d2))





namespace Tests
{
Test_crack_plate_gen::Test_crack_plate_gen()
{
    // дирректория данных теста
    dir = "./tests/crack_plate_gen/";
    // имена gnuplot файлов(+пути к файлам в опциях) и соответствующих картинок
    fn_filesList = dir + "files.txt";
    resGraph.read_gnuplot_from_file(fn_filesList);
    // инициализация графиков с невязками
    initStepsGraphs();

}
Test_base::Type Test_crack_plate_gen::get_type() const
{
    return Type::Crack_plate_gen;
}
void Test_crack_plate_gen::initTask(Solid::Task &task)
{
    gp.grid_mode = 1;
    gp.p1 = VECTOR3(-4, -4, -1);
    gp.p2 = VECTOR3(+4, +4, +1);
    gp.matetial_surface_y[0] = -0.99;
    gp.matetial_surface_y[1] = +0.99;
    gp.NU = 0.4;
    gp.E[2] = 1.e6;//1.e6;
    gp.E[1] = 2.e6;//2.e6;
    gp.E[0] = 3.e6;//3.e6;
    gp.fix_z = true; // ###не используется
    double PP = -5.e3;
    gp.Px = PP*1.e-100;
    gp.P_in_left = -PP*0;
    gp.P_in_right = -PP*1;
    int NN = 80*1;//80;//32;
    int fixGrid = 0;
    gp.fem_multipler = 4;
    gp.N = VECTOR3_uint(NN+1, NN, 2);//8
    gp.crack_y[0] = -0.85;//-0.851;
    gp.crack_y[1] = +0.85;//+0.851;
     gp.crack_x_local = -0.99999999;
     gp.crack_x_global = gp.crack_x_local * (gp.p2[0] - gp.p1[0])/gp.N[0]/2.;
     gp.XFEM_N[0] = 8;
     gp.XFEM_N[1] = 8;
     gp.XFEM_crack_Nx_left = 1;//gp.XFEM_N[0]/2 - 3;
     gp.XFEM_crack_Nx_right = gp.XFEM_N[0] - gp.XFEM_crack_Nx_left;
    gp.minSigmaEqv = -1.e100;
    gp.maxSigmaEqv = 1.6e4;//+1.e100;
    gp.vertex_tip_add = true;
    gp.tip_elastic = true;
    gp.tip_elastoplastic = !gp.tip_elastic;
    gp.pl_hardening_n = 10;
    //gp.maxSigmaEqv = 3.7e4;//+1.e100;
    //0.85: 2.8e4;
    //1.0: 3.7e4
    //bc2 0.85: 1.6e4
    //gp.crack_y[0] = -0.875;//-0.851;
    //gp.crack_y[1] = +.55;//+0.851;


    // подрубка
    task.mechTask.enabled = true;
    task.thermTask.enabled = false;
    // общие шаги и сетка
    std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
    Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
    if(gp.grid_mode == 1)
        grid->genCrack_plate_gen_XFEM(gp);
    else
        grid->genCrack_plate_gen_FEM(gp);

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
    if(gp.tip_elastoplastic && gp.grid_mode == 0)
    {
        StepsNumber = 10;
    }
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
        material = new std::vector<MechMaterialSource>(3);
        for(int i = 0; i < 3; i++)
        {
            MechMaterialSource mechMat;
            mechMat.elasticSigmaLimit = 1.e10;// ни на что не влияет
            mechMat.Ceps = 4./9.;
            mechMat.Csigma = 1;
            mechMat.elasticParameters0.Talpha = 0;
            mechMat.set_M_sigma();
            mechMat.F = {0,0,0};
            mechMat.elasticParameters0.ro = 0;
            mechMat.elasticParameters0.Talpha = 0*1.e-5; // Коэффициент линейного расширения
            mechMat.w_midPoint = w_midPoint;
            mechMat.w_project = w_project;
            mechMat.plasticityMethodType = plasticityMethodType;
            mechMat.PCUnloadingType = PCUnloadingType;
            mechMat.set_E_NU(gp.E[i], gp.NU);
            if(gp.tip_elastoplastic && gp.grid_mode == 0)
            {
                mechMat.elasticSigmaLimit = 1.e4;//###
                setPlasticMaterialCurve_Yeld_hardening(1, mechMat, 1/gp.pl_hardening_n);//##1/n??
                mechMat.plasticityMethodType = MechPlasticityMethodType::Combo_D_pl_InitialSigma;
                mechMat.PCUnloadingType = PCUnloadingType;
                mechMat.w_midPoint = w_midPoint;//1.0;//0.5//0
                mechMat.w_project = w_project;
                mechMat.cosTettaMin = 0.1;
            }
            (*material)[i] = mechMat;
        }
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



        // шаги (на каждом шаге задаются вторые краевые условия)
        std::vector<MechGlobalStepParameters> *mechStep = new std::vector<MechGlobalStepParameters>;
        // нагружение
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
        // первый шаг (1 временной слой, заведомо упругое нагружение Px)
        s = (*step)[0];
        bc2Source = new std::vector<MechBoundaryCondition2Source_base *>(4);
        setLinearBc2_el(0, gp.Px*0, s.t_start, s.t_finish,
                        (*bc2Source)[0]);
        setLinearBc2_el(0, gp.Px, s.t_start, s.t_finish,
                        (*bc2Source)[1]);
        setLinearBc2_el(0, gp.P_in_left, s.t_start, s.t_finish,
                        (*bc2Source)[2]);
        setLinearBc2_el(0, gp.P_in_right, s.t_start, s.t_finish,
                        (*bc2Source)[3]);
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
void Test_crack_plate_gen::writeResults(const Task &task, const Solid::OutData &out)
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
    // пути к файлам для отображения графиков
{
    std::ofstream f;
    if(globalStepNumber == 0)
    {
        f.open(fn_filesList, std::ofstream::out);
        resGraph.generalGraphsNomber = 0;
        resGraph.graphsPerGlobalStepNomber = 1;
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
    f.close();
}
    resGraph.read_gnuplot_from_file(fn_filesList);
    }

    FILE *f_F = fopen(fn_F.c_str(), "w");
    size_t i0, i1, i2;



    if(gp.grid_mode)
    {
        // centr_ind - индексы КЭ
        for(size_t centr_ind = 0; centr_ind < gp.centr_ind.size(); centr_ind++)
        {
            size_t feInd = gp.centr_ind[centr_ind];
            const Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<const Grid::FE_LinearHexagon_XFEM *>(task.grid->fe[feInd]);

            // вариант 1
            /*
            MechOutFeData &r = (*out.mechOut)[globalStepNumber].fe[feInd];
            POINT3_CUBE XYZ;
            double bas[16];
            // перемещение правой части трещины (левая часть правой подобласти)
            XYZ = {-1, 0, 0.5};
            feEl_XFEM->calcSubBasFuncs(1, XYZ, bas);
            double ux = 0;
            for(int i = 0; i < 16; i++)
            {
                ux += bas[i] * r.q[i][0];
            }
            fprintf(f_F, "%le %le\n", r.pd[1].p[1], ux);
            */

            // вариант 2
            /*
            // копии глобальных индексов и координат вершин 6-гранника
            POINT3 v[8];
            feEl_XFEM->getGeomVertexes(task.grid->vertex, task.grid->vertexForCurvature, v);
            // вершины подобластей
            POINT3 subv[8];
            feEl_XFEM->gn_getSubVertexes(1, subv);
            fprintf(f_F, "%le %le\n", subv[0][1], subv[0][0]);
            */

            // вариант 3
            //feEl_XFEM->gn_getSubVertexes(1, subv);
            //fprintf(f_F, "%le %le\n", subv[0][1], subv[0][0]);
            for(std::set<size_t>::iterator it = gp.subVertex_ind[centr_ind].begin(); it != gp.subVertex_ind[centr_ind].end(); ++it)
            {
                double y = feEl_XFEM->subVertex[*it][1];
                double dx = feEl_XFEM->dsubVertex[*it][0];
                //if(y >= gp.crack_y[0] && y <= gp.crack_y[1])
                    fprintf(f_F, "%le %le\n", y, dx);
            }
        }
    }
    else
    {
        // centr_ind - индексы вершин
        for(size_t centr_ind = 0; centr_ind < gp.centr_ind.size(); centr_ind++)
        {
            double ux = (*out.mechOut)[globalStepNumber].vertex[gp.centr_ind[centr_ind]].sum_du[0];
            double y = (*out.mechOut)[globalStepNumber].vertex[gp.centr_ind[centr_ind]].p[1];
            fprintf(f_F, "%le %le\n", y, ux);
        }

    }


    /*
    int FeInd = 0;  // индекс конечного элемента
    if(task.mechTask.enabled)
    for (i0 = 0; i0 < gp.N[0]; i0++)            // x
        for (i1 = 0; i1 < gp.N[1]; i1++)        // y
            for (i2 = 0; i2 < gp.N[2]; i2++)    // z
            {
                MechOutFePointData &r = (*out.mechOut)[globalStepNumber].fe[FeInd].pd[0];
                FeInd++;
            }*/
    fclose(f_F);
    // построение графиков
    makeGraphsAfterGlobalStep(globalStepNumber, task);
}
bool Test_crack_plate_gen::possibleToShow2d() const
{
    return true;
}
bool Test_crack_plate_gen::getContactSurfaceCircle(const Task &task, const int globalStepIndex, POINT2 &c, double &R) const
{
    return false;
}
void Test_crack_plate_gen::needToDrawFe(const Solid::Task &, const OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const
{
    //if(!(feInd%gp.N[2] == (size_t)z_index ||
    //        feInd >= gp.N[0] * gp.N[1] * gp.N[2]))
    if(!(out.mechOut[0][0].fe[feInd].pd[0].p[2] > 0))
    {
        faceState = {0, 0, 0, 0, 0, 0};
    }
    else
    {
        faceState = {1, 0, 0, 0, 0, 0};
    }
}
void Test_crack_plate_gen::needToDrawSubFe(const Task &task, const OutData &out, const size_t feInd, const size_t subHexInd, const int z_index, std::array<bool, 6> &faceState) const
{
    faceState = {0, 1, 0, 0, 0, 0};
    //faceState = {1, 1, 1, 1, 1, 1};
}
bool Test_crack_plate_gen::needToPaintFiniteElementSurface() const
{
    return false;
}
double Test_crack_plate_gen::min_value(const int valueIndex) const
{
    return gp.minSigmaEqv;
}
double Test_crack_plate_gen::max_value(const int valueIndex) const
{
    return gp.maxSigmaEqv;
}
}

namespace Grid
{
void Grid3D::genCrack_plate_gen_XFEM(Grid::CrackPlateGenParameters &gp)
{
    using namespace Operations;
    fe.clear();
    vertex.clear();
    vertexForCurvature.clear();
    bc1.clear();
    FESurface.clear();
    analiticalSurface.clear();
    ISurface.clear();
    gp.centr_ind.clear();
    cs.clear();


    // поверхность трещины
    CrackSurface_Analitical_vertical_x0 *cs1 = new CrackSurface_Analitical_vertical_x0;
    cs1->init(gp.crack_x_global, gp.crack_y[0], gp.crack_y[1], gp.pl_hardening_n);
    //cs1->init(0, -0.85, 0.85);
    cs.push_back(cs1);

    {
        int N0 = gp.N[0];
        int N1 = gp.N[1];
        int N2 = gp.N[2];
        int offset2 = 1;
        int offset1 = N2 + 1;
        int offset0 = (N2 + 1)*(N1 + 1);
        // построение вершин и первых краевых условий
        for (int i0 = 0; i0 <= N0; i0++)            // x
            for (int i1 = 0; i1 <= N1; i1++)        // y
                for (int i2 = 0; i2 <= N2; i2++)    // z
                {
                    POINT3 p;
                    findPointOnTheLine_1d_conc(gp.p1[0], gp.p2[0], N0, 1.- 1./64*0, i0, p[0]);
                    findPointOnTheLine_1d(gp.p1[1], gp.p2[1], N1, 1, i1, p[1]);
                    findPointOnTheLine_1d(gp.p1[2], gp.p2[2], N2, 1, i2, p[2]);
                    vertex.push_back(p);
                    // фиксируем по x слева
                    if(i0 == 0)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 0;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    // фиксируем по x справа
                    /*
                    if(i0 == N0)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 0;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    */
                    // фиксируем по y снизу
                    if(i1 == 0)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 1;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    // фиксируем по z посередине
                    if(i2 == N2/2)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 2;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    // фиксируем по z везде
                    /*
                    if(gp.fix_z)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 2;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }*/
                }

        // 1) определение множества вершин, ближайших к вершинам трещины (Вильямс)
        std::set<size_t> vertex_tip;
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
                    // XFEM
                    if(vertex[vi[0]][0] < 0 && vertex[vi[1]][0] > 0 // КЭ пересекает плоскость x = 0
                            && vertex[vi[2]][1] >= gp.crack_y[0]    // верхняя грань выше нижней вершины трещины
                            && vertex[vi[0]][1] <= gp.crack_y[1])   // нижняя грань ниже верхней вершины трещины
                    {
                        // задание обогащённых вершин
                        bool crack_top = vertex[vi[2]][1] > gp.crack_y[1]; // верхняя грань не пересекает трещину
                        bool crack_bot = vertex[vi[0]][1] < gp.crack_y[0]; // нижняя грань не пересекает трещину
                        bool crack_middle = !crack_top && !crack_bot;
                        if(crack_top)
                        {
                            for (int t = 0; t < 8; t++)
                                vertex_tip.insert(vi[t]);
                        }
                        if(crack_bot)
                        {
                            for (int t = 0; t < 8; t++)
                                vertex_tip.insert(vi[t]);
                        }
                        if(crack_middle)
                        {
                        }
                    }
                }
        // 2) определение множества вершин, ближайших к трещине и не найденных пунктом 1 (Хевисайд)
        std::set<size_t> vertex_strong;
        std::set<size_t> fe_strong;
        size_t fe_count = 0;
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
                    // XFEM
                    if(vertex[vi[0]][0] < 0 && vertex[vi[1]][0] > 0 // КЭ пересекает плоскость x = 0
                            && vertex[vi[2]][1] >= gp.crack_y[0]    // верхняя грань выше нижней вершины трещины
                            && vertex[vi[0]][1] <= gp.crack_y[1])   // нижняя грань ниже верхней вершины трещины
                    {
                        // задание обогащённых вершин
                        bool crack_top = vertex[vi[2]][1] > gp.crack_y[1]; // верхняя грань не пересекает трещину
                        bool crack_bot = vertex[vi[0]][1] < gp.crack_y[0]; // нижняя грань не пересекает трещину
                        bool crack_middle = !crack_top && !crack_bot;
                        if(crack_top)
                        {
                        }
                        if(crack_bot)
                        {
                        }
                        if(crack_middle)
                        {
                            for (int t = 0; t < 8; t++)
                            {
                                if(vertex_tip.count(vi[t]) == 0)
                                    vertex_strong.insert(vi[t]);
                            }
                            fe_strong.insert(fe_count);
                        }
                    }
                    fe_count++;
                }
        if(gp.vertex_tip_add)
        {
            // 3) дополнительный слой вершин, близких к вершинам трещины (Вильямс)
            std::set<size_t> vertex_tip_add;
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
                        bool tip_state = false;
                        for (int t = 0; t < 8; t++)
                        {
                            if(vertex_tip.count(vi[t]) != 0)
                            {
                                tip_state = true;
                                break;
                            }
                        }
                        if(tip_state)
                        {
                            for (int t = 0; t < 8; t++)
                                vertex_tip_add.insert(vi[t]);
                        }
                    }
            // 4) добавление дополнительного слоя
            for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
            {
                if(vertex_tip.count(vertexIndex) == 0 && vertex_tip_add.count(vertexIndex) != 0)
                {
                    vertex_tip.insert(vertexIndex);
                }
            }
            vertex_tip_add.clear();
        }
        //vertex_tip.clear();
        //vertex_strong.clear();

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
        // обогащающие ф-и Хевисайда
        for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
        {
            if(vertex_strong.count(vertexIndex) != 0)
            {
                Grid::FuncID funcID;
                funcID.type = Grid::FuncType::H;
                funcID.crackIndex = 0;
                DOFs->addDOFs(vertexIndex, funcID);
            }
        }
        // обогащающие ф-и для вершины трещины (Вильсона)
        //vertex_tip.clear();
        for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
        {
            if(vertex_tip.count(vertexIndex) != 0)
            {
                Grid::FuncID funcID;
                funcID.crackIndex = 0;
                if(gp.tip_elastic)
                {
                    funcID.type = Grid::FuncType::tip1;
                    DOFs->addDOFs(vertexIndex, funcID);
                    funcID.type = Grid::FuncType::tip2;
                    DOFs->addDOFs(vertexIndex, funcID);
                    funcID.type = Grid::FuncType::tip3;
                    DOFs->addDOFs(vertexIndex, funcID);
                    funcID.type = Grid::FuncType::tip4;
                    DOFs->addDOFs(vertexIndex, funcID);
                }
                if(gp.tip_elastoplastic)
                {
                    funcID.type = Grid::FuncType::tip_ep1;
                    DOFs->addDOFs(vertexIndex, funcID);
                    funcID.type = Grid::FuncType::tip_ep2;
                    DOFs->addDOFs(vertexIndex, funcID);
                    funcID.type = Grid::FuncType::tip_ep3;
                    DOFs->addDOFs(vertexIndex, funcID);
                    funcID.type = Grid::FuncType::tip_ep4;
                    DOFs->addDOFs(vertexIndex, funcID);
                    funcID.type = Grid::FuncType::tip_ep5;
                    DOFs->addDOFs(vertexIndex, funcID);
                    funcID.type = Grid::FuncType::tip_ep6;
                    DOFs->addDOFs(vertexIndex, funcID);
                }
            }
        }


        // построение конечных элементов и поверхностей
        // будет 4 поверхности для краевых
        FESurface.resize(4);
        for (size_t FEsurfaceInd = 0; FEsurfaceInd < FESurface.size(); FEsurfaceInd++)
            FESurface[FEsurfaceInd].face.clear();
        gp.subVertex_ind.resize(N0*N1*N2);
        for(size_t i = 0; i < gp.subVertex_ind.size(); i++)
            gp.subVertex_ind.clear();
        fe_count = 0;
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
                    // материал
                    bool material_bot = vertex[vi[2]][1] < gp.matetial_surface_y[0]; // нижняя часть пластины
                    bool material_top = vertex[vi[0]][1] > gp.matetial_surface_y[1]; // верхняя часть пластины
                    bool material_middle = !material_bot && !material_top;
                    int material_index;
                    if(material_bot)
                    {
                        material_index = 0;
                    }
                    if(material_middle)
                    {
                        material_index = 1;
                    }
                    if(material_top)
                    {
                        material_index = 2;
                    }
                    // Проверка наличия вершин с обогащением у КЭ (XFEM?)
                    bool XFEM_state;
                    bool XFEM_strong_state;
                    {
                        int count_strong = 0;
                        int count_tip = 0;
                        for (int t = 0; t < 8; t++)
                        {
                            if(vertex_strong.count(vi[t]) != 0)
                            {
                                count_strong++;
                            }
                            if(vertex_tip.count(vi[t]) != 0)
                            {
                                count_tip++;
                            }
                        }
                        if(count_strong == 8 || count_tip != 0)
                        {
                            XFEM_state = true;
                        }
                        else
                        {
                            XFEM_state = false;
                        }
                    }


                    //if(vertex[vi[0]][0] < 0 && vertex[vi[1]][0] > 0 // КЭ пересекает плоскость x = 0
                    //        && vertex[vi[2]][1] >= gp.crack_y[0]    // верхняя грань выше ниже нижней вершины трещины
                    //        && vertex[vi[0]][1] <= gp.crack_y[1])   // нижняя грань ниже верхней вершины трещины
                    if(XFEM_state)
                    {
                        FE_LinearHexagon_XFEM *fe_el = new FE_LinearHexagon_XFEM;
                        int Nx = 2;//8; // регулярное разбиение по x
                        int Ny = 1;//8; // регулярное разбиение по y
                        int XFEM_crack_Nx_left = 1;
                        int XFEM_crack_Nx_right = 1;
                        double crack_x_local = gp.crack_x_local;
                        fe_el->mi = material_index; // материал с индексом 0
                        for (int t = 0; t < 8; t++)
                            fe_el->vi[t] = vi[t];
                        // тип обогащения
                        int count_strong = 0;
                        int count_tip = 0;
                        for (int t = 0; t < 8; t++)
                        {
                            if(vertex_strong.count(vi[t]) == 0)
                            {
                                //fe_el->enriched_strong[t] = 0;
                            }
                            else
                            {
                                //fe_el->enriched_strong[t] = 1;
                                count_strong++;
                            }
                            if(vertex_tip.count(vi[t]) == 0)
                            {
                                //fe_el->enriched_tip[t] = 0;
                            }
                            else
                            {
                                //fe_el->enriched_tip[t] = 1;
                                count_tip++;
                            }
                        }
                        if(count_tip >= 1)
                        {
                            //fe_el->enrichment_tip = true;
                            //fe_el->enrichment_strong = true;
                            Nx = gp.XFEM_N[0];
                            Ny = gp.XFEM_N[1];
                            if(i0 == N0/2)
                            {
                                XFEM_crack_Nx_left = gp.XFEM_crack_Nx_left;
                                XFEM_crack_Nx_right = gp.XFEM_crack_Nx_right;
                                crack_x_local = gp.crack_x_local;
                            }
                            else
                            {
                                XFEM_crack_Nx_left = Nx/2;
                                XFEM_crack_Nx_right = Nx/2;
                                crack_x_local = 0;
                            }
                            //if(count_strong != 8 && count_strong != 4)
                            if(fe_strong.count(fe_count) == 0)
                            {
                                // Если трещина не проходит через КЭ, то ф-ми Хевисайда не обогащаем
                                for (int t = 0; t < 8; t++)
                                {
                                    //fe_el->enriched_strong[t] = 0;
                                }
                                XFEM_strong_state = false;
                            }
                            else
                            {
                                XFEM_strong_state = true;
                            }
                        }
                        else
                        {
                            if(count_strong == 8)
                            {
                                //fe_el->enrichment_tip = false;
                                //fe_el->enrichment_strong = true;
                            }
                            XFEM_strong_state = true;
                        }



                        /*
                        fe_el->enrichment_strong = true;
                        fe_el->enrichment_tip = false;
                        if(crack_top)
                        {
                            fe_el->enrichment_strong = true;
                            fe_el->enriched_strong[0] = 1;
                            fe_el->enriched_strong[1] = 1;
                            fe_el->enriched_strong[2] = 0;
                            fe_el->enriched_strong[3] = 0;
                            fe_el->enriched_strong[4] = 1;
                            fe_el->enriched_strong[5] = 1;
                            fe_el->enriched_strong[6] = 0;
                            fe_el->enriched_strong[7] = 0;
                            fe_el->enrichment_tip = true;
                            fe_el->enriched_tip[0] = 1;
                            fe_el->enriched_tip[1] = 1;
                            fe_el->enriched_tip[2] = 1;
                            fe_el->enriched_tip[3] = 1;
                            fe_el->enriched_tip[4] = 1;
                            fe_el->enriched_tip[5] = 1;
                            fe_el->enriched_tip[6] = 1;
                            fe_el->enriched_tip[7] = 1;
                            Nx = 8;
                            Ny = 8;
                        }
                        if(crack_bot)
                        {
                            fe_el->enrichment_strong = true;
                            fe_el->enriched_strong[0] = 0;
                            fe_el->enriched_strong[1] = 0;
                            fe_el->enriched_strong[2] = 1;
                            fe_el->enriched_strong[3] = 1;
                            fe_el->enriched_strong[4] = 0;
                            fe_el->enriched_strong[5] = 0;
                            fe_el->enriched_strong[6] = 1;
                            fe_el->enriched_strong[7] = 1;
                            fe_el->enrichment_tip = false;
                        }
                        if(crack_middle)
                        {
                            fe_el->enrichment_strong = true;
                            fe_el->enriched_strong[0] = 1;
                            fe_el->enriched_strong[1] = 1;
                            fe_el->enriched_strong[2] = 1;
                            fe_el->enriched_strong[3] = 1;
                            fe_el->enriched_strong[4] = 1;
                            fe_el->enriched_strong[5] = 1;
                            fe_el->enriched_strong[6] = 1;
                            fe_el->enriched_strong[7] = 1;
                            fe_el->enrichment_tip = false;
                        }*/
                        // значения функции Хевисайда в узлах
                        for (int t = 0; t < 8; t++)
                        {
                            /*
                            if(t%2 == 0)
                                fe_el->H_m[t] = -1;
                            else
                                fe_el->H_m[t] = +1;
                            */
                        }
                        // разбиение 6-гранника на 6-гранные подобласти
                        std::vector<POINT3_CUBE> subVertex_XYZ;
                         subVertex_XYZ.resize((Nx+2)*(Ny+1)*2);
                         fe_el->subHex.resize(Nx*Ny);
                        fe_el->subVertex.resize((Nx+2)*(Ny+1)*2);
                        fe_el->dsubVertex.resize((Nx+2)*(Ny+1)*2);
                        fe_el->subVertexData.resize((Nx+2)*(Ny+1)*2);
                         fe_el->funcsID.clear();
                        {
                            #define vi_left(ix, iy, iz)  ((((XFEM_crack_Nx_left+1)*(iy) + (ix))*2 + iz))
                            #define vi_right(ix, iy, iz) (((XFEM_crack_Nx_left+1)*(Ny+1)*2 + ((XFEM_crack_Nx_right+1)*(iy) + (ix))*2 + iz))
                            #define feIndex_left(ix, iy)  ((XFEM_crack_Nx_left)*(iy) + (ix))
                            #define feIndex_right(ix, iy)  ((XFEM_crack_Nx_left)*(Ny) + (XFEM_crack_Nx_right)*(iy) + (ix))

                            double hx_left = (crack_x_local + 1)/XFEM_crack_Nx_left; //2./Nx;
                            double hx_right = (1 - crack_x_local)/XFEM_crack_Nx_right; //2./Nx;
                            double hy = 2./Ny;
                            // локальные координаты вершин подобластей
                            // слева от трещины
                            for(int iy = 0; iy <= Ny; iy++)
                                for(int ix = 0; ix <= XFEM_crack_Nx_left; ix++)
                                {
                                    double x = -1. + hx_left*ix;
                                    double y = -1. + hy*iy;
                                    if(ix == XFEM_crack_Nx_left)
                                        x -= 1.e-10;
                                    subVertex_XYZ[vi_left(ix, iy, 0)] = POINT3_CUBE(x, y, -1);
                                    subVertex_XYZ[vi_left(ix, iy, 1)] = POINT3_CUBE(x, y, +1);
                                }
                            // справа от трещины
                            for(int iy = 0; iy <= Ny; iy++)
                                for(int ix = 0; ix <= XFEM_crack_Nx_right; ix++)
                                {
                                    double x = crack_x_local + hx_right*ix;//-1. + hx*XFEM_crack_Nx_left + hx*ix;
                                    double y = -1. + hy*iy;
                                    if(ix == 0)
                                        x += 1.e-10;
                                    subVertex_XYZ[vi_right(ix, iy, 0)] = POINT3_CUBE(x, y, -1);
                                    subVertex_XYZ[vi_right(ix, iy, 1)] = POINT3_CUBE(x, y, +1);
                                }
                            // индексы вершин подобластей, которые соответствуют вершинам КЭ
                            fe_el->FESubVertexIndex[0] = vi_left(0, 0, 0);
                            fe_el->FESubVertexIndex[1] = vi_right(XFEM_crack_Nx_right, 0, 0);
                            fe_el->FESubVertexIndex[2] = vi_left(0, Ny, 0);
                            fe_el->FESubVertexIndex[3] = vi_right(XFEM_crack_Nx_right, Ny, 0);
                            fe_el->FESubVertexIndex[4] = vi_left(0, 0, 1);
                            fe_el->FESubVertexIndex[5] = vi_right(XFEM_crack_Nx_right, 0, 1);
                            fe_el->FESubVertexIndex[6] = vi_left(0, Ny, 1);
                            fe_el->FESubVertexIndex[7] = vi_right(XFEM_crack_Nx_right, Ny, 1);
                            /*
                            fprintf(stderr, "feIndex=%d\n", (int)fe.size());
                            for(size_t i = 0; i < 8; i++)
                            {
                                POINT3_CUBE XYZ = subVertex_XYZ[fe_el->FESubVertexIndex[i]];
                                fprintf(stderr, "FESubVertexIndex[%d] = (%le, %le, %le)\n", (int)i, XYZ[0], XYZ[1], XYZ[2]);
                            }
                            fprintf(stderr, "\n");
                            for(size_t i = 0; i < subVertex_XYZ.size(); i++)
                            {
                                POINT3_CUBE XYZ = subVertex_XYZ[i];
                                fprintf(stderr, "subVertex_XYZ[%d] = (%le, %le, %le)\n", (int)i, XYZ[0], XYZ[1], XYZ[2]);
                            }
                            fprintf(stderr, "\n");
                            */
                            // индексы вершин 6-гранников
                            // слева от трещины
                            for(int iy = 0; iy < Ny; iy++)
                                for(int ix = 0; ix < XFEM_crack_Nx_left; ix++)
                                {
                                    size_t *vi = fe_el->subHex[feIndex_left(ix, iy)].vi;
                                    vi[0] = vi_left(ix, iy, 0);
                                    vi[1] = vi_left(ix + 1, iy, 0);
                                    vi[2] = vi_left(ix, iy + 1, 0);
                                    vi[3] = vi_left(ix + 1, iy + 1, 0);
                                    vi[4] = vi_left(ix, iy, 1);
                                    vi[5] = vi_left(ix + 1, iy, 1);
                                    vi[6] = vi_left(ix, iy + 1, 1);
                                    vi[7] = vi_left(ix + 1, iy + 1, 1);
                                }
                            // справа от трещины
                            for(int iy = 0; iy < Ny; iy++)
                                for(int ix = 0; ix < XFEM_crack_Nx_right; ix++)
                                {
                                    size_t *vi = fe_el->subHex[feIndex_right(ix, iy)].vi;
                                    vi[0] = vi_right(ix, iy, 0);
                                    vi[1] = vi_right(ix + 1, iy, 0);
                                    vi[2] = vi_right(ix, iy + 1, 0);
                                    vi[3] = vi_right(ix + 1, iy + 1, 0);
                                    vi[4] = vi_right(ix, iy, 1);
                                    vi[5] = vi_right(ix + 1, iy, 1);
                                    vi[6] = vi_right(ix, iy + 1, 1);
                                    vi[7] = vi_right(ix + 1, iy + 1, 1);
                                }
                        }
                        // массив базисных ф-й
                        // считаем, что нет ф-й одинаковых типов, привязанных к разным трещинам
                        for(FuncType ft = (FuncType)0; ft < FuncType::_SIZE; ft = (FuncType)((size_t)ft + 1))
                        {
                            for (size_t t = 0; t < 8; t++)
                            {
                                for (size_t i = 0; i < DOFs->DOFs[vi[t]].size(); i++)
                                {
                                    const FuncID &funcID_el = DOFs->DOFs[vi[t]][i].funcID;
                                    if(funcID_el.type == ft)
                                    {
                                        fe_el->funcsID.push_back(funcID_el);
                                        goto extft;
                                    }
                                }
                            }
                            extft:;
                        }
                        for(size_t i = 0; i < fe_el->subVertexData.size(); i++)
                        {
                            fe_el->subVertexData[i].funcValue.resize(fe_el->funcsID.size());
                        }
                        // инициализация координат вершин шестигранных подобластей
                        fe_el->gn_initSubVertexes(cs, vertex, vertexForCurvature, subVertex_XYZ);
                        // инициализация координат вершин шестигранных подобластей
                        // справа от трещины
                        if(i2 == 1 && i0 == N0/2)
                        //if(XFEM_strong_state)
                        {
                            for(int iy = 0; iy <= Ny; iy++)
                                for(int ix = 0; ix <= 0; ix++)
                                {
                                    double y = fe_el->subVertex[vi_right(ix, iy, 0)][1];
                                    if(y >= gp.crack_y[0] - 1.e-6 &&
                                       y <= gp.crack_y[1] + 1.e-6)
                                        gp.subVertex_ind[gp.centr_ind.size()].insert(vi_right(ix, iy, 0));
                                    //gp.subVertex_ind[gp.centr_ind.size()].insert(vi_right(ix, iy, 1));
                                }
                            gp.centr_ind.push_back(fe.size());
                        }
                        // вторые краевые на поверхности трещины
                        if(i0 == N0/2)
                        {
                            // слева от трещины
                            SubHexFace subHexFace_el;
                            POINT3 v[4];
                            for(int iy = 0; iy < Ny; iy++)
                                for(int ix = XFEM_crack_Nx_left - 1; ix <= XFEM_crack_Nx_left - 1; ix++)
                                {
                                    subHexFace_el.faceIndex = 3; // X = +1
                                    subHexFace_el.hexInd = feIndex_left(ix, iy);
                                    //fprintf(stderr, "%d\n", subHexFace_el.hexInd);
                                    fe_el->gn_getSubFaceVertexes(subHexFace_el.hexInd, subHexFace_el.faceIndex, v);
                                    for(int j = 0; j < 4; j++)
                                    {
                                        double y = v[j][1];
                                        if(y >= gp.crack_y[0] - 1.e-6 &&
                                           y <= gp.crack_y[1] + 1.e-6)
                                        {
                                            fe_el->subHexFace[0].push_back(subHexFace_el);
                                            break;
                                        }
                                    }
                                }
                            // справа от трещины
                            for(int iy = 0; iy < Ny; iy++)
                                for(int ix = 0; ix <= 0; ix++)
                                {
                                    subHexFace_el.faceIndex = 2; // X = -1
                                    subHexFace_el.hexInd = feIndex_right(ix, iy);
                                    fe_el->gn_getSubFaceVertexes(subHexFace_el.hexInd, subHexFace_el.faceIndex, v);
                                    for(int j = 0; j < 4; j++)
                                    {
                                        double y = v[j][1];
                                        if(y >= gp.crack_y[0] - 1.e-6 &&
                                           y <= gp.crack_y[1] + 1.e-6)
                                        {
                                            fe_el->subHexFace[1].push_back(subHexFace_el);
                                            break;
                                        }
                                    }
                                }
                            MechBoundaryCondition2Crack bc2Crack_el;
                            bc2Crack_el.FEIndex = fe.size();
                            bc2Crack_el.subSurfaceIndex = 0;
                            bc2Crack_el.bc2SourceIndex = 2;
                            bc2Crack.push_back(bc2Crack_el);
                            bc2Crack_el.subSurfaceIndex = 1;
                            bc2Crack_el.bc2SourceIndex = 3;
                            bc2Crack.push_back(bc2Crack_el);
                        }
                        fe.push_back(fe_el);
                    }
                    else
                    {
                        FE_LinearHexagon *fe_el = new FE_LinearHexagon; // шестигранник, линейное отображение
                        fe_el->mi = material_index;
                        for (int t = 0; t < 8; t++)
                            fe_el->vi[t] = vi[t];
                        fe.push_back(fe_el);
                    }
                    fe_count++;
                    // поверхности для 2-х краевых условий
                    FEFace face;
                    if (i0 == 0)
                    {
                        face.feIndex = (int)fe.size() - 1;
                        face.faceIndex = 2; // X = -1
                        FESurface[0].face.push_back(face);
                    }
                    if (i0 == N0 - 1)
                    {
                        face.feIndex = (int)fe.size() - 1;
                        face.faceIndex = 3; // X = +1
                        FESurface[1].face.push_back(face);
                    }
                }
    }
}
void Grid3D::genCrack_plate_gen_FEM(Grid::CrackPlateGenParameters &gp)
{
    using namespace Operations;
    fe.clear();
    vertex.clear();
    vertexForCurvature.clear();
    bc1.clear();
    FESurface.clear();
    analiticalSurface.clear();
    ISurface.clear();
    gp.centr_ind.clear();

    {
        int N0 = (gp.N[0]-1)*gp.fem_multipler + 1;
        int N1 = gp.N[1]*gp.fem_multipler;
        int N2 = gp.N[2];
        int offset2 = 1;
        int offset1 = N2 + 1;
        int offset0 = (N2 + 1)*(N1 + 1);
        // построение вершин и первых краевых условий
        for (int i0 = 0; i0 <= N0; i0++)            // x
            for (int i1 = 0; i1 <= N1; i1++)        // y
                for (int i2 = 0; i2 <= N2; i2++)    // z
                {
                    POINT3 p;
                    findPointOnTheLine_1d_conc(gp.p1[0], gp.p2[0], N0, 1.- 1./64*0, i0, p[0]);
                    findPointOnTheLine_1d(gp.p1[1], gp.p2[1], N1, 1, i1, p[1]);
                    findPointOnTheLine_1d(gp.p1[2], gp.p2[2], N2, 1, i2, p[2]);
                    vertex.push_back(p);
                    // фиксируем по x слева
                    if(i0 == 0)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 0;
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
                    if(i2 == N2/2)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 2;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }
                    /*
                    // фиксируем по z везде
                    if(gp.fix_z)
                    {
                        BoundaryCondition1 bc1_el;
                        bc1_el.bc1SourceIndex = 2;
                        bc1_el.vertexIndex = (int)vertex.size() - 1;
                        bc1.push_back(bc1_el);
                    }*/
                }
        // дополнительные вершины-дубли для задания разрыва
        size_t startvind_add = vertex.size();
        for (int i0 = 0; i0 <= N0; i0++)            // x
            for (int i1 = 0; i1 <= N1; i1++)        // y
                for (int i2 = 0; i2 <= N2; i2++)    // z
                {
                    POINT3 p;
                    findPointOnTheLine_1d_conc(gp.p1[0], gp.p2[0], N0, 1.- 1./64*0, i0, p[0]);
                    findPointOnTheLine_1d(gp.p1[1], gp.p2[1], N1, 1, i1, p[1]);
                    findPointOnTheLine_1d(gp.p1[2], gp.p2[2], N2, 1, i2, p[2]);
                    vertex.push_back(p);
                }
        // построение конечных элементов и поверхностей
        // будет 4 поверхности для краевых
        FESurface.resize(4);
        for (size_t FEsurfaceInd = 0; FEsurfaceInd < FESurface.size(); FEsurfaceInd++)
            FESurface[FEsurfaceInd].face.clear();
        for (int i0 = 0; i0 < N0; i0++)           // x
            for (int i1 = 0; i1 < N1; i1++)        // y
                for (int i2 = 0; i2 < N2; i2++)    // z
                {
                    FE_LinearHexagon *fe_el = new FE_LinearHexagon; // шестигранник, линейное отображение
                    bool material_bot; // нижняя часть пластины
                    bool material_top; // верхняя часть пластины
                    bool material_middle;
                    bool crack; // трещина
                    bool crack_top; // верхняя грань не пересекает трещину
                    bool crack_bot; // нижняя грань не пересекает трещину
                    bool crack_middle;

                    int vind0 = i0*offset0 + i1*offset1 + i2*offset2;       // индекс вершины, определяющей шестигранник
                    int vi[8] =
                    {
                        MOVE(vind0, 0, 0, 0), MOVE(vind0, 1, 0, 0), MOVE(vind0, 0, 1, 0), MOVE(vind0, 1, 1, 0),
                        MOVE(vind0, 0, 0, 1), MOVE(vind0, 1, 0, 1), MOVE(vind0, 0, 1, 1), MOVE(vind0, 1, 1, 1)
                    };
                    /*
                    if((int)fe.size() - 1 == 0)
                    {
                        for (int t = 0; t < 8; t++)
                        {
                            fprintf(stderr, "%le %le %le\n", vertex[fe_el->vi[t] = vi[t]][0], vertex[fe_el->vi[t] = vi[t]][1], vertex[fe_el->vi[t] = vi[t]][2]);
                        }
                    }*/
                    material_bot = vertex[vi[2]][1] < gp.matetial_surface_y[0]; // нижняя часть пластины
                    material_top = vertex[vi[0]][1] > gp.matetial_surface_y[1]; // верхняя часть пластины
                    material_middle = !material_bot && !material_top;
                    if(i0 == N0/2
                            && vertex[vi[2]][1] > gp.crack_y[0]+1.e-6    // верхняя грань выше нижней вершины трещины
                            && vertex[vi[0]][1] < gp.crack_y[1]-1.e-6)   // нижняя грань ниже верхней вершины трещины
                    {
                        crack = true;       // КЭ справа от трещины
                        crack_top = vertex[vi[2]][1] > gp.crack_y[1]-1.e-6; // верхняя грань не пересекает трещину
                        crack_bot = vertex[vi[0]][1] < gp.crack_y[0]+1.e-6; // нижняя грань не пересекает трещину
                        crack_middle = !crack_top && !crack_bot;
                        //fprintf(stderr, "FeIndex = %d\n", (int)fe.size());
                    }
                    else
                    {
                        crack = false;
                    }
                    int material_index;
                    if(material_bot)
                    {
                        material_index = 0;
                    }
                    if(material_middle)
                    {
                        material_index = 1;
                    }
                    if(material_top)
                    {
                        material_index = 2;
                    }
                    fe_el->mi = material_index;
                    //crack = false;
                    if(crack)
                    {
                        if(crack_middle)
                        {
                            for (int t = 0; t < 8; t++)
                            {
                                if(t%2 == 0)
                                    fe_el->vi[t] = startvind_add + vi[t];
                                else
                                    fe_el->vi[t] = vi[t];
                            }
                        }
                        if(crack_bot)
                        {
                            for (int t = 0; t < 8; t++)
                            {
                                if(t == 2 || t == 6)
                                    fe_el->vi[t] = startvind_add + vi[t];
                                else
                                    fe_el->vi[t] = vi[t];
                            }
                        }
                        if(crack_top)
                        {
                            for (int t = 0; t < 8; t++)
                            {
                                if(t == 0 || t == 4)
                                    fe_el->vi[t] = startvind_add + vi[t];
                                else
                                    fe_el->vi[t] = vi[t];
                            }
                        }
                        if(i2 == 1)
                        {
                            gp.centr_ind.push_back(fe_el->vi[0]);
                            gp.centr_ind.push_back(fe_el->vi[2]);
                        }
                    }
                    else
                    {
                        // нет разрыва
                        for (int t = 0; t < 8; t++)
                            fe_el->vi[t] = vi[t];
                    }
                    fe.push_back(fe_el);
                    // поверхности для 2-х краевых условий
                    FEFace face;
                    // слева
                    if (i0 == 0)
                    {
                        face.feIndex = (int)fe.size() - 1;
                        face.faceIndex = 2; // X = -1
                        FESurface[0].face.push_back(face);
                    }
                    // справа
                    if (i0 == N0 - 1)
                    {
                        face.feIndex = (int)fe.size() - 1;
                        face.faceIndex = 3; // X = +1
                        FESurface[1].face.push_back(face);
                    }
                    if(vertex[vi[2]][1] > gp.crack_y[0]+1.e-6    // верхняя грань выше нижней вершины трещины
                    && vertex[vi[0]][1] < gp.crack_y[1]-1.e-6)   // нижняя грань ниже верхней вершины трещины
                    {
                        // поверхность трещины слева
                        if (i0 == N0/2 - 1)
                        {
                            face.feIndex = (int)fe.size() - 1;
                            face.faceIndex = 3; // X = +1
                            FESurface[2].face.push_back(face);
                        }
                        // поверхность трещины справа
                        if (i0 == N0/2)
                        {
                            face.feIndex = (int)fe.size() - 1;
                            face.faceIndex = 2; // X = -1
                            FESurface[3].face.push_back(face);
                        }
                    }
                }


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
    /*
    for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
    {
        DOFs->findDOFIndex_Lagr1(vertexIndex);
    }*/
}
}
