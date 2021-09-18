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
Test_sphere_T::Test_sphere_T()
{
    // дирректория данных теста
    dir = "./tests/hollowSphere_T/";
    // имена gnuplot файлов(+пути к файлам в опциях) и соответствующих картинок
    fn_filesList = dir + "files.txt";
    resGraph.read_gnuplot_from_file(fn_filesList);
    // инициализация графиков с невязками
    initStepsGraphs();
}
Test_base::Type Test_sphere_T::get_type() const
{
    return Type::Sphere_T;
}
void Test_sphere_T::initTask(Solid::Task &task)
{
    // индекс теста решателя теплопроводности
    thermTestIndex = 4;             // 1 - t(a)=t1, t(b)=t2
                                    // 2 - t(a)=t1, конвекция на b
                                    // 3 - конвекция на r=a, конвекция на r=b
                                    // 4 - подогрев на r=a, t(b)=t2
    // подрубка
    task.thermTask.enabled = true;
    task.mechTask.enabled = false;
    // общие шаги и сетка
    std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
    Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
    // сетка
    gp.r1 = 1;
    gp.r2 = 4;
    gp.N = 16;   // 4 промежутка на четвертьокружностях    //8, 32
    gp.Nparts = 64;
    gp.q = 1;
    gp.curvilinear = 1;    // отображение (0 - линейное, 1 - квадратичное)
    gp.buildingMethod = 0; // способ построения (0 - делим дуги, 1 - отображение на куб)
    grid->genSphere(gp);

    int StepsNumber = 1;
    int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny

    double Time = 1;
    // шаги
    GlobalStep s;       // шаг
    // первый шаг
    s.t_start = 0;                          // (имеет значение только для первого глобального шага)
    s.t_finish = Time*1;
    s.dt0 = Time/StepsNumber;
    step->push_back(s);
    // второй шаг (пластичное нагружение P2)
    //s.t1 = Time*1;
    //s.t2 = Time*2;
    //s.dt0 = Time/StepsNumber2;
    //step->push_back(s);

    // сетка и шаги общие
    task.grid = grid;
    task.step = step;


    // МДТТ
    if(task.mechTask.enabled)
    {
    }
    else
    {
        using namespace Solid;
        std::vector<MechGlobalStepParameters> *mechStep = new std::vector<MechGlobalStepParameters>;
        mechStep->push_back({});
        task.mechTask.mechStep = mechStep;
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
        // Тест1
        if(thermTestIndex == 1)
            setBcSphereT_1(bc1SourceT, bc2SourceT,
                           1,   // температура на r = a
                           2    // температура на r = b
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
        // шаги
        std::vector<ThermGlobalStep> *thermStep = new std::vector<ThermGlobalStep>;
        ThermGlobalStep ts;  // шаг для теплопроводности
        ts.slausolverParameters = slausolver_parameters;
        ts.timeMode = timeMode;    // 0 - статика, 1 - динамика
        ts.material = material;
        ts.bc1SourceT = bc1SourceT;
        ts.bc2SourceT = bc2SourceT;
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
void Test_sphere_T::writeResults(const Task &task, const Solid::OutData &out)
{
    using namespace Interpolation;
    int globalStepNumber = (int)(*out.thermOut).size() - 2;

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
    fn_T_pogr_abs = subdir + "T_pogr_abs.txt";
    fn_T_pogr = subdir + "T_pogr.txt";
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
    // T
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "T.png") +
         genGnuplotCommandLineParameter("fn_T", subdir + "T.txt") +
         genGnuplotCommandLineParameter("fn_T_analit", subdir + "T_analit.txt") +
         "\" " +
         dir + "T.gnu" << "\n";
    f << subdir + "T.png" << "\n";
    f << "T.png" << "\n";
    // T_pogr
    f << std::string("-e \"") +
         genGnuplotCommandLineParameter("fn_out", subdir + "T_pogr.png") +
         genGnuplotCommandLineParameter("fn_T_pogr", subdir + "T_pogr.txt") +
         "\" " +
         dir + "T_pogr.gnu" << "\n";
    f << subdir + "T_pogr.png" << "\n";
    f << "T_pogr.png" << "\n";
    f.close();
}
    resGraph.read_gnuplot_from_file(fn_filesList);
    }




    FILE *fT_analit = fopen(fn_T_analit.c_str(), "w");
    FILE *fT = fopen(fn_T.c_str(), "w");
    FILE *fT_pogr_abs = fopen(fn_T_pogr_abs.c_str(), "w");
    FILE *fT_pogr = fopen(fn_T_pogr.c_str(), "w");



    double t0;      // момент времени, в который тело разгружено после нагрузки
    double t1;      // момент времени, в который тело нагружено
    int globalStepNumber0;
    int globalStepNumber1;  // соответствующие номера глобальных шагов
    std::vector<GlobalStep> *step = task.step;
    Grid::Grid3D *grid = task.grid;

    if(step->size() == 4)
    {
        globalStepNumber1 = 1;
        globalStepNumber0 = 3;
        t1 = (*step)[globalStepNumber1].t_finish;
        t0 = (*step)[globalStepNumber0].t_finish;
    }
    if(step->size() == 2)
    {
        globalStepNumber1 = 0;
        globalStepNumber0 = 1;
        t1 = (*step)[globalStepNumber1].t_finish;
        t0 = (*step)[globalStepNumber0].t_finish;
    }
    if(step->size() == 1)
    {
        globalStepNumber1 = 0;
        globalStepNumber0 = -123;
        t1 = (*step)[globalStepNumber1].t_finish;
        t0 = -123;
    }

    //for (int globalStepNumber = 0; globalStepNumber < (int)task.step->size(); globalStepNumber++)
    {
        double T_pogr_max = 0;
        double T_pogr_max_abs = 0;
        Heat::ThermGlobalStep &ThermGlobalStep_el = (*task.thermTask.thermStep)[globalStepNumber];

        // вывод данных задачи теплопроводности
        if(task.thermTask.enabled)
        for (int numLocalStep = 0; numLocalStep < (int)(*out.thermOut)[globalStepNumber].step.size(); numLocalStep++)
        {
            using namespace Heat;
            ThermOutStepData &thermOutStep_el = (*out.thermOut)[globalStepNumber].step[numLocalStep];
            POINT3 v1 = grid->vertex[0];
            int N = gp.N/2+1;
            int N_b = (3*N*N-3*N)*(gp.Nparts+1);
            POINT3 v2 = grid->vertex[N_b];
            //POINT3 v2 = grid->vertex[task.sphere.Nparts];
            double a = v1.abs();
            double b = v2.abs();
            if(thermOutStep_el.t == t0 || thermOutStep_el.t == t1)
            {
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
                //double T_b_ch = (*out->thermOut)[globalStepNumber].vertex[task.sphere.Nparts].T;
                double T_b_ch = (*out.thermOut)[globalStepNumber].vertex[N_b].T;
                double T_a_ch = (*out.thermOut)[globalStepNumber].vertex[0].T;
                double Q0 = (*ThermGlobalStep_el.bc2SourceT)[0].q * 4*PI*SQR(a);
                //int resSize = (int)(*out.thermOut)[globalStepNumber].vertex.size();
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
                    fprintf(fT_pogr_abs, "%le %le\n", R, T_pogr_abs*100);
                    }
                }
                // расчет аналитического решения
                int NumPoints = 1000;
                for(int j = 0; j <= NumPoints; j++)
                {
                    double R = a + (b-a)*j/NumPoints;
                    double T;
                    if(thermTestIndex == 1)
                    {
                        T = (a*T1*(b-R)+b*T2*(R-a)) / (R*(b-a));
                    }
                    if(thermTestIndex == 2)
                    {
                        T = (a*T1*(h2*SQR(b) + R*(1-h2*b))+h2*SQR(b)*Ta2*(R-a)) /
                                (R*(h2*SQR(b) + a*(1-h2*b)));
                    }
                    if(thermTestIndex == 3)
                    {
                        T = (Ta1*SQR(a)*h1*(SQR(b)*h2 - R*(b*h2-1)) + Ta2*SQR(b)*h2*(R*(a*h1+1)-SQR(a)*h1)) /
                                (R*(SQR(b)*h2*(a*h1+1)-SQR(a)*h1*(b*h2-1)));
                    }
                    if(thermTestIndex == 4)
                    {
                        T = T_b_ch + Q0/(4*PI*m0.L.m[0][0])*(1./R - 1./b);
                    }
                    if(thermTestIndex == 5)
                    {
                        T = T_a_ch + Q0/(4*PI*m0.L.m[0][0])*(1./R - 1./a);
                    }
                    fprintf(fT_analit, "%le %le\n", R, T);
                }
            }
        }
        //print("\nglobalStepNumber = %d outed\n", ARGS(globalStepNumber));
    }

    fclose(fT_analit);
    fclose(fT);
    fclose(fT_pogr);
    fclose(fT_pogr_abs);
    // построение графиков
    makeGraphsAfterGlobalStep(globalStepNumber, task);
}
}
