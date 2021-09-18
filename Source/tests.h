/* --------------------------------------------------------- */
// построение сеток для тестов
/* --------------------------------------------------------- */

#ifndef TESTS_H
#define TESTS_H

#include "elementary.h"

#include "solid.h"
#include "solving.h"


namespace Tests
{

struct GraphInf
{
    int globalStepNumber;
    std::string fn_graph;
    std::string fn_pic;
    std::string fn_pic_short;
    // построение графика с аргументом fn_graph
    void makeGraph() const;
};

typedef std::vector<GraphInf> ttt;

struct GraphInfs: public ttt
{
    // добавление информации о графике: имена файла gnuplot и файла с картинкой (до полного пути добавляется dir)
    void add_all_in_subdir(const std::string dir, const std::string fn_gnuplot, const std::string fn_picture);
    // добавление информации о графике: имена файла gnuplot и файла с картинкой
    // до полного пути добавляется dir, кроме файлов gnuplot
    void add_gnuplot_in_root(const std::string dir, const std::string fn_gnuplot, const std::string fn_picture);
    // добавление информации о графиках: имена файла gnuplot и файла с картинкой
    // все данные загружаются из файла fn_source
    void read_gnuplot_from_file(const std::string fn_source);
    // построение графиков по именам файлов
    void makeAllGraphs() const;
    // построение графиков для глобального шага globalStepNumber
    void makeglobalStepNumberGraphs(const int globalStepNumber) const;
    // построение графиков, общих для всех глобальных шагов
    void makeGeneralGraphs() const;
    int generalGraphsNomber;      // количество графиков, общих для всех глобальных шагов
    int graphsPerGlobalStepNomber;// количество графиков, которые рисуются на каждом шаге
    // первые generalGraphsNomber элементов - общие графики (для них globalStepNumber = -1)
    // затем для каждого шага по graphsPerGlobalStepNomber графиков
};

struct Test_base: public Solid::StepResultsWriter_base
{
    // типы тестов (ID)
    enum class Type
    {
        Sphere_T = 0,
        Sphere_ep = 1,
        Sphere_Te = 2,
        ContactCylinder = 3,
        ContactSphere = 4,
        Creep = 5,
        Sphere_creep = 6,
        Formovka = 7,
        Kirsch = 8,
        Crack_plate_gen = 9,
    };
    static Test_base *gen(const Type type);
    SlauSolving::SolverParameters slausolver_parameters;
protected:
    // дирректория данных теста
    std::string dir;
    // путь к файлу лога
    std::string fn_log;
    // путь к текстовому файлу с невязками
    std::string fn_steps_all_iterations;
    // путь к текстовому файлу с невязками для графика
    std::string fn_visible_steps_all_iterations;
    // путь к текстовому файлу с данными о построенных графиках
    std::string fn_filesList;
    // пути к gnuplot файлам и полученным картинкам с результатами теста
    GraphInfs resGraph;
    // пути к gnuplot файлам и полученным картинкам с невязками (выводятся данные о всех итерациях)
    GraphInfs stepsGraph_all_iterations;

    // путь к текстовому файлу с невязками (выводятся данные о результатах итераций)
    std::string fn_steps_results_of_iterations;
    // пути к gnuplot файлам и полученным картинкам с невязками (выводятся данные о результатах итераций)
    GraphInfs stepsGraph_results_of_iterations;

    // инициализация stepsGraph
    // все файлы кроме .gnu будут лежать в поддиректориях тестов
    void initStepsGraphs();
    // построение графиков после вывода результата после глобального шага globalStepNumber
    // если глобальный шаг последний, то выводятся все также общие для всех шагов графики
    void makeGraphsAfterGlobalStep(const int globalStepNumber, const Solid::Task &task);
public:
    // построение графиков с результатами
    void makeResultGraphs();
    // получение resGraph
    const GraphInfs *get_oldResultGraphs()const;
    // получение fn_log
    const std::string *get_logFileName();
    // сохранение данных о невязках в файл
    void saveLastStepInf(const Solid::OutDataForLogger &stepsInf);
    // загрузка сохранённых данных о невязках
     void loadStepsInf(std::vector<Solid::MechStepInf> &stepInf);
    // определение номера последнего шага
    int get_stepsInfLastSavedStepNomber();
    // построение данных для графиков с невязками и самих графиков
    void makeStepsGraphs(const int stepIndex1, const int stepIndex2);
    // получение stepsGraph
    void get_oldStepsGraphs(GraphInfs const *&oldStepsGraph_all_iterations, GraphInfs const *&oldStepsGraph_results_of_iterations);
    // получение ID теста
    virtual Type get_type()const = 0;
    // создание входных данных для теста
    virtual void initTask(Solid::Task &task) = 0;
    // вывод результата теста и построение графиков
    //virtual void writeResults(const Solid::Task &task, const Solid::OutData &out) = 0;
    // задача 2д и пригодна для отображения
    virtual bool possibleToShow2d()const;
    // положение поверхности (окружность). Если не задана возвращает false
    virtual bool getContactSurfaceCircle(const Solid::Task &task, const int globalStepIndex, Elementary::POINT2 &c, double &R)const;
    // положение поверхности (окружность)
    //virtual bool circle(const Solid::Task &task, Elementary::POINT2 &c, double &R);
    // индексы координат, в которых рисовать картинку
    virtual void notFixedCoordinates(int &cn1, int &cn2);
    // рисовать ли конечный элемент, возвращает статусы для всех граней
    virtual void needToDrawFe(const Solid::Task &task, const Solid::OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const;
    // рисовать ли подобласть subHexInd конечного элемента, возвращает статусы для всех граней
    virtual void needToDrawSubFe(const Solid::Task &task, const Solid::OutData &out, const size_t feInd, const size_t subHexInd, const int z_index, std::array<bool, 6> &faceState)const;
    // рисовать ли декомпозицию пространства
    virtual bool needToPaintFiniteElementSurface()const;
    virtual double min_value(const int valueIndex)const
    {
        return -1.e100;
    }
    virtual double max_value(const int valueIndex)const
    {
        return +1.e100;
    }
};




// формирует 1 из строковых параметров для помещения его в коммандную строку gnuplot
std::string genGnuplotCommandLineParameter(const std::string &var_name, const std::string &value);

// краевые
void setLinearFunction(double P1, double P2, double t1, double t2, FunParser::Function &fun);
void setLinearBc2_el(double P1, double P2, double t1, double t2, Solid::MechBoundaryCondition2Source_base *&bc2Source_el);
// инициализация кривых для упруго-пластичного материала
void setPlasticCurve(Solid::MechMaterialSource &mechMat, const double k_sigma_eps, const double k_eps_sigma);
// инициализация кривых ползучего материала
void setCreepMaterialCurve(const double A, const double n, const double m, Solid::MechMaterialSource &mechMat);
// инициализация кривых пластичного материала
void setPlasticMaterialCurve(const int plasticPlot, Solid::MechMaterialSource &mechMat, const double k_sigma_eps, const double k_eps_sigma);
// инициализация кривых пластичного материала, sigma(eps_pl)
void setPlasticMaterialCurve_Yeld(const int plasticPlot, Solid::MechMaterialSource &mechMat, const double k_sigma_eps, const double k_eps_sigma);
void setPlasticMaterialCurve_Yeld_hardening(const int plasticPlot, Solid::MechMaterialSource &mechMat, const double n);
// отсутствие 2-х краевых
void setBc2_none(std::vector<Solid::MechBoundaryCondition2Source_base *> *&bc2Source);
// внутреннее давление на полый шар линейной функцией
void setBc2Sphere(std::vector<Solid::MechBoundaryCondition2Source_base *> *&bc2Source, double P1, double P2, double t1, double t2);

void setBcSphereT_1(std::vector<Heat::ThermBoundaryCondition1Source> *&thermBc1, std::vector<Heat::ThermBoundaryCondition2Source> *&thermBc2,
                    double T1, double T2);
void setBcSphereT_2(std::vector<Heat::ThermBoundaryCondition1Source> *&thermBc1, std::vector<Heat::ThermBoundaryCondition2Source> *&thermBc2,
                    double T1, double Ta, double hi);
void setBcSphereT_3(std::vector<Heat::ThermBoundaryCondition1Source> *&thermBc1, std::vector<Heat::ThermBoundaryCondition2Source> *&thermBc2,
                    double Ta1, double hi1, double Ta2, double hi2);
void setBcSphereT_4(std::vector<Heat::ThermBoundaryCondition1Source> *&thermBc1, std::vector<Heat::ThermBoundaryCondition2Source> *&thermBc2, Heat::ThermMaterialSource &thermMat,
                    double q, double L, double T2);
void setBcSphereT_5(std::vector<Heat::ThermBoundaryCondition1Source> *&thermBc1, std::vector<Heat::ThermBoundaryCondition2Source> *&thermBc2, Heat::ThermMaterialSource &thermMat,
                    double q, double L);

// обработка результатов
struct ShperePoint
{
    double R;
    double sigma_r;
    double sigma_fi;
    double sigma_r_analit;
    double sigma_fi_analit;
    double sigma_r_pogr_abs;
    double sigma_fi_pogr_abs;
    double sigma_r_pogr;
    double sigma_fi_pogr;
};
void sortms(const Elementary::VECTOR3 &ms, int *ind);
void sortms_maxToMin(Elementary::VECTOR3 &ms);
void calcAnalit(double a,
                 double b,
                 double R,
                 double P0,
                 double elasticSigmaLimit,
                 double c,
                 double &sigma_r,
                 double &sigma_fi,
                 double &sigma_r_0,
                 double &sigma_fi_0
                 );
double calc_c(double a, double b, double P0, double elasticSigmaLimit);

// толстостенная сфера - температура
struct Test_sphere_T: public Test_base
{
    Grid::SphereParameters gp;
    int thermTestIndex;
    std::string fn_T_analit;
    std::string fn_T;
    std::string fn_T_pogr_abs;
    std::string fn_T_pogr;
    Test_sphere_T();
    virtual Type get_type()const;
    virtual void initTask(Solid::Task &task);
    virtual void writeResults(const Solid::Task &task, const Solid::OutData &out);
};

// толстостенная сфера - упруго-пластичность
struct Test_sphere_ep: public Test_base
{
    Grid::SphereParameters gp;
    int sigmaSolvingType;   // 0 - главные напряжения, 1 - проекции площадки
    std::string fn_inf;
    std::string fn_curve;
    std::string fn_dcurve;
    std::string fn_sigma_r_analit;
    std::string fn_sigma_r;
    std::string fn_sigma_fi_analit;
    std::string fn_sigma_fi;
    std::string fn_sigma_r_pogr_abs;
    std::string fn_sigma_fi_pogr_abs;
    std::string fn_sigma_r_pogr;
    std::string fn_sigma_fi_pogr;
    std::string fn_res_maxpogr;
    std::string fn_s_nonlinearStateFENumber;
    Test_sphere_ep();
    virtual Type get_type()const;
    virtual void initTask(Solid::Task &task);
    virtual void writeResults(const Solid::Task &task, const Solid::OutData &out);
};

// толстостенная сфера - ползучесть
struct Test_sphere_creep: public Test_base
{
    Grid::SphereParameters gp;
    int sigmaSolvingType;   // 0 - главные напряжения, 1 - проекции площадки
    std::string fn_sigma_r_analit;
    std::string fn_sigma_r;
    std::string fn_sigma_fi_analit;
    std::string fn_sigma_fi;
    std::string fn_sigma_eqv;
    std::string fn_eps_el_eqv;
    std::string fn_eps_pl_eqv;
    Test_sphere_creep();
    virtual Type get_type()const;
    virtual void initTask(Solid::Task &task);
    virtual void writeResults(const Solid::Task &task, const Solid::OutData &out);
};

// толстостенная сфера - термо-упругость
struct Test_sphere_Te: public Test_base
{
    Grid::SphereParameters gp;
    int thermTestIndex;
    int sigmaSolvingType;   // 0 - главные напряжения, 1 - проекции площадки
    int posledovatelnost;
    std::string fn_T_analit;
    std::string fn_T;
    std::string fn_T_pogr;
    std::string fn_u_analit;
    std::string fn_u;
    std::string fn_u_pogr;
    std::string fn_sigma_r_analit;
    std::string fn_sigma_r;
    std::string fn_sigma_r_pogr;
    std::string fn_sigma_fi_analit;
    std::string fn_sigma_fi;
    std::string fn_sigma_fi_pogr;
    std::string fn_max_pogr;
    Test_sphere_Te();
    virtual Type get_type()const;
    virtual void initTask(Solid::Task &task);
    virtual void writeResults(const Solid::Task &task, const Solid::OutData &out);
};

// формовка
struct Test_formovka: public Test_base
{
    Grid::FormovkaParameters gp;
    std::string fn_F;
    std::string fn_h;
    std::string fn_sigmaTop1;
    std::string fn_sigmaTop2;
    std::string fn_sigmaTop3;
    std::string fn_sigmaMiddle1;
    std::string fn_sigmaMiddle2;
    std::string fn_sigmaMiddle3;
    std::string fn_sigmaBottom1;
    std::string fn_sigmaBottom2;
    std::string fn_sigmaBottom3;
    std::string fn_T;
    Test_formovka();
    virtual Type get_type()const override;
    virtual void initTask(Solid::Task &task) override;
    virtual void writeResults(const Solid::Task &task, const Solid::OutData &out) override;
    virtual bool possibleToShow2d()const override;
    virtual bool getContactSurfaceCircle(const Solid::Task &task, const int globalStepIndex, Elementary::POINT2 &c, double &R)const override;
    virtual void needToDrawFe(const Solid::Task &task, const Solid::OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const override;
    virtual bool needToPaintFiniteElementSurface()const override;
};

// вдавливание цилиндра в полупространство
struct Test_contactCylinder: public Test_base
{
    Grid::ContactCylinderParameters ccp;
    // пути к текстовым файлам с данными (для каждого глобального шага в отдельной директории)
     std::string fn_inf;
    std::string fn_P;
    std::string fn_P_analit_d;
    std::string fn_P_analit_F;
    std::string fn_h;
     std::string fn_Pxyz;
     std::string fn_F;
    std::string fn_sigma1;
    std::string fn_sigma2;
    std::string fn_sigma3;
     std::string fn_stiffness;
    // пути к текстовым файлам с данными (в корневой директории, общие для всех глобальных шагов)
    std::string fn_Fn;
    std::string fn_Fn_analit;
    std::string fn_a;
    std::string fn_a_analit;
    std::string fn_Pmax;
    std::string fn_Pmax_analit;
    Test_contactCylinder();
    virtual Type get_type()const override;
    virtual void initTask(Solid::Task &task) override;
    virtual void writeResults(const Solid::Task &task, const Solid::OutData &out) override;
    virtual bool possibleToShow2d()const override;
    virtual bool getContactSurfaceCircle(const Solid::Task &task, const int globalStepIndex, Elementary::POINT2 &c, double &R)const override;
    virtual void needToDrawFe(const Solid::Task &task, const Solid::OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const override;
    virtual bool needToPaintFiniteElementSurface()const override;
    // вспомогательные методы
    double calc_r(const Elementary::POINT3 &p);
    double calc_P(const double d, const Grid::ContactCylinderParameters &ccp);
};

// вдавливание шара в полупространство
struct Test_contactSphere: public Test_base
{
    Grid::ContactSphereParameters ccp;
    // пути к текстовым файлам с данными (для каждого глобального шага в отдельной директории)
     std::string fn_input;
     std::string fn_general_inf;
     std::string fn_inf;
    // решения в контактных узлах
    std::string fn_P;
    std::string fn_P_analit;
    std::string fn_h;
     std::string fn_F;
    // численное решение в центрах приповерхностных КЭ
    std::string fn_z0_sigma1;
    std::string fn_z0_sigma2;
    std::string fn_z0_sigma3;
    std::string fn_z0_sigma_r;
    std::string fn_z0_sigma_fi;
    std::string fn_z0_sigma_z;
    // аналитическое решение во всей области
    /*
    std::string fn_2d_sigma_r_analit;
    std::string fn_2d_sigma_fi_analit;
    std::string fn_2d_sigma_z_analit;
    std::string fn_2d_sigma1_analit;
    std::string fn_2d_sigma2_analit;
    std::string fn_2d_sigma3_analit;
    */
    // численное решение в центрах всех КЭ
    std::string fn_2d_sigma_r;
    std::string fn_2d_sigma_fi;
    std::string fn_2d_sigma_z;
    std::string fn_2d_sigmaEqv;
    std::string fn_2d_sigma1;
    std::string fn_2d_sigma2;
    std::string fn_2d_sigma3;
    std::string fn_2d_sigmaEqv_unload;
    std::string fn_2d_sigma1_unload;
    std::string fn_2d_sigma2_unload;
    std::string fn_2d_sigma3_unload;
    // пути к текстовым файлам с данными (в корневой директории, общие для всех глобальных шагов)
    std::string fn_curve_diagramma;
    std::string fn_curve;
    std::string fn_dcurve;
    std::string fn_Fn;
    std::string fn_Fn_analit;
    std::string fn_a;
    std::string fn_a_analit;
    std::string fn_Pmax;
    std::string fn_Pmax_analit;
     std::string fn_bP_bd;
      std::string fn_bP_bd_analit;
     std::string fn_unload_bP_bd;
     std::string fn_unload_bP_bd_analit;




     std::string fn_bA_bd_true;
     std::string fn_bA_bd_geom;
     std::string fn_bA_bd_shamanstvo;
     std::string fn_unload_bA_bd;
     std::string fn_unload_bA_bd_analit;




    Test_contactSphere();
    virtual Type get_type()const override;
    virtual void initTask(Solid::Task &task) override;
    virtual void writeResults(const Solid::Task &task, const Solid::OutData &out) override;
    virtual bool possibleToShow2d()const override;
    virtual bool getContactSurfaceCircle(const Solid::Task &task, const int globalStepIndex, Elementary::POINT2 &c, double &R)const override;
    virtual void notFixedCoordinates(int &cn1, int &cn2) override;
    virtual void needToDrawFe(const Solid::Task &task, const Solid::OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const override;
    virtual bool needToPaintFiniteElementSurface()const override;
    // вспомогательные методы
    double calc_r(const Elementary::POINT3 &p);
    double calc_P(const double d, const Grid::ContactSphereParameters &csp);

    void calc_elastic_ms(const Grid::ContactSphereParameters &csp, const double a_analit, const double z, const double r, Elementary::VECTOR3 &r_fi_z, Elementary::VECTOR3 &ms);
    void calc_number_ms(const double a_analit, const double P_m_analit, const Solid::MechOutFePointData &fePointData,
                         double &z, double &r, Elementary::VECTOR3 &r_fi_z, Elementary::VECTOR3 &ms);
};

// ползучесть
struct Test_creep: public Test_base
{
    Grid::CreepParameters cp;
    // пути к текстовым файлам с данными
    std::string fn_eps11;
    std::string fn_eps22;
    std::string fn_eps33;
    std::string fn_eps11_analit;
    std::string fn_eps22_analit;
    std::string fn_eps33_analit;
    std::string fn_sigmaEqv;
    void setBc2Creep(double Px1, double Px2,
                     double Py1, double Py2,
                     double t1, double t2, std::vector<Solid::MechBoundaryCondition2Source_base *> *&bc2Source);
    Test_creep();
    virtual Type get_type()const override;
    virtual void initTask(Solid::Task &task) override;
    virtual void writeResults(const Solid::Task &task, const Solid::OutData &out) override;
    virtual bool possibleToShow2d()const override;
    virtual bool getContactSurfaceCircle(const Solid::Task &task, const int globalStepIndex, Elementary::POINT2 &c, double &R)const override;
    virtual void needToDrawFe(const Solid::Task &task, const Solid::OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const override;
    virtual bool needToPaintFiniteElementSurface()const override;
};

// задача Кирша (растяжение пластинки с отверстием)
struct Test_Kirsch: public Test_base
{
    Grid::KirschParameters gp;
    std::string fn_parametors;
    std::string fn_general_inf;
    std::string fn_sigma_fi;
    Test_Kirsch();
    virtual Type get_type()const override;
    virtual void initTask(Solid::Task &task) override;
    virtual void writeResults(const Solid::Task &task, const Solid::OutData &out) override;
    virtual bool possibleToShow2d()const override;
    virtual bool getContactSurfaceCircle(const Solid::Task &task, const int globalStepIndex, Elementary::POINT2 &c, double &R)const override;
    virtual void needToDrawFe(const Solid::Task &task, const Solid::OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const override;
    virtual bool needToPaintFiniteElementSurface()const override;
};

// трещина в неоднородной пластине (обобщение: базисные ф-и задаются интерполянтами)
struct Test_crack_plate_gen: public Test_base
{
    Grid::CrackPlateGenParameters gp;
    std::string fn_F;
    Test_crack_plate_gen();
    virtual Type get_type()const override;
    virtual void initTask(Solid::Task &task) override;
    virtual void writeResults(const Solid::Task &task, const Solid::OutData &out) override;
    virtual bool possibleToShow2d()const override;
    virtual bool getContactSurfaceCircle(const Solid::Task &task, const int globalStepIndex, Elementary::POINT2 &c, double &R)const override;
    virtual void needToDrawFe(const Solid::Task &task, const Solid::OutData &out, const size_t feInd, const int z_index, std::array<bool, 6> &faceState)const override;
    virtual void needToDrawSubFe(const Solid::Task &task, const Solid::OutData &out, const size_t feInd, const size_t subHexInd, const int z_index, std::array<bool, 6> &faceState)const override;
    virtual bool needToPaintFiniteElementSurface()const override;
    virtual double min_value(const int valueIndex)const override;
    virtual double max_value(const int valueIndex)const override;
};

}   // namespace Tests

#endif  // TESTS_H
