#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"

#include "heat.h"

namespace Heat
{
// 6-гранник, заданный явно координатами 8-ми вершин
struct LinearHexagon
{
    POINT3 v[8];
    inline POINT3 &operator[](const int &i)
    {
        return v[i];
    }
    inline const POINT3 &operator[](const int &i)const
    {
        return v[i];
    }
};
// 6-гранник, заданный явно координатами 27-ми точек
struct QuadraticHexagon
{
    POINT3 v[27];
    inline POINT3 &operator[](const int &i)
    {
        return v[i];
    }
    inline const POINT3 &operator[](const int &i)const
    {
        return v[i];
    }
};



// инициализация
void ThermSolver::init(ThermTask &thermTask, Threads::Logger_base *set_logger, Threads::SignalToThread_base *set_signal)
{
    initLogging(set_logger, set_signal);
    *((ThermTask *)this) = thermTask;  // копия входных данных(указателей)
    out = nullptr;
    if(!enabled) return;
    firstStep = true;              // портреты матриц не построены
    // матрицы
    matrixSizeT = (int)grid->vertex.size();
    PRINT(" size = %d\n",ARGS(matrixSizeT));
    GmatrixT.init();
    GbulderT.setMatrixSize(matrixSizeT);
    MmatrixT.init();
    MbulderT.setMatrixSize(matrixSizeT);
    // векторы
    bT.resize(matrixSizeT);
    qT.resize(matrixSizeT);
    q1T.resize(matrixSizeT);
    MqT.resize(matrixSizeT);
    // вспомогательные данные для учёта первых краевых условий
    bc1_u0T.resize(matrixSizeT);
    bc1_stateT.resize(matrixSizeT);
    // подготовка к первому шагу по времени
    for (int i = 0; i < matrixSizeT; i++)
    {
        qT[i] = (*vertex)[i].T;
        q1T[i] = 0;
    }
    // выходные данные
    out = new std::vector<ThermOutGlobalStepData>;
    GpreconditionerT = nullptr;
    GsolverT = nullptr;
}
// подготовка перед шагом по времени
void ThermSolver::prepareForStep()
{
    if(!enabled) return;
    VectorCopy(qT, q1T);
}
// шаг по времени
void ThermSolver::tryStep(const int globalStepNumber, Grid::TimeLayers &tl)
{
    if(!enabled) return;
    ThermGlobalStep &step0 = (*thermStep)[globalStepNumber];
    // подготовка к учёту краевых условий первого рода
    for (int i = 0; i < matrixSizeT; i++)
    {
        bc1_stateT[i] = false;
        bc1_u0T[i] = 0;
    }
    for (size_t bc1Ind = 0; bc1Ind < grid->bc1.size(); bc1Ind++)
    {
        Grid::BoundaryCondition1 &bc1_el = grid->bc1[bc1Ind];
        ThermBoundaryCondition1Source &bc1SourceT_el = (*step0.bc1SourceT)[bc1_el.bc1SourceIndex];
        if(bc1SourceT_el.mode != -1)
        {
            int m = bc1_el.vertexIndex;
            bc1_stateT[m] = true;
            bc1_u0T[m] = bc1SourceT_el.T0;
        }
    }
    // построение матриц
    PRINT1("T: genGMb..");
    MbulderT.start();
    GbulderT.start();
    buldGMb_T(grid, fe, step0.material, GbulderT, MbulderT, bT);
    // вторые краевые условия
    PRINT1("bc2..");
    addBc2_T(grid, step0.bc2SourceT, bT, GbulderT);
    PRINT1("bld..");
    if(firstStep)
    {
        // портрет создается только 1 раз
        GbulderT.complete(GmatrixT);
        GbulderT.fixReservedMemory();
        MbulderT.complete(MmatrixT);
        MbulderT.fixReservedMemory();
        //SSCMcopy(Gelements, tempGelements);
        GmatrixT.p->saveBMP("TG.bmp", 400);
        MmatrixT.p->saveBMP("TM.bmp", 400);
    }
    else
    {
        // сборка с учетом известного портрета
        GbulderT.completeWithPortrait(GmatrixT);
        MbulderT.completeWithPortrait(MmatrixT);
    }
    // 2-слойная неявная схема
    if(step0.timeMode == 1)
    {
        PRINT1("b..");
        SlauSolving::SSCMmulVector(MmatrixT, q1T, MqT);
        Vector1PlusCmulVector2(bT, 1./tl.dt, MqT, bT);
        PRINT1("G+M..");
        SlauSolving::SSCM1addEnclosedM2mulScalar(GmatrixT, MmatrixT, 1./tl.dt);
    }
    PRINT1("bc1..");
    SlauSolving::SSCMaddBoundaryCondition1(GmatrixT, bT, bc1_u0T, bc1_stateT);
    PRINT1("solv..");
    double timePreconditioner;
    if(GpreconditionerT == nullptr)
        GpreconditionerT = SlauSolving::SSCMPreconditioner_base::gen(step0.slausolverParameters.preconditioning);
    GpreconditionerT->bulid(GmatrixT, 0, timePreconditioner);
    if(GsolverT == nullptr)
        GsolverT = SlauSolving::SSCMSolver_base::gen(step0.slausolverParameters.type);
    if(firstStep)
        GsolverT->init(GmatrixT.getMatrixSize());
    GsolverT->solve(GmatrixT, bT,
                    GpreconditionerT, qT, step0.slausolverParameters,
                    qT, nlInf.slauResidual, nlInf.slauRelativeResidual, nlInf.slauIterations, nlInf.slauTime);
    PRINT(" iter=%6d residual=%.1le dresidual=%.1le\n",ARGS(nlInf.slauIterations, nlInf.slauResidual, nlInf.slauRelativeResidual));








    // температуры в узлах
    for (size_t vertexInd = 0; vertexInd < grid->vertex.size(); vertexInd++)
    {
        (*vertex)[vertexInd].newT = qT[vertexInd];
    }
    // температуры внутри конечных элементов
    // заполнение таблиц значений на шаблонных кубах
    using namespace Integration;
    Integrator integrationCubeSource[3];
    for(int it = (int)IntegrationType::Gauss2; it <= (int)IntegrationType::Gauss5; it++)
        integrationCubeSource[it].init3D((IntegrationType)it);
    for (size_t feInd = 0; feInd < grid->fe.size(); feInd++)
    {
        ThermFeData &feDataEl = (*fe)[feInd];
        int vi[8];
        grid->fe[feInd]->getVertexIndexes(vi);
        for (size_t pInd = 0; pInd < feDataEl.pd.size(); pInd++)    // индекс точки Гаусса внутри конечного элемента
        {
            ThermFePointData &fePD = feDataEl.pd[pInd];
            POINT3_CUBE XYZ;    // координаты точки в шаблонном кубе
            Integrator &integrationCube = integrationCubeSource[(int)feDataEl.integrationType];
            // центр конечного элемента k
            if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_homogeny)
            {
                XYZ = VECTOR3_NULL;		// центр конечного элемента в координатах шаблонного куба
            }
            if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_not_homogeny)
            {
                XYZ = integrationCube.p3[pInd]; //## метод интегрирования feDataEl.integrationType
            }
            // приращения температуры
            fePD.newT = 0;//fePD.T;
            for (int i = 0; i < 8; i++)	// i - локальный номер вершины - локальной базисной функции
            {
                int global_i = vi[i];	// глобальный номер вершины с локальным номером i
                fePD.newT += Fem::cubeLagrange1_3D(XYZ, i, Fem::dif_NULL3) * qT[global_i];
            }
        }
    }

    firstStep = false;  // теперь матрицы будут собираться с известным портретом
}
// завершение шага по времени (добавка приращений, обновление значений величин)
void ThermSolver::finalizeStep()
{
    if(!enabled) return;
    using namespace Operations;
    // конечные элементы
    for (size_t feInd = 0; feInd < grid->fe.size(); feInd++)
    {
        ThermFeData &feDataEl = (*fe)[feInd];
        for (size_t pInd = 0; pInd < feDataEl.pd.size(); pInd++)    // индекс точки Гаусса внутри конечного элемента
        {
            ThermFePointData &fePD = feDataEl.pd[pInd];
            // обновляется температура
            fePD.T = fePD.newT;
        }
    }
    // узлы сетки
    for (size_t vertexIndex = 0; vertexIndex < grid->vertex.size(); vertexIndex++)
    {
        ThermVertexData &vertexEl = (*vertex)[vertexIndex];
        vertexEl.T = vertexEl.newT;
    }
}
// сохранение найденных значений величин после шага по времени
void ThermSolver::saveResultsOfStep(const int globalStepNumber, const Grid::TimeLayers &tl)
{
    if(!enabled) return;
    using namespace Integration;
    ThermOutStepData s;
    LinearHexagon linearHexagon;
    QuadraticHexagon quadraticHexagon;
    Integrator integrationCubeSource[3];
    for(int it = (int)IntegrationType::Gauss2; it <= (int)IntegrationType::Gauss5; it++)
        integrationCubeSource[it].init3D((IntegrationType)it);
    s.t = tl.t0;
    if(out->empty())
        out->push_back({});
    out->back().step.push_back(s);
    if(tl.t0 == (*step)[globalStepNumber].t_finish)
    {
        out->back().fe.resize(grid->fe.size());
        for (size_t feInd = 0; feInd < grid->fe.size(); feInd++)
        {
            const ThermFeData &feDataEl = (*fe)[feInd];
            ThermOutFeData &outFeDataEl = (*out).back().fe[feInd];
            Integrator &integrationCube = integrationCubeSource[(int)feDataEl.integrationType];
            // копии координат вершин 6-гранника
            const Grid::FE_base *fe_el = grid->fe[feInd];
            POINT3 v[27];
            fe_el->getGeomVertexes(grid->vertex, grid->vertexForCurvature, v);
            // копируем информацию в точках элемента
            outFeDataEl.pd.resize(feDataEl.pd.size());
            for (size_t pInd = 0; pInd < feDataEl.pd.size(); pInd++)	// индекс точки Гаусса внутри конечного элемента
            {
                const ThermFePointData &fePD = feDataEl.pd[pInd];
                ThermOutFePointData &outFePD = outFeDataEl.pd[pInd];  // данные вывода для точки конечного элемента номер i
                outFePD.T = fePD.T;
                POINT3 XYZ;
                if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_homogeny)
                {
                    XYZ = POINT3(0,0,0);
                }
                if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_not_homogeny)
                {
                    XYZ = integrationCube.p3[pInd];
                }
                POINT3 xyz;
                fe_el->cubeToHexagon(v, XYZ, Fem::dif_NULL3, xyz);
                outFePD.p = xyz;
            }
            out->back().fe.push_back(outFeDataEl);
        }
        for (size_t vertexIndex = 0; vertexIndex < grid->vertex.size(); vertexIndex++)
        {
            ThermOutVertexData vertex_el;  // данные вывода для вершины vertexIndex
            vertex_el.T = (*vertex)[vertexIndex].T;
            out->back().vertex.push_back(vertex_el);
        }
        out->push_back({});
    }
}

void buldGMb_T(const Grid::Grid3D *grid, const std::vector<ThermFeData> *fe, std::vector<ThermMaterialSource> *material,
             SlauSolving::SSCMBulder &GbulderT, SlauSolving::SSCMBulder &MbulderT, Vector &bT)
{
    using namespace Integration;
    double w[Integration3D_Size_Gauss3];
    double dbas[8][Integration3D_Size_Gauss3][3];       // значения производных базисных функций(8) в каждой точке интегрирования 6-гранника

    std::vector<double> integrationwSource[3];
    std::vector<double> basCubeSource[3];
    std::vector<VECTOR3> dLinearBasCubeSource[3];
    std::vector<VECTOR3> dQuadraticBasCubeSource[3];

    int vi[8];              // глобальный индекс вершины шестигранника
    double intForG_mnij[8][8][3][3];
    double intForM_mn[8][8];
    double intForb_m[8];

    // заполнение таблиц значений на шаблонных кубах
    for(int it = (int)IntegrationType::Gauss2; it <= (int)IntegrationType::Gauss3; it++)
    {
        Fem::calcCubeLinearBasisFuncValues_ttt((IntegrationType)it,
                                            integrationwSource[it],basCubeSource[it],dLinearBasCubeSource[it]);
        Fem::calcCubeQuadraticBasisFuncValues_ttt((IntegrationType)it,
                                               dQuadraticBasCubeSource[it]);
    }

    // обнуление вектора правой части
    for (size_t i = 0; i < bT.size(); i++)
        bT[i] = 0;

    // шестигранники Solid_FeData_BasisType_1L
    for (size_t feInd = 0; feInd < grid->fe.size(); feInd++)	// индекс конечного элемента
    {
        const ThermFeData &feDataEl = (*fe)[feInd];
        ThermMaterialSource &m0 = (*material)[grid->fe[feInd]->mi];
        MATR3x3 L = m0.L;   // тензор теплопроводности
        double ro = m0.ro;  // плотность
        double c = m0.c;    // удельная теплопроводность материала
        double f = m0.f;    // мощность внутренних объёмных источников (стоков) тепла
        int it = (int)feDataEl.integrationType;
        std::vector<double> &integrationw = integrationwSource[it];
        std::vector<double> &basCube = basCubeSource[it];
        std::vector<VECTOR3> &dLinearBasCube = dLinearBasCubeSource[it];
        //if(grid->fe[feInd].type == Grid::FEType::QuadraticHexagon)
        std::vector<VECTOR3> &dQuadraticBasCube = dQuadraticBasCubeSource[it];
        int numPoints = (int)integrationw.size();
        // копии координат вершин 6-гранника
        // расчёт коэффициентов для интегрирования, модулей детерминантов
        // и производных базисных функций на шестиграннике
        grid->fe[feInd]->calcBasisFuncValues(grid->vertex, grid->vertexForCurvature, numPoints, integrationw.data(), dLinearBasCube.data(), dQuadraticBasCube.data(),
                                    w, dbas);
        grid->fe[feInd]->getVertexIndexes(vi);
        // внутренние интегралы для добавок к матрицам и к правой части
        // для матрицы G
        if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_homogeny)
            Fem::calcIntForG(numPoints, w, dbas, intForG_mnij);
        // для матрицы M
        Fem::calcIntForM_ttt(numPoints, w, basCube, intForM_mn);
        // для правой части
        Fem::calcIntForb_f(numPoints, w, basCube.data(), intForb_m);
        // I. добавка к матрице G от конечного элемента feInd
        for (int m = 0; m < 8; m++)
        {
            int global_m = vi[m];
            for (int n = 0; n < 8; n++)	// m, n - локальные номера базисных функций (вершин)
            {
                int global_n = vi[n];
                if (global_m >= global_n)
                {
                    double E = 0;
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            E += L.m[i][j] * intForG_mnij[m][n][i][j];
                    GbulderT.addElement(
                        E,
                        global_m,
                        global_n);
                }
            }
        }
        // II. добавка к матрице M от конечного элемента feInd
        for (int m = 0; m < 8; m++)
        {
            int global_m = vi[m];
            for (int n = 0; n < 8; n++)	// m, n - локальные номера базисных функций (вершин)
            {
                int global_n = vi[n];
                if(global_m >= global_n)
                {
                    MbulderT.addElement(
                        ro * c * intForM_mn[m][n],
                        global_m,
                        global_n);
                }
            }
        }
        // III добавка к вектору (внутренние объёмные источники (стоки) тепла: ro*c*...)
        for (int m = 0; m < 8; m++)	// i - локальный номер базисной функции (вершины)
            bT[vi[m]] += f * intForb_m[m];
    } // feInd
}

void addBc2_T(const Grid::Grid3D *grid, const std::vector<ThermBoundaryCondition2Source> *bc2Source,
              Vector &b, SlauSolving::SSCMBulder &GbulderT)
{
    using namespace Integration;
    Integrator integrationFoursquare;
    integrationFoursquare.init2D(IntegrationType::Gauss3);//##
    double linearBasCube[8*Integration::Integration2D_Size_Max];
    double hexagonCoef[Integration::Integration2D_Size_Max];
    VECTOR3 hexagonNormal[Integration::Integration2D_Size_Max];

    for (size_t FEsurfaceInd = 0; FEsurfaceInd < grid->FESurface.size(); FEsurfaceInd++)
    {
        int si = (int)FEsurfaceInd;  // FEsurfaceInd соответствует номеру ресурса второго краевого
        const Grid::FiniteElementSurface &s = grid->FESurface[FEsurfaceInd];
        for (size_t faceInd = 0; faceInd < s.face.size(); faceInd++)
        {
            const Grid::FEFace &face = s.face[faceInd];
            if(si != -1)
            {
                int feInd = face.feIndex;
                const Grid::FE_base *fe_el = grid->fe[feInd];
                const ThermBoundaryCondition2Source &bc2_el = (*bc2Source)[si];
                int vi[8];
                fe_el->getVertexIndexes(vi);
                if(bc2_el.mode == ThermBoundaryCondition2Type::scalar)
                {
                    fe_el->calcBc2BasisFuncValues(grid->vertex, grid->vertexForCurvature, integrationFoursquare, face.faceIndex,
                                                   linearBasCube, hexagonCoef, hexagonNormal);
                    double q = bc2_el.q;
                    double hi = bc2_el.hi;
                    double Ta = bc2_el.Ta;
                    // добавка к правой части
                    for (int m = 0; m < 8; m++)		// m - локальный номер базисной функции (вершины)
                    {
                        double E = 0;
                        for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
                        {
                            E += (q + hi*Ta) *
                                 linearBasCube[m*integrationFoursquare.size + valIndex] *
                                 hexagonCoef[valIndex];
                        }
                        b[vi[m]] += E;
                    }
                    // добавка в матрицу жесткости
                    for (int m = 0; m < 8; m++)
                        for (int n = 0; n < 8; n++)
                        {
                            int global_m = vi[m];
                            int global_n = vi[n];
                            if(global_m >= global_n)
                            {
                                double E = 0;
                                for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
                                {
                                    E += hi *
                                         linearBasCube[m*integrationFoursquare.size + valIndex] *
                                         linearBasCube[n*integrationFoursquare.size + valIndex] *
                                         hexagonCoef[valIndex];
                                }
                                GbulderT.addElement_not_null(
                                            E,
                                            global_m,
                                            global_n);
                            }
                        }
                }
            }
        }
    }
}

}   // namespace Solid
