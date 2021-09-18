#define _CRT_SECURE_NO_WARNINGS

#include <set>
#include "stdio.h"

#include "solid.h"
#include "string.h"

#define MechTimeInf_print_01(f, s) {fprintf(f, " %-40s: ", s);}
#define MechTimeInf_print_02(f, s) {fprintf(f, "  %-40s: ", s);}
#define MechTimeInf_print_11(f, s) {fprintf(f, "  \t %-40s: ", s);}
#define MechTimeInf_print_12(f, s) {fprintf(f, "  \t  %-40s: ", s);}

namespace Solid
{
void MechTimeIterInf::print(const int numIter, FILE *f)
{
    fprintf(f, "  \tnumIter = %4d\n", numIter);
    // 1
    MechTimeInf_print_11(f, "buildGMb");
    buildGMb.print(f);
     MechTimeInf_print_12(f, "buildGMb_buildNonContactLocalGMb");
     buildGMb_buildNonContactLocalGMb.print(f);
     MechTimeInf_print_12(f, "buildGMb_addNonContactLocalGMb_to_global");
     buildGMb_addNonContactLocalGMb_to_global.print(f);
     MechTimeInf_print_12(f, "buildGMb_addContact");
     buildGMb_addContact.print(f);
    // 2
    MechTimeInf_print_11(f, "SSCMaddBoundaryCondition1");
    SSCMaddBoundaryCondition1.print(f);
    // 3
    MechTimeInf_print_11(f, "solveSLAU");
    solveSLAU.print(f);
     MechTimeInf_print_12(f, "solveSLAU_Gpreconditioner_init");
     solveSLAU_Gpreconditioner_init.print(f);
     MechTimeInf_print_12(f, "solveSLAU_Gsolver_solve");
     solveSLAU_Gsolver_solve.print(f);
    // 4
    MechTimeInf_print_11(f, "solveIterResults");
    calcIterResults.print(f);
     MechTimeInf_print_12(f, "solveIterResults_fe");
     calcIterResults_fe.print(f);
     MechTimeInf_print_12(f, "solveIterResults_contact");
     calcIterResults_contact.print(f);
}

MechSlauIterInf::MechSlauIterInf()
{
    matrixG_portraitChanged = true;   // портреты матрицы G не построены
    matrixG_elementsChanged = true;   // элементы матрицы G требуется обновить
    matrixM_portraitChanged = true;   // портреты матрицы M не построены
    matrixM_elementsChanged = true;   // элементы матрицы M требуется обновить
    precChanged = true;               // предобусловливатель требуется обновить
    GFirstStrChanged = 0;             // весь предобусловливатель требуется обновить
    GsolverChanged = true;            // решатель СЛАУ требуется обновить
    debugMatrixesInformationNeedToWrite = true;   // требуется вывести отладочную информацию о матрицах(параметры и портреты)
}
void MechSlauIterInf::finalizeStep()
{
    matrixG_elementsChanged = true;   // элементы матрицы G требуется обновить
    matrixM_elementsChanged = true;   // элементы матрицы M требуется обновить
    precChanged = true;               // предобусловливатель требуется обновить
    GFirstStrChanged = 0;             // весь предобусловливатель требуется обновить
}

MechIterInf::MechIterInf()
{
    BC1Changed = true;                // первые краевые требуется обновить
    BC2Changed = true;                // вторые краевые требуется обновить
}
void MechIterInf::prepareForIter(const MechGlobalStepParameters &step_el)
{
    iterNumber = 0;
    plastic.iterNumber = 0;
    contact.iterNumber = 0;
    switch (step_el.switchIterationsMode)
    {
    case SwitchIterationsMode::Serial:
    {
        // запуск итераций контакта и пластичности по-очереди
        // начнём с контакта (пластичность по касательной)
        contact.iterationMode = IterationMode::Free;
        plastic.iterationMode = IterationMode::ReadOnly;
    }break;
    case SwitchIterationsMode::Parallel:
    {
        // запуск итераций контакта и пластичности параллельно
        contact.iterationMode = IterationMode::Free;
        plastic.iterationMode = IterationMode::Free;
    }break;
    }
}
void MechIterInf::clearCounts(const int matrixSize)
{
    plastic.clearCounts(matrixSize);
    contact.clearCounts(matrixSize);
}
void MechIterInf::calcIterResults(const MechGlobalStepParameters &step_el, const bool needToReleaseState, const size_t matrixSize, MechIterationsGeneralResult &mechIterationsGeneralResult)
{
    plastic.calcIterResults(step_el);    
    contact.calcIterResults(step_el);
    // первая изменённая строка матрицы
    slau.GFirstStrChanged = MIN(slau.GFirstStrChanged, MIN(plastic.GFirstStrChanged, contact.GFirstStrChanged));
    // нужно ли пересобирать матрицу?
    if(slau.GFirstStrChanged >= (int)matrixSize)
        slau.matrixG_elementsChanged = false;
    else
        slau.matrixG_elementsChanged = true;

    // достигнутая точность
    if(contact.accuracyAchieved && plastic.accuracyAchieved)
    {
        accuracyAchieved = true;
        improvingAccuracy = contact.improvingAccuracy || plastic.improvingAccuracy;
    }
    else
    {
        accuracyAchieved = false;
        improvingAccuracy = false;
    }
    // обновление предыдущих невязок
    slau.lastSlauResidualWithLastq = slau.slauResidualWithLastq;
    plastic.prevMaxPlasticResidual = plastic.maxPlasticResidual;
    // обновление номеров итераций
    if(plastic.iterationMode == IterationMode::Free)
        plastic.iterNumber++;
    if(contact.iterationMode == IterationMode::Free)
        contact.iterNumber++;
    iterNumber++;
    // переключение автивной нелинейности для последовательного режима
    switch (step_el.switchIterationsMode)
    {
    case SwitchIterationsMode::Serial:
    {
        // меняем активную и пассивную нелинейности, если точность для активной нелинейности достигнута
        if(contact.iterationMode == IterationMode::Free)
        {
            if(contact.accuracyAchieved)
            {
                // точность достигнута - переключаемся на пластичность
                plastic.iterationMode = IterationMode::Free;
                contact.iterationMode = IterationMode::ReadOnly;
            }
        }
        else
        {
            if(plastic.accuracyAchieved)
            {
                // точность достигнута - переключаемся на контакт
                contact.iterationMode = IterationMode::Free;
                plastic.iterationMode = IterationMode::ReadOnly;
            }
        }
    }break;
    case SwitchIterationsMode::Parallel:
    {
    }break;
    }
    // общий результат итерации
    mechIterationsGeneralResult = MechIterationsGeneralResult::Continue_NeedMoreIterations;
    if(iterNumber >= step_el.iterLimit)
    {
        mechIterationsGeneralResult = MechIterationsGeneralResult::Interrupted_OutOfIterationsLimit;
    }
    if(accuracyAchieved)
    {
        mechIterationsGeneralResult = MechIterationsGeneralResult::Interrupted_Ok;
    }
    if(needToReleaseState == true)
    {
        mechIterationsGeneralResult = MechIterationsGeneralResult::Interrupted_Terminated;
    }
}
void MechIterInf::finalizeStep()
{
    //BC1Changed = true;                // первые краевые требуется обновить
    BC2Changed = true;                // вторые краевые требуется обновить
    slau.finalizeStep();
}
void MechIterInf::save(FILE *f)
{
    fprintf(f, "%le\t%le\t%le\t%le\t%le\t%d\t%le\t%le\t%d",
            slau.slauResidual,
            slau.slauRelativeResidual,
            slau.slauResidualWithLastq,
            contact.max_deltaF_residual,
            contact.max_endPoint_residual,
            contact.contactChangedNumber,
            plastic.maxPlasticResidual.eps,
            plastic.maxPlasticResidual.sigma,
            plastic.GLocalChangedNumber);
}
void MechIterInf::load(FILE *f)
{
    fscanf(f, "%le%le%le%le%le%d%le%le%d",
           &slau.slauResidual,
           &slau.slauRelativeResidual,
           &slau.slauResidualWithLastq,
           &contact.max_deltaF_residual,
           &contact.max_endPoint_residual,
           &contact.contactChangedNumber,
           &plastic.maxPlasticResidual.eps,
           &plastic.maxPlasticResidual.sigma,
           &plastic.GLocalChangedNumber);
}

void MechStepInf::print()
{
    FILE *timeInf_f = fopen("nlInf/timeInf.txt", "a");
    fprintf(timeInf_f, "globalStepNumber = %4d, stepNumber = %4d\n", globalStepNumber, stepNumber);
    // 1
    MechTimeInf_print_01(timeInf_f, "prepareForStep");
    prepareForStep.print(timeInf_f);
    // 2
    MechTimeInf_print_01(timeInf_f, "prepareForIter");
    prepareForIter.print(timeInf_f);
    MechTimeInf_print_02(timeInf_f, "prepareForIter_prepareBC1");
    prepareForIter_prepareBC1.print(timeInf_f);
    MechTimeInf_print_02(timeInf_f, "prepareForIter_prepareBC2");
    prepareForIter_prepareBC2.print(timeInf_f);
    MechTimeInf_print_02(timeInf_f, "prepareForIter_prepareContact");
    prepareForIter_prepareContact.print(timeInf_f);
    // 3
    MechTimeInf_print_01(timeInf_f, "tryStep");
    tryStep.print(timeInf_f);
    for(int numIter = 0; numIter < (int)(iterInf).size(); numIter++)
    {
        (iterInf)[numIter].timeIterInf.print(numIter, timeInf_f);
    }
    // 4
    MechTimeInf_print_01(timeInf_f, "finalizeStep");
    finalizeStep.print(timeInf_f);
     MechTimeInf_print_02(timeInf_f, "finalizeStep_fe");
     finalizeStep_fe.print(timeInf_f);
     MechTimeInf_print_02(timeInf_f, "finalizeStep_vertex");
     finalizeStep_vertex.print(timeInf_f);
     MechTimeInf_print_02(timeInf_f, "finalizeStep_contact");
     finalizeStep_contact.print(timeInf_f);
     MechTimeInf_print_02(timeInf_f, "finalizeStep_movingGrid");
     finalizeStep_movingGrid.print(timeInf_f);
    // 5
    MechTimeInf_print_01(timeInf_f, "saveResultsOfStep");
    saveResultsOfStep.print(timeInf_f);
    // 6
    MechTimeInf_print_01(timeInf_f, "writeResults");
    writeResults.print(timeInf_f);
    fclose(timeInf_f);
}

void MechTask::initNull(const Integration::IntegrationType set_integrationType)
{
    // инициализация начальных деформаций и напряжений (=0)
    fe = new std::vector<MechFeData_base *>;
    fe->resize(grid->fe.size());
    for (size_t i = 0; i < grid->fe.size(); i++)
    {
        switch (grid->fe[i]->get_FEType())
        {
        case Grid::FEType::LinearHexagon:
        case Grid::FEType::QuadraticHexagon:
        {
            MechFeData_LinearHexagon_homogeny_E_NLE_EP *feDataEl = new MechFeData_LinearHexagon_homogeny_E_NLE_EP;
            feDataEl->init(set_integrationType, grid->fe[i]);
            (*fe)[i] = feDataEl;
        }break;
        case Grid::FEType::LinearHexagon_XFEM:
        {
            //Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<Grid::FE_LinearHexagon_XFEM *>(grid->fe[i]);
            MechFeData_LinearXFEM_gen *feDataEl = new MechFeData_LinearXFEM_gen;
            feDataEl->init(set_integrationType, grid->fe[i]);
            (*fe)[i] = feDataEl;
        }break;
        }
    }
    // инициализация начальных перемещений (=0)
    vertex = new std::vector<MechVertexData>;
    for (size_t i = 0; i < grid->vertex.size(); i++)
    {
        MechVertexData v_el;
        for (int j = 0; j < 3; j++)
        {
            v_el.sum_du[j] = 0;
            v_el.du[j] = 0; // не имеет значения
        }
        vertex->push_back(v_el);
    }
    // для криволинейных границ ###
    vertexForCurvature = new std::vector<MechVertexData>;
    for (size_t i = 0; i < grid->vertexForCurvature.size(); i++)
    {
        MechVertexData v_el;
        for (int j = 0; j < 3; j++)
        {
            v_el.sum_du[j] = 0;
            v_el.du[j] = 0; // не имеет значения
        }
        vertexForCurvature->push_back(v_el);
    }
    // инициализация начальных скоростей (=0) и ускорений (=0)
    V0 = new Vector(grid->vertex.size()*3);
    dV0 = new Vector(grid->vertex.size()*3);
    for (size_t i = 0; i < V0->size(); i++)
    {
        (*V0)[i] = 0;
        (*dV0)[i] = 0;
    }
    // подготовка к учёту реакции опоры
    contactVertex_FE_Rigid.resize(0);
    // заполнение массива контактных вершин contactVertex_FE_Rigid
    // в соответствии с набором контактов FE поверхность - жёсткая поверхность
    for(size_t CCSource_FE_Rigid_ind = 0; CCSource_FE_Rigid_ind < (*CC_FE_Rigid).size(); CCSource_FE_Rigid_ind++)
    {
        ContactCondition_FE_Rigid &CCSource_FE_Rigid_el = (*CC_FE_Rigid)[CCSource_FE_Rigid_ind];  // контакт FE поверхность - жёсткая поверхность
        std::vector<Grid::FEFace> &face = (*grid).FESurface[CCSource_FE_Rigid_el.FESurfaceInd].face;    // набор граней
        std::set<int> pushedVertexIndexes;
        pushedVertexIndexes.clear();
        for(size_t faceInd = 0; faceInd < face.size(); faceInd++)
        {
            int vi[27];   // вершины, принадлежащие данной грани
            int surfaceVertexesNumber = (*grid).fe[face[faceInd].feIndex]->getFaceVertexIndexes(face[faceInd].faceIndex, vi);
            for(int i = 0; i < surfaceVertexesNumber; i++)
            {
                if(pushedVertexIndexes.count(vi[i]) == 0)   // дважды одну и ту же вершину добавлять не будем
                {
                    //if(bc1_state[3 * vi[i] + 0] == false &&
                    //   bc1_state[3 * vi[i] + 1] == false &&
                    //   bc1_state[3 * vi[i] + 2] == false)
                    {
                        // узлы, на которых заданы первые краевые условия, не включаем
                        ContactSurfaceVertexData_FE_Rigid CVD_FE_R;
                        CVD_FE_R.init(vi[i], CCSource_FE_Rigid_ind);
                        contactVertex_FE_Rigid.push_back(CVD_FE_R);
                        pushedVertexIndexes.insert(vi[i]);
                    }
                }
            }
        }
        // инициализация карт готовых решений для поверхностей
        //(*rigidSurface)[CCSource_FE_Rigid_el.RigidSurfaceInd]->
    }
}
void MechTask::save(const std::string &subdir) const
{
    std::string fn_vertexData = subdir + "vertexData.dat";
    std::string fn_feData = subdir + "feData.dat";
    FILE *f_vertexData = fopen(fn_vertexData.c_str(), "w");
    FILE *f_feData = fopen(fn_feData.c_str(), "w");
    // узлы
    fprintf(f_vertexData, "%d\n", 6);
    for(int vertexInd = 0; vertexInd < (int)(*vertex).size(); vertexInd++)
    {
        const MechVertexData &vertex_el = (*vertex)[vertexInd];
        for(int i = 0; i < 3; i++)
            fprintf(f_vertexData, "%le ", vertex_el.du[i]);
        for(int i = 0; i < 3; i++)
            fprintf(f_vertexData, "%le ", vertex_el.sum_du[i]);
        fprintf(f_vertexData, "\n");
    }
    // КЭ
    fprintf(f_feData, "%d\n", 4);
    for(int feInd = 0; feInd < (int)(*fe).size(); feInd++)
    {
        const MechFeData_LinearHexagon_homogeny_E_NLE_EP &fe_el = *(MechFeData_LinearHexagon_homogeny_E_NLE_EP *)(*fe)[feInd];
        fprintf(f_feData, "%le ", fe_el.pd.sumSigma_eqv);
        fprintf(f_feData, "%le ", fe_el.pd.q);
        fprintf(f_feData, "%le ", fe_el.pd.it.residual.sigma);
        fprintf(f_feData, "%le ", fe_el.pd.it.residual.eps);
        fprintf(f_feData, "\n");
    }
    fclose(f_vertexData);
    fclose(f_feData);
}

void MechSolver::prepareBC1(const MechGlobalStepParameters &step_el)
{
    MechIterInf &nlInf = stepInf.back().iterInf.back();
    if(nlInf.BC1Changed == true)
    {
        // подготовка к учёту краевых условий первого рода
        for (size_t i = 0; i < matrixSize; i++)
        {
            bc1_u0[i] = 0;
            bc1_state[i] = false;
        }
        for (size_t bc1Ind = 0; bc1Ind < grid->bc1.size(); bc1Ind++)
        {
            MechBoundaryCondition1Source &bc1Source_el = (*step_el.bc1Source)[grid->bc1[bc1Ind].bc1SourceIndex];
            const VECTOR3_int &DOFindex = grid->DOFs->findDOFIndex_Lagr1(grid->bc1[bc1Ind].vertexIndex);
            //VECTOR3_int globalFuncIndex = vertexBusFuncInd->linear[grid->bc1[bc1Ind].vertexIndex];
            for(int i = 0; i < 3; i++)
            {
                if (bc1Source_el.mode[i] != -1 && DOFindex[i] != -1)
                {
                    bc1_state[DOFindex[i]] = true;
                    bc1_u0[DOFindex[i]] = bc1Source_el.u0[i];
                }
            }
        }
        nlInf.BC1Changed = false;
    }
}
void MechSolver::prepareBC2(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl)
{
    MechIterInf &nlInf = stepInf.back().iterInf.back();
    if(nlInf.BC2Changed == true)
    {
        // расчёт значений функций поверхностных сил P для данного шага по времени
        for (size_t bc2SourceIndex = 0; bc2SourceIndex < (*step_el.bc2Source).size(); bc2SourceIndex++)
            (*step_el.bc2Source)[bc2SourceIndex]->updateValue(tl.t_prev1, tl.t0);
        nlInf.BC2Changed = false;
    }
}
void MechSolver::prepareContact(const Grid::TimeLayers &tl)
{
    // обновление жёстких поверхностей
    // контакты FE поверхность - жёсткая поверхность
    for(size_t CCSource_FE_Rigid_ind = 0; CCSource_FE_Rigid_ind < (*CC_FE_Rigid).size(); CCSource_FE_Rigid_ind++)
    {
        ContactCondition_FE_Rigid &CCSource_FE_Rigid_el = (*CC_FE_Rigid)[CCSource_FE_Rigid_ind];  // контакт FE поверхность - жёсткая поверхность
        (*rigidSurface)[CCSource_FE_Rigid_el.RigidSurfaceInd]->update(tl.t0, &grid->regionIndexationParameters);
    }
    // подготовка контактных узлов
    for(size_t CV_FE_Rigid_ind = 0; CV_FE_Rigid_ind < contactVertex_FE_Rigid.size(); CV_FE_Rigid_ind++)
    {
        contactVertex_FE_Rigid[CV_FE_Rigid_ind].prepareForIter();
    }
}
void MechSolver::buildGMb(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl)
{
    MechIterInf &nlInf = stepInf.back().iterInf.back();
    {
        int next_iterNumber = nlInf.iterNumber + 1;
        int next_contact_iterNumber = nlInf.contact.iterNumber + (nlInf.contact.iterationMode==IterationMode::Free);
        int next_plastic_iterNumber = nlInf.plastic.iterNumber + (nlInf.plastic.iterationMode==IterationMode::Free);
        PRINT("%.2d ci = %.2d pi = %.2d ", ARGS(next_iterNumber, next_contact_iterNumber, next_plastic_iterNumber));
    }
    logger->sendString("genGMb..");
    nlInf.timeIterInf.buildGMb_buildNonContactLocalGMb.begin();
    buildNonContactLocalGMb(step_el, tl);
    nlInf.timeIterInf.buildGMb_buildNonContactLocalGMb.end();
    if(nlInf.slau.matrixG_elementsChanged)
        logger->sendString("bld..");
    nlInf.timeIterInf.buildGMb_addNonContactLocalGMb_to_global.begin();
    addNonContactLocalGMb_to_global(step_el);
    nlInf.timeIterInf.buildGMb_addNonContactLocalGMb_to_global.end();
    logger->sendString("c..");
    nlInf.timeIterInf.buildGMb_addContact.begin();
    addContact(step_el, tl);
    nlInf.timeIterInf.buildGMb_addContact.end();
}
void MechSolver::buildNonContactLocalGMb(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl)
{
    using namespace Integration;
    MechIterInf &nlInf = stepInf.back().iterInf.back();
    IncForsesMode incForsesMode = step_el.incForsesMode;
    bool firstIter = (nlInf.iterNumber == 0);

    // локальные матрицы и векторы
    for (size_t feInd = 0; feInd < grid->fe.size(); feInd++)   // индекс конечного элемента
    {
        MechMaterialSource &m0 = (*step_el.material)[grid->fe[feInd]->mi];
        (*fe)[feInd]->updateLocalGMb(grid, grid->DOFs, tl, grid->fe[feInd], m0, nlInf.plastic, firstIter, incForsesMode);
    } // feInd

    // вторые краевые
    Integrator integrationFoursquare;
    integrationFoursquare.init2D(IntegrationType::Gauss3);//##
    for(size_t bc2Ind = 0; bc2Ind < (*step_el.bc2).size(); bc2Ind++)
    {
        Grid::FiniteElementSurface &s = grid->FESurface[(*step_el.bc2)[bc2Ind].FEsurfaceInd];
        for (size_t faceInd = 0; faceInd < s.face.size(); faceInd++)
        {
            Grid::FEFace &face = s.face[faceInd];
            (*fe)[face.feIndex]->add_BC2_to_localb(
                        grid,
                        grid->fe[face.feIndex],
                        (*step_el.bc2Source)[(*step_el.bc2)[bc2Ind].bc2SourceIndex],
                        face.faceIndex,
                        integrationFoursquare);
        }
    }
    // 2-е кроевые для поверхностей трещин
    for(size_t bc2CrackInd = 0; bc2CrackInd < (grid->bc2Crack).size(); bc2CrackInd++)
    {
        const Grid::MechBoundaryCondition2Crack &bc2Crack_el = (grid->bc2Crack)[bc2CrackInd];
        const Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<Grid::FE_LinearHexagon_XFEM *>(grid->fe[bc2Crack_el.FEIndex]);
        MechFeData_LinearXFEM_gen *XFEM_gen = dynamic_cast<MechFeData_LinearXFEM_gen *>((*fe)[bc2Crack_el.FEIndex]);
        XFEM_gen->add_BC2_to_localb(
                    grid,
                    feEl_XFEM,
                    (*step_el.bc2Source)[bc2Crack_el.bc2SourceIndex],
                    bc2Crack_el.subSurfaceIndex,
                    integrationFoursquare);
    }
    /*for(size_t bc2CrackInd = 0; bc2CrackInd < (*step_el.bc2Crack).size(); bc2CrackInd++)
    {
        const Grid::MechBoundaryCondition2Crack &bc2Crack_el = (*step_el.bc2Crack)[bc2CrackInd];
        const Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<Grid::FE_LinearHexagon_XFEM *>(grid->fe[bc2Crack_el.FEIndex]);
        MechFeData_LinearXFEM_gen *XFEM_gen = dynamic_cast<MechFeData_LinearXFEM_gen *>((*fe)[bc2Crack_el.FEIndex]);
        XFEM_gen->add_BC2_to_localb(
                    grid,
                    feEl_XFEM,
                    (*step_el.bc2Source)[bc2Crack_el.bc2SourceIndex],
                    bc2Crack_el.subSurfaceIndex,
                    integrationFoursquare);
    }*/

    // расчёт локальных векторов невязки (в начале шага, на 0-й итерации)
    if(firstIter)
    {
        switch (incForsesMode)
        {
        case IncForsesMode::MinusIntegral:
        {
            // вычисляется в начале шага, после сдвига сетки
            buildLocalR(step_el, tl);
        }break;
        case IncForsesMode::bPlusR:
        {
            // вычисляется в конце шага, после всех итераций, до сдвига сетки
        }break;
        case IncForsesMode::IncrementP:
        {
        }break;
        }
    }
}
void MechSolver::buildLocalR(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl)
{
    using namespace Integration;
    MechIterInf &nlInf = stepInf.back().iterInf.back();
    IncForsesMode incForsesMode = step_el.incForsesMode;

    // векторы невязок
    // шестигранники Solid_FeData_BasisType_1L
    for (size_t feInd = 0; feInd < grid->fe.size(); feInd++)   // индекс конечного элемента
    {
        MechMaterialSource &m0 = (*step_el.material)[grid->fe[feInd]->mi];
        (*fe)[feInd]->calcLocalR(grid, tl, grid->fe[feInd], m0, nlInf.plastic, incForsesMode);
    } // feInd

    // векторы невязок (вторые краевые)
    Integrator integrationFoursquare;
    integrationFoursquare.init2D(IntegrationType::Gauss3);//##
    for(size_t bc2Ind = 0; bc2Ind < (*step_el.bc2).size(); bc2Ind++)
    {
        Grid::FiniteElementSurface &s = grid->FESurface[(*step_el.bc2)[bc2Ind].FEsurfaceInd];
        for (size_t faceInd = 0; faceInd < s.face.size(); faceInd++)
        {
            Grid::FEFace &face = s.face[faceInd];
            (*fe)[face.feIndex]->add_BC2_to_localR(
                        grid,
                        grid->fe[face.feIndex],
                        (*step_el.bc2Source)[(*step_el.bc2)[bc2Ind].bc2SourceIndex],
                        face.faceIndex,
                        integrationFoursquare,
                        incForsesMode);
        }
    }
}
void MechSolver::addNonContactLocalGMb_to_global(const MechGlobalStepParameters &step_el)
{
    MechIterInf &nlInf = stepInf.back().iterInf.back();
    // построение портрета и инициализация предобусловливателя
    if(nlInf.slau.matrixG_portraitChanged)
    {
        SlauSolving::SSCMPortraitBulder pb;
        if(step_el.slausolverParameters.blocks)
        {
            pb.init(matrixSize/3);
        }
        else
        {
            pb.init(matrixSize);
        }
        for (size_t feInd = 0; feInd < grid->fe.size(); feInd++)
        {
            if(step_el.slausolverParameters.blocks)
            {
                (*fe)[feInd]->addGlobalIndexesToGlobalPortrait3x3(pb);
            }
            else
            {
                (*fe)[feInd]->addGlobalIndexesToGlobalPortrait1x1(pb);
            }
        }
        if(step_el.slausolverParameters.blocks)
        {
            pb.completePortrait(*Gmatrix3x3.p);
            (*Gmatrix3x3.e).init(Gmatrix3x3.getElementsNumber(), Gmatrix3x3.getMatrixSize()); // инициализация матрицы
        }
        else
        {
            pb.completePortrait(*Gmatrix.p);
            (*Gmatrix.e).init(Gmatrix.getElementsNumber(), Gmatrix.getMatrixSize()); // инициализация матрицы
        }
        if(nlInf.slau.debugMatrixesInformationNeedToWrite)
        {
            // Единожды вывести информацию о размерностях задачи и портреты
            PRINT1("SavingPortraits G..");
            if(step_el.slausolverParameters.blocks)
            {
                Gmatrix3x3.p->saveProperties("3x3_slauInformation.txt");
                Gmatrix3x3.p->saveBMP("3x3_G.bmp", Gmatrix3x3.getMatrixSize());
            }
            else
            {
                Gmatrix.p->saveProperties("1x1_slauInformation.txt");
                Gmatrix.p->saveBMP("1x1_G.bmp", Gmatrix.getMatrixSize());
            }
        }
        // инициализация предобусловливателя
        if(step_el.slausolverParameters.blocks)
        {
            if(Gpreconditioner3x3 == nullptr)
                Gpreconditioner3x3 = SlauSolving::SSCM3x3Preconditioner_base::gen(step_el.slausolverParameters.preconditioning);
        }
        else
        {
            if(Gpreconditioner == nullptr)
                Gpreconditioner = SlauSolving::SSCMPreconditioner_base::gen(step_el.slausolverParameters.preconditioning);
        }
        double time;
        if(step_el.slausolverParameters.blocks)
        {
            Gpreconditioner3x3->initPortraitAndAllocateMemory(*Gmatrix3x3.p, time);
        }
        else
        {
            Gpreconditioner->initPortraitAndAllocateMemory(*Gmatrix.p, time);
        }
        if(nlInf.slau.debugMatrixesInformationNeedToWrite)
        {
            // Единожды вывести информацию о размерностях задачи и портреты
            PRINT1("preconditioner..");
            if(step_el.slausolverParameters.blocks)
            {
                Gpreconditioner3x3->saveProperties("3x3_slauPreconditionerInformation.txt");
                Gpreconditioner3x3->saveBMP("3x3_preconditioner.bmp", Gmatrix3x3.getMatrixSize());
            }
            else
            {
                Gpreconditioner->saveProperties("1x1_slauPreconditionerInformation.txt");
                Gpreconditioner->saveBMP("1x1_preconditioner.bmp", Gmatrix.getMatrixSize());
            }
        }
        // подготовка индексов для сборки глобальной матрицы
        for (size_t feInd = 0; feInd < grid->fe.size(); feInd++)
        {
            if(step_el.slausolverParameters.blocks)
            {
                (*fe)[feInd]->updateSSCM_indexes3x3(*Gmatrix3x3.p);
            }
            else
            {
                (*fe)[feInd]->updateSSCM_indexes1x1(*Gmatrix.p);
            }
        }
        nlInf.slau.matrixG_portraitChanged = false;
    }
    // обнуление глобальной матрицы
    if(nlInf.slau.matrixG_elementsChanged)
    {
        if(step_el.slausolverParameters.blocks)
        {
            MATR3x3 m0;
            m0.clear();
            (*Gmatrix3x3.e).fill(m0);
        }
        else
        {
            (*Gmatrix.e).fill(0);
        }
    }
    // обнуление глобального вектора правой части
    if(step_el.slausolverParameters.blocks)
    {
        VECTOR3 v0;
        v0.clear();
        for (size_t i = 0; i < b3x3.size(); i++)
        {
            b3x3[i] = v0;
        }
    }
    else
    {
        for (size_t i = 0; i < b.size(); i++)
        {
            b[i] = 0;
        }
    }
    // занесение локальных матриц и векторов в глобальные
    for (size_t feInd = 0; feInd < grid->fe.size(); feInd++) // индекс конечного элемента
    {
        // матрица
        if(nlInf.slau.matrixG_elementsChanged)
        {
            if(step_el.slausolverParameters.blocks)
            {
                (*fe)[feInd]->addLocalGM_to_global3x3(*Gmatrix3x3.e);
            }
            else
            {
                (*fe)[feInd]->addLocalGM_to_global1x1(*Gmatrix.e);
            }
        }
        // вектор
        if(step_el.slausolverParameters.blocks)
        {
            (*fe)[feInd]->addLocalb_to_global3x3(b3x3);
        }
        else
        {
            (*fe)[feInd]->addLocalb_to_global1x1(b);
        }
    }
}
void MechSolver::addContact(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &)
{
    MechIterInf &nlInf = stepInf.back().iterInf.back();
    IncForsesMode incForsesMode = step_el.incForsesMode;
    // сохранение копий диагональных элементов локальной узловой матрицы жёсткости
    if(nlInf.slau.matrixG_elementsChanged)
    {
        // (как минимум на 0-й итерации диагонали скопированы, они копируются не зависимо от статуса контакта)
        for(size_t CV_FE_RigidInd = 0; CV_FE_RigidInd < contactVertex_FE_Rigid.size(); CV_FE_RigidInd++)
        {
            if(step_el.slausolverParameters.blocks)
            {
                contactVertex_FE_Rigid[CV_FE_RigidInd].copyDiagStiffnessesFromG3x3(*Gmatrix3x3.e);
            }
            else
            {
                contactVertex_FE_Rigid[CV_FE_RigidInd].copyDiagStiffnessesFromG1x1(*Gmatrix.e);
            }
        }
    }
    // контактные локальные матрицы и векторы
    for(size_t CV_FE_RigidInd = 0; CV_FE_RigidInd < contactVertex_FE_Rigid.size(); CV_FE_RigidInd++)
    {
        contactVertex_FE_Rigid[CV_FE_RigidInd].updateLocalGMb(grid, rigidSurface, CC_FE_Rigid, bc1_state, grid->DOFs, nlInf.contact, incForsesMode);
    }
    // занесение контактных локальных матриц и векторов в глобальные
    for(size_t CV_FE_RigidInd = 0; CV_FE_RigidInd < contactVertex_FE_Rigid.size(); CV_FE_RigidInd++)
    {
        if(nlInf.slau.matrixG_elementsChanged)
        {
            if(step_el.slausolverParameters.blocks)
            {
                contactVertex_FE_Rigid[CV_FE_RigidInd].addLocalGM_to_global3x3(*Gmatrix3x3.e, *Gmatrix3x3.p);
            }
            else
            {
                contactVertex_FE_Rigid[CV_FE_RigidInd].addLocalGM_to_global1x1(*Gmatrix.e, *Gmatrix.p);
            }
        }
        if(step_el.slausolverParameters.blocks)
        {
            contactVertex_FE_Rigid[CV_FE_RigidInd].addLocalb_to_global3x3(b3x3);
        }
        else
        {
            contactVertex_FE_Rigid[CV_FE_RigidInd].addLocalb_to_global1x1(b);
        }
    }
}
void MechSolver::solveSLAU(MechGlobalStepParameters &step_el)
{
    MechIterInf &nlInf = stepInf.back().iterInf.back();
    // решение СЛАУ
    bool b_is_null = true;  // статус равенства нулю правой части (за исключением вектора невязок R)
    // проверка необходимости решения СЛАУ (если b = 0 или силы = 0)
    {
        // среди элементов правой части есть ненулевые?
        /*if(b_is_null)
        {
            if(step_el.slausolverParameters.blocks)
            {
                double norma_b = sqrt(Vector3ScalMul(b3x3, b3x3));
                if(norma_b != 0)
                    b_is_null = false;
            }
            else
            {
                double norma_b = sqrt(VectorScalMul(b, b));
                if(norma_b != 0)
                    b_is_null = false;
            }
        }*/
        // среди не контактных локальных векторов есть ненулевые?
        if(b_is_null)
        for(size_t bc2Ind = 0; bc2Ind < (*step_el.bc2).size(); bc2Ind++)
        {
            Grid::FiniteElementSurface &s = grid->FESurface[(*step_el.bc2)[bc2Ind].FEsurfaceInd];
            for (size_t faceInd = 0; faceInd < s.face.size(); faceInd++)
            {
                Grid::FEFace &face = s.face[faceInd];
                if(!(*fe)[face.feIndex]->localb_forces_is_null())
                {
                    b_is_null = false;
                    break;
                }
            }
        }
        // среди контактных сил есть ненулевые?
        if(b_is_null)
        for(size_t CV_FE_RigidInd = 0; CV_FE_RigidInd < contactVertex_FE_Rigid.size(); CV_FE_RigidInd++)
        {
            if(!contactVertex_FE_Rigid[CV_FE_RigidInd].localb_forces_is_null())
            {
                b_is_null = false;
                break;
            }
        }
    }
    //b_is_null = false;

retryWithPreconditioner:;
    if(step_el.slausolverParameters.type == SlauSolving::SolverType::Direct)
    {
        nlInf.slau.precChanged =
                nlInf.slau.precChanged ||
                nlInf.slau.matrixG_elementsChanged;
    }
    else
    {
        nlInf.slau.precChanged =
                nlInf.slau.precChanged ||
                nlInf.slau.slauTime > nlInf.slau.slauTimeSumWithoutPreconditioner/nlInf.slau.slauSolvedWithoutPreconditioner ||
                (nlInf.contact.iterationMode == IterationMode::Free && nlInf.contact.contactChangedNumber != 0)
                ;
    }

    nlInf.timeIterInf.solveSLAU_Gpreconditioner_init.begin();
    if(nlInf.slau.precChanged)
    {
        double timePreconditioner = 0;
        if(!b_is_null)
        {
            logger->sendString("Prec..");
            if(step_el.slausolverParameters.blocks)
            {
                Gpreconditioner3x3->updatePreconditioner(*Gmatrix3x3.e, nlInf.slau.GFirstStrChanged/3, timePreconditioner);
            }
            else
            {
                Gpreconditioner->updatePreconditioner(*Gmatrix.e, nlInf.slau.GFirstStrChanged, timePreconditioner);
            }
            nlInf.slau.GFirstStrChanged = matrixSize;
            nlInf.slau.precChanged = false;
        }
        nlInf.slau.slauTimeSumWithoutPreconditioner = timePreconditioner;
        nlInf.slau.slauSolvedWithoutPreconditioner = 0;
    }
    nlInf.timeIterInf.solveSLAU_Gpreconditioner_init.end();
    PRINT1("solv..");
    if(nlInf.slau.GsolverChanged)
    {
        if(step_el.slausolverParameters.blocks)
        {
            if(Gsolver3x3 == nullptr)
                Gsolver3x3 = SlauSolving::SSCM3x3Solver_base::gen(step_el.slausolverParameters.type);
            Gsolver3x3->init(Gmatrix3x3.getMatrixSize());
        }
        else
        {
            if(Gsolver == nullptr)
                Gsolver = SlauSolving::SSCMSolver_base::gen(step_el.slausolverParameters.type);
            Gsolver->init(Gmatrix.getMatrixSize());
        }
        nlInf.slau.GsolverChanged = false;
    }

    // решение СЛАУ
    //VectorMulc(dq, 0, dq);//## обнуление решения
    nlInf.timeIterInf.solveSLAU_Gsolver_solve.begin();
    if(!b_is_null)
    {
        if(step_el.slausolverParameters.blocks)
        {
            //SlauSolving::Vector1x1CopyToVector3x3(b, b3x3);
            Gsolver3x3->solve(Gmatrix3x3, b3x3,
                          Gpreconditioner3x3, dq3x3, step_el.slausolverParameters,
                          dq3x3, nlInf.slau.slauResidual, nlInf.slau.slauRelativeResidual, nlInf.slau.slauIterations, nlInf.slau.slauTime);
            // ##копия блочных векторов в не блочные, как буд то работает не блочный решатель
            SlauSolving::Vector3x3CopyToVector1x1(dq3x3, dq);
        }
        else
        {
            Gsolver->solve(Gmatrix, b,
                          Gpreconditioner, dq, step_el.slausolverParameters,
                          dq, nlInf.slau.slauResidual, nlInf.slau.slauRelativeResidual, nlInf.slau.slauIterations, nlInf.slau.slauTime);
        }
    }
    else
    {
        if(step_el.slausolverParameters.blocks)
        {
            for(size_t i = 0; i < dq3x3.size(); i++)
                dq3x3[i].clear();
            nlInf.slau.slauResidual = 0;
            nlInf.slau.slauIterations = 0;
            nlInf.slau.slauRelativeResidual = -1;
            nlInf.slau.slauResidual = 0;
            nlInf.slau.slauTime = 0;
            SlauSolving::Vector3x3CopyToVector1x1(dq3x3, dq);
        }
        else
        {
            for(size_t i = 0; i < dq.size(); i++)
                dq[i] = 0;
            nlInf.slau.slauResidual = 0;
            nlInf.slau.slauIterations = 0;
            nlInf.slau.slauRelativeResidual = -1;
            nlInf.slau.slauResidual = 0;
            nlInf.slau.slauTime = 0;
        }
    }
    nlInf.timeIterInf.solveSLAU_Gsolver_solve.end();

    PRINT(" iter=%6d residual=%.1le dresidual=%.1le\n",ARGS(nlInf.slau.slauIterations, nlInf.slau.slauResidual, nlInf.slau.slauRelativeResidual));
    nlInf.slau.slauTimeSumWithoutPreconditioner += nlInf.slau.slauTime;
    nlInf.slau.slauSolvedWithoutPreconditioner++;

    // СЛАУ не дорешивается
    if (nlInf.slau.slauIterations >= step_el.slausolverParameters.maxIter)
    {
        if(nlInf.slau.slauSolvedWithoutPreconditioner == 1 ||
           step_el.slausolverParameters.type == SlauSolving::SolverType::Direct)
        {
            // это первая итерация, предобусловливатель не помогает
            logger->sendString("WARNING: точность решения СЛАУ не достигнута!");
            if(step_el.controlMode == 0)
            {
                //logger->sendString(" Игнорирую\n");
                // => увеличиваем невязку (без суда и следствия)###
                nlInf.slau.precChanged = true;
                //step_el.slausolverParameters.eps *= 2;
                //PRINT1("\nУвеличиваю требуемую невязку\n");
                //        VectorMulc(dq, 0, dq);//## обнуление решения
                //        goto retryWithPreconditioner;
            }
            if(step_el.controlMode == 1)
            {
                // преждевременное завершение процесса
            }
        }
        else
        {
            // вызываем предобусловливатель
            nlInf.slau.precChanged = true;
            PRINT1("    resolv..");
            //        VectorMulc(dq, 0, dq);//## обнуление решения
            goto retryWithPreconditioner;
        }
    }
    // сравнение решения с решением на предыдущем шаге
    if (nlInf.iterNumber >= 1)
    {
        if(step_el.slausolverParameters.blocks)
        {
            SlauSolving::SSCM3x3calcResidual(Gmatrix3x3, b3x3, dqLastIter3x3, r3x3);
            double rr = Vector3ScalMul(r3x3, r3x3);
            if(rr == 0)
            {
                nlInf.slau.slauResidualWithLastq = 0;
            }
            else
            {
                nlInf.slau.slauResidualWithLastq = sqrt(rr / Vector3ScalMul(b3x3, b3x3));
            }
        }
        else
        {
            SlauSolving::SSCMsolveResidual(Gmatrix, b, dqLastIter, r);
            double rr = VectorScalMul(r, r);
            if(rr == 0)
            {
                nlInf.slau.slauResidualWithLastq = 0;
            }
            else
            {
                nlInf.slau.slauResidualWithLastq = sqrt(rr / VectorScalMul(b, b));
            }
        }
    }
    else
    {
        nlInf.slau.slauResidualWithLastq = -1;
    }
    // сохранение копии решения для следующей итерации
    if(step_el.slausolverParameters.blocks)
    {
        VectorCopy(dq3x3, dqLastIter3x3);
    }
    else
    {
        VectorCopy(dq, dqLastIter);
    }

    nlInf.slau.debugMatrixesInformationNeedToWrite = false;  // 1 раз вывели и хватит
}

void MechSolver::init(MechTask &mechTask, Threads::Logger_base *set_logger, Threads::SignalToThread_base *set_signal, const Grid::TimeLayers &tl)
{
    initLogging(set_logger, set_signal);
    *((MechTask *)this) = mechTask;    // копия входных данных(указателей)
    out = nullptr;
    if(!enabled) return;
    /*
    for(int i = (int)grid->vertex.size() - 1; i >= 0; i--)
    {
        vertexBusFuncInd->linear[i][1] = vertexBusFuncInd->DOF_count; vertexBusFuncInd->DOF_count++;
        vertexBusFuncInd->linear[i][0] = vertexBusFuncInd->DOF_count; vertexBusFuncInd->DOF_count++;
        vertexBusFuncInd->linear[i][2] = vertexBusFuncInd->DOF_count; vertexBusFuncInd->DOF_count++;
        int i1 = rand()%3;
        int i2 = rand()%3;
        std::swap(vertexBusFuncInd->linear[i][i1], vertexBusFuncInd->linear[i][i2]);
    }
    for(int i = (int)grid->vertex.size() - 1; i >= 0; i--)
    {
        int i1 = rand()%(int)grid->vertex.size();
        int i2 = rand()%(int)grid->vertex.size();
        std::swap(vertexBusFuncInd->linear[i1], vertexBusFuncInd->linear[i2]);
        //fprintf(stderr, "%d %d \n", i1, i2);
    }
    */
    matrixSize = grid->DOFs->DOF_count;
    // матрицы
    Gmatrix.init();
    Mmatrix.init();
    // векторы
    b.resize(matrixSize);          // правая часть
    dq.resize(matrixSize);         // решение СЛАУ на текущем шаге
    q.resize(matrixSize);          // перемещения в узлах на текущем шаге
    q1.resize(matrixSize);         // перемещения в узлах 1 шаг назад
    q2.resize(matrixSize);         // перемещения в узлах 2 шага назад
    q3.resize(matrixSize);         // перемещения в узлах 3 шага назад
    dqLastIter.resize(matrixSize); // перемещения в узлах на предыдущей итерации
    // вспомогательные векторы
    r.resize(matrixSize);          // для подсчета невязки
     Mddq.resize(matrixSize);
     qForM.resize(matrixSize);
    // вспомогательные данные для учёта первых краевых условий
    bc1_u0.resize(matrixSize);
    bc1_state.resize(matrixSize);
    // выходные данные
    out = new std::vector<MechOutGlobalStepData>;
    // подготовка к первому шагу по времени
    for (size_t i = 0; i < matrixSize; i++)
    {
        dq[i] = 0; // важно для начального приближения
        q[i] = 0;  // начальные перемещения принимаем нулевыми
        q1[i] = 0 - (*mechTask.V0)[i]*tl.dt + (*mechTask.dV0)[i]*SQR(tl.dt)/2;
        q2[i] = 0; // не имеет значения. т.к. после первого 3-точечного шага затрется
        q3[i] = 0; // не имеет значения. т.к. сразу затрется
        //R[i] = 0;
    }
    Gpreconditioner = nullptr;
    Gsolver = nullptr;


    Gpreconditioner3x3 = nullptr;
    Gsolver3x3 = nullptr;
    b3x3.resize(matrixSize/3);
    dq3x3.resize(matrixSize/3);
    dqLastIter3x3.resize(matrixSize/3); // перемещения в узлах на предыдущей итерации
    r3x3.resize(matrixSize/3);          // для подсчета невязки
    Gmatrix3x3.init();
    for (size_t i = 0; i < matrixSize/3; i++)
    {
        dq3x3[i].clear();
    }
    // заполнение таблиц значений на шаблонных кубах для всех реализаций MechFeData
    MechFeData_base::initConstantTables();
}
void MechSolver::prepareForStep(const MechGlobalStepParameters &step_el)
{
    stepInf.back().iterInf.push_back({});
    if(!enabled) return;
    // 0). Без разностной аппроксимации по времени
    if(step_el.timeMode == 0)
    {
    }
    // 1, 2). 3-слойные схемы
    if(step_el.timeMode == 1 || step_el.timeMode == 2)
    {
        VectorCopy(q1, q2);
        VectorCopy(q, q1);
    }
    // 3). 4-слойная схема
    if(step_el.timeMode == 3)
    {
        VectorCopy(q2, q3);
        VectorCopy(q1, q2);
        VectorCopy(q, q1);
    }
}
void MechSolver::prepareForIter(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl)
{
    if(!enabled) return;
    // подготовка к учёту первых краевых условий (они не будут меняться на глобальном шаге##)
    stepInf.back().prepareForIter_prepareBC1.begin();
    prepareBC1(step_el);
    stepInf.back().prepareForIter_prepareBC1.end();
    // подготовка значений функций для вторых краевых условий
    stepInf.back().prepareForIter_prepareBC2.begin();
    prepareBC2(step_el, tl);
    stepInf.back().prepareForIter_prepareBC2.end();
    // подготовка к учёту контакта
    stepInf.back().prepareForIter_prepareContact.begin();
    prepareContact(tl);
    stepInf.back().prepareForIter_prepareContact.end();
    // подготовка информации об итерациях
    MechIterInf &nlInf = stepInf.back().iterInf.back();
    nlInf.prepareForIter(step_el);
}
void MechSolver::tryStep(MechGlobalStepParameters &step_el, Grid::TimeLayers &tl,
                        MechIterationsGeneralResult &mechIterationsGeneralResult)
{
    if(!enabled) return;
    // итерации для учёта упруго-пластичности и механического контакта
    for(;;)
    {
        MechIterInf &nlInf = stepInf.back().iterInf.back();
        // построение матриц и вектора правой части
        nlInf.timeIterInf.buildGMb.begin();
        buildGMb(step_el, tl);
        nlInf.timeIterInf.buildGMb.end();
        // Построение матрицы G и вектора b из СЛАУ Gq = b
        // для разных схем аппроксимации по времени
        // (будем использовать ту же матрицу G для СЛАУ)
        if(step_el.timeMode == 0)
        {
        }
        // учёт первых краевых условий
        nlInf.timeIterInf.SSCMaddBoundaryCondition1.begin();
        //if(nlInf.slau.matrixG_elementsChanged)
        {
            logger->sendString("bc1..");
            if(step_el.slausolverParameters.blocks)
            {
                SlauSolving::SSCM3x3addBoundaryCondition1(Gmatrix3x3, b3x3, bc1_u0, bc1_state);
            }
            else
            {
                SlauSolving::SSCMaddBoundaryCondition1(Gmatrix, b, bc1_u0, bc1_state);
            }
        }
        nlInf.timeIterInf.SSCMaddBoundaryCondition1.end();
        // решение СЛАУ
        nlInf.timeIterInf.solveSLAU.begin();
        solveSLAU(step_el);
        nlInf.timeIterInf.solveSLAU.end();
        // расчет результатов итерации и общих результатов итерации
        PRINT1("res..");
        nlInf.timeIterInf.calcIterResults.begin();
        calcIterResults(step_el, tl);
        nlInf.timeIterInf.calcIterResults.end();
        // все счётчики просуммированы - расчитываем общий результат итерации
        nlInf.calcIterResults(step_el, signal->get_needToReleaseState(), matrixSize, mechIterationsGeneralResult);
        // отправка невязок логгеру
        {
            OutDataForLogger inf;
            inf.stepInf = &stepInf;
            inf.stepNumber = stepInf.back().stepNumber;
            inf.iterNumber = stepInf.back().iterInf.back().iterNumber - 1;// -1 т.к. уже увеличили на 1
            Threads::Message m;
            m.type = Threads::Message::Type::stepInf;
            m.size = sizeof inf;
            m.data = new char[m.size];
            memcpy((void *)m.data, (void *)&inf, m.size);
            logger->send(m);
        }
        // вывод на экран невязок после итерации
        PRINT("slauR = %le sigmaR = %le epsR = %le prelDChanged = %d DChanged = %d pl_nl = %d un = %d CChanged = %d CN = %d SN = %d endPointResidual = %le deltaF_residual = %le\n",
              ARGS(nlInf.slau.slauResidualWithLastq, nlInf.plastic.maxPlasticResidual.sigma, nlInf.plastic.maxPlasticResidual.eps, nlInf.plastic.preliminarily_GLocalChangedNumber, nlInf.plastic.GLocalChangedNumber, nlInf.plastic.nonlinearStateFENumber, nlInf.plastic.unloadStateFENumber,
                   nlInf.contact.contactChangedNumber, nlInf.contact.contactNumber, nlInf.contact.separationNumber, nlInf.contact.max_endPoint_residual, nlInf.contact.max_deltaF_residual));
        if(mechIterationsGeneralResult != MechIterationsGeneralResult::Continue_NeedMoreIterations)
        {
            break;
        }
        stepInf.back().iterInf.push_back(stepInf.back().iterInf.back());
    }
}
void MechSolver::finalizeStep(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl)
{
    if(!enabled) return;
    MechIterInf &nlInf = stepInf.back().iterInf.back();
    using namespace Operations;
    // добавление приращений напряжений, деформаций и т.д.
    stepInf.back().finalizeStep_fe.begin();
    for (size_t feInd = 0; feInd < grid->fe.size(); feInd++)    // индекс конечного элемента
    {
        (*fe)[feInd]->finalizeStep((*step_el.material)[grid->fe[feInd]->mi], step_el.movingGridMode);
    }
    stepInf.back().finalizeStep_fe.end();
    // добавление приращений перемещений в узлах сетки
    stepInf.back().finalizeStep_vertex.begin();
    for (size_t vertexIndex = 0; vertexIndex < grid->vertex.size(); vertexIndex++)
    {
        (*vertex)[vertexIndex].sum_du += (*vertex)[vertexIndex].du;
    }
    // добавление приращений перемещений в 6-гранных подобластях
    for (size_t feInd = 0; feInd < grid->fe.size(); feInd++)    // индекс конечного элемента
    {
        if(grid->fe[feInd]->get_FEType() == Grid::FEType::LinearHexagon_XFEM)
        {
            Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<Grid::FE_LinearHexagon_XFEM *>(grid->fe[feInd]);
            if(step_el.fixGrid == 0)
            {
                MechFeData_LinearXFEM_gen *XFEM_gen = dynamic_cast<MechFeData_LinearXFEM_gen *>((*fe)[feInd]);
                XFEM_gen->moveSubVertexes(grid, grid->DOFs, dq,
                                          grid->fe[feInd]);
            }
        }
    }
    stepInf.back().finalizeStep_vertex.end();
    // обновление сил реакции опоры
    stepInf.back().finalizeStep_contact.begin();
    for(size_t CV_FE_RigidInd = 0; CV_FE_RigidInd < contactVertex_FE_Rigid.size(); CV_FE_RigidInd++)
    {
        contactVertex_FE_Rigid[CV_FE_RigidInd].finalizeStep(vertex);
    }
    stepInf.back().finalizeStep_contact.end();

    // расчёт локальных векторов невязки (до сдвига сетки)
    {
        switch (step_el.incForsesMode)
        {
        case IncForsesMode::MinusIntegral:
        {
            // вычисляется в начале шага, после сдвига сетки
        }break;
        case IncForsesMode::bPlusR:
        {
            // вычисляется в конце шага, после всех итераций, до сдвига сетки
            buildLocalR(step_el, tl);
        }break;
        case IncForsesMode::IncrementP:
        {
        }break;
        }
    }

    // сдвиг сетки
    stepInf.back().finalizeStep_movingGrid.begin();
    if(step_el.fixGrid == 0)
    {
        // смещаются узлы сетки
        for (size_t vertexIndex = 0; vertexIndex < grid->vertex.size(); vertexIndex++)
        {
            grid->vertex[vertexIndex] += (*vertex)[vertexIndex].du;
        }
        for (size_t vertexIndex = 0; vertexIndex < grid->vertexForCurvature.size(); vertexIndex++)
        {
            grid->vertexForCurvature[vertexIndex] += (*vertexForCurvature)[vertexIndex].du;
        }
    }
    stepInf.back().finalizeStep_movingGrid.end();
    // приращение суммарных перемещений за всё время
    Vector1PlusCmulVector2(q, 1, dq, q);
    nlInf.finalizeStep();
}
void MechSolver::saveResultsOfStep(const bool saveDetailedInf, const int globalStepNumber, const int stepNumber, const Grid::TimeLayers &tl)
{
    if(!enabled) return;
    using namespace Integration;
    MechOutStepData s;
    //Integrator integrationCubeSource[3];
    //for(int it = (int)IntegrationType::Gauss2; it <= (int)IntegrationType::Gauss5; it++)
    //    integrationCubeSource[it].init3D((IntegrationType)it);
    /*
        nlInf.saveResultsOfStep(s);
        s.slauResidualWithLastq = slau.slauResidualWithLastq;
        s.maxSigmaResidual = plastic.maxSigmaResidual;
        s.maxEpsResidual = plastic.maxEpsResidual;
        s.iterationsNumber = iterNumber;
        s.slauResidual = slau.slauResidual;
        s.nonlinearStateFENumber = plastic.nonlinearStateFENumber;
        s.stepSolvingTime = (clock() - slau.stepSolvingTime) / (double)CLOCKS_PER_SEC;
    */

    if(out->empty())
        out->push_back({});
    out->back().step.push_back({});
    out->back().step.back().globalStepNumber = globalStepNumber;
    out->back().step.back().stepNumber = stepNumber;
    out->back().step.back().t0 = tl.t0;
    out->back().step.back().nlInf = stepInf.back();//#####
    //out->back().step.back().nlInf.print();


    if(saveDetailedInf)
    {
        out->back().fe.resize(grid->fe.size());
        for (int feInd = 0; feInd < (int)grid->fe.size(); feInd++)
        {
            MechOutFeData &outFeDataEl = (*out).back().fe[feInd];
            (*fe)[feInd]->saveResultsOfStep(dq, grid, grid->fe[feInd], outFeDataEl);
            out->back().fe.push_back(outFeDataEl);
        }
        for (size_t vertexIndex = 0; vertexIndex < grid->vertex.size(); vertexIndex++)
        {
            MechOutVertexData vertex_el;  // данные вывода для конечного элемента номер i
            vertex_el.p = grid->vertex[vertexIndex];
            vertex_el.sum_du = (*vertex)[vertexIndex].sum_du;
            vertex_el.contact = false;
            vertex_el.F_sum = VECTOR3_NULL;
            vertex_el.h = 0;
            out->back().vertex.push_back(vertex_el);
        }
        // контактная информация
        for(size_t CV_FE_RigidInd = 0; CV_FE_RigidInd < contactVertex_FE_Rigid.size(); CV_FE_RigidInd++)
        {
            contactVertex_FE_Rigid[CV_FE_RigidInd].saveResultsOfStep(out->back().vertex);
        }
        out->push_back({});
        // сохранение сетки и всего остального в файл

    }
}


extern FILE *f_debug;
// вычисление скоростей перемещений du - в каждой вершине
// и производных скоростей по пространству - в центре каждого конечного элемента
// после этого становятся доступны все приращения искомых величин
void MechSolver::calcIterResults(const MechGlobalStepParameters &step_el, const Grid::TimeLayers &tl)
{
    using namespace Integration;
    using namespace Elementary::Operations;
    MechIterInf &nlInf = stepInf.back().iterInf.back();
    // очистка информации о нелинейностях
    nlInf.clearCounts(matrixSize);
    // расчёт результатов на каждом конечном элементе
    PRINT1("ep..");

    fprintf(f_debug, "IterNumber = %d\n", nlInf.iterNumber);
    nlInf.timeIterInf.calcIterResults_fe.begin();
    for (int feInd = 0; feInd < (int)grid->fe.size(); feInd++)    // k - индекс конечного элемента
    {
        MechMaterialSource &m0 = (*step_el.material)[grid->fe[feInd]->mi];
        (*fe)[feInd]->calcIterResults1(dq, grid, grid->DOFs, tl, grid->fe[feInd], m0,
                                       vertex, vertexForCurvature, nlInf.plastic);
        //if(feInd >= (int)grid->fe.size() - 5)
        //    (*fe)[feInd]->debug(m0);
    }
    for (int feInd = 0; feInd < (int)grid->fe.size(); feInd++)    // k - индекс конечного элемента
    {
        MechMaterialSource &m0 = (*step_el.material)[grid->fe[feInd]->mi];
        (*fe)[feInd]->calcIterResults2(tl, m0,
                                        nlInf.plastic);
    }
    nlInf.timeIterInf.calcIterResults_fe.end();
    fprintf(f_debug, "\n");

    // контакт
    PRINT1("c..");
    nlInf.timeIterInf.calcIterResults_contact.begin();
    for(size_t CV_FE_RigidInd = 0; CV_FE_RigidInd < contactVertex_FE_Rigid.size(); CV_FE_RigidInd++)
    {
        contactVertex_FE_Rigid[CV_FE_RigidInd].calcIterResults1(grid, vertex, rigidSurface, CC_FE_Rigid, nlInf.contact.iterNumber, bc1_state);
    }
    //PRINT("[time = %le, %le, %d]", ARGS(nlInf.timeInf.solveIterResults_contact.duration, nlInf.timeInf.solveIterResults_contact.sumDurations, nlInf.timeInf.solveIterResults_contact.intervalsNumber));
    //PRINT1("c2..");
    for(size_t CV_FE_RigidInd = 0; CV_FE_RigidInd < contactVertex_FE_Rigid.size(); CV_FE_RigidInd++)
    {
        contactVertex_FE_Rigid[CV_FE_RigidInd].calcIterResults2(nlInf.contact, rigidSurface, CC_FE_Rigid);
    }
    nlInf.timeIterInf.calcIterResults_contact.end();
}



}   // namespace Solcvd_FE_Rigid
