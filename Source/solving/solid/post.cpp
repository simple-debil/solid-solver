#define _CRT_SECURE_NO_WARNINGS

#include <map>
#include <unordered_set>
#include <vector>
#include <deque>
#include "stdio.h"

#include "grid.h"
#include "solid_base.h"

namespace Post
{
using namespace Elementary;
void calcPressureByNodeForces(const Grid::FiniteElementSurface &contactSurface, const Grid::Grid3D &grid, const std::vector<VECTOR3> &vertexForce,
                               std::vector<VECTOR3> &vertexP)
{
    using namespace Integration;
    size_t matrixSize = grid.vertex.size() * 3; // размер матрицы СЛАУ
    SlauSolving::SSCM matrix;                           // матрица
    SlauSolving::SSCMBulder mBulder;                    // сборщик
    SlauSolving::SSCMPreconditioner_base *mPreconditioner;    // предобусловливатель
    SlauSolving::SSCMSolver_base *mSolver;                    // решатель СЛАУ
    Vector b;                                           // правая часть
    Vector res;                                         // коэффициенты разложения давления по базису
    SlauSolving::SolverParameters ssp;
    ssp.type = SlauSolving::SolverType::Direct;
    ssp.preconditioning = SlauSolving::Preconditioning::Profile_LLT;//SlauPreconditioning_none;
    //ssp.preconditioning = SlauSolving::Preconditioning::SSCM_LLT;//SlauPreconditioning_none;

    matrix.init();
    b.resize(matrixSize);
    res.resize(matrixSize);
    for (size_t i = 0; i < matrixSize; i++)
    {
        b[i] = 0;
        res[i] = 0;
    }
    mBulder.setMatrixSize(matrixSize);
    mBulder.start();

    // заполнение матрицы
    Integrator integrationFoursquare;
    integrationFoursquare.init2D(IntegrationType::Gauss3);

    for(size_t faceInd = 0; faceInd < contactSurface.face.size(); faceInd++)
    {
        const Grid::FEFace &faceEl = contactSurface.face[faceInd];
        const Grid::FE_base *feEl = grid.fe[faceEl.feIndex];
        int vi[8];
        feEl->getVertexIndexes(vi);

        /*int vi[4];
        feEl->getFaceVertexIndexes(faceEl.faceIndex, vi);
        if(vertexForce[vi[0]].abs() != 0 &&
           vertexForce[vi[1]].abs() != 0 &&
           vertexForce[vi[2]].abs() != 0 &&
           vertexForce[vi[3]].abs() != 0)*/
        {
        double basCube[8*Integration::Integration2D_Size_Max];
        double hexagonCoef[Integration::Integration2D_Size_Max];
        VECTOR3 hexagonNormal[Integration::Integration2D_Size_Max];

        feEl->calcBc2BasisFuncValues(grid.vertex, grid.vertexForCurvature, integrationFoursquare, faceEl.faceIndex,
                                      basCube, hexagonCoef, hexagonNormal);
        // добавка к матрице
        for (int m = 0; m < 8; m++)
            for (int n = 0; n < 8; n++)
            {
                int global_mx3 = vi[m]*3;
                int global_nx3 = vi[n]*3;
                for (int ik = 0; ik < 3; ik++)
                {
                    int &i = ik;
                    int &k = ik;
                    if(global_mx3 + i >= global_nx3 + k)
                    {
                        double E = 0;
                        for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
                        {
                            E += basCube[m*integrationFoursquare.size + valIndex] *
                                 basCube[n*integrationFoursquare.size + valIndex] *
                                 hexagonCoef[valIndex];
                        }
                        mBulder.addElement_not_null(
                                    E,
                                    global_mx3 + i,
                                    global_nx3 + k);
                    }
                }
                /*
                for (int i = 0; i < 3; i++)
                for (int k = 0; k < 3; k++)//можно обойтись одним циклом учитывая что i == k
                if(i == k && global_mx3 + i >= global_nx3 + k)
                {
                    double E = 0;
                    for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
                    {
                        E += basCube[m*integrationFoursquare.size + valIndex] *
                             basCube[n*integrationFoursquare.size + valIndex] *
                             hexagonCoef[valIndex];
                    }
                    mBulder.addElement_not_null(
                                E,
                                global_mx3 + i,
                                global_nx3 + k);
                }
                */
            }
        }
    }

    // заполнение вектора
    for(size_t global_m = 0; global_m < grid.vertex.size(); global_m++)
    {
        for(size_t i = 0; i < 3; i++)
        {
            b[global_m*3 + i] = vertexForce[global_m][i];
        }
    }

    mBulder.complete(matrix);

    // первые краевые
    for(size_t global_m = 0; global_m < grid.vertex.size(); global_m++)
    {
        for(size_t i = 0; i < 3; i++)
        {
            if(vertexForce[global_m][i] == 0)
                SlauSolving::SSCMaddBoundaryCondition1(matrix, b, global_m*3 + i, 0);
        }
    }
    // первые краевые
    /*for(size_t global_m = 0; global_m < grid->vertex.size(); global_m++)
    {
        if(vertexForce[global_m].abs() == 0)
        for(size_t i = 0; i < 3; i++)
        {
            SlauSolving::SSCMaddBoundaryCondition1(matrix, b, global_m*3 + i, 0);
        }
    }*/
    // искусственное заполнение диагонали пустых строк

    /*for(size_t i = 0; i < matrix.e->d.size(); i++)
    {
        if(matrix.e->d[i] == 0.)
            matrix.e->d[i] = 1;//1.e20;
    }*/


    // решение СЛАУ
    double residual;
    double relativeResidual;
    int iterations;
    double time;
    double timePreconditioner;
    mPreconditioner = SlauSolving::SSCMPreconditioner_base::gen(ssp.preconditioning);
    mSolver = SlauSolving::SSCMSolver_base::gen(ssp.type);
    mPreconditioner->bulid(matrix, 0, timePreconditioner);
    //mPreconditioner->saveBMP("solvePressureByNodeForces.bmp", 512);
    mSolver->init(matrix.getMatrixSize());
    mSolver->solve(matrix, b,
                  mPreconditioner, res, ssp,
                  res, residual, relativeResidual, iterations, time);

    // копия решения в виде векторов давлений для каждой вершины вместо коэффициентов разложения
    vertexP.resize(grid.vertex.size());
    for(size_t global_m = 0; global_m < grid.vertex.size(); global_m++)
    {
        for(size_t i = 0; i < 3; i++)
        {
            vertexP[global_m][i] = res[global_m*3 + i];
        }
    }

    // отладочный вывод
    /*
    matrix.p->saveInformation("___slauInformation_presure.txt");
    matrix.p->saveBMP("G_pressure.bmp", (int)matrix.getMatrixSize());
    mPreconditioner.saveInformation("___slauPreconditionerInformation_presure.txt");
    mPreconditioner.saveBMP("preconditioner_presure.bmp", (int)mPreconditioner.getMatrixSize());

    FILE *ff = fopen("0000000", "w");
    for(size_t i = 0; i < matrix.e->d.size(); i++)
    {
        fprintf(ff, "%le\t%le\t%le\n", matrix.e->d[i], res[i], b[i]);
    }
    fclose(ff);
    */

    matrix.release();
    mBulder.release();
    delete mPreconditioner;
    delete mSolver;
    b.clear();
    res.clear();
}

void calcContactArea(const Grid::FiniteElementSurface &contactSurface, const Grid::Grid3D &grid, const std::vector<Solid::MechOutVertexData> &vertexData,
                      std::vector<VECTOR3> &vertexP, double &S_true, double &S_geom, double &S_shamanstvo)
{
    std::vector<VECTOR3> vertexForce(grid.vertex.size());    // силы для каждого узла сетки
    // заполнение вектора контактных сил
    for(size_t vertexInd = 0; vertexInd < vertexData.size(); vertexInd++)
    {
        bool contacted = vertexData[vertexInd].contact;
        if(contacted)
        {
            vertexForce[vertexInd] = vertexData[vertexInd].F_sum;
        }
        else
        {
            vertexForce[vertexInd] = VECTOR3_NULL;
        }
    }
    // узловые давления по узловым силам
    calcPressureByNodeForces(contactSurface, grid, vertexForce,
                              vertexP);
    //contactSurface.solvePressureByNodeForces(vertexForce,
    //                                         vertexP);
    //return;
    std::vector<VECTOR3> feF; // суммарные контактные силы на грани каждого КЭ(интегрирование vertexP на площадках)
    std::vector<double> feS; // площади контактных граней каждого КЭ(интегрирование vertexP на площадках)
    std::vector<int> feContactVertexesCount; // количество контактных узлов для каждого КЭ
    feF.resize(grid.fe.size());
    feS.resize(grid.fe.size());
    feContactVertexesCount.resize(grid.fe.size());
    for(size_t i = 0; i < grid.fe.size(); i++)
    {
        feF[i] = VECTOR3_NULL;
        feS[i] = 0;
        feContactVertexesCount[i] = 0;
    }
    // каждый пограничный узел считаем 1 раз
    //std::vector<bool> vertexCounted;
    //vertexCounted.resize(vertexData.size());
    //for(size_t vertexInd = 0; vertexInd < vertexData.size(); vertexInd++)
    //    vertexCounted[vertexInd] = false;
    // индексы КЭ, которым принадлежит узел
    std::vector<std::set<int>> vertex_nearFeIndexes; // vertex_nearFeIndexes[i] - множество индексов КЭ, которые содержат узел i
    vertex_nearFeIndexes.resize(vertexData.size());

    // подготовка feF, feS, feContactVertexesCount
    for(size_t faceInd = 0; faceInd < contactSurface.face.size(); faceInd++)
    {
        const Grid::FEFace &faceEl = contactSurface.face[faceInd];
        const Grid::FE_base *feEl = grid.fe[faceEl.feIndex];
        int vi[4];
        int vi_local[4];
        feEl->getFaceVertexIndexes(faceEl.faceIndex, vi, vi_local);
        int local_contactPointsCount = 0;
        for(int i = 0; i < 4; i++)
        {
            if(vertexData[vi[i]].contact)
                local_contactPointsCount++;
        }
        feContactVertexesCount[faceEl.feIndex] = local_contactPointsCount;
        // подсчёт площади контактной грани КЭ и среднего давления на площадке
        //if(local_contactPointsCount == 4)
        //if(local_contactPointsCount >= 1)
        {
            double fe_S = 0; // площадь грани КЭ
            VECTOR3 fe_F = VECTOR3_NULL; // суммарная сила реакции опоры на грани КЭ
            double basCube[8*Integration::Integration2D_Size_Max];
            double hexagonCoef[Integration::Integration2D_Size_Max];
            VECTOR3 hexagonNormal[Integration::Integration2D_Size_Max];
            Integration::Integrator integrationFoursquare;
            integrationFoursquare.init2D(Integration::IntegrationType::Gauss3);
            feEl->calcBc2BasisFuncValues(grid.vertex, grid.vertexForCurvature, integrationFoursquare, faceEl.faceIndex,
                                          basCube, hexagonCoef, hexagonNormal);
            VECTOR3 fe_pressure[Integration::Integration2D_Size_Max]; // давление в точках интегрирования
            for(int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
            {
                fe_pressure[valIndex] = VECTOR3_NULL;
            }
            for(int i = 0; i < 4; i++)
            {
                int m = vi_local[i];        // локальный индекс вершины
                int vertexIndex = vi[i];    // глобальный индекс вершины
                for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
                {
                    fe_S += 1 *
                         basCube[m*integrationFoursquare.size + valIndex] *
                         hexagonCoef[valIndex];
                    fe_pressure[valIndex] += basCube[m*integrationFoursquare.size + valIndex]*vertexP[vertexIndex];
                }
            }
            for(int i = 0; i < 4; i++)
            {
                int m = vi_local[i];
                for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
                {
                    fe_F += fe_pressure[valIndex] *
                         basCube[m*integrationFoursquare.size + valIndex] *
                         hexagonCoef[valIndex];
                }
            }
            //S += fe_S*4.;
            //S += fe_S*4.*local_contactPointsCount/4.;
            feF[faceEl.feIndex] = fe_F;// суммарная сила на площадке
            feS[faceEl.feIndex] = fe_S;// площадь площадки
        }
        // заполнение списков индексов КЭ, содержащих заданную вершину
        for(int i = 0; i < 4; i++)
        {
            vertex_nearFeIndexes[vi[i]].insert(faceEl.feIndex);
        }
    }

    // подсчёт суммарной площади контакта
    S_true = 0;
    S_geom = 0;
    S_shamanstvo = 0;
    for(size_t faceInd = 0; faceInd < contactSurface.face.size(); faceInd++)
    {
        const Grid::FEFace &faceEl = contactSurface.face[faceInd];
        const Grid::FE_base *feEl = grid.fe[faceEl.feIndex];
        int vi[4];
        int vi_local[4];
        feEl->getFaceVertexIndexes(faceEl.faceIndex, vi, vi_local);
        int local_contactPointsCount = feContactVertexesCount[faceEl.feIndex];
        if(local_contactPointsCount == 4)
        {
            // полностью контактная площадка
            S_true += feS[faceEl.feIndex];
            S_geom += feS[faceEl.feIndex];
            S_shamanstvo += feS[faceEl.feIndex];
        }
        if(local_contactPointsCount >= 1 && local_contactPointsCount <= 3)
        {
            // частично контактная площадка
            // 1) S_true. Учёт количества контактных узлов на граничных ячейках
            {
                S_true += feS[faceEl.feIndex] * local_contactPointsCount / 4;
            }
            // 2) S_geom. Учёт пересечения жёсткой поверхности с площадкой
            // 3) S_shamanstvo. Шаманство с давлением: учёт величины давления на граничных ячейках
            {
                // поиск соседних граней, которые полностью контактные
                std::set<int> feind;    // множество индексов соседних КЭ, которые полностью контактные
                int nearFeCount = 0;    // количество соседних полностью контактных КЭ
                double nearFeSumF = 0;  // суммарная сила на соседних полностью контактных КЭ
                double nearFeSumS = 0;  // суммарная площадь соседних полностью контактных КЭ
                for(int i = 0; i < 4; i++)
                {
                    for (std::set<int>::iterator it=vertex_nearFeIndexes[vi[i]].begin(); it!=vertex_nearFeIndexes[vi[i]].end(); ++it)
                    {
                        // *it - индекс соседнего КЭ
                        if(feind.insert(*it).second)
                        {
                            // вставка успешна, нашли очередной соседний КЭ
                            if(feContactVertexesCount[*it] == 4)
                            if(*it != faceEl.feIndex)
                            {
                                nearFeCount++;
                                nearFeSumF += feF[*it].abs();
                                nearFeSumS += feS[*it];
                            }
                        }
                    }
                }
                // расчёт среднего давления на соседних полностью контактных КЭ
                if(nearFeCount != 0)
                {
                    double fe_near_averageP = nearFeSumF / nearFeSumS;    // среднее давление на полностью контактных соседник КЭ
                    double fe_P = feF[faceEl.feIndex].abs() / feS[faceEl.feIndex];    // среднее давление на данной площадке
                    double coef = MIN(fe_P / fe_near_averageP, 1.);
                    S_shamanstvo += feS[faceEl.feIndex] * coef;
                }
            }
        }
    }
}

}   // namespace Post
