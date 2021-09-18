#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"

#include "solidContact.h"

namespace Solid
{

void MechContactIterInf::clearCounts(const int matrixSize)
{
    contactNumber = 0;
    separationNumber = 0;
    contactChangedNumber = 0;
    max_endPoint_residual = 0;
    max_deltaF_residual = 0;
    GFirstStrChanged = matrixSize;
}
void MechContactIterInf::addContactVertexIterResultsToCounts(const ContactSurfaceVertexData_FE_Rigid &cvd, const int iterNumber, const bool constantNormal)
{
    if(cvd.it_new.contactNewChanged)
        contactChangedNumber++;
    if(cvd.it_new.contactNew)
        contactNumber++;
    if(cvd.it_new.separation)
        separationNumber++;
    if(cvd.endPoint_residual > max_endPoint_residual)
        max_endPoint_residual = cvd.endPoint_residual;
    if(cvd.deltaF_residual > max_deltaF_residual)
        max_deltaF_residual = cvd.deltaF_residual;

    //if(cvd.it_new.contactNew && cvd.endPoint_residual > max_endPoint_residual)
    //    max_endPoint_residual = cvd.endPoint_residual;
    //if(cvd.it_new.contactNew && cvd.deltaF_residual > max_deltaF_residual)
    //    max_deltaF_residual = cvd.deltaF_residual;
    if(cvd.it_new.GLocalChanged)
    {
        int global_mi = cvd.vertexIndex*3 + 0;
        if(global_mi < GFirstStrChanged)
            GFirstStrChanged = global_mi;
    }
}
void MechContactIterInf::calcIterResults(const MechGlobalStepParameters &step_el)
{
    if(max_endPoint_residual <= step_el.contactEndPointResidualLimit &&
       max_deltaF_residual <= step_el.contactDeltaFResidualLimit)
    {
        accuracyAchieved = true;
        improvingAccuracy = false;//## дожимания нет
    }
    else
    {
        accuracyAchieved = false;
        improvingAccuracy = false;
    }
}

void ContactSurfaceVertexData_FE_Rigid::init(const size_t set_vertexIndex, const size_t set_si)
{
    vertexIndex = set_vertexIndex;
    si = set_si;
    it.init();
    it_new.init();
    localNodeMatrixDiag.clear();
    ssd1.index = -1;     // данные о начальном приближении отсутствуют
    ssd2.index = -1;
    h = 0;
    endPoint_residual = 0;
    deltaF_residual = 0;
    contact = false;
    F_sum.clear();
    noContactRadius = -1;   // расстояние до контакта не определено
}
void ContactSurfaceVertexData_FE_Rigid::updateIndexes(const Grid::GlobalDOFs *DOFs)
{
    DOFindex = DOFs->findDOFIndex_Lagr1(vertexIndex);
}
void ContactSurfaceVertexData_FE_Rigid::prepareForIter()
{
}
void ContactSurfaceVertexData_FE_Rigid::updateLocalGMb(const Grid::Grid3D *grid, const std::vector<Grid::Surface_base *> *rigidSurface, const std::vector<ContactCondition_FE_Rigid> *CC_FE_Rigid, const std::vector<bool> &bc1_state, const Grid::GlobalDOFs *DOFs, const MechContactIterInf &nlInf, const IncForsesMode incForsesMode)
{
    updateIndexes(DOFs);
    if(nlInf.iterNumber == 0)
    {
        // для подвижных контактных поверхностей обновляются предполагаемые точки контакта и номали после сдвига сетки и поверхностей
        if(it.contactNew)// = contact
        {
            // нормаль в исходной точке
            POINT3 x1 = grid->vertex[vertexIndex];  // исходное положение узла
            POINT3 x2 = x1;                         // положение узла после перемещения
            POINT3 x1_nearestPoint;
            POINT3 x1_normal;
            int x1_side;
            bool x1_onBorder;
            (*rigidSurface)[(*CC_FE_Rigid)[si].RigidSurfaceInd]->findNearestPoint(x1, x2 - x1, ssd2,
                                                                                  ssd2, x1_nearestPoint, x1_normal, x1_side, x1_onBorder);
            it.updateEndPoint(x1_nearestPoint, x1);
            it.updateNormal(x1_normal, bc1_state, DOFindex, localNodeMatrixDiag, (*CC_FE_Rigid)[si].w_stiffness);
            it.GLocalChanged = true;
        }
        it.F = F_sum;
    }
    else
    {
        if(nlInf.iterationMode == IterationMode::Free)
        {
            // обновляем данные для итераций
            it = it_new;
        }
        if(nlInf.iterationMode == IterationMode::ReadOnly)
        {
            // итерации приостановлены, пока не переключимся в свободный режим
        }
    }

    // контакт
    if(it.contactNew)
    {
        if(it.GLocalChanged)
        {
            for(int k = 0; k < 3; k++)
            {
                for(int l = 0; l < 3; l++)
                {
                    G_local[k][l] = it.stiffness*it.normal[k]*it.normal[l];
                }
            }
        }

        VECTOR3 F_k;
        switch ((*CC_FE_Rigid)[si].method)
        {
        case ContactType::AugmentedLagrange:
        {
            // контактный узел стремится быть точно на поверхности
            F_k = it.F + it.stiffness*it.a;
        }break;
        case ContactType::Penalty:
        {
            // контактный узел может немного вдавливаться в поверхность
            F_k = F_sum + it.stiffness*it.a;
        }break;
        }

        for(int k = 0; k < 3; k++)
        {
            double E = 0;
            for(int l = 0; l < 3; l++)
            {
                E += F_k[l]*it.normal[k]*it.normal[l];
            }
            switch (incForsesMode)
            {
            case IncForsesMode::MinusIntegral:
            case IncForsesMode::bPlusR:
            {
                // полные силы
                b_local[k] = E;
            }break;
            case IncForsesMode::IncrementP:
            {
                // приращения сил
                b_local[k] = E - F_sum[k];
            }break;
            }
        }
    }
    // отлипание
    if(contact && !it.contactNew)
    {
        for(int k = 0; k < 3; k++)
        {
            switch (incForsesMode)
            {
            case IncForsesMode::MinusIntegral:
            case IncForsesMode::bPlusR:
            {
                // полные силы
                b_local[k] = 0;
            }break;
            case IncForsesMode::IncrementP:
            {
                // приращения сил
                b_local[k] = -F_sum[k];
            }break;
            }
        }
    }
}

void ContactSurfaceVertexData_FE_Rigid::addLocalGM_to_global1x1(SlauSolving::SSCMElements &elements, SlauSolving::SSCMPortrait &portrait)
{
    // контакт
    if(it.contactNew)
    {
        for(int k = 0; k < 3; k++)
        {
            for(int l = 0; l < 3; l++)
            {
                size_t global_mk = DOFindex[k];
                size_t global_ml = DOFindex[l];
                if(global_mk > global_ml)
                {
                    size_t SSCM_index = portrait.findSorted(global_mk, global_ml);
                    elements.a[SSCM_index] += G_local[k][l];
                }
                if(global_mk == global_ml)
                    elements.d[global_mk] += G_local[k][l];
            }
        }
    }
}
void ContactSurfaceVertexData_FE_Rigid::addLocalb_to_global1x1(Vector &b)
{
    if(it.contactNew                      // контакт
       || (contact && !it.contactNew))    // или отлипание
    {
        // в случае контакта добавляется контактный вектор
        // в случае отлипания добавляется вектор, зануляющий суммарную силу
        for(int k = 0; k < 3; k++)
        {
            size_t global_mk = DOFindex[k];
            b[global_mk] += b_local[k];
        }
    }
}
void ContactSurfaceVertexData_FE_Rigid::copyDiagStiffnessesFromG1x1(const SlauSolving::SSCMElements &Gelements)
{
    for(int k = 0; k < 3; k++)
    {
        size_t global_mk = DOFindex[k];
        localNodeMatrixDiag[k] = Gelements.d[global_mk];
    }
}

void ContactSurfaceVertexData_FE_Rigid::addLocalGM_to_global3x3(SlauSolving::SSCM3x3Elements &elements, SlauSolving::SSCMPortrait &)
{
    // контакт
    if(it.contactNew)
    {
        int i = DOFindex[0]/3;
        MATR3x3 &el = elements.d[i];
        for(int k = 0; k < 3; k++)
        {
            for(int l = k; l < 3; l++)
            {
                if(l != k)
                {
                    el.m[k][l] += G_local[k][l];
                    el.m[l][k] += G_local[k][l];
                }
                else
                    el.m[k][l] += G_local[k][l];
            }
        }
    }
}
void ContactSurfaceVertexData_FE_Rigid::addLocalb_to_global3x3(Vector3 &b)
{
    if(it.contactNew                      // контакт
       || (contact && !it.contactNew))    // или отлипание
    {
        // в случае контакта добавляется контактный вектор
        // в случае отлипания добавляется вектор, зануляющий суммарную силу
        int i = DOFindex[0]/3;
        VECTOR3 &el = b[i];
        for(int k = 0; k < 3; k++)
        {
            el[k] += b_local[k];
        }
    }
}
void ContactSurfaceVertexData_FE_Rigid::copyDiagStiffnessesFromG3x3(const SlauSolving::SSCM3x3Elements &Gelements)
{
    int i = DOFindex[0]/3;
    const MATR3x3 &el = Gelements.d[i];
    for(int k = 0; k < 3; k++)
    {
        localNodeMatrixDiag[k] = el.m[k][k];
    }
}

bool ContactSurfaceVertexData_FE_Rigid::localb_forces_is_null()
{
    if(it.contactNew                      // контакт
       || (contact && !it.contactNew))    // или отлипание
    {
        for(int k = 0; k < 3; k++)
        {
            if(b_local[k] != 0)
                return false;
        }
    }
    return true;
}


void ContactSurfaceVertexData_FE_Rigid::calcIterResults1(const Grid::Grid3D *grid, const std::vector<MechVertexData> *vertex, const std::vector<Grid::Surface_base *> *rigidSurface, const std::vector<ContactCondition_FE_Rigid> *CC_FE_Rigid, const int iterNumber, const std::vector<bool> &bc1_state)
{
    bool constantNormal = (*CC_FE_Rigid)[si].constantNormal;
    POINT3 x1 = grid->vertex[vertexIndex];          // исходное положение узла
    POINT3 x2 = x1 + (*vertex)[vertexIndex].du;     // положение узла после перемещения
    // информация о ближайшей к x2 точке поверхности
    POINT3 x2_nearestPoint;
    POINT3 x2_normal;
    int x2_side;
    bool x2_onBorder;

    it_new = it;    // будем обновлять it_new, оставив неизменным исходные значения it

    // ## случай x2_onBorder = true не рассматривается
    // проверка наличия контакта и расчёт предполагаемой точки контакта
    if(it.contactNew)
    {
        // на предыдущей итерации контакт был (т.е. контактная матрица собрана)
        // проверяем точку x2
        (*rigidSurface)[(*CC_FE_Rigid)[si].RigidSurfaceInd]->findNearestPoint(x2, x2 - x1, ssd2,
                                                                              ssd2, x2_nearestPoint, x2_normal, x2_side, x2_onBorder);
        // сила зависит от x2, it.endPoint и it.normal
        VECTOR3 d = ((x2 - it.endPoint)*it.normal)*it.normal;

        switch ((*CC_FE_Rigid)[si].method)
        {
        case ContactType::AugmentedLagrange:
        {
            // контактный узел стремится быть точно на поверхности
            it_new.F = ((it.F + it.stiffness*(it.endPoint - x2))*it.normal)*it.normal;
        }break;
        case ContactType::Penalty:
        {
            // контактный узел может немного вдавливаться в поверхность
            it_new.F = ((F_sum + it.stiffness*(it.endPoint - x2))*it.normal)*it.normal;
        }break;
        }

        if(d*it.normal >= 0)
            h = +d.abs();
        else
            h = -d.abs();
        if(x2_side == -1)
        {
            // точка x2 находится с внутренней стороны поверхности
            // сохраняется контакт
            it_new.contactNewChanged = false; // контакт был - нет изменения состояния
            it_new.contactNew = true;         // контакт есть
            it_new.separation = false;        // отлипания нет
            // используем x2
            it_new.updateEndPoint(x2_nearestPoint, x1);
            if(constantNormal)
            {
                it_new.GLocalChanged = false;
            }
            else
            {
                it_new.updateNormal(x2_normal, bc1_state, DOFindex, localNodeMatrixDiag, (*CC_FE_Rigid)[si].w_stiffness);
                it_new.GLocalChanged = true;
            }
        }
        else
        {
            // точка x2 находится с внешней стороны поверхности
            // проверка на отлипание
            if(it_new.F*it.normal < 0) // строго меньше
            {
                // отлипание
                it_new.contactNewChanged = true; // изменение состояния: контакт -> нет контакта
                it_new.contactNew = false;       // теперь нет контакта
                if(contact)
                {
                    it_new.separation = true;        // отлипание есть, т.к. в конце прошлого шага контакт был
                }
                else
                {
                    it_new.separation = false;       // отлипания нет, т.к. в конце прошлого шага контакта не было
                }
                it_new.GLocalChanged = true;   // матрица была, а станет нулевой
            }
            else
            {
                // сохраняется контакт
                it_new.contactNewChanged = false; // контакт был - нет изменения состояния
                it_new.contactNew = true;         // контакт есть
                it_new.separation = false;        // отлипания нет
                // используем x2
                it_new.updateEndPoint(x2_nearestPoint, x1);
                if(constantNormal)
                {
                    it_new.GLocalChanged = false;
                }
                else
                {
                    it_new.updateNormal(x2_normal, bc1_state, DOFindex, localNodeMatrixDiag, (*CC_FE_Rigid)[si].w_stiffness);
                    it_new.GLocalChanged = true;
                }
            }
        }
    }
    else
    {
        // на предыдущей итерации контакта не было (т.е. контактная матрица не собрана)
        it_new.F = F_sum;
        if(contact)
        {
            // в конце прошлого шага контакт был (т.е. x1 находится на поверхности)
            // проверяем точку x2
            // ##требуется только x2_side
            (*rigidSurface)[(*CC_FE_Rigid)[si].RigidSurfaceInd]->findNearestPoint(x2, x2 - x1, ssd2,
                                                                                  ssd2, x2_nearestPoint, x2_normal, x2_side, x2_onBorder);
            if(x2_side == -1)
            {
                // точка x2 находится с внутренней стороны поверхности
                // контакт
                it_new.contactNewChanged = true;  // изменение состояния: нет контакта -> контакт
                it_new.contactNew = true;         // контакт есть
                it_new.separation = false;        // отлипания нет
                // используем x1
                // информация о ближайшей к x1 точке поверхности
                POINT3 x1_nearestPoint;
                POINT3 x1_normal;
                int x1_side;
                bool x1_onBorder;
                (*rigidSurface)[(*CC_FE_Rigid)[si].RigidSurfaceInd]->findNearestPoint(x1, x1 - x1, ssd1,
                                                                                      ssd1, x1_nearestPoint, x1_normal, x1_side, x1_onBorder);
                it_new.updateEndPoint(x1_nearestPoint, x1);
                it_new.updateNormal(x1_normal, bc1_state, DOFindex, localNodeMatrixDiag, (*CC_FE_Rigid)[si].w_stiffness);
                it_new.GLocalChanged = true;
            }
            else
            {
                // точка x2 находится с внешней стороны поверхности
                // сохраняется отсутствие контакта
                // т.е. происходит отлипание
                it_new.contactNewChanged = false;  // отсутствие контакта было - нет изменения состояния
                it_new.contactNew = false;         // контакта нет
                it_new.separation = true;          // отлипание есть
                it_new.GLocalChanged = false;
            }
        }
        else
        {
            // в конце прошлого шага контакта не было
            // проверка на первый случай контакта
            // ищем точку пересечения отрезка x1-x2 с поверхностью
            // информация о пересечении x1-x2 с поверхностью
            POINT3 x1x2_intersectionPoint;
            POINT3 x1x2_normal;
            int x1x2_side;
            bool x1x2_found;

            if((*CC_FE_Rigid)[si].noContactRadiusOptimization)
            {
                if(noContactRadius <= 0)
                {
                    // расчёт радиуса отсутствия контакта
                    // информация о ближайшей к x1 точке поверхности
                    POINT3 x1_nearestPoint;
                    POINT3 x1_normal;
                    int x1_side;
                    bool x1_onBorder;
                    (*rigidSurface)[(*CC_FE_Rigid)[si].RigidSurfaceInd]->findNearestPoint(x1, x2 - x1, ssd1,
                                                                                          ssd1, x1_nearestPoint, x1_normal, x1_side, x1_onBorder);
                    noContactRadius = (x1_nearestPoint - x2).abs();
                }

                if((*vertex)[vertexIndex].du.abs() < noContactRadius)
                {
                    // если не вылезли за радиус отсутствия контакта, то его точно нет
                    x1x2_found = false;
                }
                else
                {
                    // вылезли за радиус отсутствия контакта, необходима проверка пересечения
                    (*rigidSurface)[(*CC_FE_Rigid)[si].RigidSurfaceInd]->findIntersection(x1, x2, ssd1, ssd2,
                                                                                          ssd1, ssd2, x1x2_intersectionPoint, x1x2_normal, x1x2_side, x1x2_found);
                }
            }
            else
            {
                (*rigidSurface)[(*CC_FE_Rigid)[si].RigidSurfaceInd]->findIntersection(x1, x2, ssd1, ssd2,
                                                                                      ssd1, ssd2, x1x2_intersectionPoint, x1x2_normal, x1x2_side, x1x2_found);
            }
            // проверяем точку x2
            //(*rigidSurface)[(*CC_FE_Rigid)[si].RigidSurfaceInd]->findNearestPoint(x2, x2 - x1, ssd,
            //                                                                      ssd, x2_nearestPoint, x2_normal, x2_side, x2_onBorder);

            if(x1x2_found)
            {
                // пересечение найдено
                // контакт (впервые)
                it_new.contactNewChanged = true;  // изменение состояния: нет контакта -> контакт
                it_new.contactNew = true;         // контакт есть
                it_new.separation = false;        // отлипания нет
                // используем точку пересечения с поверхностью
                it_new.updateEndPoint(x1x2_intersectionPoint, x1);
                it_new.updateNormal(x1x2_normal, bc1_state, DOFindex, localNodeMatrixDiag, (*CC_FE_Rigid)[si].w_stiffness);
                it_new.GLocalChanged = true;
            }
            else
            {
                // пересечение не найдено
                // сохраняется отсутствие контакта
                it_new.contactNewChanged = false;  // отсутствие контакта было - нет изменения состояния
                it_new.contactNew = false;         // контакта нет
                it_new.separation = false;         // отлипаня нет
                it_new.GLocalChanged = false;
            }
        }
    }

    // расчёт погрешности
    if(it_new.contactNewChanged)
    {
        // изменился статус контакта
        endPoint_residual = 1;  // погрешность не определена
        deltaF_residual = 1;    // погрешность не определена
    }
    else
    {
        if(it_new.contactNew)
        {
            // расчёт невязок
            endPoint_residual = (it_new.endPoint - it.endPoint).abs();
            VECTOR3 deltaF = it_new.F - it.F;
            double F_sum_abs = (F_sum.abs() + (it_new.F.abs() + it.F.abs())/2)/2;
            if(F_sum_abs == 0)
                deltaF_residual = 0;
            else
                deltaF_residual = deltaF.abs() / F_sum_abs;
            // а что если F_sum, deltaF и deltaF_new все маленькие одновременно?
        }
        else
        {
            // в случае отсутствия контакта невязку принимаем за 0, при отлипании тоже
            endPoint_residual = 0;
            deltaF_residual = 0;
        }
    }
}
void ContactSurfaceVertexData_FE_Rigid::calcIterResults2(MechContactIterInf &nlInf, const std::vector<Grid::Surface_base *> *rigidSurface, const std::vector<ContactCondition_FE_Rigid> *CC_FE_Rigid)
{
    nlInf.addContactVertexIterResultsToCounts(*this, nlInf.iterNumber, (*CC_FE_Rigid)[si].constantNormal);
}
void ContactSurfaceVertexData_FE_Rigid::finalizeStep(const std::vector<MechVertexData> *vertex)
{
    // контакт
    if(it.contactNew)
    {
        contact = true;
        F_sum = it_new.F;   // используем следующее приближения, т.к. для новой итерации расчитана сила с учётом влияния контактной жёсткости на текущей итерации
    }
    // отлипание
    if(contact && !it.contactNew)
    {
        contact = false;
        F_sum.clear();
    }
    // x2->x1
    ssd1 = ssd2;
    // радиус возможного контакта уменьшается (если эта оптимизация выключена, то noContactRadius всегда < 0)
    if(noContactRadius > 0)
        noContactRadius -= (*vertex)[vertexIndex].du.abs();
}
void ContactSurfaceVertexData_FE_Rigid::saveResultsOfStep(std::vector<MechOutVertexData> &outVertex)
{
    outVertex[vertexIndex].contact = contact;
    outVertex[vertexIndex].F_sum = F_sum;
    outVertex[vertexIndex].h = h;
    outVertex[vertexIndex].stiffness = it.stiffness;
}


}   // namespace Solid
