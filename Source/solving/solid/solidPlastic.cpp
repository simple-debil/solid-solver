#define _CRT_SECURE_NO_WARNINGS

#include "stdio.h"
#include <set>

#include "solidPlastic.h"
#include "fem.h"

namespace Solid
{
// константные таблицы
// Гаусс-2
double MechFeData_base::Gauss2_integrationwSource[8];
double MechFeData_base::Gauss2_basCubeSource[8*8];
VECTOR3 MechFeData_base::Gauss2_dLinearBasCubeSource[8*8];
VECTOR3 MechFeData_base::Gauss2_dQuadraticBasCubeSource[8*27];
// Гаусс-3
double MechFeData_base::Gauss3_integrationwSource[27];
double MechFeData_base::Gauss3_basCubeSource[27*8];
VECTOR3 MechFeData_base::Gauss3_dLinearBasCubeSource[27*8];
VECTOR3 MechFeData_base::Gauss3_dQuadraticBasCubeSource[27*27];

void MechFeData_base::initConstantTables()
{
    Fem::calcCubeLinearBasisFuncValues(Integration::IntegrationType::Gauss2,
                                        Gauss2_integrationwSource, Gauss2_basCubeSource, Gauss2_dLinearBasCubeSource);
    Fem::calcCubeQuadraticBasisFuncValues(Integration::IntegrationType::Gauss2,
                                           Gauss2_dQuadraticBasCubeSource);
    Fem::calcCubeLinearBasisFuncValues(Integration::IntegrationType::Gauss3,
                                        Gauss3_integrationwSource, Gauss3_basCubeSource, Gauss3_dLinearBasCubeSource);
    Fem::calcCubeQuadraticBasisFuncValues(Integration::IntegrationType::Gauss3,
                                           Gauss3_dQuadraticBasCubeSource);
}

// невязка не ухудшилась, с учётом того, что -1 соответствует неопределённой невязке
bool residualIsImproving(double prevResidual, double newResidual)
{
    return newResidual <= prevResidual || prevResidual == -1;
}

void MechPlasticIterInf::clearCounts(const int matrixSize)
{
    indexesChangedNumber = 0;
    tablesChangedNumber = 0;
     preliminarily_GLocalChangedNumber = 0;
    GLocalChangedNumber = 0;
     GGlobalChangedState = 0;
    bLocalChangedNumber = 0;
    nonlinearStateFENumber = 0;
    unloadStateFENumber = 0;
    GFirstStrChanged = matrixSize;
    maxPlasticResidual.set_null();
}
void MechPlasticIterInf::addFeIterResultsToCounts(const MechFeData_base *feEl)
{
    // необходимость итераций
    if(feEl->plasticNeedIterations)
        nonlinearStateFENumber++;
    // разгружающиеся КЭ
    if(feEl->isUnloading)
        unloadStateFENumber++;
    // невязки напряжений и деформаций
    maxPlasticResidual.compareAndSaveIfLarger(feEl->maxPlasticResidual);
    // данные об изменениях матриц, векторов, таблиц
    if(feEl->indexesChanged)
        indexesChangedNumber++;
    if(feEl->tablesChanged)
        tablesChangedNumber++;
    if(feEl->GLocalChanged)
    {
        GLocalChangedNumber++;
        GGlobalChangedState = 1;
        int minGlobalIndex = feEl->getMinGlobalIndex();//feEl->globalIndex[feEl->local_mn[0]*3 + 0];
        //int global_mk = feEl->local_mn[0]*3 + 0;//?????###
        // номера базисных ф-й отсортированы по возрастанию значит feEl->local_mn[0]*3 + 0 - минимальный номер изменённой строки
        if(minGlobalIndex < GFirstStrChanged)
            GFirstStrChanged = minGlobalIndex;
    }
    if(feEl->bLocalChanged)
        bLocalChangedNumber++;
}
void MechPlasticIterInf::calcIterResults(const MechGlobalStepParameters &step_el)
{
    // все элементы не требуют итераций
    if(nonlinearStateFENumber == 0)
    {
        accuracyAchieved = true;
        improvingAccuracy = false;    // дожать без нелинейности не получится
        maxPlasticResidual.set_null();
        return;
    }
    // точность не известна
    if(iterNumber == 0)
    {
        accuracyAchieved = false;
        improvingAccuracy = false;
        maxPlasticResidual.set_undefined();
        return;
    }
    // достигнута желаемая точность
    if(maxPlasticResidual.eps <= step_el.plasticResidualLimit.eps &&
       maxPlasticResidual.sigma <= step_el.plasticResidualLimit.sigma)
    {
        accuracyAchieved = true;
    }
    else
    {
        accuracyAchieved = false;
    }
    // проверяем, улучшаются ли невязки с итерациями (для дожимания)
    if(residualIsImproving(prevMaxPlasticResidual.eps, maxPlasticResidual.eps) &&
       residualIsImproving(prevMaxPlasticResidual.sigma, maxPlasticResidual.sigma))
    {
        // невязки не ухудшились и хотя бы одна улучшилась,
        improvingAccuracy = true;     // дожимаем
    }
    else
    {
        if(step_el.terminateIfAccuracyIsNotAchieving)
            accuracyAchieved = true;   // искусственно завершаем итерации######
        else
            improvingAccuracy = false;    // начало ухудшаться, прекращаем дожимать
    }
}

void MechFePointData::initNewIteration(MechMaterialSource &m0, const Grid::TimeLayers &tl, const bool firstIter, const IterationMode iterationMode)
{
    if(firstIter)
    {
        using namespace Elementary::Operations;
        // вычисление упругих параметров при температуре T + deltaT
        m0.calc_elastic_parameters(T + deltaT,
                                    elasticParameters_T);
        // вычисление эквивалентных напряжений и упругопластических деформаций
        sumSigma_eqv = sumSigma.eqv_SIGMA();
        // вычисление единичного девиатора напряжений и соответствующего аналога для деформаций

        if(sumSigma_eqv != 0)
        {
            z0 = (sumSigma.deviator() / sumSigma_eqv)*3./2;
            h0 = (sumSigma.deviator() / sumSigma.deviatorEuclideanNorm());
            it.z = z0;
            it.h = h0;
        }
        // температурные деформации
        if(deltaT != 0)
        {
            // изотропное тепловое расширение
            depsTerm.initDiag(elasticParameters_T.Talpha * deltaT);
            if(m0.temperatureDependence)
            {
                ElasticParameters ep_T1;
                ElasticParameters ep_T2;
                m0.calc_elastic_parameters(T, ep_T1);
                m0.calc_elastic_parameters(T + deltaT, ep_T2);
                TENSOR4 deltaA; // приращение обратного тензора упругости
                TENSOR4::setInvertDeltaIsotropicC(ep_T1, ep_T2, deltaA);
                depsTerm += deltaA * sumSigma;
            }
        }
        else
        {
            depsTerm.clear();
        }
        // вычисление начального приближения
        switch (m0.plasticityMethodType)
        {
        // линейная упругость
        case MechPlasticityMethodType::Elasticity:
        // метод начальных напряжений
        case MechPlasticityMethodType::InitialSigma:
        {
            it.A.is_elastic = true;
            it.dsigma0.clear();
        }break;
        // матрица D_pl
        case MechPlasticityMethodType::D_pl:
        // матрица D_pl + метод начальных напряжений
        case MechPlasticityMethodType::Combo_D_pl_InitialSigma:
        {
            it.dsigma0.clear();
            if(isYelding == true
                    && !it.lb.isNeutralLoading)// если на предыдущем шаге было нейтральное нагружение (легендарное), то матрицу можно оставить упругой
            {
                // касательная в точке начала шага
                it.A.is_elastic = false;
                switch (m0.PCType)
                {
                case MechPlasticityCurveType::Sigma_eps:
                {
                    it.A.value = m0.difSigma_Yeld(elasticParameters_T, q, tl.t_prev1);
                }break;
                case MechPlasticityCurveType::Eps_sigma:
                {
                    it.A.value = m0.difq_Yeld(elasticParameters_T, sumSigma_eqv, tl.t_prev1);
                }break;
                }
            }
            else
            {
                // упругость
                it.A.is_elastic = true;
            }
        }break;
        // матрица D_pl, алгоритм Соловейчика
        case MechPlasticityMethodType::D_pl_Solov:
        {
            it.dsigma0.clear();
            if(isYelding == true
                    && !it.lb.isNeutralLoading)// если на предыдущем шаге было нейтральное нагружение (легендарное), то матрицу можно оставить упругой
            {
                // касательная в точке начала шага
                it.A.is_elastic = false;
                switch (m0.PCType)
                {
                case MechPlasticityCurveType::Sigma_eps:
                {
                    it.A.value = m0.difSigma_Yeld(elasticParameters_T, q, tl.t_prev1);

                }break;
                case MechPlasticityCurveType::Eps_sigma:
                {
                    it.A.value = m0.difq_Yeld(elasticParameters_T, sumSigma_eqv, tl.t_prev1);
                }break;
                }
            }
            else
            {
                // упругость
                it.A.is_elastic = true;
            }
            TENSOR4 C_el;      // упругая матрица
            TENSOR4 Y_pl;      // матрица для нахождения компонент пластических деформаций (которая появляется из-за C_pl != C_el)
            TENSOR4::setCY_pl(elasticParameters_T, it.A, it.z,
                              C_el, it.C_pl, Y_pl);
        }break;
        }
    }
    else
    {
        switch (iterationMode)
        {
        case IterationMode::Free:
        {
            // итерации свободные
            switch (m0.plasticityMethodType)
            {
            // линейная упругость
            case MechPlasticityMethodType::Elasticity:
            {
                // нелинейностей нет
            }break;
            // матрица D_pl
            case MechPlasticityMethodType::D_pl:
            // матрица D_pl, алгоритм Соловейчика
            case MechPlasticityMethodType::D_pl_Solov:
            // метод начальных напряжений
            case MechPlasticityMethodType::InitialSigma:
            // матрица D_pl + метод начальных напряжений
            case MechPlasticityMethodType::Combo_D_pl_InitialSigma:
            {
                // обновление параметров итерации
                // предыдущие результаты итераций сохраняются в it_new (для подсчёта погрешностей)
                std::swap(it_new, it);
                //it_new = it;
                //it = it_new;
            }break;
            }
        }break;
        case IterationMode::ReadOnly:
        {
            // итерации приостановлены, пока не переключимся в свободный режим
        }break;
        }
    }
}
void MechFePointData::calcIterResults(MechMaterialSource &m0, const VECTOR3 *du, const Grid::TimeLayers &tl)
{
    using namespace Elementary::Operations;
    // обновление значений приращений
    TENSOR4 C_el;      // упругая матрица
    TENSOR4 C_pl;      // тензор, с которым строилось решение на текущей итерации
    TENSOR4 Y_pl;      // матрица для нахождения компонент пластических деформаций (которая появляется из-за C_pl != C_el)
    if(m0.plasticityMethodType == MechPlasticityMethodType::D_pl_Solov)
    {
        TENSOR4::setC_el(elasticParameters_T,
                         C_el);
        C_pl = it.C_pl;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                        Y_pl.m[i][j][k][l] = (C_el.m[i][j][k][l] - C_pl.m[i][j][k][l]) / (2*elasticParameters_T.G);
    }
    else
    {
        TENSOR4::setCY_pl(elasticParameters_T, it.A, it.z,
                          C_el, C_pl, Y_pl);
    }
    // приращение наблюдаемой деформации
    TENSOR2 deps;
    deps.m[0][0] = du[0][0];
    deps.m[1][1] = du[1][1];
    deps.m[2][2] = du[2][2];
    deps.m[1][2] = (du[1][2] + du[2][1])/2;
    deps.m[2][1] = (du[1][2] + du[2][1])/2;
    deps.m[0][2] = (du[2][0] + du[0][2])/2;
    deps.m[2][0] = (du[2][0] + du[0][2])/2;
    deps.m[0][1] = (du[0][1] + du[1][0])/2;
    deps.m[1][0] = (du[0][1] + du[1][0])/2;

    // приращение напряжения
    it.dsigma = C_pl * deps + it.dsigma0 - (C_el * depsTerm);

    // приращения деформаций
    {
        TENSOR2 depsPlastic1;
        TENSOR2 depsPlastic2;
        it.depsPlastic = (Y_pl * deps) - (it.dsigma0 / (2.*elasticParameters_T.G));
        it.depsElastic = deps - it.depsPlastic - depsTerm;
    }

    // приращение параметра Одквиста
    it.dq = it.depsPlastic.eqv_EPS();
    // пробное приращение напряжения (если считать приращение деформации упругой)
    it.dsigma_trial = C_el * (it.depsElastic + it.depsPlastic);
    // девиатор приращений упруго-пластических (не температурных) деформаций
    TENSOR2 deviator_deps = (it.depsElastic + it.depsPlastic).deviator();
    // девиатор напряжений
    TENSOR2 deviator_sigma = sumSigma.deviator();

    // работа внутренних сил
    if(sumSigma_eqv == 0)
    {
        // если sumSigmaEqv == 0, то dW = 0, но происходит нагружение
        it.dW = 1;
    }
    else
    {
        it.dW = deviator_sigma * deviator_deps;
    }

    // проверка превышения предела упругости
    bool yeld_state;
    switch (m0.plasticityMethodType)
    {
    // линейная упругость
    case MechPlasticityMethodType::Elasticity:
    {
        yeld_state = false;
    }break;
    // пластичность
    case MechPlasticityMethodType::D_pl:
    // матрица D_pl, алгоритм Соловейчика
    case MechPlasticityMethodType::D_pl_Solov:
    // метод начальных напряжений
    case MechPlasticityMethodType::InitialSigma:
    // матрица D_pl + метод начальных напряжений
    case MechPlasticityMethodType::Combo_D_pl_InitialSigma:
    {
        yeld_state = ((sumSigma + it.dsigma_trial).eqv_SIGMA() >= m0.sigma_Yeld(elasticParameters_T, q, tl.t0));
    }break;
    }
    // режим нагружения
    it.lb.clear();
    if(sumSigma_eqv == 0)
    {
        it.lb.isElasticLoading = true;
    }
    else
    {
        switch (m0.plasticityMethodType)
        {
        // линейная упругость
        case MechPlasticityMethodType::Elasticity:
        {
            if(it.dW > 0)
            {
                it.lb.isElasticLoading = true;
            }
            else
            {
                if(it.dW < 0)
                {
                    it.lb.isUnloading = true;
                }
                else
                {
                    // dW = 0
                    it.lb.isNeutralLoading = true;
                }
            }
        }break;
        // матрица D_pl
        case MechPlasticityMethodType::D_pl:
        // матрица D_pl, алгоритм Соловейчика
        case MechPlasticityMethodType::D_pl_Solov:
        {
            /*if(it.dW > 0 && yeld_state)
            {
                // активное нагружение
                it.lb.isActiveLoading = true;
            }*/
            if(yeld_state)
            {
                // предел упругости превышен
                // считаем его "активным" (хотя возможно dW <= 0)
                it.lb.isActiveLoading = true;
            }
            else
            {
                // упругое нагружение, нейтральное нагружение или разгрузка
                if(it.dW > 0)
                {
                    // упругое нагружение
                    it.lb.isElasticLoading = true;
                }
                if(it.dW < 0)
                {
                    // разгрузка
                    it.lb.isUnloading = true;
                }
                if(it.dW == 0)
                {
                    // нейтральное нагружение (легендарный случай)
                    it.lb.isNeutralLoading = true;
                }
            }
        }break;
        // метод начальных напряжений
        case MechPlasticityMethodType::InitialSigma:
        // матрица D_pl + метод начальных напряжений
        case MechPlasticityMethodType::Combo_D_pl_InitialSigma:
        {
            if(yeld_state)
            {
                // предел упругости превышен
                // считаем его "активным" (хотя возможно dW <= 0)
                it.lb.isActiveLoading = true;
                // расчёт угла между девиатором и приращением напряжений, если считать что пластичных деформаций нет
                {
                    double dsigma_trial_norm = sqrt(it.dsigma_trial * it.dsigma_trial);
                    double deviator_sigma_norm = sqrt(deviator_sigma * deviator_sigma);
                    double proj = deviator_sigma * it.dsigma_trial;
                    if(deviator_sigma_norm == 0 || dsigma_trial_norm == 0)
                        it.cosTetta = 0;
                    else
                        it.cosTetta = proj/deviator_sigma_norm/dsigma_trial_norm;
                }
                // определение it.readyToUnload
                if(it.cosTetta < m0.cosTettaMin)
                    it.lb.readyToUnload = true;
                else
                    it.lb.readyToUnload = false;
            }
            else
            {
                // предел упругости не превышен
                if(it.dW > 0)
                {
                    // упругое нагружение
                    it.lb.isElasticLoading = true;
                }
                if(it.dW < 0)
                {
                    // разгрузка
                    it.lb.isUnloading = true;
                }
                if(it.dW == 0)
                {
                    // нейтральное нагружение
                    it.lb.isNeutralLoading = true;
                }
            }
        }break;
        }
    }

    // определение неизбежности изменения матрицы с пластичной на упругую
    it.preliminarily_needUpdateMatrix = false;
    switch (m0.plasticityMethodType)
    {
    // линейная упругость
    case MechPlasticityMethodType::Elasticity:
    // пластичность
    case MechPlasticityMethodType::D_pl:
    // матрица D_pl, алгоритм Соловейчика
    case MechPlasticityMethodType::D_pl_Solov:
    // метод начальных напряжений
    case MechPlasticityMethodType::InitialSigma:
    {
    }break;
    // матрица D_pl + метод начальных напряжений
    case MechPlasticityMethodType::Combo_D_pl_InitialSigma:
    {
        if(!it.lb.isActiveLoading && !it.A.is_elastic)
        {
            it.preliminarily_needUpdateMatrix = true;
        }
    }break;
    }
}

void MechFePointData::prepareForNextIter(MechMaterialSource &m0, const int iterNumber, const Grid::TimeLayers &tl, int preliminarily_GLocalChangedNumber)
{
    // нахождение следующего приближения
    switch (m0.plasticityMethodType)
    {
    // линейная упругость
    case MechPlasticityMethodType::Elasticity:
    {
        it_new.A.is_elastic = true;
        it_new.dsigma0.clear();
        it.needIterations = false;
        it.needUpdateMatrix = false;
        it.residual.set_null();
    }break;
    // пластичность
    case MechPlasticityMethodType::D_pl:
    {
        it_new.dsigma0.clear();
        if(it.lb.isActiveLoading)
        {
            it.needIterations = true;
            it.needUpdateMatrix = false;
            switch (m0.PCType)
            {
            case MechPlasticityCurveType::Sigma_eps:
            {
                it_new.A.is_elastic = false;
                // новый направляющий тензор
                {
                    TENSOR2 sumSigma_trial_deviator = (sumSigma + it.dsigma_trial).deviator();
                    TENSOR2 h_new = (sumSigma_trial_deviator / sumSigma_trial_deviator.deviatorEuclideanNorm());
                    TENSOR2 z_new = (sumSigma_trial_deviator / sumSigma_trial_deviator.eqv_SIGMA())*3./2;

                    // средняя точка
                    it_new.z = z0 * (1. - m0.w_midPoint) + z_new * m0.w_midPoint;
                    it_new.h = h0 * (1. - m0.w_midPoint) + h_new * m0.w_midPoint;

                    // трапеция
                    //it_new.z = it.z * (1. - m0.w_midPoint) + z_new * m0.w_midPoint;
                    //it_new.h = it.h * (1. - m0.w_midPoint) + h_new * m0.w_midPoint;

                    it_new.z /= it_new.z.eqv_EPS();
                    it_new.h /= it_new.h.deviatorEuclideanNorm();
                }

                // чёткий поиск нового приближения
                {
                    // корректирующая добавка к начальному напряжению, чтобы общая пластическая деформация стала пропорциональна it_new.h
                    TENSOR2 dsigma0_corr;
                    bool search_crashed;
                    double d;
                    project(m0, tl,
                            dsigma0_corr, d, search_crashed);
                    if(search_crashed)
                    {
                        // матрица становится упругой
                        it_new.A.is_elastic = true;
                        if(!it.A.is_elastic)
                        {
                            it.needUpdateMatrix = true;
                        }
                    }
                    else
                    {
                        double a = (it_new.z * (it.depsElastic + it.depsPlastic)) * 2*elasticParameters_T.G;
                        TENSOR2 depsPlastic_want = it.depsPlastic + (dsigma0_corr + it_new.h*d*m0.w_project)/(-2.*elasticParameters_T.G);
                        double b = a * (it_new.z * it_new.h);
                        double c = depsPlastic_want * it_new.h;
                        if(b == 0 || c == 0)
                        {
                            // матрица становится упругой
                            it_new.A.is_elastic = true;
                            if(!it.A.is_elastic)
                            {
                                it.needUpdateMatrix = true;
                            }
                        }
                        else
                        {
                            it_new.A.value = b / c - 3*elasticParameters_T.G;
                            it.needUpdateMatrix = true;
                        }
                    }
                    it.residual.sigma = m0.sigma_Yeld(elasticParameters_T, q + it.dq, tl.t0) - (sumSigma + it.dsigma).eqv_SIGMA();
                    it.residual.eps = it.dq - it_new.dq;
                    it.residual.sigma /= sumSigma_eqv;
                    it.residual.eps /= (sumEpsElastic + sumEpsPlastic).eqv_EPS();
                }

                // магическая формула с умножением
                /*
                {
                    it.residual.sigma = m0.sigma_Yeld(elasticParameters_T, q + it.dq, tl.t0) - (sumSigma + it.dsigma).eqv_SIGMA();
                    it.residual.eps = it.dq - it_new.dq;
                    it.residual.sigma /= sumSigma_eqv;
                    it.residual.eps /= (sumEpsElastic + sumEpsPlastic).eqv_EPS();
                    it_new.A.is_elastic = false;


                    if(it.A.is_elastic)
                    //if(it.dq == 0)
                    {
                        it_new.A.value = m0.difSigma_Yeld(elasticParameters_T, q, tl.t_prev1);
                        //it_new.A.value = (m0.sigma_Yeld(elasticParameters_T, q + it.depsElastic.eqv_EPS(), tl.t0) - sumSigma.eqv_SIGMA()) / it.depsElastic.eqv_EPS();
                        //it_new.A.value = (m0.sigma_Yeld(elasticParameters_T, q + it.dq, tl.t0) - sumSigma.eqv_SIGMA()) / it.depsElastic.eqv_EPS();
                    }
                    else
                    {
                        it_new.A.value *= (sumSigma + it.dsigma).eqv_SIGMA() / m0.sigma_Yeld(elasticParameters_T, q + it.dq, tl.t0);
                        //if(it.residual.eps < 1.e-3)
                        //{
                        //    it_new.A.value *= (sumSigma + it.dsigma).eqv_SIGMA() / m0.sigma_Yeld(elasticParameters_T, q + it.dq, tl.t0);
                        //}
                        //else
                        //{
                        //    it_new.A.value = (m0.sigma_Yeld(elasticParameters_T, q + it.dq, tl.t0) - sumSigma.eqv_SIGMA()) / it.dq;
                        //}
                    }
                }*/
            }break;
            case MechPlasticityCurveType::Eps_sigma:
            {
            }break;
            }
        }
        else
        {
            // упругое нагружение, нейтральное нагружение или разгрузка
            if(it.A.is_elastic)
            {
                // не меняется матрица ОС
                it.needIterations = false;
                it.needUpdateMatrix = false;
                it.residual.set_null();
            }
            else
            {
                // меняется матрица ОС
                it.needIterations = true;
                it.needUpdateMatrix = true;
                it.residual.set_undefined();
            }
            it_new.A.is_elastic = true;
        }

    }break;
    // матрица D_pl, алгоритм Соловейчика
    case MechPlasticityMethodType::D_pl_Solov:
    {
        it_new.dsigma0.clear();
        if(it.lb.isActiveLoading)
        {
            it.needIterations = true;
            it.needUpdateMatrix = false;
            switch (m0.PCType)
            {
            case MechPlasticityCurveType::Sigma_eps:
            {
                it_new.A.is_elastic = false;
                // новый направляющий тензор
                {
                    TENSOR2 sumSigma_trial_deviator = (sumSigma + it.dsigma_trial).deviator();
                    TENSOR2 h_new = (sumSigma_trial_deviator / sumSigma_trial_deviator.deviatorEuclideanNorm());
                    TENSOR2 z_new = (sumSigma_trial_deviator / sumSigma_trial_deviator.eqv_SIGMA())*3./2;

                    // средняя точка
                    it_new.z = z0 * (1. - m0.w_midPoint) + z_new * m0.w_midPoint;
                    it_new.h = h0 * (1. - m0.w_midPoint) + h_new * m0.w_midPoint;

                    // трапеция
                    //it_new.z = it.z * (1. - m0.w_midPoint) + z_new * m0.w_midPoint;
                    //it_new.h = it.h * (1. - m0.w_midPoint) + h_new * m0.w_midPoint;

                    it_new.z /= it_new.z.eqv_EPS();
                    it_new.h /= it_new.h.deviatorEuclideanNorm();
                }

                // чёткий поиск нового приближения
                {
                    // корректирующая добавка к начальному напряжению, чтобы общая пластическая деформация стала пропорциональна it_new.h
                    TENSOR2 dsigma0_corr;
                    bool search_crashed;
                    double d;
                    project(m0, tl,
                            dsigma0_corr, d, search_crashed);
                    if(search_crashed)
                    {
                        // матрица становится упругой
                        it_new.A.is_elastic = true;
                        if(!it.A.is_elastic)
                        {
                            it.needUpdateMatrix = true;
                        }
                    }
                    else
                    {
                        double a = (it_new.z * (it.depsElastic + it.depsPlastic)) * 2*elasticParameters_T.G;
                        TENSOR2 depsPlastic_want = it.depsPlastic + (dsigma0_corr + it_new.h*d*m0.w_project)/(-2.*elasticParameters_T.G);
                        double b = a * (it_new.z * it_new.h);
                        double c = depsPlastic_want * it_new.h;
                        if(b == 0 || c == 0)
                        {
                            // матрица становится упругой
                            it_new.A.is_elastic = true;
                            if(!it.A.is_elastic)
                            {
                                it.needUpdateMatrix = true;
                            }
                        }
                        else
                        {
                            it_new.A.value = b / c - 3*elasticParameters_T.G;
                            it.needUpdateMatrix = true;
                        }
                    }
                    it.residual.sigma = m0.sigma_Yeld(elasticParameters_T, q + it.dq, tl.t0) - (sumSigma + it.dsigma).eqv_SIGMA();
                    it.residual.eps = it.dq - it_new.dq;
                    it.residual.sigma /= sumSigma_eqv;
                    it.residual.eps /= (sumEpsElastic + sumEpsPlastic).eqv_EPS();
                }

                // новое приближение с учётом кривой
                /*{
                    it.residual.sigma = m0.sigma_Yeld(elasticParameters_T, q + it.dq, tl.t0) - (sumSigma + it.dsigma).eqv_SIGMA();
                    it.residual.eps = it.dq - it_new.dq;
                    it.residual.sigma /= sumSigma_eqv;
                    it.residual.eps /= (sumEpsElastic + sumEpsPlastic).eqv_EPS();
                    it_new.A.is_elastic = false;
                    if(it.A.is_elastic)
                    {
                        //it_new.A.value = m0.difSigma_Yeld(elasticParameters_T, q, tl.t_prev1);
                        it_new.A.value = (m0.sigma_Yeld(elasticParameters_T, q + it.depsElastic.eqv_EPS(), tl.t0) - sumSigma.eqv_SIGMA()) / it.depsElastic.eqv_EPS();
                    }
                    else
                    {
                        it_new.A.value = (m0.sigma_Yeld(elasticParameters_T, q + it.dq, tl.t0) - sumSigma.eqv_SIGMA()) / it.dq;
                    }
                }*/

                // алгоритм Соловейчика
                {
                    TENSOR4 C_pl;      // новый касательный тензор
                    {
                        TENSOR4 C_el;      // упругий тензор
                        TENSOR4 Y_pl;      // матрица для нахождения компонент пластических деформаций (которая появляется из-за C_pl != C_el)
                        TENSOR4::setCY_pl(elasticParameters_T, it_new.A, it_new.z,
                                          C_el, C_pl, Y_pl);
                    }
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            for (int k = 0; k < 3; k++)
                                for (int l = 0; l < 3; l++)
                                    //it_new.C_pl.m[i][j][k][l] = 0.5*(it.C_pl.m[i][j][k][l] + C_pl.m[i][j][k][l]);
                                    it_new.C_pl.m[i][j][k][l] = C_pl.m[i][j][k][l];
                }
            }break;
            case MechPlasticityCurveType::Eps_sigma:
            {
            }break;
            }
        }
        else
        {
            // упругое нагружение, нейтральное нагружение или разгрузка
            if(it.A.is_elastic)
            {
                // не меняется матрица ОС
                it.needIterations = false;
                it.needUpdateMatrix = false;
                it.residual.set_null();
            }
            else
            {
                // меняется матрица ОС
                it.needIterations = true;
                it.needUpdateMatrix = true;
                it.residual.set_undefined();
            }
            it_new.A.is_elastic = true;
            TENSOR4::setC_el(elasticParameters_T,
                              it_new.C_pl);
        }

    }break;
    // метод начальных напряжений
    case MechPlasticityMethodType::InitialSigma:
    // матрица D_pl + метод начальных напряжений
    case MechPlasticityMethodType::Combo_D_pl_InitialSigma:
    {
        it_new.A = it.A;
        it_new.dsigma0.clear();
        if(it.lb.isActiveLoading)
        {
            // новый направляющий тензор
            {
                TENSOR2 sumSigma_trial_deviator = (sumSigma + it.dsigma_trial).deviator();
                TENSOR2 h_new = (sumSigma_trial_deviator / sumSigma_trial_deviator.deviatorEuclideanNorm());
                it_new.z = z0;

                // средняя точка
                it_new.h = h0 * (1. - m0.w_midPoint) + h_new * m0.w_midPoint;

                // трапеция
                //it_new.h = it.h * (1. - m0.w_midPoint) + h_new * m0.w_midPoint;

                it_new.h /= it_new.h.deviatorEuclideanNorm();
            }
            // чёткий поиск нового приближения
            {
                // корректирующая добавка к начальному напряжению, чтобы общая пластическая деформация стала пропорциональна it_new.h
                TENSOR2 dsigma0_corr;
                bool search_crashed;
                double d;
                project(m0, tl,
                        dsigma0_corr, d, search_crashed);
                if(search_crashed)
                {
                    if(!it.A.is_elastic)
                    {
                        // матрица становится упругой
                        it.needUpdateMatrix = true;
                        it_new.A.is_elastic = true;
                    }
                    // зануление приращения пластических деформаций
                    it_new.dsigma0.clear();
                }
                else
                {
                    if((preliminarily_GLocalChangedNumber >= 1 && !it.A.is_elastic && it.lb.readyToUnload)
                            || (!it.A.is_elastic && it.dW <= 0))    // при it.dW <= 0 упруго-пластическая матрица мешает итерациям
                    {
                        // если матрица пластичная, хотя бы для 1 КЭ (среди всех КЭ) точно придётся менять матрицу с пластичной на упругую, и cos угла достаточно мал, то матрица меняется с пластичной на упругую
                        it.needUpdateMatrix = true;
                        it_new.A.is_elastic = true;
                        // все пластические деформации перерасчитываются в начальное напряжение
                        TENSOR2 dsigmaPlastic = it.depsPlastic * (-2.*elasticParameters_T.G);   // начальное напряжение, соответствующее пластической деформации it.depsPlastic
                        it_new.dsigma0 = dsigmaPlastic + dsigma0_corr + it_new.h*d*m0.w_project;
                    }
                    else
                    {
                        // матрица остаётся прежней
                        it.needUpdateMatrix = false;
                        it_new.dsigma0 = it.dsigma0 + dsigma0_corr + it_new.h*d*m0.w_project;
                    }
                }
                it.needIterations = true;
                it.residual.sigma = m0.sigma_Yeld(elasticParameters_T, q + it.dq, tl.t0) - (sumSigma + it.dsigma).eqv_SIGMA();
                it.residual.eps = it.dq - it_new.dq;
                it.residual.sigma /= sumSigma_eqv;
                it.residual.eps /= (sumEpsElastic + sumEpsPlastic).eqv_EPS();
            }
        }
        else
        {
            if(it.A.is_elastic)
            {
                // не меняется матрица ОС
                it.needIterations = false;
                it.needUpdateMatrix = false;
                it.residual.set_null();
            }
            else
            {
                // меняется матрица ОС
                it.needIterations = true;
                it.needUpdateMatrix = true;
                it.residual.set_undefined();
            }
            it_new.A.is_elastic = true;
        }
    }break;
    }
}
void MechFePointData::finalizeStep(const MechMaterialSource &m0, const int)
{
    using namespace Operations;
    // добавление приращений
    q += it.dq;
    sumEpsTerm += depsTerm;
    sumEpsElastic += it.depsElastic;
    sumEpsPlastic += it.depsPlastic;
    sumSigma += it.dsigma;
    // обновление статуса состояния текучести
    if(isYelding == false && it.lb.isActiveLoading)
    {
        isYelding = true;
    }
    if(isYelding == true && it.lb.isUnloading)
    {
        isYelding = false;
    }
}
void MechFePointData::saveResultsOfStep(MechOutFePointData &outFePD)
{
    outFePD.q = q;
    outFePD.sumEpsElastic = sumEpsElastic;
    outFePD.sumEpsPlastic = sumEpsPlastic;
    outFePD.sumEpsTerm = sumEpsTerm;
    outFePD.sumSigma = sumSigma;
    outFePD.isYelding = isYelding;
    outFePD.residual = it.residual;
}

void MechFePointData::project(MechMaterialSource &m0, const Grid::TimeLayers &tl, TENSOR2 &dsigma0_corr, double &d, bool &search_crashed)
{
    // чёткий поиск нового приближения
    // начальное напряжение, соответствующее пластической деформации it.depsPlastic
    TENSOR2 dsigmaPlastic = it.depsPlastic * (-2.*elasticParameters_T.G);
    // корректирующая добавка к начальному напряжению, чтобы общая пластическая деформация стала пропорциональна it_new.h
    dsigma0_corr = it_new.h * (dsigmaPlastic * it_new.h) - dsigmaPlastic;
    search_crashed = false;
    d = 0; // искомая добавка начальных напряжений
    // (it_new.dsigma0 = it.dsigma0 + dsigma0_corr + it_new.h*d*m0.w_initialSigma)
    double d_residual;
    double d_residual_prev = 1.e100;
    // если !it.A.is_elastic и it.depsPlastic < 0, то корректируем d
    /*
    {
        TENSOR2 dsigmaPlastic_0 = dsigmaPlastic + dsigma0_corr;
        double l = dsigmaPlastic_0*it_new.h;
        if(l > 0)
        {
            d = -l;
        }
    }*/
    for(int i = 0; i < 10000; i++)
    {
        TENSOR2 dsigma0_add_i = dsigma0_corr + it_new.h*d;
        TENSOR2 depsPlastic_i = it.depsPlastic + dsigma0_add_i/(-2.*elasticParameters_T.G);
        // проверка знака общей пластической деформации (не должна быть < 0)
        //if(i >= 1)
        {
            double sign = depsPlastic_i * it_new.h;
            if(sign < 0)
            {
                // поиск крашнулся, т.к. получилась отрицательная пластическая деформация
                search_crashed = true;
                return;
            }
        }
        double q_i = q + depsPlastic_i.eqv_EPS();
        double sumSigmaEqv_i = (sumSigma + it.dsigma + dsigma0_add_i).eqv_SIGMA();
        d_residual = m0.sigma_Yeld(elasticParameters_T, q_i, tl.t0) - sumSigmaEqv_i;
        d += d_residual*0.1;
        if(fabs(d_residual/sumSigma_eqv) < 1.e-14 && fabs(d_residual) >= fabs(d_residual_prev))
            break;
        d_residual_prev = d_residual;
    }
}

void MechFeData_LinearHexagon::G_local_clear()
{
    for(int i = 0; i < 300; i++)
        GLocalL[i] = 0;
}
void MechFeData_LinearHexagon::b_local_clear()
{
    for(int i = 0; i < 24; i++)
        b_local_el(i) = 0;
}
void MechFeData_LinearHexagon::R_local_clear()
{
    for(int i = 0; i < 24; i++)
        R_local_el(i) = 0;
}
void MechFeData_LinearHexagon::updateIndexes(const Grid::FE_base *feEl, const Grid::GlobalDOFs *DOFs)
{
    // локальные и глобальные индексы
    if(indexesChanged)
    {
        // копии глобальных индексов и координат вершин 6-гранника
        int vi[8];
        feEl->getVertexIndexes(vi);
        // заполнение исходными значениями по умолчанию
        int local_ind = 0;
        for (int mn = 0; mn < 8; mn++)
        {
            const VECTOR3_int &DOFIndex = DOFs->findDOFIndex_Lagr1(vi[mn]);
            for (int ik = 0; ik < 3; ik++)
            {
                local_mn[local_ind] = mn;
                local_ik[local_ind] = ik;
                globalIndex[local_ind] = DOFIndex[ik];
                //globalIndex[local_ind] = vertexBusFuncInd->linear[vi[mn]][ik];//vi[mn] * 3 + ik;
                local_ind++;
            }
        }
        // сортировка по возрастанию значений глобальных индексов
        for (int local_ind1 = 0; local_ind1 < 24; local_ind1++)
            for (int local_ind2 = local_ind1 + 1; local_ind2 < 24; local_ind2++)
            {
                if(globalIndex[local_ind1] > globalIndex[local_ind2])
                {
                    std::swap(local_mn[local_ind1], local_mn[local_ind2]);
                    std::swap(local_ik[local_ind1], local_ik[local_ind2]);
                    std::swap(globalIndex[local_ind1], globalIndex[local_ind2]);
                }
            }
        // индексы обновлены, значит остальное тоже требуется обновить
        indexesChanged = false;
        tablesChanged = true;
        GLocalChanged = true;
        bLocalChanged = true;
    }
}

void MechFeData_LinearHexagon::add_BC2_to_localb(const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechBoundaryCondition2Source_base *bc2SourceEl, const int faceIndex, const Integration::Integrator &integrationFoursquare)
{
    if(bLocalChanged)
    {
        double basCube[8*Integration::Integration2D_Size_Max];
        double hexagonCoef[Integration::Integration2D_Size_Max];
        VECTOR3 hexagonNormal[Integration::Integration2D_Size_Max];
        feEl->calcBc2BasisFuncValues(grid->vertex, grid->vertexForCurvature, integrationFoursquare, faceIndex,
                                      basCube, hexagonCoef, hexagonNormal);
        // добавка к локальному вектору
        for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
        {
            int m = local_mn[local_mi_ind];
            int i = local_ik[local_mi_ind];
            double E = 0;
            for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
            {
                VECTOR3 vectorP =
                        bc2SourceEl->calcVectorValueInPoint2(hexagonNormal[valIndex]) -
                        bc2SourceEl->calcVectorValueInPoint1(hexagonNormal[valIndex]);
                E += vectorP[i] *
                     basCube[m*integrationFoursquare.size + valIndex]*
                     hexagonCoef[valIndex];
            }
            b_local_el(local_mi_ind) += E;
        }
    }
};
void MechFeData_LinearHexagon::add_BC2_to_localR(const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechBoundaryCondition2Source_base *bc2SourceEl, const int faceIndex, const Integration::Integrator &integrationFoursquare, const IncForsesMode incForsesMode)
{
    if(incForsesMode == IncForsesMode::IncrementP)
        return;
    double basCube[8*Integration::Integration2D_Size_Max];
    double hexagonCoef[Integration::Integration2D_Size_Max];
    VECTOR3 hexagonNormal[Integration::Integration2D_Size_Max];
    feEl->calcBc2BasisFuncValues(grid->vertex, grid->vertexForCurvature, integrationFoursquare, faceIndex,
                                  basCube, hexagonCoef, hexagonNormal);
    // добавка к локальному вектору невязки
    for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
    {
        int m = local_mn[local_mi_ind];
        int i = local_ik[local_mi_ind];
        double E = 0;
        for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
        {
            VECTOR3 vectorP;
            switch (incForsesMode)
            {
            case IncForsesMode::MinusIntegral:
            {
                // вычисляется в начале шага, после сдвига сетки
                // со знаком +
                vectorP = bc2SourceEl->calcVectorValueInPoint1(hexagonNormal[valIndex]);
                E += vectorP[i] *
                     basCube[m*integrationFoursquare.size + valIndex]*
                     hexagonCoef[valIndex];
            }break;
            case IncForsesMode::bPlusR:
            {
                // вычисляется в конце шага, после всех итераций, до сдвига сетки
                // со знаком +
                vectorP = bc2SourceEl->calcVectorValueInPoint2(hexagonNormal[valIndex]);
                E += vectorP[i] *
                     basCube[m*integrationFoursquare.size + valIndex]*
                     hexagonCoef[valIndex];
            }break;
            case IncForsesMode::IncrementP:
            {
            }break;
            }
        }
        R_local_el(local_mi_ind) += E;
    }
};

void MechFeData_LinearHexagon::addGlobalIndexesToGlobalPortrait1x1(SlauSolving::SSCMPortraitBulder &pb)
{
    for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
    {
        for (int local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
        {
            int global_mi = globalIndex[local_mi_ind];
            int global_nk = globalIndex[local_nk_ind]; // local_nk_ind < local_mi_ind
            pb.addElement(global_mi, global_nk); // global_nk < global_mi
        }
    }
}
void MechFeData_LinearHexagon::updateSSCM_indexes1x1(const SlauSolving::SSCMPortrait &portrait)
{
    int G_local_ind = 0;
    for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
    {
        for (int local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
        {
            int global_mi = globalIndex[local_mi_ind];
            int global_nk = globalIndex[local_nk_ind];
            if(global_mi == global_nk)
                SSCM_index[G_local_ind] = -1;//######
            else
                SSCM_index[G_local_ind] = portrait.findSorted(global_mi, global_nk);
                        //portrait.findUnsorted(global_mi, global_nk);
            G_local_ind++;
        }
    }
}
void MechFeData_LinearHexagon::addLocalGM_to_global1x1(SlauSolving::SSCMElements &elements)
{
    int G_local_ind = 0;
    for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
    {
        for (int local_nk_ind = 0; local_nk_ind < local_mi_ind; local_nk_ind++)
        {
            //int global_mi = feDataEl.globalIndex[local_mi_ind];
            //int global_nk = feDataEl.globalIndex[local_nk_ind];
            //(*Gmatrix.e).a[Gmatrix.p->findUnsorted(global_mi, global_nk)] += feDataEl.G_local_el(G_local_ind);
            //(*Gmatrix.e).a[Gmatrix.p->findSorted(global_mi, global_nk)] += feDataEl.G_local_el(G_local_ind);
            elements.a[SSCM_index[G_local_ind]] += G_local_el(G_local_ind);
            G_local_ind++;
        }
        elements.d[globalIndex[local_mi_ind]] += G_local_el(G_local_ind);
        G_local_ind++;
    }
}
void MechFeData_LinearHexagon::addLocalb_to_global1x1(Vector &b)
{
    for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
    {
        b[globalIndex[local_mi_ind]] += b_local_el(local_mi_ind) + R_local_el(local_mi_ind);
    }
}

void MechFeData_LinearHexagon::addGlobalIndexesToGlobalPortrait3x3(SlauSolving::SSCMPortraitBulder &pb)
{
    for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
    {
        for (int local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
        {
            int global_m = globalIndex[local_mi_ind]/3;
            int global_n = globalIndex[local_nk_ind]/3; // local_nk_ind <= local_mi_ind
            pb.addElement(global_m, global_n); // global_n <= global_m
        }
    }
}
void MechFeData_LinearHexagon::updateSSCM_indexes3x3(const SlauSolving::SSCMPortrait &portrait)
{
    int G_local_ind = 0;
    for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
    {
        for (int local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
        {
            int global_m = globalIndex[local_mi_ind]/3;
            int global_n = globalIndex[local_nk_ind]/3;
            if(global_m == global_n)
                SSCM_index[G_local_ind] = -1;//######
            else
                SSCM_index[G_local_ind] = portrait.findSorted(global_m, global_n);
            G_local_ind++;
        }
    }
}
void MechFeData_LinearHexagon::addLocalGM_to_global3x3(SlauSolving::SSCM3x3Elements &elements)
{
    int G_local_ind = 0;
    for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
    {
        for (int local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
        {
            size_t bi = local_mi_ind%3;
            size_t bj = local_nk_ind%3;
            if(SSCM_index[G_local_ind] != (size_t)-1)
            {
                // ниже диагонали
                MATR3x3 &el = elements.a[SSCM_index[G_local_ind]];
                el.m[bi][bj] += G_local_el(G_local_ind);
            }
            else
            {
                // диагональ
                int global_m = globalIndex[local_mi_ind]/3;
                MATR3x3 &el = elements.d[global_m];
                if(bi == bj)
                {
                    // диагональ блока
                    el.m[bi][bj] += G_local_el(G_local_ind);
                }
                else
                {
                    // вне диагонали
                    el.m[bi][bj] += G_local_el(G_local_ind);
                    el.m[bj][bi] += G_local_el(G_local_ind);
                }
            }
            G_local_ind++;
        }
    }
}
void MechFeData_LinearHexagon::addLocalb_to_global3x3(Vector3 &b)
{
    for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
    {
        int global_m = globalIndex[local_mi_ind]/3;
        size_t bi = local_mi_ind%3;
        b[global_m][bi] += b_local_el(local_mi_ind) + R_local_el(local_mi_ind);
    }
}

bool MechFeData_LinearHexagon::localb_forces_is_null()
{
    for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
    {
        if(b_local_el(local_mi_ind) != 0)
            return false;
    }
    return true;
}
size_t MechFeData_LinearHexagon::getMinGlobalIndex() const
{
    return globalIndex[0];
}

void MechFeData_LinearHexagon_homogeny::updateTables(const Grid::Grid3D *grid, const Grid::FE_base *feEl, const int numPoints, const double *integrationw, const double *basCube, const VECTOR3 *dLinearBasCube, const VECTOR3 *dQuadraticBasCube)
{
    // Таблицы коэффииентов для расчёта элементов локальных матриц и векторов
    // Для однородного материала константы материала можно вынести за знак интеграла
    if(tablesChanged)
    {
        using namespace Integration;
        double w[Integration3D_Size_Gauss3];
        double dbas[8][Integration3D_Size_Gauss3][3];       // значения производных базисных функций(8) в каждой точке интегрирования 6-гранника
        // расчёт коэффициентов для интегрирования, модулей детерминантов
        // и производных базисных функций на шестиграннике
        feEl->calcBasisFuncValues(grid->vertex, grid->vertexForCurvature, numPoints, integrationw, dLinearBasCube, dQuadraticBasCube,
                                    w, dbas);
        // внутренние интегралы для добавок к матрицам и к правой части
        // для матрицы G
        Fem::calcIntForG(numPoints, w, dbas, intForG_mnjl);
        // для матрицы M
        //if(ro != 0)
        //    Fem::solveIntForM(numPoints, w, basCube, intForM_mn);
        // для 1-го слагаемого правой части
        Fem::calcIntForb_f(numPoints, w, basCube, intForb1_m);
        // для 2-го и 3-го слагаемого правой части
        Fem::calcIntForb_df(numPoints, w, dbas, intForb23_mj);
        tablesChanged = false;
    }
}
void MechFeData_LinearHexagon_homogeny::init(Integration::IntegrationType set_integrationType, const Grid::FE_base *feEl)
{
    indexesChanged = true;
    tablesChanged = true;
    GLocalChanged = true;
    bLocalChanged = true;
    needUpdateMatrix = false;
    plasticNeedIterations = false;
    isUnloading = false;
    integrationType = set_integrationType;
    pd.init();
    R_local_clear();
}
void MechFeData_LinearHexagon_homogeny::setT(const double T, const double deltaT)
{
    pd.T = T;
    pd.deltaT = deltaT;
}
void MechFeData_LinearHexagon_homogeny::finalizeStep(const MechMaterialSource &m0, const int movingGridMode)
{
    pd.finalizeStep(m0, movingGridMode);
    // на нулевой итерации следующего шага требуется обновить матрицу
    tablesChanged = true;
    GLocalChanged = true;
    bLocalChanged = true;
}
void MechFeData_LinearHexagon_homogeny::saveResultsOfStep(const Vector &q, const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechOutFeData &outFeDataEl)
{
    // копии глобальных индексов и координат вершин 6-гранника
    POINT3 v[27];
    feEl->getGeomVertexes(grid->vertex, grid->vertexForCurvature, v);
    // копируем информацию в точках элемента
    outFeDataEl.pd.resize(1);
    MechOutFePointData &outFePD = outFeDataEl.pd[0];  // данные вывода для точки конечного элемента номер i
    pd.saveResultsOfStep(outFePD);
    POINT3 XYZ = POINT3(0,0,0);
    POINT3 xyz;
    feEl->cubeToHexagon(v, XYZ, Fem::dif_NULL3, xyz);
    outFePD.p = xyz;
}

void MechFeData_LinearHexagon_homogeny_E_NLE_EP::updateLocalGMb(const Grid::Grid3D *grid, const Grid::GlobalDOFs *DOFs, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0, const MechPlasticIterInf &nlInfPlastic, const bool firstIter, const IncForsesMode incForsesMode)
{
    // подготовка параметров для итерации
    pd.initNewIteration(m0, tl, firstIter, nlInfPlastic.iterationMode);
    // расчёт матрицы определяющих соотношений
    TENSOR4 C_pl0;  // матрица D с волной, записанная в виде тензора
    if(m0.plasticityMethodType == MechPlasticityMethodType::D_pl_Solov)
    {
        C_pl0 = pd.it.C_pl;
    }
    else
    {
        TENSOR4 C_el0;
        TENSOR4 Y_pl0;
        TENSOR4::setCY_pl(pd.elasticParameters_T, pd.it.A, pd.it.z,
                          C_el0, C_pl0, Y_pl0);
    }
    //VECTOR3 F = m0.F;     // объемная сила
    //double ro = m0.ro;    // плотность
    //const double Talpha = pd.elasticParameters_T.Talpha;
    const double deltaT = pd.deltaT;
    // локальные и глобальные индексы
    updateIndexes(feEl, DOFs);
    // таблицы с коэффициентами для текущей формы шестигранника
    if(integrationType == Integration::IntegrationType::Gauss2)
        updateTables(grid, feEl, 8, Gauss2_integrationwSource, Gauss2_basCubeSource, Gauss2_dLinearBasCubeSource, Gauss2_dQuadraticBasCubeSource);
    if(integrationType == Integration::IntegrationType::Gauss3)
        updateTables(grid, feEl, 27, Gauss3_integrationwSource, Gauss3_basCubeSource, Gauss3_dLinearBasCubeSource, Gauss3_dQuadraticBasCubeSource);

    // I. добавка к матрице G
    if(GLocalChanged)
    {
        // очистка локальной матрицы G
        G_local_clear();
        // заполнение локальной матрицы G (без XFEM)
        int G_local_ind = 0;
        for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
        {
            int m = local_mn[local_mi_ind];
            int i = local_ik[local_mi_ind];
            for (int local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
            {
                int n = local_mn[local_nk_ind];
                int k = local_ik[local_nk_ind];
                double E = 0;
                for (int j = 0; j < 3; j++)
                    for (int l = 0; l < 3; l++)
                        E += C_pl0.m[i][j][k][l] * intForG_mnjl[m][n][j][l];
                G_local_el(G_local_ind) += E;
                G_local_ind++;
            }
        }
        GLocalChanged = false;
    }
    // III добавка к локальному вектору
    if(bLocalChanged)
    {
        MATR3x3x3x3 C_el;  // упругая матрица D, записанная в виде тензора
        {
            MATR6x6 D_el;
            m0.calcD(pd.elasticParameters_T, D_el);
            Elementary::Operations::convert6x6to3x3x3x3(D_el, C_el);
        }
        // очистка локального вектора b
        b_local_clear();
        // заполнение локального вектора b
        // III.3 добавка к локальному вектору (температурные слагаемые: alpha*deltaT*... и sigmakl*dAijkl...)
        if(deltaT != 0)
        {
            // добавка в правую часть узлавых усилий, соответствующих температурным деформациям
            for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
            {
                int m = local_mn[local_mi_ind];
                int i = local_ik[local_mi_ind];
                double E = 0;
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            E += C_el.m[i][j][k][l] * intForb23_mj[m][j] *
                                    pd.depsTerm.m[k][l];
                b_local_el(local_mi_ind) += E;
            }
        }
        // III.5 добавка к локальному вектору (начальные напряжения)
        if(pd.it.dsigma0.eqv_SIGMA() != 0)
        {
            for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
            {
                int m = local_mn[local_mi_ind];
                int i = local_ik[local_mi_ind];
                double E = 0;
                for (int j = 0; j < 3; j++)
                    E += pd.it.dsigma0.m[i][j] * intForb23_mj[m][j];
                b_local_el(local_mi_ind) -= E;
            }
        }
    }
}
void MechFeData_LinearHexagon_homogeny_E_NLE_EP::calcLocalR(const Grid::Grid3D *grid, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0, const MechPlasticIterInf &nlInfPlastic, const IncForsesMode incForsesMode)
{
    R_local_clear();
    if(incForsesMode == IncForsesMode::IncrementP)
        return;
    //const int map3x3to6[3][3] =
    //{{0,5,4},
    // {5,1,3},
    // {4,3,2},};
    for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
    {
        int m = local_mn[local_mi_ind];
        int i = local_ik[local_mi_ind];
        double E = 0;
        for (int j = 0; j < 3; j++)
        {
            switch (incForsesMode)
            {
            case IncForsesMode::MinusIntegral:
            {
                // вычисляется в начале шага, после сдвига сетки
                // со знаком -
                E -= pd.sumSigma.m[i][j] * intForb23_mj[m][j];
                //E -= pd.sumSigma[map3x3to6[i][j]] * intForb23_mj[m][j];
            }break;
            case IncForsesMode::bPlusR:
            {
                // вычисляется в конце шага, после всех итераций, после приращения напряжений, до сдвига сетки
                // со знаком -
                //E -= (pd.sumSigma + pd.dsigma)[map3x3to6[i][j]] * intForb23_mj[m][j];
                E -= pd.sumSigma.m[i][j] * intForb23_mj[m][j];
                //E -= pd.sumSigma[map3x3to6[i][j]] * intForb23_mj[m][j];
            }break;
            case IncForsesMode::IncrementP:
            {
            }break;
            }
        }
        R_local_el(local_mi_ind) += E;
    }
}
void MechFeData_LinearHexagon_homogeny_E_NLE_EP::calcIterResults1(const Vector &q, const Grid::Grid3D *grid, const Grid::GlobalDOFs *DOFs, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0,
                                                                  std::vector<MechVertexData> *vertex, std::vector<MechVertexData> *vertexForCurvature, MechPlasticIterInf &nlInfPlastic)
{
    // расчёт результатов на каждом конечном элементе
    POINT3_CUBE XYZ;    // координаты точки в шаблонном кубе
    VECTOR3 du[3];  // производные перемещений для центра конечного элемента
    int vi[8];  // глобальные индексы вершин
    feEl->getVertexIndexes(vi);

    // вычисление перемещений в каждой вершине
    if(feEl->get_FEType() == Grid::FEType::LinearHexagon)
    {
        // обход вершин
        for(int t = 0; t < 8; t++)
        {
            // XYZ - координаты вершины шаблонного куба (-1 или +1 каждая координата)
            XYZ[0] = (t%2)*2-1;
            XYZ[1] = ((t/2)%2)*2-1;
            XYZ[2] = (t/4)*2-1;
            // в вершине XYZ вычисляем сумму значений базисных функций
            double dux = 0;
            double duy = 0;
            double duz = 0;
            for(int localFuncIndex = 0; localFuncIndex < 8; localFuncIndex++)
            {
                double basisFuncValue = Fem::cubeLagrange1_3D(XYZ, localFuncIndex, Fem::dif_NULL3);
                VECTOR3_int DOFindex = DOFs->findDOFIndex_Lagr1(vi[localFuncIndex]);
                dux += basisFuncValue * q[DOFindex[0]];
                duy += basisFuncValue * q[DOFindex[1]];
                duz += basisFuncValue * q[DOFindex[2]];
            }
            // индекс вершины XYZ
            int vertexInd = vi[t];
            // сохраняем результат
            (*vertex)[vertexInd].du[0] = dux;
            (*vertex)[vertexInd].du[1] = duy;
            (*vertex)[vertexInd].du[2] = duz;
        }
    }
    if(feEl->get_FEType() == Grid::FEType::QuadraticHexagon)
    {
        int vi_geom[27];
        feEl->getGeomVertexIndexes(vi_geom);
        // обход вершин
        for(int t = 0; t < 27; t++)
        {
            // XYZ - координаты вершины шаблонного куба (-1 или +1 каждая координата)
            XYZ[0] = (t%3)-1;
            XYZ[1] = ((t/3)%3)-1;
            XYZ[2] = (t/9)-1;
            // в вершине XYZ вычисляем сумму значений базисных функций
            VECTOR3 du = VECTOR3_NULL;
            for(int localFuncIndex = 0; localFuncIndex < 8; localFuncIndex++)
            {
                double basisFuncValue = Fem::cubeLagrange1_3D(XYZ, localFuncIndex, Fem::dif_NULL3);
                VECTOR3_int DOFindex = DOFs->findDOFIndex_Lagr1(vi[localFuncIndex]);
                du[0] += basisFuncValue * q[DOFindex[0]];
                du[1] += basisFuncValue * q[DOFindex[1]];
                du[2] += basisFuncValue * q[DOFindex[2]];

            }
            int vertexInd = vi_geom[t];
            if(t == 0 ||
               t == 2 ||
               t == 6 ||
               t == 8 ||
               t == 18 ||
               t == 20 ||
               t == 24 ||
               t == 26)
            {
                // 1) данная точка - одна из 8-ми вершин куба
                // глобальный индекс вершины XYZ, он же - индекс соответствующей базисной функции
                // сохраняем результат
                (*vertex)[vertexInd].du = du;
            }
            else
            {
                // 2) данная точка - одна из (27-8)-ми точек, не вершин куба
                // в этом случае индекс вершины - индекс массива vertexForCurvature
                // сохраняем результат
                (*vertexForCurvature)[vertexInd].du = du;
            }
        }
    }

    // напряжения в центре конечного элемента k
    XYZ = VECTOR3_NULL;     // центр конечного элемента в координатах шаблонного куба
    for(int i = 0; i < 3; i++)
        du[i] = VECTOR3_NULL;

    // копии координат вершин 6-гранника
    POINT3 v[27];
    feEl->getGeomVertexes(grid->vertex, grid->vertexForCurvature, v);

    // вычисление производных
    for (int i = 0; i < 8; i++) // i - локальный номер вершины
    {
        VECTOR3 q0;
        VECTOR3 dbas;
        /*
        int global_i = vi[i];
        for(int j = 0; j < 3; j++)
            q0[j] = q[3 * global_i + j];
        */
        const VECTOR3_int &DOFIndex = DOFs->findDOFIndex_Lagr1(vi[i]);
        for(int j = 0; j < 3; j++)
            q0[j] = q[DOFIndex[j]];
            //q0[j] = q[vertexBusFuncInd->fs[vi[i]][0].DOF_index[j]];
            //q0[j] = q[vertexBusFuncInd->linear[vi[i]][j]];
        // производные вектора перемещений по x,y,z в точке XYZ
        // линейное отображение
        if(feEl->get_FEType() == Grid::FEType::LinearHexagon)
            Fem::lagrange1Hexagon_difLagrange1(v, XYZ, i, dbas);
        // квадратичное отображение
        if(feEl->get_FEType() == Grid::FEType::QuadraticHexagon)
            Fem::lagrange2Hexagon_difLagrange1(v, XYZ, i, dbas);
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                du[i][j] += dbas[j] * q0[i];
    }

    // вычисление результатов итерации
    pd.calcIterResults(m0, du, tl);

    // расчёт количества КЭ, для которых неизбежно придётся менять пластичную матрицу на упругую
    if(pd.it.preliminarily_needUpdateMatrix)
    {
        nlInfPlastic.preliminarily_GLocalChangedNumber++;
    }
}

void MechFeData_LinearHexagon_homogeny_E_NLE_EP::calcIterResults2(const Grid::TimeLayers &tl, MechMaterialSource &m0,
                                                                   MechPlasticIterInf &nlInfPlastic)
{
    // расчёт нового нового приближения для следующей итерации
    pd.prepareForNextIter(m0, nlInfPlastic.iterNumber, tl, nlInfPlastic.preliminarily_GLocalChangedNumber);

    switch (m0.plasticityMethodType)
    {
    // линейная упругость
    case MechPlasticityMethodType::Elasticity:
    // пластичность
    case MechPlasticityMethodType::D_pl:
    // матрица D_pl, алгоритм Соловейчика
    case MechPlasticityMethodType::D_pl_Solov:
    // метод начальных напряжений
    case MechPlasticityMethodType::InitialSigma:
    // матрица D_pl + метод начальных напряжений
    case MechPlasticityMethodType::Combo_D_pl_InitialSigma:
    {
        maxPlasticResidual.set_abs(pd.it.residual);
        needUpdateMatrix = pd.it.needUpdateMatrix;
        plasticNeedIterations = pd.it.needIterations;
        GLocalChanged = pd.it.needUpdateMatrix;
        isUnloading = (pd.it.lb.isUnloading && pd.q > 0);
        if(plasticNeedIterations == true)
        {
            // на следующей итерации требуется обновить вектор
            bLocalChanged = true;
        }
        else
        {
            bLocalChanged = false;
        }
    }break;
    }
    nlInfPlastic.addFeIterResultsToCounts(this);
}

FILE *f_debug = fopen("__plasticDebug.txt", "w");
void MechFePointData::debug(MechMaterialSource &m0) const
{
}
void MechFeData_LinearHexagon_homogeny_E_NLE_EP::debug(MechMaterialSource &m0)
{
    pd.debug(m0);
}








void MechFeData_LinearXFEM_gen::G_local_clear()
{
    for(size_t i = 0; i < GLocalL.size(); i++)
        GLocalL[i] = 0;
}
void MechFeData_LinearXFEM_gen::b_local_clear()
{
    for(size_t i = 0; i < bLocal.size(); i++)
        b_local_el(i) = 0;
}
void MechFeData_LinearXFEM_gen::R_local_clear()
{
    for(size_t i = 0; i < RLocal.size(); i++)
        R_local_el(i) = 0;
}
void MechFeData_LinearXFEM_gen::updateIndexes(const Grid::FE_base *feEl, const Grid::GlobalDOFs *DOFs)
{
    // локальные и глобальные индексы
    if(indexesChanged)
    {
        const Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<const Grid::FE_LinearHexagon_XFEM *>(feEl);
        // копии глобальных индексов и координат вершин 6-гранника
        int vi[8];
        feEl->getVertexIndexes(vi);
        // выделение памяти для local_DOF_ID
        {
            size_t local_DOF_ID_count = 0;
            for (size_t m = 0; m < 8; m++) // индекс вершины КЭ
            {
                local_DOF_ID_count += DOFs->DOFs[vi[m]].size();
            }
            local_DOF_ID.resize(local_DOF_ID_count);
        }
        // заполнение local_DOF_ID
        {
            size_t local_DOF_ID_index = 0;
            for (size_t m = 0; m < 8; m++) // индекс вершины КЭ
            {
                for(size_t DOFs_index = 0; DOFs_index < DOFs->DOFs[vi[m]].size(); DOFs_index++)    // индекс степени свободы узла vi[m]
                {
                    const Grid::VertexDOF &vd = DOFs->DOFs[vi[m]][DOFs_index];
                    local_DOF_ID[local_DOF_ID_index].m = m;
                    local_DOF_ID[local_DOF_ID_index].DOFs_index = DOFs_index;
                    local_DOF_ID[local_DOF_ID_index].funcID_index = 0;
                    for(size_t funcID_index = 0; funcID_index < feEl_XFEM->funcsID.size(); funcID_index++)
                    {
                        if(feEl_XFEM->funcsID[funcID_index] == vd.funcID)
                        {
                            local_DOF_ID[local_DOF_ID_index].funcID_index = funcID_index;
                            break;
                        }
                    }
                    local_DOF_ID_index++;
                }
            }
        }
        // заполнение local_mn, local_ik, globalIndex
        {
            local_mn.resize(local_DOF_ID.size()*3);
            local_ik.resize(local_DOF_ID.size()*3);
            globalIndex.resize(local_DOF_ID.size()*3);
            size_t local_ind = 0;
            for (size_t mn = 0; mn < local_DOF_ID.size(); mn++) // индекс базисной ф-и
            {
                const Grid::VertexDOF &vd = DOFs->DOFs[vi[local_DOF_ID[mn].m]][local_DOF_ID[mn].DOFs_index];
                for (int ik = 0; ik < 3; ik++)
                {
                    local_mn[local_ind] = mn;
                    local_ik[local_ind] = ik;
                    globalIndex[local_ind] = vd.DOF_index[ik];
                    local_ind++;
                }
            }
        }
        // сортировка по возрастанию значений глобальных индексов
        for (size_t local_ind1 = 0; local_ind1 < local_mn.size(); local_ind1++)
            for (size_t local_ind2 = local_ind1 + 1; local_ind2 < local_mn.size(); local_ind2++)
            {
                if(globalIndex[local_ind1] > globalIndex[local_ind2])
                {
                    std::swap(local_mn[local_ind1], local_mn[local_ind2]);
                    std::swap(local_ik[local_ind1], local_ik[local_ind2]);
                    std::swap(globalIndex[local_ind1], globalIndex[local_ind2]);
                }
            }
        // выделение памяти для локальных матриц и векторов
        GLocalL.resize(((globalIndex.size() + 1)*globalIndex.size()/2));
        SSCM_index.resize(((globalIndex.size() + 1)*globalIndex.size()/2));
        bLocal.resize(globalIndex.size());
        RLocal.resize(globalIndex.size());
        R_local_clear();
        // индексы обновлены, значит остальное тоже требуется обновить
        indexesChanged = false;
        tablesChanged = true;
        GLocalChanged = true;
        bLocalChanged = true;
    }
}

void MechFeData_LinearXFEM_gen::add_BC2_to_localb(const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechBoundaryCondition2Source_base *bc2SourceEl, const int subSurfaceIndex, const Integration::Integrator &integrationFoursquare)
{
    if(bLocalChanged)
    {
        double basCube[Grid::FE_LinearHexagon_XFEM::DOF_max*Integration::Integration2D_Size_Max];
        double hexagonCoef[Integration::Integration2D_Size_Max];
        VECTOR3 hexagonNormal[Integration::Integration2D_Size_Max];
        const Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<const Grid::FE_LinearHexagon_XFEM *>(feEl);
        const std::vector<Grid::SubHexFace> &subHexFace = feEl_XFEM->subHexFace[subSurfaceIndex];
        for(size_t subHexFaceInd = 0; subHexFaceInd < subHexFace.size(); subHexFaceInd++)
        {
            feEl_XFEM->gn_calcBc2BasisFuncValues(local_DOF_ID, integrationFoursquare, subHexFace[subHexFaceInd].hexInd, subHexFace[subHexFaceInd].faceIndex,
                                                 basCube, hexagonCoef, hexagonNormal);
            // добавка к локальному вектору
            for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
            {
                int m = local_mn[local_mi_ind];
                int i = local_ik[local_mi_ind];
                double E = 0;
                for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
                {
                    VECTOR3 vectorP =
                            bc2SourceEl->calcVectorValueInPoint2(hexagonNormal[valIndex]) -
                            bc2SourceEl->calcVectorValueInPoint1(hexagonNormal[valIndex]);
                    E += vectorP[i] *
                         basCube[m*integrationFoursquare.size + valIndex]*
                         hexagonCoef[valIndex];
                }
                b_local_el(local_mi_ind) += E;
            }
        }



        /*double basCube[8*Integration::Integration2D_Size_Max];
        double hexagonCoef[Integration::Integration2D_Size_Max];
        VECTOR3 hexagonNormal[Integration::Integration2D_Size_Max];
        feEl->calcBc2BasisFuncValues(grid->vertex, grid->vertexForCurvature, integrationFoursquare, faceIndex,
                                      basCube, hexagonCoef, hexagonNormal);
        // добавка к локальному вектору
        // базисные ф-и 1-го порядка
        for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
        {
            int m = local_mn[local_mi_ind];
            int i = local_ik[local_mi_ind];
            double E = 0;
            for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
            {
                VECTOR3 vectorP =
                        bc2SourceEl->calcVectorValueInPoint2(hexagonNormal[valIndex]) -
                        bc2SourceEl->calcVectorValueInPoint1(hexagonNormal[valIndex]);
                E += vectorP[i] *
                     basCube[m*integrationFoursquare.size + valIndex]*
                     hexagonCoef[valIndex];
            }
            b_local_el(local_mi_ind) += E;
        }
        // XFEM
        for (int local_mi_ind = 24; local_mi_ind < 144; local_mi_ind++)
        {
            b_local_el(local_mi_ind) += 0;
        }*/
    }
};
void MechFeData_LinearXFEM_gen::add_BC2_to_localR(const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechBoundaryCondition2Source_base *bc2SourceEl, const int faceIndex, const Integration::Integrator &integrationFoursquare, const IncForsesMode incForsesMode)
{
    if(incForsesMode == IncForsesMode::IncrementP)
        return;
    double basCube[8*Integration::Integration2D_Size_Max];
    double hexagonCoef[Integration::Integration2D_Size_Max];
    VECTOR3 hexagonNormal[Integration::Integration2D_Size_Max];
    feEl->calcBc2BasisFuncValues(grid->vertex, grid->vertexForCurvature, integrationFoursquare, faceIndex,
                                  basCube, hexagonCoef, hexagonNormal);
    // добавка к локальному вектору невязки
    // базисные ф-и 1-го порядка
    for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
    {
        int m = local_mn[local_mi_ind];
        int i = local_ik[local_mi_ind];
        double E = 0;
        for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
        {
            VECTOR3 vectorP;
            switch (incForsesMode)
            {
            case IncForsesMode::MinusIntegral:
            {
                // вычисляется в начале шага, после сдвига сетки
                // со знаком +
                vectorP = bc2SourceEl->calcVectorValueInPoint1(hexagonNormal[valIndex]);
                E += vectorP[i] *
                     basCube[m*integrationFoursquare.size + valIndex]*
                     hexagonCoef[valIndex];
            }break;
            case IncForsesMode::bPlusR:
            {
                // вычисляется в конце шага, после всех итераций, до сдвига сетки
                // со знаком +
                vectorP = bc2SourceEl->calcVectorValueInPoint2(hexagonNormal[valIndex]);
                E += vectorP[i] *
                     basCube[m*integrationFoursquare.size + valIndex]*
                     hexagonCoef[valIndex];
            }break;
            case IncForsesMode::IncrementP:
            {
            }break;
            }
        }
        R_local_el(local_mi_ind) += E;
    }
    // XFEM
    /*
    for (int local_mi_ind = 24; local_mi_ind < 144; local_mi_ind++)
    {
        R_local_el(local_mi_ind) += 0;
    }*/

};
void MechFeData_LinearXFEM_gen::addGlobalIndexesToGlobalPortrait1x1(SlauSolving::SSCMPortraitBulder &pb)
{
    for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
    {
        for (size_t local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
        {
            int global_mi = globalIndex[local_mi_ind];
            int global_nk = globalIndex[local_nk_ind]; // local_nk_ind < local_mi_ind
            pb.addElement(global_mi, global_nk); // global_nk < global_mi
        }
    }
}
void MechFeData_LinearXFEM_gen::updateSSCM_indexes1x1(const SlauSolving::SSCMPortrait &portrait)
{
    int G_local_ind = 0;
    for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
    {
        for (size_t local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
        {
            int global_mi = globalIndex[local_mi_ind];
            int global_nk = globalIndex[local_nk_ind];
            if(global_mi == global_nk)
                SSCM_index[G_local_ind] = -1;//######
            else
                SSCM_index[G_local_ind] = portrait.findSorted(global_mi, global_nk);
                        //portrait.findUnsorted(global_mi, global_nk);
            G_local_ind++;
        }
    }
}
void MechFeData_LinearXFEM_gen::addLocalGM_to_global1x1(SlauSolving::SSCMElements &elements)
{
    int G_local_ind = 0;
    for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
    {
        for (size_t local_nk_ind = 0; local_nk_ind < local_mi_ind; local_nk_ind++)
        {
            elements.a[SSCM_index[G_local_ind]] += G_local_el(G_local_ind);
            G_local_ind++;
        }
        elements.d[globalIndex[local_mi_ind]] += G_local_el(G_local_ind);
        G_local_ind++;
    }
}
void MechFeData_LinearXFEM_gen::addLocalb_to_global1x1(Vector &b)
{
    for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
    {
        b[globalIndex[local_mi_ind]] += b_local_el(local_mi_ind) + R_local_el(local_mi_ind);
    }
}

void MechFeData_LinearXFEM_gen::addGlobalIndexesToGlobalPortrait3x3(SlauSolving::SSCMPortraitBulder &pb)
{
    for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
    {
        for (size_t local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
        {
            int global_m = globalIndex[local_mi_ind]/3;
            int global_n = globalIndex[local_nk_ind]/3; // local_nk_ind <= local_mi_ind
            pb.addElement(global_m, global_n); // global_n <= global_m
        }
    }
}
void MechFeData_LinearXFEM_gen::updateSSCM_indexes3x3(const SlauSolving::SSCMPortrait &portrait)
{
    int G_local_ind = 0;
    for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
    {
        for (size_t local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
        {
            int global_m = globalIndex[local_mi_ind]/3;
            int global_n = globalIndex[local_nk_ind]/3;
            if(global_m == global_n)
                SSCM_index[G_local_ind] = -1;//######
            else
                SSCM_index[G_local_ind] = portrait.findSorted(global_m, global_n);
            G_local_ind++;
        }
    }
}
void MechFeData_LinearXFEM_gen::addLocalGM_to_global3x3(SlauSolving::SSCM3x3Elements &elements)
{
    int G_local_ind = 0;
    for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
    {
        for (size_t local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
        {
            size_t bi = globalIndex[local_mi_ind]%3;
            size_t bj = globalIndex[local_nk_ind]%3;
            if(SSCM_index[G_local_ind] != (size_t)-1)
            {
                // ниже диагонали
                MATR3x3 &el = elements.a[SSCM_index[G_local_ind]];
                el.m[bi][bj] += G_local_el(G_local_ind);
            }
            else
            {
                // диагональ
                int global_m = globalIndex[local_mi_ind]/3;
                MATR3x3 &el = elements.d[global_m];
                if(bi == bj)
                {
                    // диагональ блока
                    el.m[bi][bj] += G_local_el(G_local_ind);
                }
                else
                {
                    // вне диагонали
                    el.m[bi][bj] += G_local_el(G_local_ind);
                    el.m[bj][bi] += G_local_el(G_local_ind);
                }
            }
            G_local_ind++;
        }
    }
}
void MechFeData_LinearXFEM_gen::addLocalb_to_global3x3(Vector3 &b)
{
    for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
    {
        int global_m = globalIndex[local_mi_ind]/3;
        size_t bi = globalIndex[local_mi_ind]%3;
        b[global_m][bi] += b_local_el(local_mi_ind) + R_local_el(local_mi_ind);
    }
}

bool MechFeData_LinearXFEM_gen::localb_forces_is_null()
{
    for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
    {
        if(b_local_el(local_mi_ind) != 0)
            return false;
    }
    return true;
}
size_t MechFeData_LinearXFEM_gen::getMinGlobalIndex() const
{
    return globalIndex[0];
}

void MechFeData_LinearXFEM_gen::updateTables(const Grid::FE_base *feEl)
{
    // Таблицы коэффииентов для расчёта элементов локальных матриц и векторов
    // Для однородного материала константы материала можно вынести за знак интеграла
    if(tablesChanged)
    {
        // внутренние интегралы для добавок к матрицам и к правой части
        // для матрицы G
        //const Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<const Grid::FE_LinearHexagon_XFEM *>(feEl);
        //feEl_XFEM->gn_calcIntegralsTables(Gauss3_integrationwSource, Gauss3_dLinearBasCubeSource, local_DOF_ID,
        //                                  hexagonIntTable);
        tablesChanged = false;
    }
}
void MechFeData_LinearXFEM_gen::init(Integration::IntegrationType set_integrationType, const Grid::FE_base *feEl)
{
    indexesChanged = true;
    tablesChanged = true;
    GLocalChanged = true;
    bLocalChanged = true;
    needUpdateMatrix = false;
    plasticNeedIterations = false;
    isUnloading = false;
    integrationType = set_integrationType;
    const Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<const Grid::FE_LinearHexagon_XFEM *>(feEl);
    pd.resize(feEl_XFEM->subHex.size());
    //hexagonIntTable.resize(feEl_XFEM->subHex.size());
    for(size_t hexInd = 0; hexInd < pd.size(); hexInd++)
        pd[hexInd].init();
}
void MechFeData_LinearXFEM_gen::setT(const double T, const double deltaT)
{
    for(size_t hexInd = 0; hexInd < pd.size(); hexInd++)
    {
        pd[hexInd].T = T;
        pd[hexInd].deltaT = deltaT;
    }
}
void MechFeData_LinearXFEM_gen::finalizeStep(const MechMaterialSource &m0, const int movingGridMode)
{
    for(size_t hexInd = 0; hexInd < pd.size(); hexInd++)
        pd[hexInd].finalizeStep(m0, movingGridMode);
    // на нулевой итерации следующего шага требуется обновить матрицу
    tablesChanged = true;
    GLocalChanged = true;
    bLocalChanged = true;
}
void MechFeData_LinearXFEM_gen::saveResultsOfStep(const Vector &q, const Grid::Grid3D *grid, const Grid::FE_base *feEl, MechOutFeData &outFeDataEl)
{
    const Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<const Grid::FE_LinearHexagon_XFEM *>(feEl);
    // копии глобальных индексов и координат вершин 6-гранника
    //POINT3 v[8];
    //feEl_XFEM->getGeomVertexes(grid->vertex, grid->vertexForCurvature, v);
    // копируем информацию в центрах 6-гранных подобластей
    outFeDataEl.pd.resize(pd.size());
    for(size_t i = 0; i < pd.size(); i++)
    {
        MechOutFePointData &outFePD = outFeDataEl.pd[i];  // данные вывода для точки конечного элемента номер i
        pd[i].saveResultsOfStep(outFePD);
        POINT3 XYZ = POINT3(0,0,0);
        POINT3 xyz;
        POINT3 subv[8];
        feEl_XFEM->gn_getSubVertexes(i, subv);
        feEl->cubeToHexagon(subv, XYZ, Fem::dif_NULL3, xyz);
        outFePD.p = xyz;
    }
    // множители перед базисными функциями
    outFeDataEl.q.resize(local_DOF_ID.size());
    for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
    {
        int m = local_mn[local_mi_ind];
        int i = local_ik[local_mi_ind];
        outFeDataEl.q[m][i] = q[globalIndex[local_mi_ind]];
    }
}

void MechFeData_LinearXFEM_gen::updateLocalGMb(const Grid::Grid3D *grid, const Grid::GlobalDOFs *DOFs, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0, const MechPlasticIterInf &nlInfPlastic, const bool firstIter, const IncForsesMode incForsesMode)
{
    // локальные и глобальные индексы
    updateIndexes(feEl, DOFs);
    // таблицы с коэффициентами для текущей формы шестигранника
    updateTables(feEl);
    if(GLocalChanged)
    {
        // очистка локальной матрицы G
        G_local_clear();
    }
    // III добавка к локальному вектору
    if(bLocalChanged)
    {
        // очистка локального вектора b
        b_local_clear();
        R_local_clear();
    }
    for(size_t hexInd = 0; hexInd < pd.size(); hexInd++) // индекс 6-гранной подобласти КЭ
    {
        MechFePointData &pd_el = pd[hexInd];
        // подготовка параметров для итерации
        pd_el.initNewIteration(m0, tl, firstIter, nlInfPlastic.iterationMode);
        // расчёт матрицы определяющих соотношений
        TENSOR4 C_pl0;  // матрица D с волной, записанная в виде тензора
        if(m0.plasticityMethodType == MechPlasticityMethodType::D_pl_Solov)
        {
            C_pl0 = pd_el.it.C_pl;
        }
        else
        {
            TENSOR4 C_el0;
            TENSOR4 Y_pl0;
            TENSOR4::setCY_pl(pd_el.elasticParameters_T, pd_el.it.A, pd_el.it.z,
                              C_el0, C_pl0, Y_pl0);
        }
        const double deltaT = pd_el.deltaT;

        // внутренние интегралы для добавок к матрицам и к правой части
        // для матрицы G
        const Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<const Grid::FE_LinearHexagon_XFEM *>(feEl);
        Grid::HexagonIntTable hexagonIntTable;
        feEl_XFEM->gn_calcIntegralsTable(Gauss3_integrationwSource, Gauss3_dLinearBasCubeSource, local_DOF_ID, hexInd,
                                          hexagonIntTable);

        // I. добавка к матрице G
        if(GLocalChanged)
        {
            using namespace Integration;
            int G_local_ind = 0;
            for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
            {
                int m = local_mn[local_mi_ind];
                int i = local_ik[local_mi_ind];
                for (size_t local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
                {
                    int n = local_mn[local_nk_ind];
                    int k = local_ik[local_nk_ind];
                    double E = 0;
                    for (int j = 0; j < 3; j++)
                        for (int l = 0; l < 3; l++)
                            E += C_pl0.m[i][j][k][l] * hexagonIntTable.intForG_mnjl[m][n][j][l];
                    G_local_el(G_local_ind) += E;
                    G_local_ind++;
                }
            }
        }
        // III добавка к локальному вектору
        if(bLocalChanged)
        {
            MATR3x3x3x3 C_el;  // упругая матрица D, записанная в виде тензора
            {
                MATR6x6 D_el;
                m0.calcD(pd_el.elasticParameters_T, D_el);
                Elementary::Operations::convert6x6to3x3x3x3(D_el, C_el);
            }
            // заполнение локального вектора b
            // III.3 добавка к локальному вектору (температурные слагаемые: alpha*deltaT*... и sigmakl*dAijkl...)
            if(deltaT != 0)
            {
                // добавка в правую часть узлавых усилий, соответствующих температурным деформациям
                for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
                {
                    int m = local_mn[local_mi_ind];
                    int i = local_ik[local_mi_ind];
                    double E = 0;
                    for (int j = 0; j < 3; j++)
                        for (int k = 0; k < 3; k++)
                            for (int l = 0; l < 3; l++)
                                E += C_el.m[i][j][k][l] * hexagonIntTable.intForb23_mj[m][j] *
                                        pd_el.depsTerm.m[k][l];
                    b_local_el(local_mi_ind) += E;
                }
            }
            // III.5 добавка к локальному вектору (начальные напряжения)
            if(pd_el.it.dsigma0.eqv_SIGMA() != 0)
            {
                for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
                {
                    int m = local_mn[local_mi_ind];
                    int i = local_ik[local_mi_ind];
                    double E = 0;
                    for (int j = 0; j < 3; j++)
                        E += pd_el.it.dsigma0.m[i][j] * hexagonIntTable.intForb23_mj[m][j];
                    b_local_el(local_mi_ind) -= E;
                }
            }
            // невязка R
            if(incForsesMode == IncForsesMode::MinusIntegral)
            {
                //MechFePointData &pd_el = pd[hexInd];
                for (size_t local_mi_ind = 0; local_mi_ind < globalIndex.size(); local_mi_ind++)
                {
                    int m = local_mn[local_mi_ind];
                    int i = local_ik[local_mi_ind];
                    double E = 0;
                    for (int j = 0; j < 3; j++)
                    {
                        switch (incForsesMode)
                        {
                        case IncForsesMode::MinusIntegral:
                        {
                            // вычисляется в начале шага, после сдвига сетки
                            // со знаком -
                            E -= pd_el.sumSigma.m[i][j] * hexagonIntTable.intForb23_mj[m][j];
                        }break;
                        case IncForsesMode::bPlusR:
                        {
                            // вычисляется в конце шага, после всех итераций, после приращения напряжений, до сдвига сетки
                            // со знаком -
                            E -= pd_el.sumSigma.m[i][j] * hexagonIntTable.intForb23_mj[m][j];
                        }break;
                        case IncForsesMode::IncrementP:
                        {
                        }break;
                        }
                    }
                    R_local_el(local_mi_ind) += E;
                }
            }

        }
    }
    GLocalChanged = false;
}
void MechFeData_LinearXFEM_gen::calcLocalR(const Grid::Grid3D *grid, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0, const MechPlasticIterInf &nlInfPlastic, const IncForsesMode incForsesMode)
{
}
void MechFeData_LinearXFEM_gen::calcIterResults1(const Vector &q, const Grid::Grid3D *grid, const Grid::GlobalDOFs *DOFs, const Grid::TimeLayers &tl, const Grid::FE_base *feEl, MechMaterialSource &m0,
                                             std::vector<MechVertexData> *vertex, std::vector<MechVertexData> *vertexForCurvature, MechPlasticIterInf &nlInfPlastic)
{
    // расчёт результатов на каждом конечном элементе
    POINT3_CUBE XYZ;    // координаты точки в шаблонном кубе
    int vi[8];  // глобальные индексы вершин
    feEl->getVertexIndexes(vi);

    // вычисление перемещений в каждой вершине
    for(int t = 0; t < 8; t++) // t - локальный индекс вершины
    {
        // XYZ - координаты вершины шаблонного куба (-1 или +1 каждая координата)
        XYZ[0] = (t%2)*2-1;
        XYZ[1] = ((t/2)%2)*2-1;
        XYZ[2] = (t/4)*2-1;
        // в вершине XYZ вычисляем сумму значений базисных функций
        double dux = 0;
        double duy = 0;
        double duz = 0;
        for(int localFuncIndex = 0; localFuncIndex < 8; localFuncIndex++)
        {
            double basisFuncValue = Fem::cubeLagrange1_3D(XYZ, localFuncIndex, Fem::dif_NULL3);
            const VECTOR3_int &DOFindex = DOFs->findDOFIndex_Lagr1(vi[localFuncIndex]);
            //VECTOR3_int globalFuncIndex = vertexBusFuncInd->fs[vi[localFuncIndex]][0].DOF_index;//vertexBusFuncInd->linear[vi[localFuncIndex]];
            dux += basisFuncValue * q[DOFindex[0]];
            duy += basisFuncValue * q[DOFindex[1]];
            duz += basisFuncValue * q[DOFindex[2]];
        }
        // глобальный индекс вершины XYZ
        int vertexInd = vi[t];
        // сохраняем результат
        (*vertex)[vertexInd].du[0] = dux;
        (*vertex)[vertexInd].du[1] = duy;
        (*vertex)[vertexInd].du[2] = duz;
    }

    // центр 6-гранной подобласти КЭ в координатах шаблонного куба
    XYZ = VECTOR3_NULL;
    // вычисление результатов итерации в центрах 6-гранных подобластей
    const Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<const Grid::FE_LinearHexagon_XFEM *>(feEl);
    for(size_t hexInd = 0; hexInd < pd.size(); hexInd++) // индекс 6-гранной подобласти КЭ
    {
        MechFePointData &pd_el = pd[hexInd];
        // производные базисных функций в центре 6-гранной подобласти
        double dbas[Grid::FE_LinearHexagon_XFEM::DOF_max][3];
        feEl_XFEM->gn_calcSubDifBasFuncs(local_DOF_ID, hexInd, XYZ,
                                         dbas);
        // производные перемещений в центре 6-гранной подобласти
        VECTOR3 du[3];
        for(int i = 0; i < 3; i++)
            du[i] = VECTOR3_NULL;
        for(size_t mn = 0; mn < local_DOF_ID.size(); mn++)
        {
            const Grid::VertexDOF &vd = DOFs->DOFs[vi[local_DOF_ID[mn].m]][local_DOF_ID[mn].DOFs_index];
            for(int i = 0; i < 3; i++) // номер координаты перемещения
            {
                for(int j = 0; j < 3; j++) // номер координаты, по которой происходит диффиренцирование
                {
                    du[i][j] += dbas[mn][j] * q[vd.DOF_index[i]];
                }
            }
        }
        // вычисление результатов итерации
        pd_el.calcIterResults(m0, du, tl);
        // расчёт количества КЭ, для которых неизбежно придётся менять пластичную матрицу на упругую
        if(pd_el.it.preliminarily_needUpdateMatrix)
        {
            nlInfPlastic.preliminarily_GLocalChangedNumber++;
        }
    }
}
void MechFeData_LinearXFEM_gen::calcIterResults2(const Grid::TimeLayers &tl, MechMaterialSource &m0,
                                                                   MechPlasticIterInf &nlInfPlastic)
{
    for(size_t hexInd = 0; hexInd < pd.size(); hexInd++) // индекс 6-гранной подобласти КЭ
    {
        MechFePointData &pd_el = pd[hexInd];
        // расчёт нового нового приближения для следующей итерации
        pd_el.prepareForNextIter(m0, nlInfPlastic.iterNumber, tl, nlInfPlastic.preliminarily_GLocalChangedNumber);

        switch (m0.plasticityMethodType)
        {
        // линейная упругость
        case MechPlasticityMethodType::Elasticity:
        // пластичность
        case MechPlasticityMethodType::D_pl:
        // матрица D_pl, алгоритм Соловейчика
        case MechPlasticityMethodType::D_pl_Solov:
        // метод начальных напряжений
        case MechPlasticityMethodType::InitialSigma:
        // матрица D_pl + метод начальных напряжений
        case MechPlasticityMethodType::Combo_D_pl_InitialSigma:
        {
            maxPlasticResidual.set_abs(pd_el.it.residual);
            needUpdateMatrix = pd_el.it.needUpdateMatrix;
            plasticNeedIterations = pd_el.it.needIterations;
            GLocalChanged = pd_el.it.needUpdateMatrix;
            isUnloading = (pd_el.it.lb.isUnloading && pd_el.q > 0);
            if(plasticNeedIterations == true)
            {
                // на следующей итерации требуется обновить вектор
                bLocalChanged = true;
            }
            else
            {
                bLocalChanged = false;
            }
        }break;
        }
        nlInfPlastic.addFeIterResultsToCounts(this);
    }
}

void MechFeData_LinearXFEM_gen::moveSubVertexes(const Grid::Grid3D *grid, const Grid::GlobalDOFs *DOFs, const Vector &dq, Grid::FE_base *feEl)
{
    Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<Grid::FE_LinearHexagon_XFEM *>(feEl);
    int vi[8];  // глобальные индексы вершин
    feEl->getVertexIndexes(vi);
    VECTOR3 q_local[Grid::FE_LinearHexagon_XFEM::DOF_max];
    for (size_t mn = 0; mn < local_DOF_ID.size(); mn++) // индекс базисной ф-и
    {
        const Grid::VertexDOF &vd = DOFs->DOFs[vi[local_DOF_ID[mn].m]][local_DOF_ID[mn].DOFs_index];
        for(int i = 0; i < 3; i++) // номер координаты перемещения
        {
            q_local[mn][i] = dq[vd.DOF_index[i]];
        }
    }
    feEl_XFEM->gn_moveSubVertexes(local_DOF_ID, q_local);
}

}   // namespace Solid
