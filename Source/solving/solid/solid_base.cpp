#include "solid_base.h"

namespace Solid
{

const MATR6x6 M_sigma_XYZ =
{
    {
    {1,     -0.5,  -0.5,  0,  0,  0,},
    {-0.5,  1,     -0.5,  0,  0,  0,},
    {-0.5,  -0.5,  1,     0,  0,  0,},
    {0,     0,     0,     3,  0,  0,},
    {0,     0,     0,     0,  3,  0,},
    {0,     0,     0,     0,  0,  3,},
    }
};
const MATR6x6 M_sigma_XY =
{
    {
    {1,     -0.5,  0,     0,  0,  0,},
    {-0.5,  1,     0,     0,  0,  0,},
    {0,     0,     0,     0,  0,  0,},
    {0,     0,     0,     0,  0,  0,},
    {0,     0,     0,     0,  0,  0,},
    {0,     0,     0,     0,  0,  3,},
    }
};

void MaxSigmaEpsResidual::set_null()
{
    sigma = 0;
    eps = 0;
}
void MaxSigmaEpsResidual::set_undefined()
{
    sigma = -1;
    eps = -1;
}
void MaxSigmaEpsResidual::set_abs(const MaxSigmaEpsResidual &r)
{
    eps = fabs(r.eps);
    sigma = fabs(r.sigma);
}
void MaxSigmaEpsResidual::compareAndSaveIfLarger(const MaxSigmaEpsResidual &r)
{
    if(r.eps > eps)
        eps = r.eps;
    if(r.sigma > sigma)
        sigma = r.sigma;
}


void TanAngle_Yeld::setValue(const double tan_elastic, const double dsigma_want, const double deps_want)
{
    // разбор случаев
    if(deps_want == 0)
    {
        if(dsigma_want == 0)
        {
            // отсутствует движение по кривой - пассивное нагружение (упругая матрицы)
            is_elastic = true;
        }
        else
        {
            // tan = бесконечность
            is_elastic = true;    // делаем упругость##
        }
    }
    else
    {
        // обычное движение по кривой
        double tan = dsigma_want / deps_want;
        if(tan == tan_elastic)
        //if(fabs((tan - tan_elastic) / tan_elastic) < 1.e-5)
        {
            is_elastic = true;
        }
        else
        {
            is_elastic = false;
            value = (tan_elastic*tan)/(tan_elastic - tan);
        }
    }
}

/*double A = 1*1.e6, epsEqv00 = 0.0005, k00 = +0.5;//-0.001;
double df_min = 1.e-10;//1
double C = 1*1.e-8;

double c_A = 1.e-14;
double c_n = 1.5;*/
//double c_n = 2;

MechMaterialSource::MechMaterialSource()
{
    plasticityMethodType = MechPlasticityMethodType::Elasticity;
    elasticParameters0.clear();
    temperatureDependence = false;
    is2Dxy = false;
}
double MechMaterialSource::eps_dep(const double sigmaEqv, const double t)
{
    switch (PCDependenceType)
    {
    case MechPlasticityCurveDependenceType::None:
    {
        return epsFun.calcValue(sigmaEqv);
    }break;
    case MechPlasticityCurveDependenceType::Time:
    {
        return epsFun.calcValue(sigmaEqv, t);
    }break;
    }
}
double MechMaterialSource::difEps_dep(const double sigmaEqv, const double t)
{
    switch (PCDependenceType)
    {
    case MechPlasticityCurveDependenceType::None:
    {
        return difEpsFun.calcValue(sigmaEqv);
    }break;
    case MechPlasticityCurveDependenceType::Time:
    {
        return difEpsFun.calcValue(sigmaEqv, t);
    }break;
    }
}
double MechMaterialSource::sigma_dep(const double epsEqv, const double t)
{
    switch (PCDependenceType)
    {
    case MechPlasticityCurveDependenceType::None:
    {
        return sigmaFun.calcValue(epsEqv);
    }break;
    case MechPlasticityCurveDependenceType::Time:
    {
        return sigmaFun.calcValue(epsEqv, t);
    }break;
    }
}
double MechMaterialSource::difSigma_dep(const double epsEqv, const double t)
{
    switch (PCDependenceType)
    {
    case MechPlasticityCurveDependenceType::None:
    {
        return difSigmaFun.calcValue(epsEqv);
    }break;
    case MechPlasticityCurveDependenceType::Time:
    {
        return difSigmaFun.calcValue(epsEqv, t);
    }break;
    }
}

double MechMaterialSource::eps(const ElasticParameters &ep_T, const double sigmaT, const double sigmaEqv, const double t)
{
    switch (PCUnloadingType)
    {
    // разгрузка по кривой
    case MechPlasticityCurveUnloadingType::Curve:
    {
        return eps_dep(sigmaEqv, t);
    }break;
    // разгрузка упругая
    case MechPlasticityCurveUnloadingType::Elastic:
    {
        if(sigmaEqv >= sigmaT)
            return eps_dep(sigmaEqv, t);
        else
            return eps_dep(sigmaT, t) + (sigmaEqv - sigmaT)/(tan_el(ep_T));
    }break;
    }
    //return (sigma_eqv / G) * (1 + g2*(sigma_eqv / G)*(sigma_eqv / G));
}
double MechMaterialSource::difEps(const ElasticParameters &ep_T, const double sigmaT, const double sigmaEqv, const double t)
{
    switch (PCUnloadingType)
    {
    // разгрузка по кривой
    case MechPlasticityCurveUnloadingType::Curve:
    {
        return difEps_dep(sigmaEqv, t);
    }break;
    // разгрузка упругая
    case MechPlasticityCurveUnloadingType::Elastic:
    {
        if(sigmaEqv >= sigmaT)
            return difEps_dep(sigmaEqv, t);
        else
            return tan_el(ep_T);
    }break;
    }
    //return 1. / G + (3. / G)*g2*(sigma_eqv / G)*(sigma_eqv / G);
}
double MechMaterialSource::sigma(const ElasticParameters &ep_T, const double epsT, const double epsEqv, const double t)
{
    switch (PCUnloadingType)
    {
    // разгрузка по кривой
    case MechPlasticityCurveUnloadingType::Curve:
    {
        return sigma_dep(epsEqv, t);
    }break;
    // разгрузка упругая
    case MechPlasticityCurveUnloadingType::Elastic:
    {
        if(epsEqv >= epsT)
            return sigma_dep(epsEqv, t);
        else
            return sigma_dep(epsT, t) + (epsEqv - epsT)*(tan_el(ep_T));
    }break;
    }
}
double MechMaterialSource::difSigma(const ElasticParameters &ep_T, const double epsT, const double epsEqv, const double t)
{
    switch (PCUnloadingType)
    {
    // разгрузка по кривой
    case MechPlasticityCurveUnloadingType::Curve:
    {
        return difSigma_dep(epsEqv, t);
    }break;
    // разгрузка упругая
    case MechPlasticityCurveUnloadingType::Elastic:
    {
        if(epsEqv >= epsT)
            return difSigma_dep(epsEqv, t);
        else
            return tan_el(ep_T);
    }break;
    }
    //return 1/sigmaFun.solveDif(epsEqv);
}

double MechMaterialSource::sigma_Yeld(const ElasticParameters &ep_T, const double q, const double t)
{
    return sigma_dep(q, t);
}

double MechMaterialSource::difSigma_Yeld(const ElasticParameters &ep_T, const double q, const double t)
{
    return difSigma_dep(q, t);
}

double MechMaterialSource::q_Yeld(const ElasticParameters &ep_T, const double sigmaEqv, const double t)
{
    return eps_dep(sigmaEqv, t);
}

double MechMaterialSource::difq_Yeld(const ElasticParameters &ep_T, const double sigmaEqv, const double t)
{
    return difEps_dep(sigmaEqv, t);
}

double MechMaterialSource::tan_el(const ElasticParameters &ep_T)
{
    return sqrt(Csigma/Ceps)*2*ep_T.G;//3*ep_T.G;
}

void MechMaterialSource::calc_elastic_parameters(const double T, ElasticParameters &ep_T)
{
    ep_T.ro = elasticParameters0.ro;
    ep_T.Talpha = elasticParameters0.Talpha;
    if(temperatureDependence)
    {
        ep_T.E = E_Fun.calcValue(T);
        ep_T.NU = NU_Fun.calcValue(T);
        ep_T.K = ep_T.E / (3 * (1 - 2 * ep_T.NU));
        ep_T.LAMBDA = ep_T.E*ep_T.NU / (1. + ep_T.NU) / (1. - 2.*ep_T.NU);
        ep_T.G = ep_T.E / 2. / (1 + ep_T.NU);
    }
    else
    {
        ep_T.E = elasticParameters0.E;
        ep_T.NU = elasticParameters0.NU;
        ep_T.K = elasticParameters0.K;
        ep_T.LAMBDA = elasticParameters0.LAMBDA;
        ep_T.G = elasticParameters0.G;
    }
}
void MechMaterialSource::calcD(const ElasticParameters &ep_T, MATR6x6 &D)
{
    for (int k = 0; k < 6; k++)
        for (int l = 0; l < 6; l++)
            D.m[k][l] = 0;
    for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
            D.m[k][l] = ep_T.LAMBDA;
    for (int k = 0; k < 3; k++)
        D.m[k][k] += 2 * ep_T.G;
    for (int k = 3; k < 6; k++)
        D.m[k][k] = ep_T.G;
    if(is2Dxy)
        set_matrix_to_2D(D);
}
void MechMaterialSource::calcInvertIsotropicD(const ElasticParameters &ep_T, MATR3x3x3x3 &A)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    A.m[i][j][k][l] = 1/ep_T.E*((1 + ep_T.NU)/2*((i==k)*(j==l) + (i==l)*(j==k)) - ep_T.NU*(i==j)*(k==l));
}
void MechMaterialSource::calcInvertDeltaIsotropicD(const ElasticParameters &ep_T1, const ElasticParameters &ep_T2, MATR3x3x3x3 &deltaA)
{
// приращение
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                {
                    double a1 = 1/ep_T1.E*((1 + ep_T1.NU)/2*((i==k)*(j==l) + (i==l)*(j==k)) - ep_T1.NU*(i==j)*(k==l));
                    double a2 = 1/ep_T2.E*((1 + ep_T2.NU)/2*((i==k)*(j==l) + (i==l)*(j==k)) - ep_T2.NU*(i==j)*(k==l));
                    deltaA.m[i][j][k][l] = a2 - a1;
                }
}
void MechMaterialSource::set_E_NU(const double E, const double NU)
{
    elasticParameters0.K = E / (3 * (1 - 2 * NU));
    elasticParameters0.E = E;
    elasticParameters0.LAMBDA = E*NU / (1. + NU) / (1. - 2.*NU);
    elasticParameters0.G = E / 2. / (1 + NU);
    elasticParameters0.NU = NU;
}
void MechMaterialSource::set_K_G(const double K, const double G)
{
    elasticParameters0.K = K;
    elasticParameters0.E = 9 * K*G / (3 * K + G);
    elasticParameters0.LAMBDA = K - 2 * G / 3;
    elasticParameters0.G = G;
    elasticParameters0.NU = (3 * K - 2 * G) / 2 / (3 * K + G);
}
void MechMaterialSource::set_M_sigma()
{
    M_sigma = M_sigma_XYZ;
}
void MechMaterialSource::set_2Dxy()
{
    is2Dxy = true;
    M_sigma = M_sigma_XY;
}
void MechMaterialSource::set_matrix_to_2D(MATR6x6 &D)
{
    for (int k = 2; k <= 4; k++)
        for (int l = 0; l < 6; l++)
        {
            D.m[k][l] = 0;
            D.m[l][k] = 0;
        }
}

void MechBoundaryCondition2Source_None::updateValue(const double, const double)
{
};
VECTOR3 MechBoundaryCondition2Source_None::calcVectorValueInPoint1(const VECTOR3 &)
{
    return VECTOR3_NULL;
};
VECTOR3 MechBoundaryCondition2Source_None::calcVectorValueInPoint2(const VECTOR3 &)
{
    return VECTOR3_NULL;
};

void MechBoundaryCondition2Source_VectorConstant::updateValue(const double, const double)
{
};
VECTOR3 MechBoundaryCondition2Source_VectorConstant::calcVectorValueInPoint1(const VECTOR3 &)
{
    return P1;
}
VECTOR3 MechBoundaryCondition2Source_VectorConstant::calcVectorValueInPoint2(const VECTOR3 &)
{
    return P2;
}
void MechBoundaryCondition2Source_VectorConstant::init(const VECTOR3 set_P)
{
    P1 = set_P;
    P2 = set_P;
}

void MechBoundaryCondition2Source_ScalarConstant::updateValue(const double, const double)
{
};
VECTOR3 MechBoundaryCondition2Source_ScalarConstant::calcVectorValueInPoint1(const VECTOR3 &normal)
{
    return -normal * P1;    // давление действует против внешней нормали
};
VECTOR3 MechBoundaryCondition2Source_ScalarConstant::calcVectorValueInPoint2(const VECTOR3 &normal)
{
    return -normal * P2;    // давление действует против внешней нормали
};
void MechBoundaryCondition2Source_ScalarConstant::init(const double set_P)
{
    P1 = set_P;
    P2 = set_P;
}

void MechBoundaryCondition2Source_VectorFunction::updateValue(const double t_prev1, const double t0)
{
    P1 = value(t_prev1);
    P2 = value(t0);
}
VECTOR3 MechBoundaryCondition2Source_VectorFunction::value(const double t)
{
    VECTOR3 v;
    for(int i = 0; i < 3; i++)
    {
        (*Pf[i]).setArgumentValue(0, t);
        v[i] = (*Pf[i]).calcValue();
    }
    return v;
}
void MechBoundaryCondition2Source_VectorFunction::init(FunParser::Function *set_Pf1, FunParser::Function *set_Pf2, FunParser::Function *set_Pf3)
{
    Pf[0] = set_Pf1;
    Pf[1] = set_Pf2;
    Pf[2] = set_Pf3;
}

void MechBoundaryCondition2Source_ScalarFunction::updateValue(const double t_prev1, const double t0)
{
    P1 = value(t_prev1);
    P2 = value(t0);
}
double MechBoundaryCondition2Source_ScalarFunction::value(const double t)
{
    double v;
    (*Pf).setArgumentValue(0, t);
    v = (*Pf).calcValue();
    return v;
}
void MechBoundaryCondition2Source_ScalarFunction::init(FunParser::Function *set_Pf)
{
    Pf = set_Pf;
}

void MechOutFePointData::mainStressesXY(VECTOR3 &ms)
{
    double tauMax = 0.5*sqrt(SQR(sumSigma.m[0][0]-sumSigma.m[1][1]) + 4*SQR(sumSigma.m[0][1]));
    double sigma0 = (sumSigma.m[0][0] + sumSigma.m[1][1]) / 2;
    ms[0] = sigma0 + tauMax;
    ms[1] = sigma0;
    ms[2] = sigma0 - tauMax;
}
void MechOutFePointData::mainStresses(VECTOR3 &ms)const
{
    VECTOR3 I;
    sumSigma.calcInvariants(I);
    Operations::solvePolynom3(-I[0], I[1], -I[2], ms[0], ms[1], ms[2]);
}
void MechOutFeData::init(const int setSize)
{
    pd.resize(setSize);
}

}   // namespace Solid
