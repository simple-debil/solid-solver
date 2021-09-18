#ifndef _TRASH_CODE_H
#define _TRASH_CODE_H



#if false
/*
double Solid_Material::g(const double t0)
{
    //return eps_sigma(G * sqrt(t0)) / sqrt(t0);
    g_fun.setArgumentValue(0, t0);
    return g_fun.solve();
}

double Solid_Material::k(const double)
{
    return 1;
}*/
/*double Solver::solveAijPrism3(const int i, const int j, const int N, const int M, const int valIndex)const
{
//##можно интегрировать аналитически
    int nu_i = i % 3;
    int nu_j = j % 3;
    int Znomi = i / 3;
    int Znomj = j / 3;
    double E = 0, EE;
    POINT1 p = integration_interval.p1[valIndex];
    double fi_0 = Fem::dif_rectangle_linear_1d(op.z1, op.z2, p, Znomi, 0);
    double fi_1 = Fem::dif_rectangle_linear_1d(op.z1, op.z2, p, Znomi, 1);
    double fj_0 = Fem::dif_rectangle_linear_1d(op.z1, op.z2, p, Znomj, 0);
    double fj_1 = Fem::dif_rectangle_linear_1d(op.z1, op.z2, p, Znomj, 1);

    for (int k = 0; k < 3; k++)
    {
        for (int l = 0; l < 3; l++)
        {
            // обе производные не по z
            if (l <= 1 && k <= 1)
            {
                EE =
                    op.inverse_D.m[nu_j][l + 1] *
                    op.inverse_D.m[nu_i][k + 1] * op.Abs_det_D / 2 *
                    fi_0*fj_0;
            }
            else
            // одна производная не по z
            if (l <= 1 && k == 2)
            {
                EE =
                    op.inverse_D.m[nu_j][l + 1] * op.Abs_det_D / 6 *
                    fi_1*fj_0;
            }else
            if (l == 2 && k <= 1)
            {
                EE =
                    op.inverse_D.m[nu_i][k + 1] * op.Abs_det_D / 6 *
                    fi_0*fj_1;

            }else
            // обе производные по z
            {
                EE =
                    op.Abs_det_D / 24 *
                    fi_1*fj_1;
            }
            E += EE*C_pl.m[N][k][M][l];
        }
    }
    return E;
}*/
/*switch (g->fe[k].type)
{
case GridFEType::Prism:
{
    // напряжения в центре конечного элемента k
    dudx = VECTOR3_NULL;
    dudy = VECTOR3_NULL;
    dudz = VECTOR3_NULL;	// обнуление
    MATR3x3 inverse_D;
    double det;
    double x1, x2, x3, y1, y2, y3, z1, z2;
    POINT3 p;
    x1 = g->v[g->fe[k].vi[0]][0];
    x2 = g->v[g->fe[k].vi[1]][0];
    x3 = g->v[g->fe[k].vi[2]][0];
    y1 = g->v[g->fe[k].vi[0]][1];
    y2 = g->v[g->fe[k].vi[1]][1];
    y3 = g->v[g->fe[k].vi[2]][1];
    z1 = g->v[g->fe[k].vi[0]][2];
    z2 = g->v[g->fe[k].vi[3]][2];

    p[0] = (x1 + x2 + x3) / 3;
    p[1] = (y1 + y2 + y3) / 3;
    p[2] = (z1 + z2) / 2;
    solve_inverse_D(x1, x2, x3, y1, y2, y3, inverse_D, det);

    for (int i = 0; i < 6; i++)	// i - локальный номер вершины - локальной базисной функции
    {
        int global_i = globalVertexIndex(k, i);	// глобальный номер вершины с локальным номером i
        double q1, q2, q3;
        VECTOR3 dbas;
        q1 = x[3 * global_i + 0];
        q2 = x[3 * global_i + 1];
        q3 = x[3 * global_i + 2];
        // производные вектора перемещений по x,y,z
        dbas[0] = dif_prism3_linear(inverse_D, z1, z2, p, i, dif_XYZ[0]);
        dbas[1] = dif_prism3_linear(inverse_D, z1, z2, p, i, dif_XYZ[1]);
        dbas[2] = dif_prism3_linear(inverse_D, z1, z2, p, i, dif_XYZ[2]);
        // производные u[0] по x,y,z
        dudx[0] += dbas[0] * q1;
        dudy[0] += dbas[1] * q1;
        dudz[0] += dbas[2] * q1;
        // производные u[1] по x,y,z
        dudx[1] += dbas[0] * q2;
        dudy[1] += dbas[1] * q2;
        dudz[1] += dbas[2] * q2;
        // производные u[2] по x,y,z
        dudx[2] += dbas[0] * q3;
        dudy[2] += dbas[1] * q3;
        dudz[2] += dbas[2] * q3;
    }

    (*fe)[k].deps[0] = dudx[0];
    (*fe)[k].deps[1] = dudy[1];
    (*fe)[k].deps[2] = dudz[2];
    (*fe)[k].deps[3] = dudz[1] + dudy[2];
    (*fe)[k].deps[4] = dudx[2] + dudz[0];
    (*fe)[k].deps[5] = dudy[0] + dudx[1];
    solveCD(k);
    // D*eps -> sigma
    MmulV6(D_pl, (*fe)[k].deps, (*fe)[k].dsigma);
}
break;
case GridFEType::Hexagon:
{
}
break;
default:
{

}
break;
}
*/
/*
// напечатать матрицу A
void SlauSolver::print_A(const double *a, const double *d)
{
    int i, j, t;
    double **m = new double*[N];
    for (i = 0; i < N; i++)
        m[i] = new double[N];
    for (j = 0; j < N; j++)
        for (i = 0; i < N; i++)
            m[j][i] = 0;
    // вне диагонали
    for (j = 0; j < N; j++)	// строка номер j
    {
        //printf("ind[j,j+1] = %d, %d\n", ind[j], ind[j+1]);
        for (t = ind[j]; t < ind[j + 1]; t++)	// t - индекс элемента j-й строки
        {
            //printf("i = ai[%d] = %d\n", t, ai[t]);

            i = ai[t];		// столбец номер i
            // a[t] - элемент A(j, i)
            m[j][i] = a[t];
            m[i][j] = a[t];

        }
    }
    // диагональ
    for (i = 0; i < N; i++)
        m[i][i] = d[i];
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < N; i++)
        {
            printf("%.1le	", m[j][i]);
        }printf("%.1le\n", f[j]);
    }
    for (i = 0; i < N; i++)
        delete m[i];
    delete m;
}*/

// Деформация балки
void genTestBeam(TaskTest &/*task*/)
{/*
    Beam_parameters bp;

    // первые краевые условия
    if (task.beamTask.type == BeamTask_x || task.beamTask.type == BeamTask_xyz || task.beamTask.type == BeamTask_x_y)
    {
        bp.bc1_x.u = { { 0., BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
        bp.bc1_y.u = { { BOLSHOE_CHISLO, 0., BOLSHOE_CHISLO } };
        bp.bc1_z.u = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, 0. } };
        bp.bc1_xyz.u = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
    }
    if (task.beamTask.type == BeamTask_around_z)
    {
        bp.bc1_x.u = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
        bp.bc1_y.u = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
        bp.bc1_z.u = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
        bp.bc1_xyz.u = { { 0, 0, 0 } };
    }

    // вторые краевые условия
    double P = task.beamTask.P;
    if (task.beamTask.type == BeamTask_x)
    {
        bp.bc2_x1 = { Solid_BoundaryCondition2_source_vector, { { -P, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_x2 = { Solid_BoundaryCondition2_source_vector, { { +P, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_y1 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_y2 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_z1 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_z2 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
    }
    if (task.beamTask.type == BeamTask_xyz)
    {
        bp.bc2_x1 = { Solid_BoundaryCondition2_source_vector, { { -P, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_x2 = { Solid_BoundaryCondition2_source_vector, { { +P, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_y1 = { Solid_BoundaryCondition2_source_vector, { { 0, -P, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_y2 = { Solid_BoundaryCondition2_source_vector, { { 0, +P, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_z1 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, -P } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_z2 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, +P } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
    }
    if (task.beamTask.type == BeamTask_x_y)
    {
        bp.bc2_x1 = { Solid_BoundaryCondition2_source_vector, { { -P, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_x2 = { Solid_BoundaryCondition2_source_vector, { { +P, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_y1 = { Solid_BoundaryCondition2_source_vector, { { 0, +P, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_y2 = { Solid_BoundaryCondition2_source_vector, { { 0, -P, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_z1 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_z2 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
    }
    if (task.beamTask.type == BeamTask_around_z)
    {
        bp.bc2_x1 = { Solid_BoundaryCondition2_source_vector, { { 0, -P, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };	//tau_yx
        bp.bc2_x2 = { Solid_BoundaryCondition2_source_vector, { { 0, +P, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };	//tau_xy
        bp.bc2_y1 = { Solid_BoundaryCondition2_source_vector, { { -P, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_y2 = { Solid_BoundaryCondition2_source_vector, { { +P, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_z1 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
        bp.bc2_z2 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
    }
    // построение сетки
    //////BuldGridRegularParallelepiped(task.bt.rpp, grid);
    // материал
    (*task.materials)[0] = *(m[task.beamTask.materialNumber].gen());
    // построение краевых условий
    set_beam_sources(bp, task.bc1Source, task.bc2Source);
    // определение точек наблюдения
    //task.
    //points.clear();
    */
}

void setBc2Rotation(std::vector<Solid::MechBoundaryCondition2Source> *&bc2Source, double P1, double P2, double t1, double t2, double a1, double a2)
{
    using namespace Solid;
    bc2Source = new std::vector<MechBoundaryCondition2Source>(7);
    for(int i = 0; i < 7; i++)
    {
       (*bc2Source)[i].mode = Solid_BoundaryCondition2_source_none;
    }
    // линейная функция
    FunParser::Function *fun[3];
    for(int i = 1; i <= 1; i++)
    {
        (*bc2Source)[i].mode = Solid_BoundaryCondition2_source_vector;
        for(int j = 0; j < 3; j++)
            fun[j] = new FunParser::Function;
        char Px[1000];
        char Py[1000];
        if(i == 1)
        {
            double aA, aB, aC;
            double PA, PB, PC;
            solveCoef(t1, t2, a1, a2, aA, aB, aC);
            solveCoef(t1, t2, P1, P2, PA, PB, PC);
            sprintf(Px, "0+((t-%.16le)*%.16le+%.16le)*cos(((t-%.16le)*%.16le+%.16le))\0", abs(PA), PB, PC, abs(aA), aB, aC);
            sprintf(Py, "0+((t-%.16le)*%.16le+%.16le)*sin(((t-%.16le)*%.16le+%.16le))\0", abs(PA), PB, PC, abs(aA), aB, aC);
            //FILE *ttt = fopen("ttt", "w");
            //fprintf(ttt, "%s\n", Px);
            //fprintf(ttt, "%s\n", Py);
            //fclose(ttt);
        }
        (*fun[0]).setExpression(std::string(Px));
        (*fun[1]).setExpression(std::string(Py));
        (*fun[2]).setExpression("0\0");
        // переменная
        for(int j = 0; j < 3; j++)
        {
            (*fun[j]).args.clear();
            (*fun[j]).addArgument("t\0");
            (*fun[j]).parse();
            (*bc2Source)[i].x[j] = fun[j];
        }
    }
}

// Вращение или изгиб балки
void genTestRotation(TaskTest &task)
{
    // индекс теста
    task.testIndex = 4;
#if false
    // общие параметры нагружения и дискретизации
    int StepsNumber1 = 32;    // шаги 128
    int StepsNumber2 = 64;    // шаги 256
    int StepsNumber3 = 32;    // шаги 128
    int sigma0Mode = 1;         // 0 - dSigma = Gq, 1 - dSigma = integral(...)
    int movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
    int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny
    //int HomogenyMode = Solid_FeData_HomogenyMode_not_homogeny;
    double P = 10.e7;//16.e7;

    // геометрия
    GridRegularParallelepiped_parameters rpp;
    rpp.N[0] = 16;
    rpp.N[1] = 4;
    rpp.N[2] = 4;
    rpp.p1 = POINT3(0,-1,-1);
    rpp.p2 = POINT3(8, 1, 1);
    //rpp.p2 = POINT3(rpp.N[0],rpp.N[1]/2,rpp.N[2]/2);
    rpp.bang_x1 = 0;rpp.bang_x2 = 0;rpp.bang_y1 = 0;rpp.bang_y2 = 0;
    // построение сетки
    BuldGridRegularParallelepiped(rpp, task.grid);
    task.grid->buldOpenSCADModel("setka.scad");

    // материал
    task.materials = new std::vector<MechMaterialSource>(1);
    MechMaterialSource &mat = (*task.materials)[0];
    mat.elasticSigmaLimit = 200000000.e7;
    mat.Ceps = 4./9.;
    mat.Csigma = 1;
    mat.set_E_NU(1.e10,0.3);
    mat.set_D_isotropic();
    mat.set_M_sigma_3D();
    mat.F = {0,0,0};
    mat.ro = 0;//1e10;
    // I. Идеально-пластичный материал
    {
        mat.mode = MechMaterialType::PerfectPlasticity;
    }
    // II. Сглаженная кривая
    {
        /*
        mat.mode = MaterialMode::Plasticity;
        mat.usingSecantMethod = true;   // метод секущих
        // строковые представления функций
        double a = 2;  // степень
        int b = 1e4;   // на это число делим для гладкости
        char ttt[1000];
        sprintf(ttt, "(x/(3*G))*(1+sign(s0-x))/2+(1-sign(s0-x))/2*(x+pow(abs(x-s0),%.1lf)/%d)/(3*G)\0", a, b);
        //sprintf(ttt, "(x/(3*G))*(1+sign(s0-x))/2+(1-sign(s0-x))/2*(x+pow(abs(x-s0),2.1)/%d)/(3*G)\0", b);
        printf("%s\n", ttt);
        std::string t(ttt);
        mat.epsFun.setExpression(t);
        char tttt[1000];
        sprintf(tttt, "(1/(3*G))*(1+sign(s0-x))/2+(1-sign(s0-x))/2*(1+pow(abs(x-s0),%.1lf)*%.1lf/%d)/(3*G)\0", a-1., a, b);
        printf("%s\n", tttt);
        std::string tt(tttt);
        mat.difEpsFun.setExpression(tt);
        // переменная и аргументы
        mat.epsFun.args.clear();
        mat.difEpsFun.args.clear();
        mat.epsFun.addArgument("x");
        mat.difEpsFun.addArgument("x");
        mat.epsFun.addArgument("K", mat.K);
        mat.difEpsFun.addArgument("K", mat.K);
        mat.epsFun.addArgument("G", mat.G);
        mat.difEpsFun.addArgument("G", mat.G);
        mat.epsFun.addArgument("s0", mat.elasticSigmaLimit);
        mat.difEpsFun.addArgument("s0", mat.elasticSigmaLimit);
        // разбор
        mat.epsFun.parse();
        mat.difEpsFun.parse();
        mat.epsFun.setName("m");
        FILE *f = fopen("___plastic.txt", "w");
        double p1 = 0;
        double p2 = mat.elasticSigmaLimit + 0.1*1e7; // strToDouble(m0.vals[1]) - предел упругости
        int N = 1000;
        for(int i = 0; i < N; i++)
        {
            double p = p1+(p2-p1)*i/N;
            double e = mat.eps(p, 0);
            fprintf(f, "%le %le\n", e, p);
        }
        fclose(f);
        */
    }

    // первые краевые условия
    task.bc1Source = new std::vector<MechBoundaryCondition1Source>(10);
    for(int i = 0; i < 10; i++)
    {
        (*task.bc1Source)[i].mode = {{-1, -1, -1}};
        (*task.bc1Source)[i].u0 = {{-1, -1, -1}};
    }
    // середина по z
    (*task.bc1Source)[2].mode = {{-1, -1, 0}};
    (*task.bc1Source)[2].u0 = {{-1, -1, 0}};

    // x1
    (*task.bc1Source)[4].mode = {{0, 0, 0}};
    (*task.bc1Source)[4].u0 = {{0, 0, 0}};



    //(*task.bc1Source)[0] = {{0, 0, 0}};   // середина по x

    //(*task.bc1Source)[0] = {{0, BOLSHOE_CHISLO, BOLSHOE_CHISLO}};   // середина по x
    //(*task.bc1Source)[2] = {{0, 0, 0}}; // середина по z
    //(*task.bc1Source)[3] = {{0, 0, 0}}; // центр

    //(*task.bc1Source)[0] = {{0, BOLSHOE_CHISLO, BOLSHOE_CHISLO}};

    // нагружение
    // глобальные шаги
    task.step = new std::vector<MechGlobalStep>;
    std::vector<MechBoundaryCondition2Source> *bc2Source;
    MechGlobalStep s0;
    s0.slausolverParameters = task.slausolver_parameters;
    s0.timeMode = 0;    // 0 - квазистатическая задача
    s0.sigma0Mode = sigma0Mode;  // 0 - dSigma = Gq, 1 - dSigma = integral(...)
    s0.fixGrid = 0;     // 1 - зафиксировать сетку
    s0.movingGridMode = movingGridMode;  // 0 - простые приращения, 1 - приращения Яумана
    s0.controlMode = 0; // 0 - игнорировать не достижение невязок(выда.тся предупреждения в консоль)
    s0.nonlinearIterLimit = 32;  // максимум итераций на шаг
    s0.slauResidualLimit = 1;
    s0.epsResidualLimit = 1.e-15;
    s0.sigmaResidualLimit = 1;
    s0.perfectPlasticitySigmaTResidualLimit = 0.01;
    double Time = 1;
    // первый шаг - нагружение
    s0.t1 = 0;                          // (имеет значение только для первого глобального шага)
    s0.t2 = s0.t1 + Time;
    s0.dt0 = (s0.t2 - s0.t1)/StepsNumber1;
    setBc2Rotation(bc2Source, 0, P, s0.t1, s0.t2, 0, 0);
    s0.bc2Source = bc2Source;
    task.step->push_back(s0);
    // вращение
    for(int i = 0; i < StepsNumber2; i++)
    {
        s0.t1 = s0.t2;                   // (имеет значение только для первого глобального шага)
        s0.t2 = s0.t1 + Time/StepsNumber2;
        s0.dt0 = (s0.t2 - s0.t1)/StepsNumber3;
        double a = PI/2.*(s0.t2-1);
        setBc2Rotation(bc2Source, P, P, s0.t1, s0.t2, a, a);
        s0.bc2Source = bc2Source;
        task.step->push_back(s0);
    }
    /*
    // первый шаг - нагружение
    s0.t1 = 0;                          // (имеет значение только для первого глобального шага)
    s0.t2 = s0.t1 + Time;
    s0.dt0 = (s0.t2 - s0.t1)/StepsNumber1;
    setBc2Rotation(bc2Source, 0, P, s0.t1, s0.t2, 0, 0);
    s0.bc2Source = bc2Source;
    task.step->push_back(s0);
    // второй шаг - вращение
    s0.t1 = s0.t2;
    s0.t2 = s0.t1 + Time;
    s0.dt0 = (s0.t2 - s0.t1)/StepsNumber2;
    setBc2Rotation(bc2Source, P, P, s0.t1, s0.t2, 0, PI/2);
    s0.bc2Source = bc2Source;
    task.step->push_back(s0);
    // третий шаг - закрепление
    s0.t1 = s0.t2;
    s0.t2 = s0.t1 + Time;
    s0.dt0 = (s0.t2 - s0.t1)/StepsNumber1;
    setBc2Rotation(bc2Source, P, P, s0.t1, s0.t2, PI/2, PI/2);
    s0.bc2Source = bc2Source;
    task.step->push_back(s0);
    */
    // инициализация начальных деформаций и напряжений (=0)
    task.fe = new std::vector<MechFeData>;
    task.fe->resize(task.grid->fe.size());
    for (size_t i = 0; i < task.grid->fe.size(); i++)
    {
        (*task.fe)[i].inBasisType_1L_1L,
                           HomogenyMode,
                           Integration::IntegrationType::Gauss3);
    }
    // инициализация начальных перемещений (=0)
    task.vertex = new std::vector<MechVertexData>;
    for (size_t i = 0; i < task.grid->vertex.size(); i++)
    {
        MechVertexData v_el;
        for (int j = 0; j < 3; j++)
        {
            v_el.sum_du[j] = 0;
            v_el.du[j] = 0; // не имеет значения
        }
        task.vertex->push_back(v_el);
    }
    // инициализация начальных скоростей (=0) и ускорений (=0)
    task.V0 = new Vector(task.grid->vertex.size()*3);
    task.dV0 = new Vector(task.grid->vertex.size()*3);
    for (size_t i = 0; i < task.V0->size(); i++)
    {
        (*task.V0)[i] = 0;
        (*task.dV0)[i] = 0;
    }
#endif
}

// Удар по пластине
void genTestBang(TaskTest &task)
{
    // индекс теста
    task.testIndex = 2;
#if false
    // количество снимков
    int N = 10;
    // выбор дробления сетки
    int K_xyz = 0;
    int K_T = 2;
    // начальная сетка
    int Nx_0 = 16;
    int Ny_0 = 16;
    int Nz_0 = 8;
    //int NT1_0 = 40;     //~1
    int NT2_0 = 40;     //~2

    int mnoj_xy = 1;//pow(2,K_xyz);
    int mnoj_z = pow(2,K_xyz);
    int mnoj_T = pow(4,K_T);
    // дробленая сетка
    int Nx = Nx_0;//*mnoj_xyz;
    int Ny = Ny_0;//*mnoj_xyz;
    int Nz = Nz_0*mnoj_z;
    //int NT1 = NT1_0*mnoj_T;
    int NT2 = NT2_0*mnoj_T;

    double t1 = 0;
    double t2 = 0;//0.0002;//0.00001;
    double t3 = t2 + 0.0002;
    //double dt_1 = (t2-t1)/NT1;
    double dt_2 = (t3-t2)/NT2;
    double t2_inc = (t3-t2)/N;
    //double P_bang = 1e10;
    double ro = 3000;

    task.bang.K_T = K_T;
    task.bang.K_xyz = K_xyz;
    task.bang.mnoj_xy = mnoj_xy;
    task.bang.mnoj_z = mnoj_z;
    task.bang.mnoj_T = mnoj_T;
    task.bang.Nx = Nx;
    task.bang.Ny = Ny;
    task.bang.Nz = Nz;


    task.beamTask.rpp.N[0] = Nx;
    task.beamTask.rpp.N[1] = Ny;
    task.beamTask.rpp.N[2] = Nz;
    task.beamTask.rpp.p1 = {-1, -1, -0.1};
    task.beamTask.rpp.p2 = {+1, +1, +0.1};
    task.beamTask.rpp.bang_x1 = Nx/2 - mnoj_xy;
    task.beamTask.rpp.bang_x2 = Nx/2 + mnoj_xy;
    task.beamTask.rpp.bang_y1 = Ny/2 - mnoj_xy;
    task.beamTask.rpp.bang_y2 = Ny/2 + mnoj_xy;
    //task.beamTask.rpp.bang_x1 = 0;
    //task.beamTask.rpp.bang_x2 = Nx;
    //task.beamTask.rpp.bang_y1 = 0;
    //task.beamTask.rpp.bang_y2 = Ny;


    // материал
    MechMaterialSource &mat = (*task.materials)[0];
        mat.elasticSigmaLimit = 2.e7;
    mat.mode = MechMaterialType::Elasticity;
    mat.set_E_NU(1.e10,0.4);
    mat.set_D_isotropic();
    mat.F = {0,0,0};
    mat.ro = ro;

    // сетка
    Beam_parameters bp;
    // первые краевые условия
    bp.bc1_x.u0 = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
    bp.bc1_y.u0 = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
    bp.bc1_z.u0 = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
    bp.bc1_xyz.u0 = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
    bp.bc1_x1.u0 = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
    bp.bc1_x2.u0 = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
    bp.bc1_y1.u0 = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
    bp.bc1_y2.u0 = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
    bp.bc1_z1.u0 = { { BOLSHOE_CHISLO, BOLSHOE_CHISLO, BOLSHOE_CHISLO } };
    bp.bc1_z2.u0 = { { 0, 0, 0 } };

    // вторые краевые условия
    bp.bc2_x1 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
    bp.bc2_x2 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
    bp.bc2_y1 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
    bp.bc2_y2 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
    bp.bc2_z1 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
    bp.bc2_z2 = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
    bp.bc2_bang = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
    //bp.bc2_bang = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, P_bang } } };
    //bp.bc2_bang = { Solid_BoundaryCondition2_source_scalar, { { 0, 0, 0 } }, { { P_bang, 0, 0 } } };


    // построение краевых условий
    std::vector<MechBoundaryCondition2Source> *bc2Source;
    set_beam_sources(bp, task.bc1Source, bc2Source);
    // построение сетки
    BuldGridRegularParallelepiped(task.beamTask.rpp, task.grid);
    //task.grid->buldOpenSCADModel("setka.scad");

    // шаги
    // первый шаг - удар по плите
    task.step = new std::vector<MechGlobalStep>;
    MechGlobalStep s0;
    // первый шаг - удар
    s0.t1 = t1;
    s0.t2 = t2;
    //s0.dt0 = dt_1;
    s0.slausolverParameters = task.slausolver_parameters;
    s0.timeMode = 3;    // 0 - квазистатическая задача
    s0.fixGrid = 1;     // 1 - зафиксировать сетку
    s0.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
    s0.controlMode = 0; // 0 - игнорировать не достижение невязок(выда.тся предупреждения в консоль)
    s0.nonlinearIterLimit = 128;
    s0.slauResidualLimit = 1;
    s0.epsResidualLimit = 1.e-14;
    s0.sigmaResidualLimit = 1;
    s0.perfectPlasticitySigmaTResidualLimit = 0.01;
    s0.bc2Source = bc2Source;
    //task.step->push_back(s0);
    // второй шаг - наблюдение за последствиями 1

    bp.bc2_bang = { Solid_BoundaryCondition2_source_vector, { { 0, 0, 0 } }, { { 0, 0, 0 } } };
    set_beam_sources(bp, task.bc1Source, bc2Source);
    //s0.t1 = s0.t2;
    s0.t1 = 0;
    s0.t2 += t2_inc;
    s0.dt0 = dt_2;
    s0.bc2Source = bc2Source;
    task.step->push_back(s0);

    for(int i = 1; i < N; i++)
    {
        // наблюдение за последствиями номер i
        s0.t1 = s0.t2;
        s0.t2 += t2_inc;
        s0.dt0 = dt_2;
        task.step->push_back(s0);
    }

    // инициализация начальных деформаций и напряжений (=0)
    task.fe = new std::vector<MechFeData>;
    task.fe->resize(task.grid->fe.size());
    for (size_t i = 0; i < task.grid->fe.size(); i++)
    {
    BasisType_1Llid_FeData_BasisType_1L,
                           Solid_FeData_HomogenyMode_homogeny,
                           Integration::IntegrationType::Gauss3);
    }
    task.vertex = new std::vector<MechVertexData>;
    // инициализация перемещений (=0)
    for (size_t i = 0; i < task.grid->vertex.size(); i++)
    {
        MechVertexData v_el;
        for (int j = 0; j < 3; j++)
        {
            v_el.sum_du.x[j] = 0;
            v_el.du[j] = 0; // не имеет значения
        }
        task.vertex->push_back(v_el);
    }

    // инициализация начальных скоростей и ускорений
    task.V0 = new Vector(task.grid->vertex.size()*3);
    task.dV0 = new Vector(task.grid->vertex.size()*3);
    for (size_t i = 0; i < task.V0->size(); i++)
    {
        (*task.V0)[i] = 0;
        (*task.dV0)[i] = 0;
    }

    // задать начальные скорости
    // дробление сетки все портит

    for (size_t i = 0; i < task.grid->bc2.size(); i++)
    {
        BoundaryCondition2_4gonal &bc2_el = task.grid->bc2[i];
        if(bc2_el.si == 6)  //z1
        {
            for(int k = 0; k < 4; k++)
            {
                int vInd = bc2_el.vi[k];
                //(*task.dV0)[vInd*3 + 2] = 10; // начальная скорость вместо краевого условия
                (*task.V0)[vInd*3 + 2] = 10; // начальная скорость вместо краевого условия

            }
        }
    }
#endif
}

void set_beam_sources(const Beam_parameters &bp,
    std::vector<Solid::MechBoundaryCondition1Source> *&bc1_sources,
    std::vector<Solid::MechBoundaryCondition2Source> *&bc2_sources)
{
    using namespace Solid;
    // первые краевые условия
    bc1_sources = new std::vector<MechBoundaryCondition1Source>(10);
    (*bc1_sources)[0] = bp.bc1_x;
    (*bc1_sources)[1] = bp.bc1_y;
    (*bc1_sources)[2] = bp.bc1_z;
    (*bc1_sources)[3] = bp.bc1_xyz;
    (*bc1_sources)[4] = bp.bc1_x1;
    (*bc1_sources)[5] = bp.bc1_x2;
    (*bc1_sources)[6] = bp.bc1_y1;
    (*bc1_sources)[7] = bp.bc1_y2;
    (*bc1_sources)[8] = bp.bc1_z1;
    (*bc1_sources)[9] = bp.bc1_z2;

    // вторые краевые условия: векторы P
    bc2_sources = new std::vector<MechBoundaryCondition2Source>(7);
    (*bc2_sources)[0] = bp.bc2_x1;
    (*bc2_sources)[1] = bp.bc2_x2;
    (*bc2_sources)[2] = bp.bc2_y1;
    (*bc2_sources)[3] = bp.bc2_y2;
    (*bc2_sources)[4] = bp.bc2_z1;
    (*bc2_sources)[5] = bp.bc2_z2;
    (*bc2_sources)[6] = bp.bc2_bang;
}

// учет вторых краевых условий (второе слагаемое в правой части)
// (старый вариант)
/*void addBc2(const Grid::Grid3D *grid, const std::vector<MechBoundaryCondition2Source> *bc2Source, Vector &b)
{
    using namespace Integration;
    LinearHexagon hexagon;
    POINT3 dxyz[2];
    POINT3_CUBE cubePoint;
    VECTOR3 P[Integration::Integration2D_Size_Max];
    double norm_P;
    VECTOR3 normal;
    Integrator integrationFoursquare;
    integrationFoursquare.init2D(IntegrationType::Gauss3);//##
    Integrator integrationCube;
    integrationCube.init3D(IntegrationType::Gauss3);
    for (size_t bc2Ind = 0; bc2Ind < grid->bc2.size(); bc2Ind++)
    {
        int si = grid->bc2[bc2Ind].si;
        if(si != -1)
        {
            const MechBoundaryCondition2Source &bc2_el = (*bc2Source)[si];
            if(bc2_el.mode != Solid_BoundaryCondition2_source_none)
            {
                // копируем заданные 4 вершины
                for (int t = 0; t < 4; t++)	// t - локальный номер вершины
                    hexagon[t] = grid->vertex[grid->bc2[bc2Ind].vi[t]];
                // нагрузка задана скаляром и направлена по нормали к грани
                if (bc2_el.mode == Solid_BoundaryCondition2_source_scalar)
                {
                    VECTOR3 normal1, normal2, normal3, normal4;
                    // вычисляем нормаль к грани
                    Elementary::Operations::solveNormal(hexagon[0], hexagon[1], hexagon[2], normal1);
                    Elementary::Operations::solveNormal(hexagon[2], hexagon[1], hexagon[3], normal2);
                    Elementary::Operations::solveNormal(hexagon[0], hexagon[1], hexagon[3], normal3);
                    Elementary::Operations::solveNormal(hexagon[0], hexagon[3], hexagon[2], normal4);
                    normal = (normal1 + normal2 + normal3 + normal4)/4;
                    normal = normal/normal.abs();
                    // в этой модификации задан скаляр (в P[0])
                    norm_P = (*bc2Source)[grid->bc2[bc2Ind].si].P[0];
                }
                // нагрузка задана скаляром и направлена к заданной точке
                if (bc2_el.mode == Solid_BoundaryCondition2_source_scalar_to_point)
                {
                    POINT3 c = (hexagon[0] + hexagon[1] + hexagon[2] + hexagon[3])/4;
                    normal = c - bc2_el.c0;
                    normal = normal/normal.abs();
                    // в этой модификации задан скаляр (в P[0])
                    norm_P = (*bc2Source)[grid->bc2[bc2Ind].si].P[0];
                }
                // достраиваем грань до 6-гранника
                for (int t = 4; t < 8; t++)	// t - локальный номер вершины
                    hexagon[t] = hexagon[t - 4] - normal*(hexagon[0]-hexagon[1]).abs();




        //        // вычисление нормали с помощью отображения с трилинейными базисными функциями
        //        POINT3 c, p1, p2;
        //       for(int j = 0; j < 3; j++)
        //        {
        //            for (int i = 0; i < 8; i++)
        //            {
        //                c[j] += Fem::scalarStandartLinear3D(POINT3(-1,0,0), i, Fem::dif_NULL) * hexagon[i][j];
        //            }
        //        }




                // задаем на грани вектор нагрузки P
                for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
                {
                    if (bc2_el.mode == Solid_BoundaryCondition2_source_vector)
                    {	// просто копируем вектор нагрузки P
                        P[valIndex] = (*bc2Source)[grid->bc2[bc2Ind].si].P;
                    }
                    if (bc2_el.mode == Solid_BoundaryCondition2_source_scalar)
                    {
                        P[valIndex] = (normal * norm_P);
                    }
                    if (bc2_el.mode == Solid_BoundaryCondition2_source_scalar_to_point)
                    {
                        P[valIndex] = (normal * norm_P);
                    }
                }
                for (int i = 0; i < 4; i++)		// i - локальный номер базисной функции (вершины)
                    for (int N = 0; N < 3; N++)
                    {
                        for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
                        {
                            cubePoint[0] = integrationFoursquare.p2[valIndex].x[0];
                            cubePoint[1] = integrationFoursquare.p2[valIndex].x[1];
                            cubePoint[2] = -1;

//                            // к центру###
//                            POINT3 xyz;
//                            Fem::cubeToHexagonLinear(hexagon, cubePoint, Fem::dif_NULL, xyz);
//                            normal = xyz - bc2_el.c0;
//                            normal = normal*norm_P/normal.abs();
//                            P[valIndex][N] = normal[N];

                            Fem::cubeToLinearHexagon(hexagon, cubePoint, Fem::dif_XYZ[0], dxyz[0]);
                            Fem::cubeToLinearHexagon(hexagon, cubePoint, Fem::dif_XYZ[1], dxyz[1]);
                            integrationFoursquare.value[valIndex] =
                                    P[valIndex][N] *
                                    Fem::linear3D(cubePoint, i, Fem::dif_NULL) *
                                    sqrt(
                                    (SQR(dxyz[0][0]) + SQR(dxyz[0][1]) + SQR(dxyz[0][2]))*
                                    (SQR(dxyz[1][0]) + SQR(dxyz[1][1]) + SQR(dxyz[1][2])) -
                                    SQR(dxyz[0][0] * dxyz[1][0] + dxyz[0][1] * dxyz[1][1] + dxyz[0][2] * dxyz[1][2]));
                        }
                        b[grid->bc2[bc2Ind].vi[i] * 3 + N] += integrationFoursquare.integrate();
                    }
            }
        }
    }
}*/

/*
            POINT3 p0 = {};
            VECTOR6 sumEps0 = {};		// деформации
            VECTOR6 sumSigma0 = {};	    // напряжения

            // считаем средние значения по шестиграннику
            for (size_t pInd = 0; pInd < 27; pInd++)	// индекс точки Гаусса внутри конечного элемента
            {
                const FePointData &fePD = (*fe)[feIndex].pd[pInd];
                sumEps0 += fePD.sumEps;
                sumSigma0 += fePD.sumSigma;
                POINT3 p(0,0,0);
                Fem::cubeToHexagonLinear(hexagon, integrationCube.p3[pInd], Fem::dif_NULL, p);
                p0 += p;
            }
            p0 /= 27;
            sumEps0 /= 27;
            sumSigma0 /= 27;
            // заполняем
            for (size_t pInd = 0; pInd < 27; pInd++)	// индекс точки Гаусса внутри конечного элемента
            {
                OutFePointData &outFePD = out->back().fe[feIndex].pd[pInd];  // данные вывода для точки конечного элемента номер i
                const FePointData &fePD = (*fe)[feIndex].pd[pInd];
                outFePD.sumEps = sumEps0;
                outFePD.sumSigma = sumSigma0;
                    outFePD.sigmaT = fePD.sigmaT;
                outFePD.p = p0;
            }*/
            //outFePD.dp = c;
            //out0.pres[size].p = g->v[g->v.getCount() - 1];        // положение крайнего узла//##
            //out0.pres[size].u = (*vertex)[g->v.getCount() - 1].u0;    // суммерные перемещения крайнего узла//##
            //out.res[size].p = move(out.p0)... // перемещение точки измерения вместе с телом
            //out.res[size].u = интерполяция...
            // перемещения = ?

// задача с ударом
void SolvingThread::outRes2()
{
#if false
    //FILE *f = fopen("___points.txt", "w");
    using namespace Interpolation;
    GridRectangleRegular3D gr;
    gr.init(task.beamTask.rpp.p1, task.beamTask.rpp.p2,
            task.beamTask.rpp.N[0]/2,
            task.beamTask.rpp.N[1]/2,
            task.beamTask.rpp.N[2]/2);
    Interpolator ip;
    for (int nnn = 0; nnn < (int)out->size() - 1; nnn++)
    {
        MechOutGlobalStepData &o = (*out)[nnn];
        // сохранение скоростей и суммарных перемещений
        if(nnn == (int)out->size() - 2)
        {
            char fn[1000];
            sprintf(fn, "%d%d_u.datu", task.bang.K_xyz, task.bang.K_T);
            FILE *fu = fopen(fn, "w");
            sprintf(fn, "%d%d_q.datq", task.bang.K_xyz, task.bang.K_T);
            FILE *fq = fopen(fn, "w");
            for(int z = 0; z <= task.bang.Nz; z += task.bang.mnoj_z)
                for(int y = 0; y <= task.bang.Ny; y += task.bang.mnoj_xy)
                    for(int x = 0; x <= task.bang.Nx; x += task.bang.mnoj_xy)
                    {
                        int k = z*(task.bang.Ny+1)*(task.bang.Nx+1) + y*(task.bang.Nx+1) + x;
                        MechOutVertexData &vertex_el = o.vertex[k];
                        //fprintf(fq, "%le\t%le\t%le\n", vertex_el.q[0], vertex_el.q[1], vertex_el.q[2]);
                        fprintf(fu, "%le\t%le\t%le\n", vertex_el.sum_du[0], vertex_el.sum_du[1], vertex_el.sum_du[2]);
                    }
            //for(size_t k = 0; k < o.vertex.size(); k++);
            fclose(fu);
            fclose(fq);
        }


        //ip.init(&gr, InterpolatorType::Lagrange3D1, 0.0001, 0, o.fe.size());
        ip.init(&gr, InterpolatorType::Hermite2D3, 0.0001, 0, o.fe.size()*27);
        for(size_t k = 0; k < o.fe.size(); k++)
        {
            for (size_t pInd = 0; pInd < o.fe[k].pd.size(); pInd++)	// индекс точки Гаусса внутри конечного элемента
            {
                MechOutFePointData &fe_el = o.fe[k].pd[pInd];
                POINT3 p = fe_el.p;
                double max;
                // главные напряжения
                VECTOR3 ms;
                fe_el.sumSigma.solveMain(ms);
                max = abs(ms[0]);
                if(abs(ms[1]) > max) max = abs(ms[1]);
                if(abs(ms[2]) > max) max = abs(ms[2]);
                if(!(max == max))
                    max = 0;
                //max = fe_el.sigma[2];
                //max = sin(fe_el.p[0]);
                /*double x1 = task.beamTask.rpp.p1[0];
                double y1 = task.beamTask.rpp.p1[1];
                double lx = task.beamTask.rpp.p2[0] - task.beamTask.rpp.p1[0];
                double ly = task.beamTask.rpp.p2[1] - task.beamTask.rpp.p1[1];
                int nx = task.beamTask.rpp.N[0];
                int ny = task.beamTask.rpp.N[1];

                if(p[0] >= x1 + lx*task.beamTask.rpp.bang_x1/nx &&
                   p[0] <= x1 + lx*task.beamTask.rpp.bang_x2/nx &&
                   p[1] >= y1 + ly*task.beamTask.rpp.bang_y1/ny &&
                   p[1] <= y1 + ly*task.beamTask.rpp.bang_y2/ny)
                {
                    ip.add_point(p, 0, UNKNOWN_FE_INDEX);
                }
                else
                {
                    ip.add_point(p, max, UNKNOWN_FE_INDEX);
                }*/
                //max = log(fe_el.dp.abs());
                if(abs(p[2]) <= 0.10001/task.beamTask.rpp.N[2]&&p[2]>=0)
                {
                    ip.addPoint(p, log(max + 1), UNKNOWN_FE_INDEX);
                    //fprintf(f, "max(%le, %le, %le) = %le\n", p[0], p[1], p[2], max);
                }
            }
        }
        ip.makeInterpolant();
        char fn[1000];
        sprintf(fn, "out%d_.bmp", nnn);
        std::vector<Circle> circle(0);
        // ч/б, с градациями
        ip.saveInterpolant(task.beamTask.rpp.p1[2]*0.0, fn, 320, 320, 1, 0, circle);
        sprintf(fn, "outc%d_.bmp", nnn);
        // цветной, с градациями
        ip.saveInterpolant(task.beamTask.rpp.p1[2]*0.0, fn, 320, 320, 3, 0, circle);
        //OutStepData &s = o.step[nn];
        ip.release();
        print("\nglobalStepNumber = %d outed\n", ARGS(nnn));
    }
    //fclose(f);
#endif
}


// вращение балки
void SolvingThread::outRes4()
{
    task.mechTask.grid->buldOpenSCADModel("setka_.scad");
}

/*{
    for (size_t bc2SourceIndex = 0; bc2SourceIndex < 1; bc2SourceIndex++)
    {
        PRINT("\nfun = %s\n", ARGS((*bc2Source)[bc2SourceIndex].x[0]->getExpression().c_str()));
    }
}*/

/*
for(int i = 0; i < Gportrait.matrixSize; i++)
{
    if(bc1_state[i])
        SlauSolving::SSCMaddBoundaryCondition1(Gportrait, Gelements, b, i, bc1_u0[i]);
}
*/


// построение матриц G, M и вектора b
void buildGMb(const int nonlinearIter, const Grid::Grid3D *grid, const std::vector<MechFeData> *fe, std::vector<MechMaterialSource> *material,
             SlauSolving::SSCMBulder &Gbulder,
             SlauSolving::SSCMBulder &/*Mbulder*/,
             Vector &b, const int sigma0Mode)
{
    using namespace Integration;
    MATR6x6 D_pl;			// матрица D с волной
    LinearHexagon linearHexagon;        // вершины шестигранника
    QuadraticHexagon quadraticHexagon;  // вершины шестигранника с криволинейными границами
    double w[Integration3D_Size_Gauss3];
    double dbas[8][Integration3D_Size_Gauss3][3];       // значения производных базисных функций(8) в каждой точке интегрирования 6-гранника
    double absDet[Integration3D_Size_Gauss3];			// детерминанты отображений шестигранников в шаблонный куб
    MATR3x3x3x3 C_pl0;		// матрица D с волной, записанная в виде тензора
    MATR3x3x3x3 C_pl[Integration3D_Size_Gauss3];		// матрица D с волной, записанная в виде тензора

    std::vector<double> integrationwSource[3];
    std::vector<double> basCubeSource[3];
    std::vector<VECTOR3> dLinearBasCubeSource[3];
    std::vector<VECTOR3> dQuadraticBasCubeSource[3];

    int vi[8];              // глобальные индексы вершин шестигранника
    double intForG_mnjl[8][8][3][3];
    double intForM_mn[8][8];
    double intForb1_m[8];
    double intForb23_mj[8][3];
    int map3x3to6[3][3] =
    {{0,5,4},
     {5,1,3},
     {4,3,2},};
    //double Gloc[24][24];

    // заполнение таблиц значений на шаблонных кубах
    for(int it = (int)IntegrationType::Gauss2; it <= (int)IntegrationType::Gauss3; it++)
    {
        Fem::solveCubeLinearBasisFuncValues((IntegrationType)it,
                                            integrationwSource[it],basCubeSource[it],dLinearBasCubeSource[it]);
        Fem::solveCubeQuadraticBasisFuncValues((IntegrationType)it,
                                               dQuadraticBasCubeSource[it]);
    }

    // обнуление вектора правой части
    for (size_t i = 0; i < b.size(); i++)
        b[i] = 0;

    // шестигранники Solid_FeData_BasisType_1L
    for (size_t feInd = 0; feInd < grid->fe.size(); feInd++)	// индекс конечного элемента
    {
        const MechFeData &feDataEl = (*fe)[feInd];
        MechMaterialSource &m0 = (*material)[grid->fe[feInd].mi];
        VECTOR3 F = m0.F;     // объемная сила
        double ro = m0.ro;    // плотность
        //if(feDataEl.basisType == Solid_FeData_BasisType_1L)
        int it = (int)feDataEl.integrationType;
        std::vector<double> &integrationw = integrationwSource[it];
        std::vector<double> &basCube = basCubeSource[it];
        std::vector<VECTOR3> &dLinearBasCube = dLinearBasCubeSource[it];
        //if(grid->fe[feInd].type == Grid::FEType::QuadraticHexagon)
        std::vector<VECTOR3> &dQuadraticBasCube = dQuadraticBasCubeSource[it];
        int numPoints = (int)integrationw.size();
        if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_homogeny)
            m0.solveCD(nonlinearIter, feDataEl.pd[0], D_pl, C_pl0);
        if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_not_homogeny)
            for (int pInd = 0; pInd < numPoints; pInd++)	// индекс точки Гаусса внутри конечного элемента
                m0.solveCD(nonlinearIter, feDataEl.pd[pInd], D_pl, C_pl[pInd]);
        // копии глобальных индексов и координат вершин 6-гранника
        grid->getFeVertexIndexes((int)feInd, vi);
        if(grid->fe[feInd].type == Grid::FEType::LinearHexagon)
            grid->getFeVertexes((int)feInd, &linearHexagon[0]);
        if(grid->fe[feInd].type == Grid::FEType::QuadraticHexagon)
            grid->getFeVertexes((int)feInd, &quadraticHexagon[0]);
        // расчёт коэффициентов для интегрирования, модулей детерминантов
        // и производных базисных функций на шестиграннике
        if(grid->fe[feInd].type == Grid::FEType::LinearHexagon)
            Fem::solveLinearHexagonLinearBasisFuncValues(linearHexagon, numPoints, integrationw, dLinearBasCube,
                                 w, absDet, dbas);
        if(grid->fe[feInd].type == Grid::FEType::QuadraticHexagon)
            Fem::solveQuadraticHexagonLinearBasisFuncValues(quadraticHexagon, numPoints, integrationw, dLinearBasCube, dQuadraticBasCube,
                                 w, absDet, dbas);
        // внутренние интегралы для добавок к матрицам и к правой части
        // для матрицы G
        if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_homogeny)
            Fem::solveIntForG(numPoints, w, dbas, intForG_mnjl);
        // для матрицы M
        if(ro != 0)
            Fem::solveIntForM(numPoints, w, basCube, intForM_mn);
        // для 1-го слагаемого правой части
        Fem::solveIntForb_f(numPoints, w, basCube, intForb1_m);
        // для 2-го и 3-го слагаемого правой части
        Fem::solveIntForb_df(numPoints, w, dbas, intForb23_mj);
        // I. добавка к матрице G от конечного элемента feInd
        if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_not_homogeny)
        {
            for (int m = 0; m < 8; m++)
                for (int n = 0; n < 8; n++)	// i, j - локальные номера базисных функций (вершин)
                    for (int i = 0; i < 3; i++)
                        for (int k = 0; k < 3; k++)
                        {
                            int ind_mi = vi[m] * 3 + i;
                            int ind_nk = vi[n] * 3 + k;
                            if (ind_mi >= ind_nk)
                            {
                                double E = 0;
                                for (int valIndex = 0; valIndex < numPoints; valIndex++)
                                {
                                    double EE = 0;
                                    for (int j = 0; j < 3; j++)
                                    {
                                        double EEE = 0;
                                        for (int l = 0; l < 3; l++)
                                            EEE += C_pl[valIndex].m[i][j][k][l] * dbas[n][valIndex][l];
                                        EE += dbas[m][valIndex][j] * EEE;
                                    }
                                    E += w[valIndex] * EE;
                                }
                                Gbulder.addElement_not_null(
                                       E,
                                       ind_mi,
                                       ind_nk);
                            }
                        }
        }
        if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_homogeny)
        {
            for (int m = 0; m < 8; m++)
            {
                int global_mx3 = vi[m] * 3;
                for (int n = 0; n < 8; n++)	// m, n - локальные номера базисных функций (вершин)
                {
                    int global_nx3 = vi[n] * 3;
                    for (int i = 0; i < 3; i++)
                        for (int k = 0; k < 3; k++)
                            if (global_mx3 + i >= global_nx3 + k)
                            {
                                double E = 0;
                                for (int j = 0; j < 3; j++)
                                    for (int l = 0; l < 3; l++)
                                        E += C_pl0.m[i][j][k][l] * intForG_mnjl[m][n][j][l];
                                Gbulder.addElement_not_null(
                                    E,
                                    global_mx3 + i,
                                    global_nx3 + k);
                            }
                }
            }
        }
        // II. добавка к матрице M от конечного элемента feInd
        /*
        if(ro != 0)
        {
            for (int m = 0; m < 8; m++)
            {
                int global_mx3 = vi[m] * 3;
                for (int n = 0; n < 8; n++)	// m, n - локальные номера базисных функций (вершин)
                {
                    int global_nx3 = vi[n] * 3;
                    for (int i = 0; i < 3; i++)
                        for (int k = 0; k < 3; k++)//можно обойтись одним циклом учитывая что i == k
                            if(i == k && global_mx3 + i >= global_nx3 + k)
                            {
                                Mbulder.addElement(
                                    ro*intForM_mn[m][n],
                                    global_mx3 + i,
                                    global_nx3 + k);
                            }
                }
            }
        }
        */
        // III.1 добавка к вектору (объемные силы: ro*Fi*...)
        /*
        for (int m = 0; m < 8; m++)	// i - локальный номер базисной функции (вершины)
            for (int i = 0; i < 3; i++)
                b[vi[m] * 3 + i] += F[i] * intForb1_m[m];
                */
        // III.2 добавка к вектору (накопленные напряжения: -sumj(sigma0ij*...))
        if(sigma0Mode == 1)
        {
            if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_homogeny)
            {

                for (int m = 0; m < 8; m++)	// i - локальный номер базисной функции (вершины)
                    for (int i = 0; i < 3; i++)
                    {
                        double E = 0;
                        for (int j = 0; j < 3; j++)
                            E += feDataEl.pd[0].sumSigma[map3x3to6[i][j]] * intForb23_mj[m][j];
                        b[vi[m] * 3 + i] -= E;
                    }
                /*
                for (int m = 0; m < 8; m++)	// i - локальный номер базисной функции (вершины)
                    for (int i = 0; i < 3; i++)
                    {
                        double E = 0;
                        for (int valIndex = 0; valIndex < numPoints; valIndex++)
                        {
                            double EE = 0;
                            for (int j = 0; j < 3; j++)
                                EE += feDataEl.pd[0].sumSigma[map3x3to6[i][j]] * dbas[m][valIndex][j];
                            E += w[valIndex] * EE;
                        }
                        b[vi[m] * 3 + i] -= E;
                    }*/
            }
            if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_not_homogeny)
            {
                for (int m = 0; m < 8; m++)	// i - локальный номер базисной функции (вершины)
                    for (int i = 0; i < 3; i++)
                    {
                        double E = 0;
                        for (int valIndex = 0; valIndex < numPoints; valIndex++)
                        {
                            double EE = 0;
                            for (int j = 0; j < 3; j++)
                                EE += feDataEl.pd[valIndex].sumSigma[map3x3to6[i][j]] * dbas[m][valIndex][j];
                            E += w[valIndex] * EE;
                        }
                        b[vi[m] * 3 + i] -= E;
                    }
            }
        }
        // III.3 добавка к вектору (температурное слагаемое: alpha*deltaT*...)
        if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_homogeny)
        {
            double Talpha = m0.Talpha;
            double deltaT = feDataEl.pd[0].deltaT;
            for (int m = 0; m < 8; m++)	// i - локальный номер базисной функции (вершины)
                for (int i = 0; i < 3; i++)
                {
                    double E = 0;
                    if(m0.is_XY)
                    {
                        // плоская задача
                        for (int j = 0; j < 3; j++)
                            for (int l = 0; l < 2; l++)
                                E += C_pl0.m[i][j][l][l] * intForb23_mj[m][j];
                    }
                    else
                    {
                        for (int j = 0; j < 3; j++)
                            for (int l = 0; l < 3; l++)
                                E += C_pl0.m[i][j][l][l] * intForb23_mj[m][j];
                    }
                    b[vi[m] * 3 + i] += Talpha * deltaT * E;
                }
        }
        if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_not_homogeny)
        {
            //##
        }
    } // feInd
}

// итерации по нелинейности
nlInf.plasticIterNumber = 0;
for (;;)
{
    if(nlInf.plasticIterNumber == step_el.plasticIterLimit - 1)
    {
        // на последней итерации выбираем лучшее решение
    }

    // построение матриц
    if(contactElasticType != 1)
    PRINT("%.2d genGMb..", nlInf.plasticIterNumber);
    if(nlInf.plasticIterNumber == 0 || (!(contactElasticType==1)))
    {
        buildGM_b_T_sigma0_Bc2(step_el);
        if(contactElasticType == 1)
        {
            solveContactCoefficients(bc1_state, Gmatrix,
                                vertexContactForces);
        }
    }
    // сборка суммы слагаемых в вектор правой части
    VectorCopy(b_P_T_const, b);
    if(contactElasticType == 1)
    {
        // добавление сил реакций опоры в правую часть
        addContactForces(vertexContactForces, b, nlInf, step_el);
        // ###занулить всё-таки может надо?
        for (int i = 0; i < Gmatrix.getMatrixSize(); i++)
        {
            if(bc1_state[i])
            {
                b[i] = bc1_u0[i];
            }
        }
    }
    // Построение матрицы G и вектора b из СЛАУ Gq = b
    // для разных схем аппроксимации по времени
    // (будем использовать ту же матрицу G для СЛАУ)
    int timeMode = step_el.timeMode;
    // 0. Без разностной аппроксимации по времени
    // Квазистатическая задача, заведомо отсутствуют инертные слагаемые
    if(timeMode == 0)
    {
    }
    /*if(timeMode == 3 && firstStep)
    {
        timeMode = 2;
    }
    // (1. Неявная 3-слойная схема с постоянным шагом во времени)
    if(timeMode == 1)
    {
        PRINT1("b..");
        C1mulVector1PlusC2mulVector2(-1/tl.dt/tl.dt, q1, 1/tl.dt/tl.dt, q2, qForM);
        SlauSolving::SSCMmulVector(Mmatrix, qForM, Mddq);
        Vector1PlusCmulVector2(b, -1, Mddq, b);
        PRINT1("G+M..");
        SlauSolving::SSCM1addEnclosedM2mulScalar(Gmatrix, Mmatrix, 1./tl.dt/tl.dt);
    }
    // 2. Неявная 3-слойная схема
    if(timeMode == 2)
    {
        double k0 = 2./(tl.t-tl.t2)/(tl.t-tl.t1);
        double k1 = 2./(tl.t1-tl.t2)/(tl.t1-tl.t);
        double k2 = 2./(tl.t2-tl.t1)/(tl.t2-tl.t);
        PRINT1("b..");
        // qForM = k1*dq1+k2*dq2
        C1mulVector1PlusC2mulVector2(k0 + k1, q1, k2, q2, qForM);
        // Mddq = M*qForM
        SlauSolving::SSCMmulVector(Mmatrix, qForM, Mddq);
        // b += -Mddq
        Vector1PlusCmulVector2(b, -1, Mddq, b);
        PRINT1("G+M..");
        // G += M*k0
        SlauSolving::SSCM1addEnclosedM2mulScalar(Gmatrix, Mmatrix, k0);
    }
    // 3. Неявная 4-слойная схема
    if(timeMode == 3)
    {
        double k0 = 2*(3*tl.t-tl.t1-tl.t2-tl.t3)/(tl.t -tl.t3)/(tl.t -tl.t2)/(tl.t -tl.t1);
        double k1 = 2*(2*tl.t-tl.t2-tl.t3)   /(tl.t1-tl.t3)/(tl.t1-tl.t2)/(tl.t1-tl.t );
        double k2 = 2*(2*tl.t-tl.t1-tl.t3)   /(tl.t2-tl.t3)/(tl.t2-tl.t1)/(tl.t2-tl.t );
        double k3 = 2*(2*tl.t-tl.t1-tl.t2)   /(tl.t3-tl.t2)/(tl.t3-tl.t1)/(tl.t3-tl.t );
        PRINT1("b..");
        // qForM = k1*dq1+k2*dq2
        C1mulVector1PlusC2mulVector2(k0 + k1, q1, k2, q2, qForM);
        // qForM = k1*dq1+k2*dq2+k3*dq3
        Vector1PlusCmulVector2(qForM, k3, q3, qForM);
        // Mddq = M*qForM
        SlauSolving::SSCMmulVector(Mmatrix, qForM, Mddq);
        // b += -Mddq
        Vector1PlusCmulVector2(b, -1, Mddq, b);
        PRINT1("G+M..");
        // G += M*k0
        SlauSolving::SSCM1addEnclosedM2mulScalar(Gmatrix, Mmatrix, k0);
    }*/

    // учёт первых краевых условий
    if(contactElasticType != 1)
    PRINT1("bc1..");

    if(nlInf.plasticIterNumber == 0 || (!(contactElasticType==1)))
    {
        SlauSolving::SSCMaddBoundaryCondition1(Gmatrix, b, bc1_u0, bc1_state);
    }
    // решение СЛАУ
    solveSLAU(step_el);

    if(contactElasticType != 1)
    PRINT1("result.. ");
    // расчет результатов итерации
    solveIterResults(dq, grid, material,
                     vertex, vertexForCurvature, fe, nlInf);
    // вывод на экран невязок после итерации
    if(contactElasticType != 1)
    PRINT("slauR = %le sigmaR = %le epsR = %le NCount = %d\n",
          ARGS(nlInf.slauResidualWithLastq, nlInf.maxSigmaResidual, nlInf.maxEpsResidual,
               nlInf.nonlinearStateFENumber));
    if(contactElasticType == 1)
    {
        // результаты итераций для учёта реакции опоры
        //solveContactIterResults(grid, vertex, bc1_state, Gportrait, Gelements_copy,//Gelements,//Gelements_copy,
        //                    vertexContactForces, nlInf, step_el);
        solveContactIterResults(grid, vertex,
                            vertexContactForces, nlInf, step_el);
        // вывод на экран невязок после итерации
        PRINT("iter=%.3d maxh=%le J=%le J_last=%le dF=%le dh=%le maxDeltah=%le maxDeltaP=%le zalezNomber=%d/%d\n",
              ARGS(nlInf.plasticIterNumber, nlInf.maxh, nlInf.ContactJ, nlInf.ContactJ_last, nlInf.F0_F1, nlInf.h0_h1, nlInf.maxDeltah, nlInf.maxDeltaF, nlInf.zalezNomberIncrement, nlInf.zalezNomber));
    }
    // условия прекращения итераций
    if(checkTerminationConditions(step_el) && checkContactTerminationConditions())
        break;
    else
    {
        nlInf.matrixG_elementsChanged = true;
    }
    nlInf.lastMaxEpsResidual = nlInf.maxEpsResidual;
    nlInf.lastSlauResidualWithLastq = nlInf.slauResidualWithLastq;
    nlInf.lastMaxSigmaResidual = nlInf.maxSigmaResidual;
    nlInf.plasticIterNumber++;
    if(nlInf.plasticIterNumber == step_el.plasticIterLimit)
        break;
}


    MATR3x3x3x3 C_pl[Integration3D_Size_Gauss3];		// матрица D с волной, записанная в виде тензора
    if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_not_homogeny)
        for (int pInd = 0; pInd < numPoints; pInd++)	// индекс точки Гаусса внутри конечного элемента
        {
            m0.solveAddingPlasticEps(firstIter, feDataEl.pd[pInd]);
            m0.solveCD(firstIter, feDataEl.pd[pInd], D_pl, C_pl[pInd]);
        }
    // Неоднородный материал (в точках интегрирования Гаусса различные материалы)
    if(feDataEl.homogenyMode == Solid_FeData_HomogenyMode_not_homogeny)
    {
        // I. добавка к матрице G от конечного элемента feInd
        for (int m = 0; m < 8; m++)
        {
            for (int n = 0; n < 8; n++)	// i, j - локальные номера базисных функций (вершин)
                for (int i = 0; i < 3; i++)
                    for (int k = 0; k < 3; k++)
                    {
                        int ind_mi = vi[m] * 3 + i;
                        int ind_nk = vi[n] * 3 + k;
                        if (ind_mi >= ind_nk)
                        {
                            double E = 0;
                            for (int valIndex = 0; valIndex < numPoints; valIndex++)
                            {
                                double EE = 0;
                                for (int j = 0; j < 3; j++)
                                {
                                    double EEE = 0;
                                    for (int l = 0; l < 3; l++)
                                        EEE += C_pl[valIndex].m[i][j][k][l] * dbas[n][valIndex][l];
                                    EE += dbas[m][valIndex][j] * EEE;
                                }
                                E += w[valIndex] * EE;
                            }
                            Gbulder.addElement_not_null(
                                   E,
                                   ind_mi,
                                   ind_nk);
                        }
                    }
        }
        // II. добавка к матрице M от конечного элемента feInd
        // III.1 добавка к вектору (объемные силы: ro*Fi*...)
        // III.2 добавка к вектору (накопленные напряжения: -sumj(sigma0ij*...))
        if(sigma0Mode == 1)
        {
            for (int m = 0; m < 8; m++)	// i - локальный номер базисной функции (вершины)
                for (int i = 0; i < 3; i++)
                {
                    double E = 0;
                    for (int valIndex = 0; valIndex < numPoints; valIndex++)
                    {
                        double EE = 0;
                        for (int j = 0; j < 3; j++)
                            EE += feDataEl.pd[valIndex].sumSigma[map3x3to6[i][j]] * dbas[m][valIndex][j];
                        E += w[valIndex] * EE;
                    }
                    b_P_T_const[vi[m] * 3 + i] -= E;
                }
        }
        // III.3 добавка к вектору (температурное слагаемое: alpha*deltaT*...)
        // III.4 добавка к вектору (добавочные пластические деформации)
    }
    /*
    // III.2 добавка к вектору (накопленные напряжения: -sumj(sigma0ij*...))
    for (int m = 0; m < 8; m++)	// i - локальный номер базисной функции (вершины)
        for (int i = 0; i < 3; i++)
        {
            double E = 0;
            for (int valIndex = 0; valIndex < numPoints; valIndex++)
            {
                double EE = 0;
                for (int j = 0; j < 3; j++)
                    EE += feDataEl.pd[0].sumSigma[map3x3to6[i][j]] * dbas[m][valIndex][j];
                E += w[valIndex] * EE;
            }
            b[vi[m] * 3 + i] -= E;
        }*/

    // //(II. добавка к матрице M от конечного элемента feInd)
    /*
    if(ro != 0)
    {
        for (int m = 0; m < 8; m++)
        {
            int global_mx3 = vi[m] * 3;
            for (int n = 0; n < 8; n++)	// m, n - локальные номера базисных функций (вершин)
            {
                int global_nx3 = vi[n] * 3;
                for (int i = 0; i < 3; i++)
                    for (int k = 0; k < 3; k++)//можно обойтись одним циклом учитывая что i == k
                        if(i == k && global_mx3 + i >= global_nx3 + k)
                        {
                            Mbulder.addElement(
                                ro*intForM_mn[m][n],
                                global_mx3 + i,
                                global_nx3 + k);
                        }
            }
        }
    }
    */
    // //(III.1 добавка к вектору (объемные силы: ro*Fi*...))
    /*
    for (int m = 0; m < 8; m++)	// i - локальный номер базисной функции (вершины)
        for (int i = 0; i < 3; i++)
            b[vi[m] * 3 + i] += F[i] * intForb1_m[m];
            */
    // III.3 добавка к вектору (температурное слагаемое: alpha*deltaT*...)
    for (int m = 0; m < 8; m++)	// i - локальный номер базисной функции (вершины)
    {
        for (int i = 0; i < 3; i++)
        {
            double E = 0;
            if(m0.is_XY)
            {
                // плоская задача
                for (int j = 0; j < 3; j++)
                    for (int l = 0; l < 2; l++)
                        E += C_pl0.m[i][j][l][l] * feDataEl.intForb23_mj[m][j];
            }
            else
            {
                for (int j = 0; j < 3; j++)
                    for (int l = 0; l < 3; l++)
                        E += C_pl0.m[i][j][l][l] * feDataEl.intForb23_mj[m][j];
            }
            b_P_T_const[vi[m] * 3 + i] += Talpha * deltaT * E;
        }
    }

    // занесение локальных матриц и векторов в глобальные
    /*for (size_t feInd = 0; feInd < grid->fe.size(); feInd++)	// индекс конечного элемента
    {
        MechFeData &feDataEl = (*fe)[feInd];
        // матрица
        if(nlInf.matrixG_elementsChanged)
        {
            int G_local_ind = 0;
            for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
            {
                for (int local_nk_ind = 0; local_nk_ind <= local_mi_ind; local_nk_ind++)
                {
                    int global_mi = feDataEl.globalIndex[local_mi_ind];
                    int global_nk = feDataEl.globalIndex[local_nk_ind];
                    //if(global_mi >= global_nk)
                    {
                        Gbulder.addElement(
                            feDataEl.G_local_el(G_local_ind),
                            global_mi,
                            global_nk);
                        G_local_ind++;
                    }
                }
            }
        }
        // занесение локальных векторов в глобальный вектор
        for (int local_mi_ind = 0; local_mi_ind < 24; local_mi_ind++)
            b_P_T_const[feDataEl.globalIndex[local_mi_ind]] += feDataEl.b_constant_local_el(local_mi_ind);
        // матрица и вектор построены, нужно ли что-то обновлять выяснится после расчёта результата итерации
        feDataEl.indexes_changed = false;
        feDataEl.tables_changed = false;
        feDataEl.CD_changed = false;
        feDataEl.b_constant_changed = false;
    }*/









    // сборщик матрицы
    class SSCMBulder
    {
    public:
        SSCMBulder(const int matrixSize = 0);
        ~SSCMBulder();
        // освободить память
        void release();
        // задать размер матрицы
        void setMatrixSize(const int matrixSize);
        // вызвать перед началом сборки
        void start();
        // добавление элемента el в ячейку (i, j) во временный массив, если el не нулевой
        inline void addElement_not_null(const double el, const int i, const int j)
        {
            if(el != 0 && i >= j)   //####
                str[i].push_back(SlauElement{el, j});
        }
        // добавление элемента el в ячейку (i, j) во временный массив
        inline void addElement(const double el, const int i, const int j)
        {
            if(i >= j)   //####
                str[i].push_back(SlauElement{el, j});
        }
        // сборка (суммирование) всех добавленных элементов матрицы в разреженую матрицу
        void complete(SSCM &matrix);
        void completeWithPortrait(SSCM &matrix);
        void fixReservedMemory();
    private:
        // сравнение 2-х элементов типа SlauElement для сортировки
        static int elementsCmp(const void* x1, const void* x2);
        std::vector<std::vector<SlauElement>> str;  // массив строк из элементов
    };
    SSCMBulder::SSCMBulder(const int matrixSize)
    {
        str.resize(matrixSize);
    }
    SSCMBulder::~SSCMBulder()
    {
        release();
    }
    void SSCMBulder::release()
    {
        str.clear();
    }
    void SSCMBulder::setMatrixSize(const int matrixSize)
    {
        str.resize(matrixSize);
    }
    void SSCMBulder::start()
    {
        srand(0);	//##для отладки
        for(size_t i = 0; i < str.size(); i++)
        {
            str[i].resize(0);
        }
    }
    void SSCMBulder::complete(SSCM &matrix)
    {
        SSCMPortrait &p = *matrix.p;
        SSCMElements &e = *matrix.e;
        int matrixSize = (int)str.size();
        int elementsNumber = 0;
        for(int i = 0; i < matrixSize; i++)
        {
            // сортировка строки
            if(str[i].size() != 0)
            {
                qsort(&str[i][0], str[i].size(), sizeof(SlauElement), elementsCmp);
                // подсчет различных внедиагональных ячеек матрицы
                int jLast = -1;
                for (int j = 0; j < str[i].size(); j++)
                {
                    if(str[i][j].j == i) break;
                    if(str[i][j].j != jLast)
                    {
                        elementsNumber++;
                        jLast = str[i][j].j;
                    }
                }
            }
        }
        p.init(elementsNumber, matrixSize);
        e.init(elementsNumber, matrixSize);
        // очистка диагонали
        for (int i = 0; i < matrixSize; i++)
            e.d[i] = 0;
        // теперь можно представить матрицу в разреженном формате
        // a[ind[i]] - первый элемент строки i. (k = ind[i]...(ind[i+1]-1) для строки i)
        // ai[k] - номер столбца элемента a[k]
        // элементы нижнего треугольника (доступ по индексу k)
        // заполнение ai[] и a[]
        int count_a = 0;	// индекс в массивах ai[] и a[]
        p.ind[0] = 0;
        int i = 0;
        int str_ind = 0;
        for (;;)
        {
            // ищем первый элемент перед началом суммирования
            for(;;)
            {
                if(str_ind >= str[i].size())
                {
                    i++;
                    str_ind = 0;
                    p.ind[i] = count_a;
                    if(i >= matrixSize) break;
                }
                else
                    break;
            }
            if(i >= matrixSize) break;
            int key_i = i;
            int key_j = str[i][str_ind].j;    // координаты очередной ячейки
            double E = 0;                     // суммарное значение в этой ячейке
            for (;;)
            {
                E += str[i][str_ind].a;
                str_ind++;
                // суммируем элементы, принадлежащие одной ячейке
                if(str_ind >= str[i].size() || str[i][str_ind].j != key_j) break;
            }
            if(key_i != key_j)
            {
                e.a[count_a] = E;
                p.ai[count_a] = key_j;
                count_a++;
            }
            else
                e.d[key_i] = E;
        }
    }
    void SSCMBulder::completeWithPortrait(SSCM &matrix)
    {
        SSCMPortrait &p = *matrix.p;
        SSCMElements &e = *matrix.e;
        int matrixSize = (int)str.size();
        for(int i = 0; i < matrixSize; i++)
        {
            // сортировка строки
            if(str[i].size() != 0)
                qsort(&str[i][0], str[i].size(), sizeof(SlauElement), elementsCmp);
        }
        e.init(p.elementsNumber, p.matrixSize);
        // очистка диагонали
        for (int i = 0; i < p.matrixSize; i++)
            e.d[i] = 0;
        int count_a = 0;	// индекс в массивах ai[] и a[]
        int i = 0;
        int str_ind = 0;
        for (;;)
        {
            // ищем первый элемент перед началом суммирования
            for(;;)
            {
                if(str_ind >= str[i].size())
                {
                    i++;
                    str_ind = 0;
                    if(i >= matrixSize) break;
                }
                else
                    break;
            }
            if(i >= matrixSize) break;
            int key_i = (int)i;
            int key_j = str[i][str_ind].j;    // координаты очередной ячейки
            double E = 0;                     // суммарное значение в этой ячейке
            for (;;)
            {
                E += str[i][str_ind].a;
                str_ind++;
                // суммируем элементы, принадлежащие одной ячейке
                if(str_ind >= str[i].size() || str[i][str_ind].j != key_j) break;
            }
            if(key_i != key_j)
            {
                e.a[count_a] = E;
                count_a++;
            }
            else
                e.d[key_i] = E;
        }
    }
    void SSCMBulder::fixReservedMemory()
    {
        // временный массив строк
        std::vector<std::vector<SlauElement>> t(str.size());
        for(size_t i = 0; i < str.size(); i++)
        {
            t[i].resize(str[i].size());
            for(size_t j = 0; j < str[i].size(); j++)
                t[i][j] = str[i][j];
        }
        str.clear();
        str.resize(t.size());
        for(size_t i = 0; i < t.size(); i++)
        {
            str[i].resize(t[i].size());
            str[i].reserve(t[i].size());
            for(size_t j = 0; j < t[i].size(); j++)
                str[i][j] = t[i][j];
        }
        t.clear();
        /*
        // резервируем память для следующего сбора матрицы
        for(size_t i = 0; i < str.size(); i++)
            str[i].reserve(str[i].size());
            */
    }
    int SSCMBulder::elementsCmp(const void* x1, const void* x2)
    {
        SlauElement &px1 = *(SlauElement *)x1;
        SlauElement &px2 = *(SlauElement *)x2;
        if (px1.j < px2.j) return -1;
        if (px1.j > px2.j) return 1;
        //if (px1.a > px2.a) return -1;
        //if (px1.a < px2.a) return 1;	// значения тоже сортируются, по убыванию
        double a1 = abs(px1.a);
        double a2 = abs(px2.a);
        if (a1 > a2) return -1;
        if (a1 < a2) return 1;	// значения тоже сортируются, по убыванию модуля
        return 0;
    }













    void solveContactUpdateState(const Grid::Grid3D *grid, const std::vector<MechVertexData> *vertex,
                                 std::vector<MechVertexContactForce> &vcf, MechNonlinearInf &nlInf)
    {
        nlInf.maxh = 0;
        nlInf.zalezNomber = 0;
        nlInf.ContactJ = 0;
        nlInf.F_average = 0;
        int numberOfNeBeznadejnieVertexes = 0;
        // вычисление функционала, h(dF)
        for(int vertexInd = 0; vertexInd < vcf.size(); vertexInd++)
        {
            // вычисление глубины
            POINT3 v = grid->vertex[vertexInd] + (*vertex)[vertexInd].du;   // положение узла сетки после перемещения
            double h0 = R - sqrt(SQR(v[0] - C[0]) + SQR(v[1] - C[1]));
            if(/*false && */(vertexInd/2)%6 == 0)    //###
            {
                // новое значение h(dF)
                vcf[vertexInd].h = h0;
                if(nlInf.plasticIterNumber == 0)
                    vcf[vertexInd].h_start = h0;
                if(vcf[vertexInd].zalez)
                {
                    // подсчёт количества
                    nlInf.zalezNomber++;
                    if(vcf[vertexInd].sumF + vcf[vertexInd].dF == 0 && vcf[vertexInd].h <= 0)
                    {
                        // не залезла и суммарная реакция опоры нулевая - тут сделать ничего нельзя - скорее всего эта вершина отлипнет
                    }
                    else
                    {
                        numberOfNeBeznadejnieVertexes++;
                        // средний прирост реакции опоры
                        nlInf.F_average += abs(vcf[vertexInd].dF1);
                        // функционал
                        nlInf.ContactJ += h0*h0;
                        // невязка
                        if(abs(h0) > nlInf.maxh)
                            nlInf.maxh = abs(h0);
                    }
                }
            }
        }
        // значение функционала
        if(numberOfNeBeznadejnieVertexes == 0)
        {
            nlInf.ContactJ = 0;
            nlInf.F_average = 10;    // ## начальное приближение
        }
        else
        {
            nlInf.ContactJ = sqrt(nlInf.ContactJ/numberOfNeBeznadejnieVertexes);
            nlInf.F_average = nlInf.F_average/numberOfNeBeznadejnieVertexes;
            nlInf.F_average = 10;
            //if(nlInf.F_average < 100)
            //    nlInf.F_average = 100;
        }
        //###
        nlInf.ContactJ = nlInf.maxh;
    }

    void solveContactUpdateRegularization(std::vector<MechVertexContactForce> &vcf, MechNonlinearInf &nlInf)
    {
        // ###
        //nlInf.ContactBetta *= 1;
        nlInf.ContactBetta *= 0.9;
        // вычисляем среднее h
        double average_h = 0;  // среднее |h|
        double average_dh = 0;  // среднее |h-h1|
        double average_vliyanie = 0;  // среднее влияние
        double max_vliyanie = 0;  // среднее влияние
        int numberOfNeBeznadejnieVertexes = 0;
        for(int vertexInd = 0; vertexInd < vcf.size(); vertexInd++)
        {
            if(vcf[vertexInd].zalez)
            {
                //average_h += abs(vcf[vertexInd].h)*abs(vcf[vertexInd].dF_new - vcf[vertexInd].dF1)*vcf[vertexInd].w;
                //average_h += abs(vcf[vertexInd].h)*vcf[vertexInd].w;
                if(vcf[vertexInd].sumF + vcf[vertexInd].dF == 0 && vcf[vertexInd].h_start <= 0 && vcf[vertexInd].h <= 0)
                {
                    // не залезла и суммарная реакция опоры нулевая - тут сделать ничего нельзя
                }
                else
                {
                    numberOfNeBeznadejnieVertexes++;
                    average_h += abs(vcf[vertexInd].h);
                    average_dh += abs(vcf[vertexInd].h - vcf[vertexInd].h1);
                    double vliyanie = abs(vcf[vertexInd].dF_new - vcf[vertexInd].dF1)*vcf[vertexInd].w;
                    average_vliyanie += vliyanie;
                    if(vliyanie > max_vliyanie)
                        max_vliyanie = vliyanie;
                }
            }
        }
        average_h /= numberOfNeBeznadejnieVertexes;
        average_dh /= numberOfNeBeznadejnieVertexes;
        average_vliyanie /= numberOfNeBeznadejnieVertexes;


        nlInf.F0_F1 = 0;
        nlInf.h0_h1 = 0;
        for(int vertexInd = 0; vertexInd < vcf.size(); vertexInd++)
        {
            if(vcf[vertexInd].zalez)
            {
                if(vcf[vertexInd].sumF + vcf[vertexInd].dF == 0 && vcf[vertexInd].h_start <= 0 && vcf[vertexInd].h <= 0)
                {
                    // не залезла и суммарная реакция опоры нулевая - тут сделать ничего нельзя
                }
                else
                {
                    // тест сходимости F1+(F2-F1)w -> F1
                    nlInf.h0_h1 += abs(vcf[vertexInd].h - vcf[vertexInd].h1);
                    nlInf.F0_F1 += abs(vcf[vertexInd].dF - vcf[vertexInd].dF1);

                    //if(abs(vcf[vertexInd].h) >= abs(vcf[vertexInd].h1) && abs(vcf[vertexInd].h) >= average_h*0.99)
                    //if(abs(vcf[vertexInd].h) >= abs(vcf[vertexInd].h1) && abs(vcf[vertexInd].h)*abs(vcf[vertexInd].dF_new - vcf[vertexInd].dF1)*vcf[vertexInd].w >= average_h*0.5)
                    //if(abs(vcf[vertexInd].h) >= abs(vcf[vertexInd].h1))
                    //if(abs(vcf[vertexInd].h) >= abs(vcf[vertexInd].h1) || abs(vcf[vertexInd].h)*abs(vcf[vertexInd].dF_new - vcf[vertexInd].dF1)*vcf[vertexInd].w >= average_h*1)
                    //if(abs(vcf[vertexInd].h) >= abs(vcf[vertexInd].h1) && abs(vcf[vertexInd].h)*vcf[vertexInd].w >= average_h)
                    //if(abs(vcf[vertexInd].h)*vcf[vertexInd].w >= average_h)
                    bool lezet_dalshe = abs(vcf[vertexInd].h) >= abs(vcf[vertexInd].h1);                // стало хуже
                    //bool silno_lezet_dalshe = abs(vcf[vertexInd].h - vcf[vertexInd].h1) >= average_h;   // сильно лезет дальше
                    bool silno_zalez = abs(vcf[vertexInd].h) >= average_h;                              // залез сильнее чем в среднем
                    double vliyanie = abs(vcf[vertexInd].dF_new - vcf[vertexInd].dF1)*vcf[vertexInd].w;
                    bool silno_vliyaet = vliyanie >= average_vliyanie + (max_vliyanie - average_vliyanie)*0.9;  // большой вклад этого узла
                    //if((lezet_dalshe && silno_vliyaet) || silno_zalez)
                    if((lezet_dalshe && silno_zalez) || silno_vliyaet)
                    //if(lezet_dalshe)
                    {
                        // тут плохо
                        //double h1 = vcf[vertexInd].h1;
                        double h = vcf[vertexInd].h;
                        double vectorF = vcf[vertexInd].dF_new-vcf[vertexInd].dF1;
                        double mul_inc = 2;//1.01;//1.01;//1.5
                        double mul_dec = 0.5;//0.2;//1./4;//1./10;
                        double mul_minimum_inc = 1.2;//1.01
                        //double mul_minimum_dec = 0.9;//1.01
                        //double mul_inc = 1.01;//1.01;//1.01;//1.5
                        //double mul_dec = 0.5;//0.2;//1./4;//1./10;
                        //double mul_minimum_inc = 1.01;//1.01
                        //if(lezet_dalshe && (silno_zalez || silno_vliyaet))
                        //if(lezet_dalshe && silno_zalez)
                        //if(lezet_dalshe && !silno_vliyaet)
                        if(lezet_dalshe)
                        {
                            /*
                            if(vectorF*h1 <= 0)
                            {
                                vcf[vertexInd].w = 0;//mul_dec;
                                //vcf[vertexInd].w *= 0;
                            }*/
                            if(vectorF*h > 0)
                            {
                                // слишком слабо надавили
                                // if(vcf[vertexInd].w < 1.e10)
                                 if(!silno_vliyaet)
                                     vcf[vertexInd].w *= mul_inc;
                                 else
                                     vcf[vertexInd].w *= mul_minimum_inc;
                            }
                            if(vectorF*h < 0)
                            {
                                // слишком сильно надавили и перегнули палку
                                vcf[vertexInd].w *= mul_dec;
                                /*
                                if(!silno_vliyaet)
                                    vcf[vertexInd].w *= mul_dec;
                                else
                                    vcf[vertexInd].w *= mul_minimum_dec;*/
                            }
                        }
                        /*
                        if(!lezet_dalshe)//&& !silno_zalez
                        {
                            vcf[vertexInd].w *= 0.99;
                        }*/
                    }
                }
            }
        }
    }


    void solveContactDirection(std::vector<MechVertexContactForce> &vcf, MechNonlinearInf &nlInf)
    {
        // сохраняем значение нового функционала, которое меньше предыдущего
        nlInf.ContactJ_last = nlInf.ContactJ;
        // возвращаем коэффициент регуляризации к начальному значению
        nlInf.ContactBetta = 1;
        // ищем новое направление (следующая итерация)
        for(int vertexInd = 0; vertexInd < vcf.size(); vertexInd++)
        {
            // вершина залезла
            if(vcf[vertexInd].zalez)
            {
                double dF0 = vcf[vertexInd].dF;    // эту силу приложили на текущем шаге
                double h0 = vcf[vertexInd].h;      // такую глубину в итоге получили
                //double mnoj = 0.0001;     // 1%
                //double dF_start = 1;
                //double mnoj = 1;//0.0001;     // 1%
                double dF_start = nlInf.F_average;//1.e7;

                double dF_new;      // искомая новая сила (направление)
                double dF_inc;      // небольшое изменение силы - первое приближение ###
                double dF_inc_minimum = dF_start;//10;
                double deltaF = 0;
                double deltah = 0;
                if(dF0 == 0)
                    dF_inc = dF_start;
                else
                    //dF_inc = abs(dF0*mnoj);
                    dF_inc = dF_start;
                if(vcf[vertexInd].zalezPerviyRaz)
                {
                    // первое приближение
                    dF_new = dF0 + dF_inc;
                }
                else
                {
                    double dF1 = vcf[vertexInd].dF1;
                    double h1 = vcf[vertexInd].h1; // это данные с прошлого раза
                    deltaF = abs(dF1 - dF0);
                    deltah = abs(h1 - h0);
                    dF_new = dF0 - (dF1 - dF0) * h0 / (h1 - h0);   // метод секущих
                    if(abs(h1 - h0) < 1.e-16 || abs(dF1 - dF0) < dF_inc || nlInf.plasticIterNumber == 0)
                    //if(h1 == h0 || nlInf.nonlinearIter == 0)
                    {
                        // секущая не работает или первая итерация => используем первое приближение
                        if(h0 < 0)
                            dF_new = dF0 - dF_inc;
                        if(h0 > 0)
                            dF_new = dF0 + dF_inc;
                        if(h0 == 0)
                            dF_new = dF0 + dF_inc_minimum;
                    }
                    // проверка направления на адекватность
                    if(h0 < 0 && dF_new >= dF0)
                    {
                        // вершина не залезла, но силу мы увеличили - это не адекватно!
                        dF_new = dF0 - dF_inc;
                    }
                    if(h0 > 0 && dF_new <= dF0)
                    {
                        // вершина залезла, но силу мы уменьшили - это не адекватно!
                        dF_new = dF0 + dF_inc;
                    }
                    if(h0 == 0)
                    {
                        // вершина на границе - делаем силу чуть больше
                        dF_new = dF0 + dF_inc_minimum;
                    }
                    /*
                    if(vcf[vertexInd].sumF + dF_new < 0)
                    {
                        // реакция опоры стала < 0, так нельзя!
                        // просто уменьшаем реакцию опоры
                        dF_new = -vcf[vertexInd].sumF;
                    }
                    */
                    /*if(dF0 != 0)
                    {
                        double delta = abs(dF_new / dF0);
                        if(delta > 10)
                        {
                            if(h <= 0)
                                dF_new = -abs(dF0 * 10);
                            if(h > 0)
                                dF_new = +abs(dF0 * 10);
                        }
                        if(delta < 0.1)
                        {
                            if(h <= 0)
                                dF_new = -abs(dF0 * 0.1);
                            if(h > 0)
                                dF_new = +abs(dF0 * 0.1);
                        }
                    }*/
                }
                // возвращаем коэффициент регуляризации к начальному значению
                vcf[vertexInd].w = 1;
                // обновление направления
                vcf[vertexInd].dF_new = dF_new;
                // сохраняем dF0, h0 для следующих итераций
                vcf[vertexInd].dF1 = dF0;
                vcf[vertexInd].h1 = h0;
                // обновление невязок
                if(deltaF > nlInf.maxDeltaF)    //##
                    nlInf.maxDeltaF = deltaF;
                if(deltah > nlInf.maxDeltah)    //##
                    nlInf.maxDeltah = deltah;
                /*if(h > 0)
                {
                    // вершина залезла => всегда можно увеличить реакцию опоры
                    if(abs(h) > nlInf.maxh)
                        nlInf.maxh = abs(h0);
                    if(deltaF > nlInf.maxDeltaF)
                        nlInf.maxDeltaF = deltaF;
                    if(deltah > nlInf.maxDeltah)
                        nlInf.maxDeltah = deltah;
                }
                if(h < 0)
                {
                    // вершина залезала ранее и реакция опоры > 0 => можно уменьшить реакцию опоры
                    if(abs(h) > nlInf.maxh)
                        nlInf.maxh = abs(h0);
                    if(deltaF > nlInf.maxDeltaF)
                        nlInf.maxDeltaF = deltaF;
                    if(deltah > nlInf.maxDeltah)
                        nlInf.maxDeltah = deltah;
                }*/
                // условие исключения вершины из итерационного процесса
                // #####
            }
        }
    }

    SlauSolving::Vector du, db;
    void solveContactDirection_Au0(const Grid::Grid3D *grid, const std::vector<MechVertexData> *vertex, const std::vector<bool> bc1_state, const SlauSolving::SSCM &Gmatrix,
                                   std::vector<MechVertexContactForce> &vcf, MechNonlinearInf &nlInf)
    {
        // сохраняем значение нового функционала, которое меньше предыдущего
        nlInf.ContactJ_last = nlInf.ContactJ;
        // возвращаем коэффициент регуляризации к начальному значению
        nlInf.ContactBetta = 1;
        // находим такие добавочные перемещения(по нормали к границе), чтобы узлы попали строго на границу
        du.resize(vcf.size()*3);
        db.resize(vcf.size()*3);
        for(int vertexInd = 0; vertexInd < vcf.size(); vertexInd++)
        {
            // обнуление
            for(int k = 0; k < 3; k++)
                du[vertexInd*3 + k] = 0;
            // вершина залезла
            if(vcf[vertexInd].zalez)
            {
                POINT3 C_v;     // центр в плоскости вершины по Z
                C_v = C;
                C_v[2] = grid->vertex[vertexInd][2];
                POINT3 r1, r2, r2_new, delta_r;     // положения вершины: до и после перемещения
                r1 = grid->vertex[vertexInd];          // положение узла сетки до перемещения
                r2 = r1 + (*vertex)[vertexInd].du;     // положение узла сетки после перемещения
                double R0 = (r2 - C_v).abs();  // на таком уравне находимся, надо попасть на уровень R
                {
                    VECTOR3 n = vcf[vertexInd].normal;
                    // |r2 + x*n| = R
                    // |r2 + x*n| - R = 0
                    // функция растёт с ростом t
                    double t = 0;
                    double dt = R - R0;
                    for(;;)
                    {
                        double R_old = (r2 + n*t - C_v).abs();
                        t += dt;
                        double R_new = (r2 + n*t - C_v).abs();
                        if((R_old - R)*(R_new - R) <= 0)
                        {
                            // перешагнули окружность - меняем направление и уменьшаем
                            dt = -dt/2;
                        }
                        if(dt < 1.e-20)
                            break;
                    }
                    delta_r = n*t;   // добавочное перемещение, по которому будут расчитаны силы
                }
                //r2_new = C_v + (r2 - C_v)/R0 * R;
                //r2_new = C_v + (r2 - C_v)/R0 * (R0 + (R-R0)*1);
                //delta_r = r2_new - r2;   // добавочное перемещение, по которому будут расчитаны силы
                for(int k = 0; k < 3; k++)
                if(!bc1_state[vertexInd*3 + k]) // если для этого перемещения есть первое краевое, то не трогаем
                    du[vertexInd*3 + k] = delta_r[k];
                //POINT3 v = grid->vertex[vertexInd] + (*vertex)[vertexInd].du;   // положение узла сетки после перемещения
                //double h0 = R - sqrt(SQR(v[0] - C[0]) + SQR(v[1] - C[1]));
            }
        }
        // по добавочным перемещениям находим силы
        SlauSolving::SSCMmulVector(Gmatrix, du, db);              // A*x -> r
        // проецируем найденные силы на нормаль
        for(int vertexInd = 0; vertexInd < vcf.size(); vertexInd++)
        {
            // вершина залезла
            if(vcf[vertexInd].zalez)
            {
                double dF0 = vcf[vertexInd].dF;    // эту силу приложили на текущем шаге
                double dF_new;                     // искомая новая сила (направление)
                VECTOR3 dF_new_vector;
                for(int k = 0; k < 3; k++)
                    dF_new_vector[k] = db[vertexInd*3 + k];    // если краевое, то и так 0
                // проецируем
                dF_new_vector = (dF_new_vector*vcf[vertexInd].normal)*vcf[vertexInd].normal;
                dF_new = dF0 + (dF_new_vector*vcf[vertexInd].normal)*1;
                // возвращаем коэффициент регуляризации к начальному значению
                vcf[vertexInd].w = 1;
                // обновление направления
                vcf[vertexInd].dF_new = dF_new;
                // сохраняем dF0, h0 для следующих итераций
                vcf[vertexInd].dF1 = dF0;
                vcf[vertexInd].h1 = vcf[vertexInd].h;
            }
        }
    }

    void solveContactDirection_AxNormal(const Grid::Grid3D *grid, const std::vector<MechVertexData> *vertex, std::vector<MechVertexContactForce> &vcf, MechNonlinearInf &nlInf)
    {
        // возвращаем коэффициент регуляризации к начальному значению
        nlInf.ContactBetta = 1;
        for(int vertexInd = 0; vertexInd < vcf.size(); vertexInd++)
        {
            // вершина залезла
            if(vcf[vertexInd].zalez)
            {
                POINT3 C_v;     // центр в плоскости вершины по Z
                C_v = C;
                C_v[2] = grid->vertex[vertexInd][2];
                POINT3 r1, r2, delta_r;     // положения вершины: до и после перемещения
                r1 = grid->vertex[vertexInd];          // положение узла сетки до перемещения
                r2 = r1 + (*vertex)[vertexInd].du;     // положение узла сетки после перемещения
                double R0 = (r2 - C_v).abs();  // на таком уравне находимся, надо попасть на уровень R
                double t = 0;
                {
                    VECTOR3 n = vcf[vertexInd].normal;
                    // |r2 + x*n| = R
                    // |r2 + x*n| - R = 0
                    // функция растёт с ростом t
                    double dt = R - R0;
                    for(;;)
                    {
                        double R_old = (r2 + n*t - C_v).abs();
                        t += dt;
                        double R_new = (r2 + n*t - C_v).abs();
                        if((R_old - R)*(R_new - R) <= 0)
                        {
                            // перешагнули окружность - меняем направление и уменьшаем
                            dt = -dt/2;
                        }
                        if(dt < 1.e-20)
                            break;
                    }
                    delta_r = n*t;   // добавочное перемещение, по которому будут расчитаны силы
                }
                // по добавочным перемещениям находим обновлённое направление
                double dF0 = vcf[vertexInd].dF;    // эту силу приложили на текущем шаге
                double dF_new;                     // искомая новая сила (направление)
                dF_new = dF0 + t*vcf[vertexInd].forseCoeff.abs();
                vcf[vertexInd].dF_new = dF_new;
                // возвращаем коэффициент регуляризации к начальному значению
                vcf[vertexInd].w = 1;
                // сохраняем dF0, h0 для следующих итераций
                vcf[vertexInd].dF1 = dF0;
                vcf[vertexInd].h1 = vcf[vertexInd].h;
            }
        }
    }

    void solveContactDirection_Kh(const Grid::Grid3D * /*grid*/, const std::vector<MechVertexData> * /*vertex*/, std::vector<MechVertexContactForce> &vcf, MechNonlinearInf &nlInf)
    {
        double K = 1.e10;
        // возвращаем коэффициент регуляризации к начальному значению
        nlInf.ContactBetta = 1;
        for(int vertexInd = 0; vertexInd < vcf.size(); vertexInd++)
        {
            // вершина залезла
            if(vcf[vertexInd].zalez)
            {
                // находим обновлённое направление
                double dF0 = vcf[vertexInd].dF;    // эту силу приложили на текущем шаге
                double dF_new;                     // искомая новая сила (направление)
                double h = vcf[vertexInd].h;
                //dF_new = dF0 + sqrt(abs(h)*K)*SIGN(h);
                dF_new = dF0 + h*K;
                vcf[vertexInd].dF_new = dF_new;
                // возвращаем коэффициент регуляризации к начальному значению
                vcf[vertexInd].w = 1;
                // сохраняем dF0, h0 для следующих итераций
                vcf[vertexInd].dF1 = dF0;
                vcf[vertexInd].h1 = vcf[vertexInd].h;
            }
        }
    }

/*
// параметры теста - деформация балки
struct Beam_parameters
{
    Solid::MechBoundaryCondition1Source
        bc1_x,      // (0) q.x = const на плоскости x = center
        bc1_y,      // (1) q.y = const на плоскости y = center
        bc1_z,      // (2) q.z = const на плоскости z = center
        bc1_xyz,    // (3) q = const в центре
        bc1_x1,     // (4) q.x = const на плоскости x = min
        bc1_x2,     // (5) q.x = const на плоскости x = max
        bc1_y1,     // (6) q.y = const на плоскости y = min
        bc1_y2,     // (7) q.y = const на плоскости y = max
        bc1_z1,     // (8) q.z = const на плоскости z = min
        bc1_z2;     // (9) q.z = const на плоскости z = max
    Solid::MechBoundaryCondition2Source
        bc2_x1,     // (0)
        bc2_x2,     // (1)
        bc2_y1,     // (2)
        bc2_y2,     // (3)
        bc2_z1,     // (4)
        bc2_z2,     // (5) // поверхностные силы, действующие на балку с разных сторон
        bc2_bang;   // (6) // поверхностная сила действующая на квадратную область в центре балки
};

enum BeamTaskType
{
    BeamTask_x,
    BeamTask_xyz,
    BeamTask_x_y,
    BeamTask_around_z
};

struct BeamTask
{
    GridRegularParallelepiped_parameters rpp;   // сетка
      BeamTaskType type;                          // режим нагружения
      double P;                                   // поверхностные силы, которые требуется приложить
      int materialNumber;                         // номер материала
};
*/


/*
void MainWindow::on_res_graph_clicked()
{
    /*
    //P/3*(1/(3*K))
    if(out == nullptr) return;

    QChart *chart = new QChart();
    // положение легенды
    //chart->legend()->hide();
    chart->legend()->setAlignment(Qt::AlignBottom);
    // общее название графиков
    //chart->setTitle("привет");
    QValueAxis *axisX;
    QValueAxis *axisY;
    double min = +1.e200;
    double max = -1.e200;
    //int count = ui->res_functions->selectionModel()->selectedRows().count();
    // находим список строк, среди яцеек которых есть выделенные
    int count = ui->res_functions->selectionModel()->selectedIndexes().count();
    std::vector<int> rows;
    for( int i = 0; i < count; i++)
    {
        bool finded = false;
        int row = ui->res_functions->selectionModel()->selectedIndexes().at(i).row();
        for(size_t j = 0; j < rows.size() && !finded; j++)
        {
            if(rows[j] == row)
                finded = true;
        }
        if(!finded)
            rows.insert(rows.end(), row);
    }
    //return;
    for(size_t i = 0; i < rows.size(); i++)
    {
        // выбрана функция номер row
        //QList<QTableWidgetItem *> selectedList = ui->res_functions->selected->selectedItems();
        //int row = ui->res_functions->selectionModel()->selectedRows().at(i).row();
        int row = rows[i];
        if(row < 0 || row >= (int)outFun.size()) return;
        Ui_Solid_OutFunction &f = outFun[row];
        Function x, y;
        genResFun(f.name_x, f.expression_x, x);
        genResFun(f.name_y, f.expression_y, y);
        // подготовка данных
        QLineSeries *series = new QLineSeries();
        // значения аргументов для каждого шага
        double P = 0;
        for (size_t j = 0; j <= (size_t)sp.PN; j++)  // цикл по шагам//(*out)[0].pres.size()
        {
            double Padd;
            if (j == 0)
                Padd = task.bt.P * sp.proportion_of_P_elastic;
            else
                Padd = task.bt.P * (1 - sp.proportion_of_P_elastic) / sp.PN;
            P += Padd;

            genResFunArgs(x, j, P);
            //return;
            genResFunArgs(y, j, P);

            double valueX = x.solve();
            double valueY = y.solve();
            series->append(valueX, valueY);
            if(valueY < min) min = valueY;
            if(valueY > max) max = valueY;
        }

        // название графика
        series->setName(f.name.c_str());

        // добавление графика
        chart->addSeries(series);

        if(i == 0)
        {
            // оси
            axisX = new QValueAxis;
            axisX->setLabelFormat("%.3g");
            axisX->setTitleText(x.name.c_str());
            axisX->setMinorTickCount(1);
            axisX->setTickCount(5);
            axisY = new QValueAxis;
            axisY->setLabelFormat("%.3g");
            axisY->setTitleText(y.name.c_str());
            axisY->setMinorTickCount(1);
            axisY->setTickCount(5);
            QFont font; //фонт
            font.setPixelSize(18);
            axisY->setLabelsFont(font);
            axisX->setLabelsFont(font);
            chart->addAxis(axisX, Qt::AlignBottom);
            chart->addAxis(axisY, Qt::AlignLeft);

            //axisY->setTitleVisible(false);

        }
        // оси
        series->attachAxis(axisX);
        series->attachAxis(axisY);
        if(i == rows.size() - 1)
        {
            axisY->setMax(max);
            axisY->setMin(min);
        }
    }
    chart->legend()->hide();    // скрывает названия графиков
    // Customize chart title
        //QFont font;
        //font.setPixelSize(18);
        //chart->setFont(font);
        //chart->setTitleBrush(QBrush(Qt::white));
        //chart->setTitle("Customchart example");

    outChartView->setChart(chart);
    delete outChart;
    outChart = chart;
    */
}

*/


/*
        static POINT3 Gauss2_pointSource[8];
        static POINT3 Gauss3_pointSource[27];

        Integration::Integrator integrationCubeSource;
        integrationCubeSource.init3D(Integration::IntegrationType::Gauss2);
        for (size_t pInd = 0; pInd < integrationCubeSource.size; pInd++)    // индекс точки Гаусса
            Gauss2_pointSource[pInd] = integrationCubeSource.p3[pInd];
        integrationCubeSource.init3D(Integration::IntegrationType::Gauss3);
        for (size_t pInd = 0; pInd < integrationCubeSource.size; pInd++)    // индекс точки Гаусса
            Gauss3_pointSource[pInd] = integrationCubeSource.p3[pInd];

*/





    // отладка
    /*
    // сохраниение локальных матриц и векторов
    char fn_G[1000];
    char fn_b[1000];
    char fn_x[1000];

    sprintf(fn_G, "G_%d", nlInf.numStep);
    sprintf(fn_b, "b_%d", nlInf.numStep);
    sprintf(fn_x, "x_%d", nlInf.numStep);
    FILE *fb = fopen(fn_b, "w");
    for(int feInd = 0; feInd < (*fe).size(); feInd++)
    {
        MechFeData_LinearHexagon *feEl = (MechFeData_LinearHexagon *)(*fe)[feInd];
        fprintf(fb, "%d\n", feInd);
        for(int i = 0; i < 24; i++)
        {
            fprintf(fb, "%le\n", feEl->bLocal[i]);
        }
        fprintf(fb, "\n");
    }
    fclose(fb);
    FILE *fG = fopen(fn_G, "w");
    for(int feInd = 0; feInd < (*fe).size(); feInd++)
    {
        MechFeData_LinearHexagon *feEl = (MechFeData_LinearHexagon *)(*fe)[feInd];
        fprintf(fG, "%d\n", feInd);
        for(int i = 0; i < 300; i++)
        {
            fprintf(fG, "%le\n", feEl->GLocalL[i]);
        }
        fprintf(fG, "\n");
    }
    fclose(fG);
    */

    /*
    FILE *fx = fopen(fn_x, "w");
    for(int i = 0; i < dq.size(); i++)
        fprintf(fx, "x[%d] = %le\n", i, dq[i]);
    fclose(fx);
    */


    /*
    fprintf(ff, "faceIndex:\n");
    fprintf(ff, "%d", faceIndex);
    fprintf(ff, "\n");
    fprintf(ff, "integrationFoursquare_detJ:\n");
    fprintf(ff, "%le", integrationFoursquare.detJ);
    fprintf(ff, "\n");
    fprintf(ff, "integrationFoursquare_p2:\n");
    for (int i = 0; i < integrationFoursquare.size; i++)
    {
        fprintf(ff, "(%le, %le) ", integrationFoursquare.p2[i].x[0], integrationFoursquare.p2[i].x[1]);
    }
    fprintf(ff, "\n");
    fprintf(ff, "integrationFoursquare_w:\n");
    for (int i = 0; i < integrationFoursquare.size; i++)
    {
        fprintf(ff, "%le ", integrationFoursquare.w[i]);
    }
    fprintf(ff, "\n");


    POINT3 v[8];
    // копируем вершины шестигранника
    feEl->getVertexes(grid->vertex, grid->vertexForCurvature, v);
    fprintf(ff, "hexagon:\n");
    for (int vInd = 0; vInd < 4; vInd++)
    {
        fprintf(ff, "(%le, %le, %le) ", v[vInd][0], v[vInd][1], v[vInd][2]);
    }
    fprintf(ff, "\n");
    for (int vInd = 4; vInd < 8; vInd++)
    {
        fprintf(ff, "(%le, %le, %le) ", v[vInd][0], v[vInd][1], v[vInd][2]);
    }
    fprintf(ff, "\n");



    fprintf(ff, "hexagonNormal:\n");
    for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
    {
        for(int i = 0; i < 3; i++)
            fprintf(ff, "%le ", hexagonNormal[valIndex][i]);
        fprintf(ff, "\n");
    }
    fprintf(ff, "\n");

    fprintf(ff, "vectorP:\n");
    for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
    {
        VECTOR3 vectorP = bc2SourceEl->solveVectorValueInPoint(hexagonNormal[valIndex]);
        for(int i = 0; i < 3; i++)
            fprintf(ff, "%le ", vectorP[i]);
        fprintf(ff, "\n");
    }
    fprintf(ff, "\n");
    fprintf(ff, "basCube:\n");
    for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
    {
        for(int m = 0; m < 8; m++)
            fprintf(ff, "%le ", basCube[m*integrationFoursquare.size + valIndex]);
        fprintf(ff, "\n");
    }
    fprintf(ff, "\n");
    fprintf(ff, "hexagonCoef:\n");
    for (int valIndex = 0; valIndex < integrationFoursquare.size; valIndex++)
    {
        fprintf(ff, "%le\n", hexagonCoef[valIndex]);
    }
    fprintf(ff, "\n");
    */
//fprintf(ff, "%le %le %le\n", vectorP[i], basCube[m*integrationFoursquare.size + valIndex], hexagonCoef[valIndex]);


    // толстостенная сфера
    void genTestSphere(TaskTest &task)
    {
        FILE *f = fopen("___gridInformation.txt", "w");
        // индекс теста
        task.testIndex = 1;
        // индекс теста решателя теплопроводности
        task.thermTestIndex = 0;        // 1 - t(a)=t1, t(b)=t2
                                        // 2 - t(a)=t1, конвекция на b
                                        // 3 - конвекция на r=a, конвекция на r=b
                                        // 4 - подогрев на r=a, t(b)=t2
        //bool fidesys = true;
        bool fidesys = false;
        // подрубка
        if(task.thermTestIndex != 0)
        {
            task.thermTask.enabled = true;
            task.mechTask.enabled = false;
        }
        else
        {
            task.thermTask.enabled = false;
            task.mechTask.enabled = true;
        }
        // общие шаги и сетка
        std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
        Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
        // сетка
        SphereParameters &gp = task.sphere;
        if(fidesys)
        {
            gp.r1 = 2.5;
            gp.r2 = 5;
        }
        else
        {
            gp.r1 = 1;
            gp.r2 = 4;
        }
        gp.N = 8;   // 4 промежутка на четвертьокружностях    //8, 32
        gp.Nparts = 32;
        if(task.thermTestIndex != 0)
        {
            gp.q = 1;
        }
        else
        {
            if(fidesys)
                gp.q = 1. + 0*1./(gp.Nparts/2);//1.0625;//1.0625;//1.03125;
            else
                gp.q = 1. + 1*1./(gp.Nparts/2);//1.0625;//1.0625;//1.03125;
        }
        gp.curvilinear = 1;    // отображение (0 - линейное, 1 - квадратичное)
        gp.buildingMethod = 0; // способ построения (0 - делим дуги, 1 - отображение на куб)
        grid->genSphere(gp);
        //grid->buldOpenSCADModel("setka.scad");
        double NU = 0.3; //0.4

        // параметры
        //task.slausolver_parameters.preconditioning = SlauSolving::Preconditioning::SlauPreconditioning_LLT;

        task.sigmaSolvingType = 0;  // способ расчёта решения (0 - главные напряжения, 1 - проекции площадки)
        int fixGrid;                // 0 - подвижная сетка, 1 - зафиксировать сетку
        int sigma0Mode = 2;         // 1 - dSigma = integral(...), 2 - в приращениях (P->dP)
        if(fidesys)
        {
            fixGrid = 0;
        }
        else
        {
            fixGrid = 0;
        }
        int plasticPlot = 2;            // 0 - упругость, 1 - сплайн, 2 - Безье
        int plasticityCurveMode = 1;    // 0 - eps(sigma), 1 - sigma(eps)
        bool usingSecantMethod = false;//true;
        double epsResidualLimit = 1.e-13;    // желаемая невязка по эквиволентным деформациям
        int nonlinearIterLimit = 100;        // ограничение на количество итераций
        int StepsNumber1 = 1;//8;//8;//128/2;      // шаги
        int StepsNumber2 = 4;//300;//512;//40;//128/2;
        int StepsNumber3 = 1;
        int StepsNumber4 = 1;
        double P = 30;//2.0e7;//2.5e7;//1000;//2.5e7;
        if(fidesys)
            P = 30;
        else
            P = 2.0e7;
        if(task.thermTestIndex != 0)
            P = 0;
        double P1 = 0.1*P;//P/2;//0.4*P;
        double P2 = P;
        double P3 = P*0.9;//P/2.;//P-((P2-P1)/StepsNumber2)/100;//-P*1/2;//-(P2/StepsNumber);
        double P4 = 0;//-P*1/2;//-(P - P2/StepsNumber);
        /*
        double P1 = P*1/2;
        double P2 = P;
        double P3 = P*1/2;
        double P4 = 0;
        */
        int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny

        double Time = 1;
        // шаги
        GlobalStep s;       // шаг
        // первый шаг (1 временной слой, заведомо упругое нагружение P1)
        s.t1 = 0;                          // (имеет значение только для первого глобального шага)
        s.t2 = Time*1;
        s.dt0 = Time/StepsNumber1;
        step->push_back(s);
        // второй шаг (пластичное нагружение P2)
        s.t1 = Time*1;
        s.t2 = Time*2;
        s.dt0 = Time/StepsNumber2;
        step->push_back(s);
        // третий шаг (1 временной слой, начало разгрузки)
        s.t1 = Time*2;
        s.t2 = Time*3;
        s.dt0 = Time/StepsNumber3;
        step->push_back(s);
        // четвертый шаг (1 временной слой, заведомо упругая разгрузка до 0)
        s.t1 = Time*3;
        s.t2 = Time*4;
        s.dt0 = Time/StepsNumber4;
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
            // материал
            task.mechTask.material = new std::vector<MechMaterialSource>(1);
            MechMaterialSource &mechMat = (*task.mechTask.material)[0];
            if(task.thermTestIndex != 0)
            {
                mechMat.elasticSigmaLimit = 1.e10;//2.e7;
            }
            else
            {
                if(fidesys)
                    mechMat.elasticSigmaLimit = 24;//2.e7;
                else
                    mechMat.elasticSigmaLimit = 2.e7;
            }
            mechMat.Ceps = 4./9.;
            mechMat.Csigma = 1;

            if(task.thermTestIndex == 1)
            {
                mechMat.set_E_NU(200*1.e9,NU);        //1.e10;
                mechMat.Talpha = 0.0001;
            }
            else
            {
                if(fidesys)
                {
                    mechMat.set_E_NU(21000., NU);        //1.e10;
                    mechMat.Talpha = 0;
                }
                else
                {
                    mechMat.set_E_NU(1.e10, NU);
                    mechMat.Talpha = 0;
                }
            }
            //mechMat.set_K_G(70.0e9, 26.0e9);
            mechMat.set_D_isotropic();
            mechMat.set_M_sigma();
            mechMat.F = {0,0,0};
            mechMat.ro = 0;//1e10;
            if(task.thermTestIndex != 0)
                setPlasticMaterialCurve(0, plasticityCurveMode, usingSecantMethod, mechMat);
            else
            {
                // setPlasticMaterialCurve(plasticPlot, plasticityCurveMode, Solid::MechStressStrainMethod::initialStrain, usingSecantMethod, mechMat);
                // setPlasticMaterialCurve(plasticPlot, plasticityCurveMode, Solid::MechStressStrainMethod::initialStrain, usingSecantMethod, mechMat);
                 setPlasticMaterialCurve(plasticPlot, plasticityCurveMode, usingSecantMethod, mechMat);
            }
            //mechMat.mode = MechMaterialType::Elasticity;


            // шаги (на каждом шаге задаются вторые краевые условия)
            std::vector<MechGlobalStep> *mechStep = new std::vector<MechGlobalStep>;
            // нагружение - разгрузка
            MechGlobalStep ms;  // шаг для МДТТ
            GlobalStep s;       // шаг
            std::vector<MechBoundaryCondition2Source_base *> *bc2Source;
            ms.slausolverParameters = task.slausolver_parameters;
            ms.timeMode = 0;    // 0 - квазистатическая задача
            ms.sigma0Mode = sigma0Mode;  // 1 - dSigma = integral(...), 2 - приращения
            ms.fixGrid = fixGrid;     // 1 - зафиксировать сетку
            ms.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
            ms.controlMode = 0; // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
            ms.iterLimit = nonlinearIterLimit;
            ms.slauResidualLimit = 1000;
            ms.epsResidualLimit = epsResidualLimit;
            ms.sigmaResidualLimit = 1000;
            // первый шаг (1 временной слой, заведомо упругое нагружение P1)
            s = (*step)[0];
            setBc2Sphere(bc2Source, 0, P1, s.t1, s.t2);
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
            // второй шаг (пластичное нагружение P2)
            s = (*step)[1];
            setBc2Sphere(bc2Source, P1, P2, s.t1, s.t2);
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
            // третий шаг (1 временной слой, начало разгрузки)
            s = (*step)[2];
            setBc2Sphere(bc2Source, P2, P3, s.t1, s.t2);
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
            // четвертый шаг (1 временной слой, заведомо упругая разгрузка до 0)
            s = (*step)[3];
            setBc2Sphere(bc2Source, P3, P4, s.t1, s.t2);
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
            task.mechTask.mechStep = mechStep;
            // первые краевые условия
            std::vector<MechBoundaryCondition1Source> *mechBc1 = new std::vector<MechBoundaryCondition1Source>(5);
            (*mechBc1)[0].mode = {{ 0, -1, -1}};
            (*mechBc1)[0].u0 =   {{ 0, -1, -1}};
            (*mechBc1)[1].mode = {{-1,  0, -1}};
            (*mechBc1)[1].u0 =   {{-1,  0, -1}};
            (*mechBc1)[2].mode = {{-1, -1,  0}};
            (*mechBc1)[2].u0 =   {{-1, -1,  0}};
            (*mechBc1)[3].mode = {{-1, -1, -1}};
            (*mechBc1)[3].u0 =   {{-1, -1, -1}};
            (*mechBc1)[4].mode = {{-1, -1, -1}};
            (*mechBc1)[4].u0 =   {{-1, -1, -1}};
            //(*mechBc1)[4].mode = {{0, 0, 0}};
            //(*mechBc1)[4].u0 =   {{0, 0, 0}};
            task.mechTask.bc1Source = mechBc1;
            // инициализация начальных деформаций и напряжений (=0)
            task.mechTask.fe = new std::vector<MechFeData_base *>;
            task.mechTask.fe->resize(grid->fe.size());
            for (size_t i = 0; i < grid->fe.size(); i++)
            {
                MechFeData_LinearHexagon_homogeny_E_NLE_EP *feDataEl = new MechFeData_LinearHexagon_homogeny_E_NLE_EP;
                feDataEl->init(Integration::IntegrationType::Gauss3);
                (*task.mechTask.fe)[i] = feDataEl;
            }
            // инициализация начальных перемещений (=0)
            task.mechTask.vertex = new std::vector<MechVertexData>;
            for (size_t i = 0; i < grid->vertex.size(); i++)
            {
                MechVertexData v_el;
                for (int j = 0; j < 3; j++)
                {
                    v_el.sum_du[j] = 0;
                    v_el.du[j] = 0; // не имеет значения
                }
                task.mechTask.vertex->push_back(v_el);
            }

            // для криволинейных границ ###
            task.mechTask.vertexForCurvature = new std::vector<MechVertexData>;
            for (size_t i = 0; i < grid->vertexForCurvature.size(); i++)
            {
                MechVertexData v_el;
                for (int j = 0; j < 3; j++)
                {
                    v_el.sum_du[j] = 0;
                    v_el.du[j] = 0; // не имеет значения
                }
                task.mechTask.vertexForCurvature->push_back(v_el);
            }

            // инициализация начальных скоростей (=0) и ускорений (=0)
            task.mechTask.V0 = new Vector(grid->vertex.size()*3);
            task.mechTask.dV0 = new Vector(grid->vertex.size()*3);
            for (size_t i = 0; i < task.mechTask.V0->size(); i++)
            {
                (*task.mechTask.V0)[i] = 0;
                (*task.mechTask.dV0)[i] = 0;
            }

            // вторые краевые
            task.mechTask.bc2 = new std::vector<MechBoundaryCondition2>;
            MechBoundaryCondition2 bc2_el;
            bc2_el.FEsurfaceInd = 0;
            bc2_el.bc2SourceIndex = 0;
            (*task.mechTask.bc2).push_back(bc2_el);
            bc2_el.FEsurfaceInd = 1;
            bc2_el.bc2SourceIndex = 1;
            (*task.mechTask.bc2).push_back(bc2_el);

            // Контакт (отсутствует)
            // поверхность: жёсткий неподвижный цилиндр
            task.mechTask.rigidSurface = new std::vector<Surface_base *>;
            // контакты КЭ-поверхность - аналитическая поверхность
            task.mechTask.CC_FE_Rigid = new std::vector<ContactCondition_FE_Rigid>;
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
            // материал
            task.thermTask.material = new std::vector<ThermMaterialSource>(1);
            ThermMaterialSource &thermMat = (*task.thermTask.material)[0];
            thermMat.ro = 1;     // плотность
            thermMat.c = 5;      // удельная теплоёмкость
            thermMat.f = 0;      // мощность внутренних объёмных источников (стоков) тепла
            double L0 = 10;//10;
            thermMat.L =         // тензор теплопроводности
            {{
                {L0, 0, 0},
                {0, L0, 0},
                {0, 0, L0},
            }};
            int timeMode = 0;
            // Тест1
            if(task.thermTestIndex == 1)
                setBcSphereT_1(task.thermTask.bc1SourceT, task.thermTask.bc2SourceT,
                               1,   // температура на r = a
                               2    // температура на r = b
                               );
            // Тест2
            if(task.thermTestIndex == 2)
                setBcSphereT_2(task.thermTask.bc1SourceT, task.thermTask.bc2SourceT,
                               1,   // температура на r = a
                               2,   // температура окружающей среды r > b
                               1    // коэффициент конвективного теплообмена
                               );
            // Тест3
            if(task.thermTestIndex == 3)
                setBcSphereT_3(task.thermTask.bc1SourceT, task.thermTask.bc2SourceT,
                               1,   // температура окружающей среды r < a
                               10,  // коэффициент конвективного теплообмена
                               2,   // температура окружающей среды r > b
                               20   // коэффициент конвективного теплообмена
                               );
            // Тест4
            if(task.thermTestIndex == 4)
            {
                setBcSphereT_4(task.thermTask.bc1SourceT, task.thermTask.bc2SourceT, thermMat,
                               100, // плотность набегающего теплового потока на r = a
                               15,  // коэффициент теплопроводности
                               1    // температура на r = b
                               );
            }
            // Тест5
            if(task.thermTestIndex == 5)
            {
                setBcSphereT_5(task.thermTask.bc1SourceT, task.thermTask.bc2SourceT, thermMat,
                               10, // плотность набегающего теплового потока на r = a
                               1   // коэффициент теплопроводности
                               );
                timeMode = 1;
            }
            // шаги
            std::vector<ThermGlobalStep> *thermStep = new std::vector<ThermGlobalStep>;
            ThermGlobalStep ts;  // шаг для теплопроводности
            ts.slausolverParameters = task.slausolver_parameters;
            ts.timeMode = timeMode;    // статика
            thermStep->push_back(ts);
            thermStep->push_back(ts);
            thermStep->push_back(ts);
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
                if(task.thermTestIndex == 5)
                {
                    double a = task.grid->vertex[0].abs();
                    double R = grid->vertex[i].abs();
                    double T_a_ch = 100;
                    double Q0 = (*task.thermTask.bc2SourceT)[0].q * 4*PI*SQR(a);
                    T0 = T_a_ch + Q0/(4*PI*(*task.thermTask.material)[0].L.m[0][0])*(1./R - 1./a);
                }
                v_el.T = T0;
                v_el.newT = T0;
                task.thermTask.vertex->push_back(v_el);
            }
        }





        // сохренение информации о сетке в файл
        fprintf(f, "r1 = %lf, r2 = %lf\n", gp.r1, gp.r2);
        fprintf(f, "Grid: %dx%dx%d, q = %le\n", gp.N, gp.N, gp.Nparts, gp.q);
        if(gp.curvilinear == 0)
            fprintf(f, "Linear\n");
        if(gp.curvilinear == 1)
            fprintf(f, "Quadratic\n");
        if(gp.buildingMethod == 0)
            fprintf(f, "buildingMethod = Divide arcs\n");
        if(gp.buildingMethod == 1)
            fprintf(f, "Quadratic = Cube to sphere\n");
        fprintf(f, "FENumber = %d\n", (int)grid->fe.size());
        fprintf(f, "VertexesNumber = %d\n", (int)grid->vertex.size());
        fprintf(f, "VertexesForCurvatureNumber = %d\n", (int)grid->vertexForCurvature.size());
        fprintf(f, "bc1Number = %d\n", (int)grid->bc1.size());
        fprintf(f, "FEsurface = %d\n", (int)grid->FESurface.size());
        for (size_t FEsurfaceInd = 0; FEsurfaceInd < grid->FESurface.size(); FEsurfaceInd++)
        {
            fprintf(f, "FEsurface[%d]Size = %d\n", (int)FEsurfaceInd, (int)grid->FESurface[FEsurfaceInd].face.size());
        }
        //grid->buldOpenSCADModel("setka.scad");
        // сохренение информации о параметрах в файл
        fprintf(f, "fixGrid = %d\n", fixGrid);
        fprintf(f, "plasticPlot = %d\n", plasticPlot);
        fprintf(f, "epsResidualLimit = %le\n", epsResidualLimit);
        fprintf(f, "nonlinearIterLimit = %d\n", nonlinearIterLimit);
        fprintf(f, "StepsNumber1 = %d\n", StepsNumber1);
        fprintf(f, "StepsNumber2 = %d\n", StepsNumber2);
        fprintf(f, "StepsNumber3 = %d\n", StepsNumber3);
        fprintf(f, "StepsNumber4 = %d\n", StepsNumber4);
        fprintf(f, "P = %lf\n", P);
        fprintf(f, "P1 = %lf\n", P1);
        fprintf(f, "P2 = %lf\n", P2);
        fprintf(f, "P3 = %lf\n", P3);
        fprintf(f, "P4 = %lf\n", P4);
        fprintf(f, "sigma0Mode = %d\n", sigma0Mode);
        fprintf(f, "HomogenyMode = %d\n", HomogenyMode);

        fclose(f);
    }

    // Формовка
    void genTestFormovka(TaskTest &task)
    {
        // индекс теста
        task.testIndex = 5;
        // подрубка
        task.mechTask.enabled = true;
        task.thermTask.enabled = false;
        // общие шаги и сетка
        std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
        Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
        // сетка
        FormovkaParameters &gp = task.formovka;
        // способ представления поверхности
        //gp.surfaceType = SurfaceType::none;
        gp.surfaceType = SurfaceType::AnaliticalSurface_Cylinder;
        //gp.surfaceType = SurfaceType::InterpolantSurface_Hermite3;
        //gp.surfaceType = SurfaceType::InterpolantSurface_Lagrange3;
            //gp.contact_mode = SurfaceType::FiniteElementSurface;//нет поиска пересечения

        // способ декомпозиции поверхности
        gp.decompositionType = DecompositionType::none;
        //gp.decompositionType = DecompositionType::cubeVertexIndexation;
        //gp.decompositionType = DecompositionType::surfaceAsSpheres;

        // подрубка универсальной оптимизации
        //bool noContactRadiusOptimization = false;
        bool noContactRadiusOptimization = true;


        int NN = 5; //9//5
        gp.p1 = VECTOR3(-10, -1, -1.0);
        gp.p2 = VECTOR3(10, -0.5, -0.5);
        //gp.p1 = VECTOR3(-10, -1, -1);
        //gp.p2 = VECTOR3(10, -0.5, -0.5);
        gp.N = VECTOR3_uint(NN*40, NN, 2);//8
        gp.N_react = 0;//(NN+1)/2;
        gp.N_p = (NN+1)/2;
        /*
        gp.p1 = VECTOR3(-10, -1, -1);
        gp.p2 = VECTOR3(10, 1, 2);
        gp.N = VECTOR3_int(90, 9, 2);
        gp.N_react = 5;
        gp.N_p = 5;
        */

        grid->genFormovka(gp);
        //for(;;);
        //grid->buldOpenSCADModel("setka.scad");
        //for(;;);
        // параметры
        int fixGrid = 0;                // 0 - подвижная сетка, 1 - зафиксировать сетку
        int plasticPlot = 0;            // 0 - упругость, 1 - сплайн, 2 - Безье
        int plasticityCurveMode = 1;    // 0 - eps(sigma), 1 - sigma(eps)
        bool usingSecantMethod = false;
        double epsResidualLimit = 1.e-10;            // желаемая относительная погрешность по эквиволентным деформациям
        double contactDeltaFResidualLimit = 1.e-7;   // желаемая относительная погрешность по силе реакции опоры
        int nonlinearIterLimit = 200;        // ограничение на количество итераций
        int StepsNumber1 = 1;        //
        int StepsNumber2 = 40;    // механическое нагружение
        int StepsNumber3 = 1;        // механическая разгрузка
      //int StepsNumber4 = 1;        // температурное нагружение
        double P1 = 0;
        double P2 = 0.20e8/2;//-0.20e8/2;//-0.20e8/3.34;//-0.20e8/3.7(если 400 шагов);//-0.20e8/3.36;(если 40 шагов)  //-0.05e8;//-0.08e7;//-0.77e7;//-0.075e8;//-0.075e8;//-1.5e8;//-1.1e8;//-2.e8;//-1.e7*1;//-5.e7*1; //-5.e7/20*1; - для пластичности
        // 9, 100, -0.20e8/3.34, 2.e7*10
        double elasticSigmaLimit = 2.e7*10;
        //double P2 = -0.20e8;
        //double elasticSigmaLimit = 2.e7*16;
      //double P3 = P2;
        int sigma0Mode = 2;         // 1 - dSigma = integral(...), 2 - в приращениях (P->dP)
        int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny

        double Time = 100;
        // шаги
        GlobalStep s;       // шаг
        // 1-й шаг (нагревание)
        s.t1 = 0;           // (имеет значение только для первого глобального шага)
        s.t2 = Time*1;
        s.dt0 = Time/StepsNumber1;
        step->push_back(s);
        // 2-й шаг (нагружение P)
        for(int i = 0; i < StepsNumber2; i++)
        {
            s.t1 = Time*1 + Time*i/StepsNumber2;
            s.t2 = Time*1 + Time*(i + 1)/StepsNumber2;
            s.dt0 = Time/StepsNumber2;
            step->push_back(s);
        }
        /*
        // 3-й шаг (разгрузка)
        for(int i = 0; i < StepsNumber3; i++)
        {
            s.t1 = Time*2 + Time*i/StepsNumber3;
            s.t2 = Time*2 + Time*(i + 1)/StepsNumber3;
            s.dt0 = Time/StepsNumber3;
            step->push_back(s);
        }
        // 4-й шаг (остужение)
        s.t1 = Time*3;
        s.t2 = Time*4;
        s.dt0 = Time/StepsNumber4;
        step->push_back(s);
        */
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
            // материал
            task.mechTask.material = new std::vector<MechMaterialSource>(1);
            MechMaterialSource &mechMat = (*task.mechTask.material)[0];
            mechMat.elasticSigmaLimit = elasticSigmaLimit;
            mechMat.Ceps = 4./9.;
            mechMat.Csigma = 1;
            mechMat.set_E_NU(1.e10,0.3);
            //mechMat.set_E_NU(1.e10,0.4);
            //mechMat.set_D_isotropic_XY();
            //mechMat.set_M_sigma_XY();
            //mechMat.set_D_isotropic();
            //mechMat.set_M_sigma();
            mechMat.set_D_isotropic();      // обычная матрица D
            mechMat.set_M_sigma();       //mechMat.set_M_sigma_XY();       // эквивалентные напряжения плоские
            mechMat.F = {0,0,0};
            mechMat.ro = 0;     // Плотность
            mechMat.Talpha = 1.e-5; // Коэффициент линейного расширения
            setPlasticMaterialCurve(plasticPlot, plasticityCurveMode, usingSecantMethod, mechMat);
            // шаги (на каждом шаге задаются вторые краевые условия)
            std::vector<MechGlobalStep> *mechStep = new std::vector<MechGlobalStep>;
            // нагревание - нагружение - разгрузка
            MechGlobalStep ms;  // шаг для МДТТ
            GlobalStep s;       // шаг
            std::vector<MechBoundaryCondition2Source_base *> *bc2Source;
            ms.slausolverParameters = task.slausolver_parameters;
            ms.timeMode = 0;    // 0 - квазистатическая задача
            ms.sigma0Mode = sigma0Mode; // 1 - dSigma = integral(...), 2 - приращения
            ms.fixGrid = fixGrid;     // 1 - зафиксировать сетку
            ms.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
            ms.controlMode = 0;     // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
            ms.iterLimit = nonlinearIterLimit;
            ms.slauResidualLimit = 10000;
            ms.epsResidualLimit = epsResidualLimit;
            ms.sigmaResidualLimit = 10000;
            ms.contactEndPointResidualLimit = 1000;//1.e-10;//1.e-13;
            ms.contactDeltaFResidualLimit = contactDeltaFResidualLimit;//1.e-13;
            // 1-й шаг (нагревание)
            s = (*step)[0];
            setBc2Sphere(bc2Source, 0, 0, s.t1, s.t2);
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
            // 2-й шаг (нагружение)
            for(int i = 0; i < StepsNumber2; i++)
            {
                s.t1 = Time*1 + Time*i/StepsNumber2;
                s.t2 = Time*1 + Time*(i + 1)/StepsNumber2;
                setBc2Sphere(bc2Source, P1 + (P2 - P1)*i/StepsNumber2, P1 + (P2 - P1)*(i + 1)/StepsNumber2, s.t1, s.t2);
                ms.bc2Source = bc2Source;
                mechStep->push_back(ms);
            }
            /*
            // 3-й шаг (разгрузка)
            for(int i = 0; i < StepsNumber3; i++)
            {
                s.t1 = Time*2 + Time*i/StepsNumber3;
                s.t2 = Time*2 + Time*(i + 1)/StepsNumber3;
                //setBc2Sphere(bc2Source, P2, P2, s.t1, s.t2);
                setBc2Sphere(bc2Source, P2 + (P3 - P2)*i/StepsNumber3, P2 + (P3 - P2)*(i + 1)/StepsNumber3, s.t1, s.t2);
                ms.bc2Source = bc2Source;
                mechStep->push_back(ms);
            }
            // 4-й шаг (остужение)
            s = (*step)[1 + StepsNumber2 + StepsNumber3];
            setBc2Sphere(bc2Source, 0, 0, s.t1, s.t2);
            ms.bc2Source = bc2Source;
            mechStep->push_back(ms);
            */
            task.mechTask.mechStep = mechStep;
            // первые краевые условия
            std::vector<MechBoundaryCondition1Source> *mechBc1 = new std::vector<MechBoundaryCondition1Source>(5);
            (*mechBc1)[0].mode = {{ 0, -1, -1}};
            (*mechBc1)[0].u0 =   {{ 0, -1, -1}};
            //(*mechBc1)[1].mode = {{0,  0, 0}};
            //(*mechBc1)[1].u0 =   {{0,  0, 0}};
            (*mechBc1)[1].mode = {{-1,  0, -1}};
            (*mechBc1)[1].u0 =   {{-1,  0, -1}};
            //(*mechBc1)[2].mode = {{-1, -1,  -1}};
            //(*mechBc1)[2].u0 =   {{-1, -1,  -1}};
            (*mechBc1)[2].mode = {{-1, -1,  0}};
            (*mechBc1)[2].u0 =   {{-1, -1,  0}};

            (*mechBc1)[3].mode = {{-1, -1,  -1}};
            (*mechBc1)[3].u0 =   {{-1, -1,  -1}};
            (*mechBc1)[4].mode = {{-1, -1,  -1}};
            (*mechBc1)[4].u0 =   {{-1, -1,  -1}};

            task.mechTask.bc1Source = mechBc1;
            // инициализация начальных деформаций и напряжений (=0)
            task.mechTask.fe = new std::vector<MechFeData_base *>;
            task.mechTask.fe->resize(grid->fe.size());
            for (size_t i = 0; i < grid->fe.size(); i++)
            {
                MechFeData_LinearHexagon_homogeny_E_NLE_EP *feDataEl = new MechFeData_LinearHexagon_homogeny_E_NLE_EP;
                feDataEl->init(Integration::IntegrationType::Gauss3);
                (*task.mechTask.fe)[i] = feDataEl;
            }
            // инициализация начальных перемещений (=0)
            task.mechTask.vertex = new std::vector<MechVertexData>;
            for (size_t i = 0; i < grid->vertex.size(); i++)
            {
                MechVertexData v_el;
                for (int j = 0; j < 3; j++)
                {
                    v_el.sum_du[j] = 0;
                    v_el.du[j] = 0; // не имеет значения
                }
                task.mechTask.vertex->push_back(v_el);
            }
            // для криволинейных границ ###
            task.mechTask.vertexForCurvature = new std::vector<MechVertexData>;
            for (size_t i = 0; i < grid->vertexForCurvature.size(); i++)
            {
                MechVertexData v_el;
                for (int j = 0; j < 3; j++)
                {
                    v_el.sum_du[j] = 0;
                    v_el.du[j] = 0; // не имеет значения
                }
                task.mechTask.vertexForCurvature->push_back(v_el);
            }
            // инициализация начальных скоростей (=0) и ускорений (=0)
            task.mechTask.V0 = new Vector(grid->vertex.size()*3);
            task.mechTask.dV0 = new Vector(grid->vertex.size()*3);
            for (size_t i = 0; i < task.mechTask.V0->size(); i++)
            {
                (*task.mechTask.V0)[i] = 0;
                (*task.mechTask.dV0)[i] = 0;
            }

            // вторые краевые
            task.mechTask.bc2 = new std::vector<MechBoundaryCondition2>;
            MechBoundaryCondition2 bc2_el;
            bc2_el.FEsurfaceInd = 0;
            bc2_el.bc2SourceIndex = 0;
            (*task.mechTask.bc2).push_back(bc2_el);
            bc2_el.FEsurfaceInd = 1;
            bc2_el.bc2SourceIndex = 1;
            (*task.mechTask.bc2).push_back(bc2_el);

            // Контакт

            // поверхность: жёсткий цилиндр
            task.mechTask.rigidSurface = new std::vector<Surface_base *>;
            (*task.mechTask.rigidSurface).resize(3);
            (*task.mechTask.rigidSurface)[0] = grid->analiticalSurface[0];   // поверхность заданная аналитически
            (*task.mechTask.rigidSurface)[1] = grid->ISurface[0];            // поверхность заданная интерполянтом
            (*task.mechTask.rigidSurface)[2] = &(grid->FESurface[3]);        // поверхность заданная сеткой (если contact_mode != 3 то эта дополнительная сетка не строится)



            // 1)цилиндр задан аналитически
            //AnaliticalSurface_Cylinder *cylinder = new AnaliticalSurface_Cylinder;
            //cylinder->initCylinderZ({0, -1 - 5, 0}, 5);
            //(*task.mechTask.rigidSurface).push_back(cylinder);

            // 2) цилиндр задан сеткой
            //grid->FEsurface[3].init(true);  // неподвижна - обноврять декомпозицию пространства не нужно
            //(*task.mechTask.rigidSurface).push_back(&grid->FESurface[3]);


            // контакты КЭ-поверхность - жёсткая поверхность
            task.mechTask.CC_FE_Rigid = new std::vector<ContactCondition_FE_Rigid>;
            if(gp.surfaceType != SurfaceType::none)
            {
                ContactCondition_FE_Rigid CC_FE_Rigid_el;
                CC_FE_Rigid_el.stiffness = 1.e11;
                CC_FE_Rigid_el.FEsurfaceInd = 2;
                if(gp.surfaceType == SurfaceType::AnaliticalSurface_Cylinder)
                    CC_FE_Rigid_el.RigidSurfaceInd = 0;
                if(gp.surfaceType == SurfaceType::InterpolantSurface_Hermite3 ||
                   gp.surfaceType == SurfaceType::InterpolantSurface_Lagrange3)
                    CC_FE_Rigid_el.RigidSurfaceInd = 1;
                if(gp.surfaceType == SurfaceType::FiniteElementSurface)
                    CC_FE_Rigid_el.RigidSurfaceInd = 2;
                CC_FE_Rigid_el.noContactRadiusOptimization = noContactRadiusOptimization;
                (*task.mechTask.CC_FE_Rigid).push_back(CC_FE_Rigid_el);
            }

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
            // материал
            task.thermTask.material = new std::vector<ThermMaterialSource>(1);
            ThermMaterialSource &thermMat = (*task.thermTask.material)[0];
            thermMat.ro = 2200;//1;//2200;     // плотность
            thermMat.c = 1;//5;//1000;      // удельная теплоёмкость
            thermMat.f = 0;        // мощность внутренних объёмных источников (стоков) тепла
            double L0 = 10;//0.7;//10;
            thermMat.L =         // тензор теплопроводности
            {{
                {L0, 0, 0},
                {0, L0, 0},
                {0, 0, L0},
            }};
            int timeMode = 1;       // 0 - статика, 1 - динамика
            setBcFormovkaT(task.thermTask.bc1SourceT, task.thermTask.bc2SourceT,
                           10000000*98 // плотность набегающего теплового потока
                           );
            // шаги
            std::vector<ThermGlobalStep> *thermStep = new std::vector<ThermGlobalStep>;
            ThermGlobalStep ts;         // шаг для теплопроводности
            ts.slausolverParameters = task.slausolver_parameters;
            ts.timeMode = timeMode;     // статика/динамика
            thermStep->push_back(ts);
            for(int i = 0; i < StepsNumber2; i++)
                thermStep->push_back(ts);
            for(int i = 0; i < StepsNumber3; i++)
                thermStep->push_back(ts);
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
                v_el.T = T0;
                v_el.newT = T0;
                task.thermTask.vertex->push_back(v_el);
            }

            // ## 2-е краевые нужно переделать
        }
    }

    // Вдавливание цилиндра
    void genTestContactCylinder(TaskTest &task)
    {
        // индекс теста
        task.testIndex = 6;
        // подрубка
        task.mechTask.enabled = true;
        task.thermTask.enabled = false;
        // общие шаги и сетка
        std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
        Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
        // сетка
        ContactCylinderParameters &ccp = task.contactCylinder;
        // способ представления поверхности
        //gp.surfaceType = SurfaceType::none;
        ccp.surfaceType = SurfaceType::AnaliticalSurface_Cylinder;
        //gp.surfaceType = SurfaceType::InterpolantSurface_Hermite3;
        //gp.surfaceType = SurfaceType::InterpolantSurface_Lagrange3;
            //gp.contact_mode = SurfaceType::FiniteElementSurface;//нет поиска пересечения

        // способ декомпозиции поверхности
        ccp.decompositionType = DecompositionType::none;
        //gp.decompositionType = DecompositionType::cubeVertexIndexation;
        //gp.decompositionType = DecompositionType::surfaceAsSpheres;

        // подрубка универсальной оптимизации
        bool noContactRadiusOptimization = false;
        //bool noContactRadiusOptimization = true;

        double Time = 100;

        int NN = 20; //9//5
        int N0 = NN*2;
        int N1 = NN/2;
        int N2 = NN;
        int N3 = NN;
        ccp.Nz = 1;

        ccp.N[0][0] = N0;
        ccp.N[0][1] = N1;
        ccp.N[0][2] = ccp.Nz;

        ccp.N[1][0] = N0;
        ccp.N[1][1] = N2;
        ccp.N[1][2] = ccp.Nz;

        ccp.N[2][0] = N2;
        ccp.N[2][1] = N3;
        ccp.N[2][2] = ccp.Nz;

        ccp.N[3][0] = N2;
        ccp.N[3][1] = N3;
        ccp.N[3][2] = ccp.Nz;

        ccp.N[4][0] = N0;
        ccp.N[4][1] = N3;
        ccp.N[4][2] = ccp.Nz;

        POINT2 p[10];

        p[0] = POINT2(-10, -20);
        p[1] = POINT2(+10, -20);
        p[2] = POINT2(-10, -10);
        p[3] = POINT2(+10, -10);
        p[4] = POINT2(-3, -3);
        p[5] = POINT2(+3, -3);
        p[6] = POINT2(-10, 0);
        p[7] = POINT2(-3, 0);
        p[8] = POINT2(+3, 0);
        p[9] = POINT2(+10, 0);

        ccp.a[0][0] = p[0];
        ccp.a[0][1] = p[1];
        ccp.a[0][2] = p[2];
        ccp.a[0][3] = p[3];

        ccp.a[1][0] = p[2];
        ccp.a[1][1] = p[3];
        ccp.a[1][2] = p[4];
        ccp.a[1][3] = p[5];

        ccp.a[2][0] = p[2];
        ccp.a[2][1] = p[4];
        ccp.a[2][2] = p[6];
        ccp.a[2][3] = p[7];

        ccp.a[3][0] = p[5];
        ccp.a[3][1] = p[3];
        ccp.a[3][2] = p[8];
        ccp.a[3][3] = p[9];

        ccp.a[4][0] = p[4];
        ccp.a[4][1] = p[5];
        ccp.a[4][2] = p[7];
        ccp.a[4][3] = p[8];

        ccp.z0 = -0.5;
        ccp.z1 = +0.5;

        ccp.R = 5;
        ccp.C = {0, 5, 0};
        ccp.d = +0.5;
        ccp.V = -ccp.d / Time;
        ccp.stiffness = 1.e14;//1.e13
        ccp.E = 1.e10;
        ccp.Nu = 0.3;



        /*
        int NN = 5; //9//5
        ccp.p1 = VECTOR3(-10, -1, -1.0);
        ccp.p2 = VECTOR3(10, -0.5, -0.5);
        //gp.p1 = VECTOR3(-10, -1, -1);
        //gp.p2 = VECTOR3(10, -0.5, -0.5);
        ccp.N = VECTOR3_uint(NN*40, NN, 4);
        ccp.k = 1. - 1./64;
        */

        grid->genContactCylinder(ccp);
        //grid->buldOpenSCADModel("setka.scad");
        // параметры
        int fixGrid = 0;                // 0 - подвижная сетка, 1 - зафиксировать сетку
        int plasticPlot = 0;            // 0 - упругость, 1 - сплайн, 2 - Безье
        int plasticityCurveMode = 1;    // 0 - eps(sigma), 1 - sigma(eps)
        bool usingSecantMethod = false;
        double epsResidualLimit = 1.e-10;            // желаемая относительная погрешность по эквиволентным деформациям
        double contactDeltaFResidualLimit = 1.e-7;  // желаемая относительная погрешность по силе реакции опоры
        int nonlinearIterLimit = 200;                // ограничение на количество итераций
        int stepsNumber = 10;//64*5;//40/10;                        // количество шагов
        double elasticSigmaLimit = 3.e8;
        int sigma0Mode = 2;                         // 1 - dSigma = integral(...), 2 - в приращениях (P->dP)
        //int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny

        // шаги
        GlobalStep s;       // шаг
        // вдавливание
        for(int i = 0; i < stepsNumber; i++)
        {
            s.t1 = Time*0 + Time*i/stepsNumber;
            s.t2 = Time*0 + Time*(i + 1)/stepsNumber;
            s.dt0 = Time/stepsNumber;
            step->push_back(s);
        }
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
            // материал
            task.mechTask.material = new std::vector<MechMaterialSource>(1);
            MechMaterialSource &mechMat = (*task.mechTask.material)[0];
            mechMat.elasticSigmaLimit = elasticSigmaLimit;
            mechMat.Ceps = 4./9.;
            mechMat.Csigma = 1;
            mechMat.set_E_NU(ccp.E, ccp.Nu);
            mechMat.set_D_isotropic();      // обычная матрица D
            mechMat.set_M_sigma();
            mechMat.F = {0,0,0};// объёмные силы
            mechMat.ro = 0;     // плотность
                mechMat.Talpha = 0*1.e-5; // Коэффициент линейного расширения
            setPlasticMaterialCurve(plasticPlot, plasticityCurveMode, usingSecantMethod, mechMat);
            // шаги (на каждом шаге задаются вторые краевые условия)
            std::vector<MechGlobalStep> *mechStep = new std::vector<MechGlobalStep>;
            MechGlobalStep ms;  // шаг для МДТТ
            GlobalStep s;       // шаг
            std::vector<MechBoundaryCondition2Source_base *> *bc2Source;
            ms.slausolverParameters = task.slausolver_parameters;
            ms.timeMode = 0;    // 0 - квазистатическая задача
            ms.sigma0Mode = sigma0Mode; // 1 - dSigma = integral(...), 2 - приращения
            ms.fixGrid = fixGrid;     // 1 - зафиксировать сетку
            ms.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
            ms.controlMode = 0;     // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
            ms.iterLimit = nonlinearIterLimit;
            ms.slauResidualLimit = 10000;
            ms.epsResidualLimit = epsResidualLimit;
            ms.sigmaResidualLimit = 10000;
            ms.contactEndPointResidualLimit = 1000;//1.e-10;//1.e-13;
            ms.contactDeltaFResidualLimit = contactDeltaFResidualLimit;//1.e-13;
            // вдавливание цилиндра
            s = (*step)[0];
            for(int i = 0; i < stepsNumber; i++)
            {
                s.t1 = Time*0 + Time*i/stepsNumber;
                s.t2 = Time*0 + Time*(i + 1)/stepsNumber;
                setBc2ContactCylinder(bc2Source);
                ms.bc2Source = bc2Source;
                mechStep->push_back(ms);
            }
            task.mechTask.mechStep = mechStep;
            // первые краевые условия
            std::vector<MechBoundaryCondition1Source> *mechBc1 = new std::vector<MechBoundaryCondition1Source>(4);
            (*mechBc1)[0].mode = {{ 0, -1, -1}};
            (*mechBc1)[0].u0 =   {{ 0, -1, -1}};
            (*mechBc1)[1].mode = {{-1,  0, -1}};
            (*mechBc1)[1].u0 =   {{-1,  0, -1}};
            (*mechBc1)[2].mode = {{-1, -1,  0}};
            (*mechBc1)[2].u0 =   {{-1, -1,  0}};
            (*mechBc1)[3].mode = {{-1, -1,  -1}};
            (*mechBc1)[3].u0 =   {{-1, -1,  -1}};
            task.mechTask.bc1Source = mechBc1;
            // инициализация начальных деформаций и напряжений (=0)
            task.mechTask.fe = new std::vector<MechFeData_base *>;
            task.mechTask.fe->resize(grid->fe.size());
            for (size_t i = 0; i < grid->fe.size(); i++)
            {
                MechFeData_LinearHexagon_homogeny_E_NLE_EP *feDataEl = new MechFeData_LinearHexagon_homogeny_E_NLE_EP;
                feDataEl->init(Integration::IntegrationType::Gauss3);
                (*task.mechTask.fe)[i] = feDataEl;
            }
            // инициализация начальных перемещений (=0)
            task.mechTask.vertex = new std::vector<MechVertexData>;
            for (size_t i = 0; i < grid->vertex.size(); i++)
            {
                MechVertexData v_el;
                for (int j = 0; j < 3; j++)
                {
                    v_el.sum_du[j] = 0;
                    v_el.du[j] = 0; // не имеет значения
                }
                task.mechTask.vertex->push_back(v_el);
            }
            // для криволинейных границ ###
            task.mechTask.vertexForCurvature = new std::vector<MechVertexData>;
            for (size_t i = 0; i < grid->vertexForCurvature.size(); i++)
            {
                MechVertexData v_el;
                for (int j = 0; j < 3; j++)
                {
                    v_el.sum_du[j] = 0;
                    v_el.du[j] = 0; // не имеет значения
                }
                task.mechTask.vertexForCurvature->push_back(v_el);
            }
            // инициализация начальных скоростей (=0) и ускорений (=0)
            task.mechTask.V0 = new Vector(grid->vertex.size()*3);
            task.mechTask.dV0 = new Vector(grid->vertex.size()*3);
            for (size_t i = 0; i < task.mechTask.V0->size(); i++)
            {
                (*task.mechTask.V0)[i] = 0;
                (*task.mechTask.dV0)[i] = 0;
            }

            // вторые краевые
            task.mechTask.bc2 = new std::vector<MechBoundaryCondition2>;
            task.mechTask.bc2->clear();

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
                CC_FE_Rigid_el.stiffness = ccp.stiffness;
                CC_FE_Rigid_el.FEsurfaceInd = 0;
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
        }

        // Температура
        if(task.thermTask.enabled)
        {
        }

    }

    // Ползучесть
    void genTestCreep(TaskTest &task)
    {
        // индекс теста
        task.testIndex = 7;
        // подрубка
        task.mechTask.enabled = true;
        task.thermTask.enabled = false;
        // общие шаги и сетка
        std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
        Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
        // сетка
        CreepParameters &cp = task.creep;

        // X
        cp.cube.i[0][0] = -10;
        cp.cube.i[0][1] = +10;
        // Y
        cp.cube.i[1][0] = -10;
        cp.cube.i[1][1] = +10;
        // Z
        cp.cube.i[2][0] = -1;
        cp.cube.i[2][1] = +1;

        cp.N[0] = 20;
        cp.N[1] = 20;
        cp.N[2] = 2;

        cp.E = 200000;//1.e5;//200*1.e3*100*100;
        cp.Nu = 0.3;

        cp.sigma0 = 200;    //double P = 100;//200*100*100 * 1.e0;     // давление
        cp.alpha = -0.1;      //double alpha = -1;
        cp.time = 1000/100;
        cp.A = 3.125*1.e-14;
        cp.n = 5;
        cp.m = 0.5;

        grid->genCreep(cp);
        //grid->buldOpenSCADModel("setka.scad");
        // параметры
        int fixGrid = 0;                // 0 - подвижная сетка, 1 - зафиксировать сетку
        int plasticPlot = 0;            // 0 - упругость, 1 - сплайн, 2 - Безье
        int plasticityCurveMode = 1;    // 0 - eps(sigma), 1 - sigma(eps)
        bool usingSecantMethod = false;
        double epsResidualLimit = 1.e-10;            // желаемая относительная погрешность по эквиволентным деформациям
        double creepResidualLimit = 1.e+10;           // желаемая относительная погрешность по деформациям ползучести
        int nonlinearIterLimit = 20;                // ограничение на количество итераций
        int stepsNumber0 = 1;                        // пустой шаг
        int stepsNumber1 = 1;                        // количество шагов нагружения
        int stepsNumber2 = 10;                       // количество шагов ползучести
        double Time0 = 0.0000001;           // время пустого шага
        double Time1 = 0.0000001;           // время нагружения
        double Time2 = cp.time;               // время ползучести
        double Time_start = -Time0 - Time1;


        double Px1 = 0;
        double Px2 = -cp.sigma0;
        double Py1 = 0;
        double Py2 = -cp.alpha*cp.sigma0;
        double elasticSigmaLimit = 3.e8;
        int sigma0Mode = 2;                         // 1 - dSigma = integral(...), 2 - в приращениях (P->dP)

        // шаги
        GlobalStep s;       // шаг
        // пустой шаг
        for(int i = 0; i < stepsNumber0; i++)
        {
            s.t1 = Time_start + Time0*i/stepsNumber0;
            s.t2 = Time_start + Time0*(i + 1)/stepsNumber0;
            s.dt0 = Time0/stepsNumber0;
            step->push_back(s);
        }
        // нагружение
        for(int i = 0; i < stepsNumber1; i++)
        {
            s.t1 = Time_start + Time0 + Time1*i/stepsNumber1;
            s.t2 = Time_start + Time0 + Time1*(i + 1)/stepsNumber1;
            s.dt0 = Time1/stepsNumber1;
            step->push_back(s);
        }
        // ползучесть
        for(int i = 0; i < stepsNumber2; i++)
        {
            s.t1 = Time_start + Time0 + Time1 + Time2*i/stepsNumber2;
            s.t2 = Time_start + Time0 + Time1 + Time2*(i + 1)/stepsNumber2;
            s.dt0 = Time2/stepsNumber2;
            step->push_back(s);
        }
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
            // материал
            task.mechTask.material = new std::vector<MechMaterialSource>(1);
            MechMaterialSource &mechMat = (*task.mechTask.material)[0];
            mechMat.elasticSigmaLimit = elasticSigmaLimit;
            mechMat.Ceps = 4./9.;
            mechMat.Csigma = 1;
            mechMat.set_E_NU(cp.E, cp.Nu);
            mechMat.set_D_isotropic();      // обычная матрица D
            mechMat.set_M_sigma();
            mechMat.F = {0,0,0};// объёмные силы
            mechMat.ro = 0;     // плотность
                mechMat.Talpha = 0*1.e-5; // Коэффициент линейного расширения
            setPlasticMaterialCurve(plasticPlot, plasticityCurveMode, usingSecantMethod, mechMat);
            setCreepMaterialCurve(cp.A, cp.n, cp.m, mechMat);
            // шаги (на каждом шаге задаются вторые краевые условия)
            std::vector<MechGlobalStep> *mechStep = new std::vector<MechGlobalStep>;
            MechGlobalStep ms;  // шаг для МДТТ
            GlobalStep s;       // шаг
            std::vector<MechBoundaryCondition2Source_base *> *bc2Source;
            ms.slausolverParameters = task.slausolver_parameters;
            ms.timeMode = 0;    // 0 - квазистатическая задача
            ms.sigma0Mode = sigma0Mode; // 1 - dSigma = integral(...), 2 - приращения
            ms.fixGrid = fixGrid;     // 1 - зафиксировать сетку
            ms.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
            ms.controlMode = 0;     // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
            ms.iterLimit = nonlinearIterLimit;
            ms.slauResidualLimit = 10000;
            ms.epsResidualLimit = epsResidualLimit;
            ms.creepResidualLimit = creepResidualLimit;
            ms.sigmaResidualLimit = 10000;
            ms.contactEndPointResidualLimit = 1000;//1.e-10;//1.e-13;
            // пустой шаг
            for(int i = 0; i < stepsNumber0; i++)
            {
                s = (*step)[0 + i];
                setBc2Creep(0, 0,
                            0, 0,
                            s.t1, s.t2, bc2Source);
                ms.bc2Source = bc2Source;
                mechStep->push_back(ms);
            }
            // нагружение
            for(int i = 0; i < stepsNumber1; i++)
            {
                s = (*step)[stepsNumber0 + i];
                //s.t1 = Time*0 + Time*i/stepsNumber;
                //  s.t2 = Time*0 + Time*(i + 1)/stepsNumber;
                setBc2Creep(Px1 + (Px2 - Px1)*i/stepsNumber1, Px1 + (Px2 - Px1)*(i + 1)/stepsNumber1,
                            Py1 + (Py2 - Py1)*i/stepsNumber1, Py1 + (Py2 - Py1)*(i + 1)/stepsNumber1,
                            s.t1, s.t2, bc2Source);
                ms.bc2Source = bc2Source;
                mechStep->push_back(ms);
            }
            // ползучесть
            for(int i = 0; i < stepsNumber2; i++)
            {
                s = (*step)[stepsNumber0 + stepsNumber1 + i];
                setBc2Creep(Px2, Px2,
                            Py2, Py2,
                            s.t1, s.t2, bc2Source);
                ms.bc2Source = bc2Source;
                mechStep->push_back(ms);
            }
            task.mechTask.mechStep = mechStep;
            // первые краевые условия
            std::vector<MechBoundaryCondition1Source> *mechBc1 = new std::vector<MechBoundaryCondition1Source>(4);
            (*mechBc1)[0].mode = {{ 0, -1, -1}};
            (*mechBc1)[0].u0 =   {{ 0, -1, -1}};
            (*mechBc1)[1].mode = {{-1,  0, -1}};
            (*mechBc1)[1].u0 =   {{-1,  0, -1}};
            (*mechBc1)[2].mode = {{-1, -1,  0}};
            (*mechBc1)[2].u0 =   {{-1, -1,  0}};
            (*mechBc1)[3].mode = {{-1, -1, -1}};
            (*mechBc1)[3].u0 =   {{-1, -1, -1}};
            task.mechTask.bc1Source = mechBc1;
            // инициализация начальных деформаций и напряжений (=0)
            task.mechTask.fe = new std::vector<MechFeData_base *>;
            task.mechTask.fe->resize(grid->fe.size());
            for (size_t i = 0; i < grid->fe.size(); i++)
            {
                MechFeData_LinearHexagon_homogeny_E_NLE_EP *feDataEl = new MechFeData_LinearHexagon_homogeny_E_NLE_EP;
                feDataEl->init(Integration::IntegrationType::Gauss3);
                (*task.mechTask.fe)[i] = feDataEl;
            }
            // инициализация начальных перемещений (=0)
            task.mechTask.vertex = new std::vector<MechVertexData>;
            for (size_t i = 0; i < grid->vertex.size(); i++)
            {
                MechVertexData v_el;
                for (int j = 0; j < 3; j++)
                {
                    v_el.sum_du[j] = 0;
                    v_el.du[j] = 0; // не имеет значения
                }
                task.mechTask.vertex->push_back(v_el);
            }
            // для криволинейных границ ###
            task.mechTask.vertexForCurvature = new std::vector<MechVertexData>;
            for (size_t i = 0; i < grid->vertexForCurvature.size(); i++)
            {
                MechVertexData v_el;
                for (int j = 0; j < 3; j++)
                {
                    v_el.sum_du[j] = 0;
                    v_el.du[j] = 0; // не имеет значения
                }
                task.mechTask.vertexForCurvature->push_back(v_el);
            }
            // инициализация начальных скоростей (=0) и ускорений (=0)
            task.mechTask.V0 = new Vector(grid->vertex.size()*3);
            task.mechTask.dV0 = new Vector(grid->vertex.size()*3);
            for (size_t i = 0; i < task.mechTask.V0->size(); i++)
            {
                (*task.mechTask.V0)[i] = 0;
                (*task.mechTask.dV0)[i] = 0;
            }

            // вторые краевые
            task.mechTask.bc2 = new std::vector<MechBoundaryCondition2>;
            MechBoundaryCondition2 bc2_el;
            bc2_el.FEsurfaceInd = 0;
            bc2_el.bc2SourceIndex = 0;
            (*task.mechTask.bc2).push_back(bc2_el);
            bc2_el.FEsurfaceInd = 1;
            bc2_el.bc2SourceIndex = 1;
            (*task.mechTask.bc2).push_back(bc2_el);

            // Контакт

            // поверхности
            task.mechTask.rigidSurface = new std::vector<Surface_base *>;
            task.mechTask.rigidSurface->clear();

            // контакты КЭ-поверхность - жёсткая поверхность
            task.mechTask.CC_FE_Rigid = new std::vector<ContactCondition_FE_Rigid>;
            task.mechTask.CC_FE_Rigid->clear();
        }

        // Температура
        if(task.thermTask.enabled)
        {
        }

    }

    // материал анизотропного пластичного тела для мдтт
    struct Ui_Solid_Material
    {
        std::string name = "";       // название материала
        std::string K = "";
        std::string G = "";
        std::string Csigma = "";
        std::string Ceps = "";
        std::string eps_sigma_expression = "";
        std::string deps_dsigma_expression = "";
        std::vector<std::string> args;
        std::vector<std::string> vals;
        MATR6x6 M_sigma;
        VECTOR3 F;
        Solid::MechMaterialSource *gen()
        {
            using namespace Solid;
            MechMaterialSource *m = new MechMaterialSource;
            m->set_E_NU(1.e10,0.3);
            //m->set_K_G(strToDouble(K), strToDouble(G));
            //m->K = strToDouble(K);
            //m->G = strToDouble(G);
            m->Csigma = FunParser::strToDouble(Csigma);
            m->Ceps = FunParser::strToDouble(Ceps);
            // строковые представления функций
            m->epsFun.setExpression(eps_sigma_expression);
            m->difEpsFun.setExpression(deps_dsigma_expression);
            // аргументы
            m->epsFun.args.clear();
            m->difEpsFun.args.clear();
            m->epsFun.addArgument("x");
            m->difEpsFun.addArgument("x");
            m->epsFun.addArgument("K", m->K);
            m->difEpsFun.addArgument("K", m->K);
            m->epsFun.addArgument("G", m->G);
            m->difEpsFun.addArgument("G", m->G);
            for(size_t i = 0; i < args.size(); i++)
            {
                m->epsFun.addArgument(args[i], FunParser::strToDouble(vals[i]));
                m->difEpsFun.addArgument(args[i], FunParser::strToDouble(vals[i]));
            }
            // разбор
            m->epsFun.parse();
            m->difEpsFun.parse();
            m->epsFun.setName("m");


            m->set_K_G(m->K, m->G);
            m->set_D_isotropic();
            m->M_sigma = M_sigma;
            m->F = F;
            return m;
        }
    };

    // функция для построения графиков
    struct Ui_Solid_OutFunction
    {
        std::string name;             // название графика
        std::string name_y;           // название функции
        std::string name_x;           // название аргумента
        std::string expression_y;     // функция
        std::string expression_x;     // аргумент
    };
    //std::vector<Ui_Solid_Material> m;               // материалы
    //std::vector<Ui_Solid_OutFunction> outFun;       // функции для построения графиков
    /*
    string tstr;
    size_t size;
    // материалы
    ifstream in;
    in.open("_materials.txt");
    in >> size;
    m.resize(size);
    getline(in, tstr);
    for(size_t k = 0; k < m.size(); k++)
    {
        Ui_Solid_Material &m0 = m[k];
        getline(in, m0.name);
        getline(in, m0.Csigma);
        getline(in, m0.Ceps);
        for(int i = 0; i < 6; i++)
            for(int j = 0; j < 6; j++)
                in >> m0.M_sigma.m[i][j];
        for(int i = 0; i < 3; i++)
            in >> m0.F.x[i];
        getline(in, tstr);
        getline(in, m0.K);
        getline(in, m0.G);
        getline(in, m0.eps_sigma_expression);
        getline(in, m0.deps_dsigma_expression);
        in >> size;
        m0.args.resize(size);
        m0.vals.resize(size);
        getline(in, tstr);
        for(size_t i = 0; i < m0.args.size(); i++)
        {
            getline(in, m0.args[i]);
            getline(in, m0.vals[i]);
        }
    }
    in.close();
    // функции для построения графиков
    in.open("_functions.txt");
    in >> size;
    outFun.resize(size);
    getline(in, tstr);
    for(size_t k = 0; k < outFun.size(); k++)
    {
        getline(in, outFun[k].name);
        getline(in, outFun[k].name_y);
        getline(in, outFun[k].name_x);
        getline(in, outFun[k].expression_y);
        getline(in, outFun[k].expression_x);
    }
    in.close();
*/

    /*
    // материалы
    ofstream fout("_materials.txt");
    fout << m.size() << endl;
    for(size_t k = 0; k < m.size(); k++)
    {
        Ui_Solid_Material &m0 = m[k];
        fout << m0.name << endl;
        fout << m0.Csigma << endl;
        fout << m0.Ceps << endl;
        for(int i = 0; i < 6; i++)
        {
            for(int j = 0; j < 6; j++)
                fout << m0.M_sigma.m[i][j] << "\t";
            fout << endl;
        }
        for(int i = 0; i < 3; i++)
            fout << m0.F.x[i] << "\t";
        fout << endl;
        fout << m0.K << endl;
        fout << m0.G << endl;
        fout << m0.eps_sigma_expression << endl;
        fout << m0.deps_dsigma_expression << endl;
        fout << m0.args.size() << endl;
        for(size_t i = 0; i < m0.args.size(); i++)
        {
            fout << m0.args[i] << endl;
            fout << m0.vals[i] << endl;
        }
    }
    fout.close();*/
    /*
    // функции для построения графиков
    fout.open("_functions.txt");
    fout << outFun.size() << endl;
    for(size_t k = 0; k < outFun.size(); k++)
    {
        fout << outFun[k].name << endl;
        fout << outFun[k].name_y << endl;
        fout << outFun[k].name_x << endl;
        fout << outFun[k].expression_y << endl;
        fout << outFun[k].expression_x << endl;
    }
    fout.close();*/

#if false
    {
        FILE *fx = fopen("fx", "w");
        FILE *fy = fopen("fy", "w");
        FILE *fz = fopen("fz", "w");
        FILE *fpr = fopen("fpr", "w");
        GridRectangleRegular2D it_grid;
        double R = 2;
        int N[2] = {100, 100};
        int color = 3;
        //int NN[2] = {N[0]*3, N[1]*3};   // без регуляризации
        double CoefAlpha = 0.00;
        //int NN[2] = {N[0], N[1]};     // с регуляризацией
        //double CoefAlpha = 0.001;
        it_grid.init(0, PI, 0, 5., N[0], N[1]);
        //Interpolation::interpolant2D_Hermite3 itx, ity, itz;
        Interpolation::Interpolant2D_Lagrange3 itx, ity, itz;
        itx.init(it_grid, CoefAlpha);
        ity.init(it_grid, CoefAlpha);
        itz.init(it_grid, CoefAlpha);

        // 1) добавление точек
/*
        for(int a = 0; a <= NN[0]; a++)    // угол (от 0 до 2*PI, N[0] штук)
        {
            for(int z = 0; z <= NN[1]; z++)    // координата z (от 0 до 5, N[1] штук)
            {
                double alpha = 0.999999*PI*a/NN[0];
                double zz = 4.999999*z/NN[1];
                itx.addPoint({alpha, zz}, R*cos(alpha));
                ity.addPoint({alpha, zz}, R*sin(alpha));
                itz.addPoint({alpha, zz}, zz);
            }
        }
        itx.buildInterpolant();
        ity.buildInterpolant();
        itz.buildInterpolant();
        //for(;;);

        itx.save("x.bmp", 800, 800, color, 1);
        ity.save("y.bmp", 800, 800, color, 1);
        itz.save("z.bmp", 800, 800, color, 1);

        for(int a = 0; a <= NN[0]; a++)    // угол (от 0 до 2*PI, N[0] штук)
        {
            for(int z = 0; z <= NN[1]; z++)    // координата z (от 0 до 5, N[1] штук)
            {
                double alpha = 0.999999*PI*a/NN[0];
                double zz = 4.999999*z/NN[1];
                fprintf(fx, "alpha = %lf, zz = %lf, Fx = %lf, ip = %lf\n", alpha, zz, R*cos(alpha), itx.fun({alpha, zz}));
                fprintf(fy, "alpha = %lf, zz = %lf, Fy = %lf, ip = %lf\n", alpha, zz, R*sin(alpha), ity.fun({alpha, zz}));
                fprintf(fz, "alpha = %lf, zz = %lf, Fz = %lf, ip = %lf\n", alpha, zz, zz, itz.fun({alpha, zz}));
            }
        }
*/

        // 2) определение ф-и в наборе точек

        std::vector<POINT2> coordinates;
        itx.getNodesCoordinates(coordinates);
        Vector nodeValuex, nodeValuey, nodeValuez;
        nodeValuex.resize(coordinates.size());
        nodeValuey.resize(coordinates.size());
        nodeValuez.resize(coordinates.size());
        for(int i = 0; i < coordinates.size(); i++)
        {
            double alpha = coordinates[i][0];
            double zz = coordinates[i][1];
            nodeValuex[i] = R*cos(alpha);
            nodeValuey[i] = R*sin(alpha);
            nodeValuez[i] = zz;
            //itx.addPoint(coordinates[i], nodeValuex[i]);
            //ity.addPoint(coordinates[i], nodeValuey[i]);
            //itz.addPoint(coordinates[i], nodeValuez[i]);

        }
        itx.buildInterpolantByAllNodes(nodeValuex);
        ity.buildInterpolantByAllNodes(nodeValuey);
        itz.buildInterpolantByAllNodes(nodeValuez);

        itx.save("x.bmp", 400, 400, color, 1);
        ity.save("y.bmp", 400, 400, color, 1);
        itz.save("z.bmp", 400, 400, color, 1);


        for(int i = 0; i < coordinates.size(); i++)
        {
            double alpha = coordinates[i][0];
            double zz = coordinates[i][1];
            fprintf(fx, "alpha = %lf, zz = %lf, Fx = %lf, ip = %lf\n", alpha, zz, nodeValuex[i], itx.fun({alpha, zz}));
            fprintf(fy, "alpha = %lf, zz = %lf, Fy = %lf, ip = %lf\n", alpha, zz, nodeValuey[i], ity.fun({alpha, zz}));
            fprintf(fz, "alpha = %lf, zz = %lf, Fz = %lf, ip = %lf\n", alpha, zz, nodeValuez[i], itz.fun({alpha, zz}));
        }


        Interpolation::InterpolantSurface is;
        is.init(&itx, &ity, &itz);
        // проекции на поверхность
        for(int i = -1; i <= 101; i++)
        {
            double alpha = PI*i/100;
            double zz = 1;  // 0..5
            POINT3 p0;
            //p0 = {0, R + 1, 0};
            p0[0] = (R + 1)*cos(alpha);
            p0[1] = (R + 1)*sin(alpha);
            p0[2] = zz;
            POINT3 nearestPointAn;
            nearestPointAn[0] = (R)*cos(alpha);
            nearestPointAn[1] = (R)*sin(alpha);
            nearestPointAn[2] = zz;
            POINT3 normalAn = (p0 - nearestPointAn) / (p0 - nearestPointAn).abs();

            POINT3 nearestPoint;
            POINT3 normal;
            int side;
            bool onBorder;
            is.findNearestPoint(p0, {0,0,0},
                                nearestPoint, normal, side, onBorder);
            fprintf(fpr, "p0 = (%lf, %lf, %lf):\n\tnearestPointAn = (%lf, %lf, %lf)\n\tnearestPoint   = (%lf, %lf, %lf)\n\tnormalAn = (%lf, %lf, %lf)\n\tnormal   = (%lf, %lf, %lf)\n\tside = %d, onBorder = %d\n",
                    p0[0], p0[1], p0[2],
                    nearestPointAn[0], nearestPointAn[1], nearestPointAn[2],
                    nearestPoint[0], nearestPoint[1], nearestPoint[2],
                    normalAn[0], normalAn[1], normalAn[2],
                    normal[0], normal[1], normal[2],
                    side, (int)onBorder);
        }


        fclose(fx);
        fclose(fy);
        fclose(fz);
        fclose(fpr);
    }
#endif
#if false
void myWidget_out::paintEvent000(QPaintEvent *)
{
    using namespace Elementary;
    TaskTest &task = *outTask;
    OutData &out = *outData;
    QPainter p(this);

    // поиск области вывода в координатах сетки
    POINT3 pMin = {1.e100, 1.e100, 1.e100};
    POINT3 pMax = {-1.e100, -1.e100, -1.e100};
    for(int feInd = 0; feInd < task.grid->vertex.size(); feInd++)
    {
        POINT3 vertex_el = task.grid->vertex[feInd];
        // минимальные координаты узлов сетки
        if(vertex_el[0] < pMin[0]) pMin[0] = vertex_el[0];
        if(-vertex_el[1] < pMin[1]) pMin[1] = -vertex_el[1];
        if(vertex_el[2] < pMin[2]) pMin[2] = vertex_el[2];
        // максимальные координаты узлов сетки
        if(vertex_el[0] > pMax[0]) pMax[0] = vertex_el[0];
        if(-vertex_el[1] > pMax[1]) pMax[1] = -vertex_el[1];
        if(vertex_el[2] > pMax[2]) pMax[2] = vertex_el[2];
    }
    // настройка области вывода
    int w = 1280;
    int w_add = 100;
    int h_add = 100;
    POINT3 sizes1 = pMax - pMin;
    //POINT3 centre1 = (pMax + pMin) / 2;
    POINT3 sizes2;
    double scale = w/sizes1[0];
    sizes2[0] = (double)w;
    sizes2[1] = (double)w*sizes1[1]/sizes1[0];
    resize(sizes2[0] + w_add, sizes2[1] + h_add);
    p.translate(sizes2[0] * (0 - pMin[0])/sizes1[0] + w_add, sizes2[1] * (0 - pMin[1])/sizes1[1] + h_add);
/*
    double scale = 100;
    int w_add = 100*0;
    int h_add = 100*0;
    POINT3 sizes1 = pMax - pMin;
    POINT3 sizes2;
    sizes2[0] = sizes1[0]*scale + w_add;
    sizes2[1] = sizes1[1]*scale + h_add;
    resize(sizes2[0], sizes2[1]);
    p.translate(width() / 2 + w_add, height() / 2 + h_add);
    p.scale(scale, scale);
*/
/*
    POINT3 sizes = pMax - pMin;
    p.translate(width() / 2, height() / 2);
    double w1 = sizes[0]*1.2;
    double w2 = width();
    double scale = w2 / w1;
    //p.scale(scale, scale);
*/
    //QRect viewport;
    //p.shear(pMax[0] - pMin[0], pMax[1] - pMin[1]);
    //p.translate(sizes[0] / 2, sizes[1] / 2);


    // поиск минимального и максимального значений отображаемого поля
    double valueMin = 1.e100;
    double valueMax = -1.e100;
    for(int feInd = 0; feInd < task.grid->fe.size(); feInd++)
    {
        double value = getValue(out, stepIndex, feInd, valueIndex);
        // минимальное значение
        if(value < valueMin) valueMin = value;
        // максимальное значение
        if(value > valueMax) valueMax = value;
    }

    // вывод четырёхугольников
    for(int feInd = 0; feInd < task.grid->fe.size(); feInd++)
    {
        int reindex[4] = {0, 1, 3, 2};
        QPointF points[4];
        //FiniteElement fe_el = task.grid->fe[feInd];
        for(int i = 0; i < 4; i++)
        {
            int vertexIndex = task.grid->globalVertexIndex(feInd, reindex[i]);
            POINT3 p = task.grid->vertex[vertexIndex];
            points[i].setX(p[0] * scale);
            points[i].setY(-p[1] * scale);
        }
        double value = getValue(out, stepIndex, feInd, valueIndex);
        int c = 255 - (int)((value - valueMin) / (valueMax - valueMin) * 255.);
        QColor color(c, c, c, 255);
        p.setBrush(QBrush(color, Qt::SolidPattern));
        p.drawPolygon(points, 4);
    }


    /*p.drawRect(-4, -4, 4, 4);
    p.setPen(QPen({100,100,100}));*/
    //p.drawText(0, 0, "Hello!");

}
#endif










#if false

// сфера
struct Test_sphere: public Test_base
{
    Grid::SphereParameters gp;
    int thermTestIndex;
    int sigmaSolvingType;   // 0 - главные напряжения, 1 - проекции площадки
    std::string fn_inf;
    std::string fn_T;
    std::string fn_T_ch;
    std::string fn_T_pogr_abs;
    std::string fn_T_pogr;
    std::string fn_u;
    std::string fn_u_ch;
    std::string fn_u_pogr_abs;
    std::string fn_u_pogr;
    std::string fn_sigma_r;
    std::string fn_sigma_r_ch;
    std::string fn_sigma_fi;
    std::string fn_sigma_fi_ch;
    std::string fn_sigma_r_0;
    std::string fn_sigma_r_ch_0;
    std::string fn_sigma_fi_0;
    std::string fn_sigma_fi_ch_0;
    std::string fn_sigma_r_pogr_abs;
    std::string fn_sigma_r_pogr_0_abs;
    std::string fn_sigma_fi_pogr_abs;
    std::string fn_sigma_fi_pogr_0_abs;
    std::string fn_sigma_r_pogr;
    std::string fn_sigma_r_pogr_0;
    std::string fn_sigma_fi_pogr;
    std::string fn_sigma_fi_pogr_0;
    std::string fn_res_maxpogr;
    std::string fn_s_iterationsNumber;
    std::string fn_s_slauResidual;
    std::string fn_s_stepSolvingTime;
    std::string fn_s_nonlinearStateFENumber;
    std::string fn_s_slauNev;
    std::string fn_s_epsNev;
    std::string fn_s_sigmaNev;
    Test_sphere();
    virtual Type get_type()const;
    virtual void initTask(Solid::Task &task);
    virtual void writeResults(const Solid::Task &task, const Solid::OutData &out);
};

Test_sphere::Test_sphere()
{
    // дирректория данных теста
    dir = "./tests/hollowSphere/";
    // пути к текстовым файлам с данными
    fn_inf = dir + "_inf.txt";
    fn_T = dir + "__res_T.txt";
    fn_T_ch = dir + "__res_T_ch.txt";
    fn_T_pogr_abs = dir + "__res_T_pogr_abs.txt";
    fn_T_pogr = dir + "__res_T_pogr.txt";
    fn_u = dir + "__res_u.txt";
    fn_u_ch = dir + "__res_u_ch.txt";
    fn_u_pogr_abs = dir + "__res_u_pogr_abs.txt";
    fn_u_pogr = dir + "__res_u_pogr.txt";

    fn_sigma_r = dir + "load_sigma_r.txt";
    fn_sigma_r_ch = dir + "load_sigma_r_ch.txt";
    fn_sigma_fi = dir + "load_sigma_fi.txt";
    fn_sigma_fi_ch = dir + "load_sigma_fi_ch.txt";

    fn_sigma_r_0 = dir + "unload_sigma_r.txt";
    fn_sigma_r_ch_0 = dir + "unload_sigma_r_ch.txt";
    fn_sigma_fi_0 = dir + "unload_sigma_fi.txt";
    fn_sigma_fi_ch_0 = dir + "unload_sigma_fi_ch.txt";

    fn_sigma_r_pogr_abs = dir + "load_sigma_r_pogr_abs.txt";
    fn_sigma_r_pogr_0_abs = dir + "unload_sigma_r_pogr_abs.txt";
    fn_sigma_fi_pogr_abs = dir + "load_sigma_fi_pogr_abs.txt";
    fn_sigma_fi_pogr_0_abs = dir + "unload_sigma_fi_pogr_abs.txt";

    fn_sigma_r_pogr = dir + "load_sigma_r_pogr.txt";
    fn_sigma_r_pogr_0 = dir + "unload_sigma_r_pogr.txt";
    fn_sigma_fi_pogr = dir + "load_sigma_fi_pogr.txt";
    fn_sigma_fi_pogr_0 = dir + "unload_sigma_fi_pogr.txt";

    fn_res_maxpogr = dir + "___maxpogr.txt";
    fn_s_iterationsNumber = dir + "__res_s_iterationsNumber.txt";
    fn_s_slauResidual = dir + "__res_s_slauResidual.txt";
    fn_s_stepSolvingTime = dir + "__res_s_stepSolvingTime.txt";
    fn_s_nonlinearStateFENumber = dir + "__res_s_nonlinearStateFENumber.txt";
    fn_s_slauNev = dir + "__res_s_slauNev.txt";
    fn_s_epsNev = dir + "__res_s_epsNev.txt";
    fn_s_sigmaNev = dir + "__res_s_sigmaNev.txt";

    // имена gnuplot файлов и соответствующих картинок
    resGraph.add_all_in_subdir(dir, "sigma_fi.gnu", "sigma_fi.png");
    resGraph.add_all_in_subdir(dir, "sigma_fi_pogr_load.gnu", "sigma_fi_pogr_load.png");
    resGraph.add_all_in_subdir(dir, "sigma_fi_pogr_unload.gnu", "sigma_fi_pogr_unload.png");
    resGraph.add_all_in_subdir(dir, "sigma_r.gnu", "sigma_r.png");
    resGraph.add_all_in_subdir(dir, "sigma_r_pogr_load.gnu", "sigma_r_pogr_load.png");
    resGraph.add_all_in_subdir(dir, "sigma_r_pogr_unload.gnu", "sigma_r_pogr_unload.png");
    // инициализация графиков с невязками
    initStepsGraphs();
}
Test_base::Type Test_sphere::get_type() const
{
    return Type::Sphere;
}
void Test_sphere::initTask(Solid::Task &task)
{
    FILE *f_inf = fopen(fn_inf.c_str(), "w");
    // индекс теста решателя теплопроводности
    thermTestIndex = 0;             // 1 - t(a)=t1, t(b)=t2
                                    // 2 - t(a)=t1, конвекция на b
                                    // 3 - конвекция на r=a, конвекция на r=b
                                    // 4 - подогрев на r=a, t(b)=t2
    //bool fidesys = true;
    bool fidesys = false;
    // подрубка
    if(thermTestIndex != 0)
    {
        task.thermTask.enabled = true;
        task.mechTask.enabled = false;
    }
    else
    {
        task.thermTask.enabled = false;
        task.mechTask.enabled = true;
    }
    // общие шаги и сетка
    std::vector<GlobalStep> *step = new std::vector<GlobalStep>;   // общие параметры глобальных шагов
    Grid::Grid3D *grid = new Grid::Grid3D;  //сетка
    // сетка
    if(fidesys)
    {
        gp.r1 = 2.5;
        gp.r2 = 5;
    }
    else
    {
        gp.r1 = 1;
        gp.r2 = 4;
    }
    gp.N = 8;   // 4 промежутка на четвертьокружностях    //8, 32
    gp.Nparts = 32;
    if(thermTestIndex != 0)
    {
        gp.q = 1;
    }
    else
    {
        if(fidesys)
            gp.q = 1. + 0*1./(gp.Nparts/2);//1.0625;//1.0625;//1.03125;
        else
            gp.q = 1. + 1*1./(gp.Nparts/2);//1.0625;//1.0625;//1.03125;
    }
    gp.curvilinear = 1;    // отображение (0 - линейное, 1 - квадратичное)
    gp.buildingMethod = 0; // способ построения (0 - делим дуги, 1 - отображение на куб)
    grid->genSphere(gp);
    //grid->buldOpenSCADModel("setka.scad");
    double NU = 0.3; //0.4

    // параметры
    //task.slausolver_parameters.preconditioning = SlauSolving::Preconditioning::SlauPreconditioning_LLT;

    sigmaSolvingType = 0;  // способ расчёта решения (0 - главные напряжения, 1 - проекции площадки)
    int fixGrid;                // 0 - подвижная сетка, 1 - зафиксировать сетку
    int sigma0Mode = 2;         // 1 - dSigma = integral(...), 2 - в приращениях (P->dP)
    if(fidesys)
    {
        fixGrid = 0;
    }
    else
    {
        fixGrid = 0;
    }
    int plasticPlot = 2;            // 0 - упругость, 1 - сплайн, 2 - Безье
    int plasticityCurveMode = 1;    // 0 - eps(sigma), 1 - sigma(eps)
    bool usingSecantMethod = false;//true;
    double epsResidualLimit = 1.e-13;    // желаемая невязка по эквиволентным деформациям
    int nonlinearIterLimit = 100;        // ограничение на количество итераций
    int StepsNumber1 = 1;//8;//8;//128/2;      // шаги
    int StepsNumber2 = 4;//300;//512;//40;//128/2;
    int StepsNumber3 = 1;
    int StepsNumber4 = 1;
    double P = 30;//2.0e7;//2.5e7;//1000;//2.5e7;
    if(fidesys)
        P = 30;
    else
        P = 2.0e7;
    if(thermTestIndex != 0)
        P = 0;
    double P1 = 0.1*P;//P/2;//0.4*P;
    double P2 = P;
    double P3 = P*0.9;//P/2.;//P-((P2-P1)/StepsNumber2)/100;//-P*1/2;//-(P2/StepsNumber);
    double P4 = 0;//-P*1/2;//-(P - P2/StepsNumber);
    /*
    double P1 = P*1/2;
    double P2 = P;
    double P3 = P*1/2;
    double P4 = 0;
    */
    int HomogenyMode = Solid_FeData_HomogenyMode_homogeny;  //Solid_FeData_HomogenyMode_not_homogeny

    double Time = 1;
    // шаги
    GlobalStep s;       // шаг
    // первый шаг (1 временной слой, заведомо упругое нагружение P1)
    s.t1 = 0;                          // (имеет значение только для первого глобального шага)
    s.t2 = Time*1;
    s.dt0 = Time/StepsNumber1;
    step->push_back(s);
    // второй шаг (пластичное нагружение P2)
    s.t1 = Time*1;
    s.t2 = Time*2;
    s.dt0 = Time/StepsNumber2;
    step->push_back(s);
    // третий шаг (1 временной слой, начало разгрузки)
    s.t1 = Time*2;
    s.t2 = Time*3;
    s.dt0 = Time/StepsNumber3;
    step->push_back(s);
    // четвертый шаг (1 временной слой, заведомо упругая разгрузка до 0)
    s.t1 = Time*3;
    s.t2 = Time*4;
    s.dt0 = Time/StepsNumber4;
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
        // материал
        task.mechTask.material = new std::vector<MechMaterialSource>(1);
        MechMaterialSource &mechMat = (*task.mechTask.material)[0];
        if(thermTestIndex != 0)
        {
            mechMat.elasticSigmaLimit = 1.e10;//2.e7;
        }
        else
        {
            if(fidesys)
                mechMat.elasticSigmaLimit = 24;//2.e7;
            else
                mechMat.elasticSigmaLimit = 2.e7;
        }
        mechMat.Ceps = 4./9.;
        mechMat.Csigma = 1;

        if(thermTestIndex == 1)
        {
            mechMat.set_E_NU(200*1.e9,NU);        //1.e10;
            mechMat.Talpha = 0.0001;
        }
        else
        {
            if(fidesys)
            {
                mechMat.set_E_NU(21000., NU);        //1.e10;
                mechMat.Talpha = 0;
            }
            else
            {
                mechMat.set_E_NU(1.e10, NU);
                mechMat.Talpha = 0;
            }
        }
        //mechMat.set_K_G(70.0e9, 26.0e9);
        mechMat.set_D_isotropic();
        mechMat.set_M_sigma();
        mechMat.F = {0,0,0};
        mechMat.ro = 0;//1e10;
        if(thermTestIndex != 0)
            setPlasticMaterialCurve(0, plasticityCurveMode, usingSecantMethod, mechMat);
        else
        {
            // setPlasticMaterialCurve(plasticPlot, plasticityCurveMode, Solid::MechStressStrainMethod::initialStrain, usingSecantMethod, mechMat);
            // setPlasticMaterialCurve(plasticPlot, plasticityCurveMode, Solid::MechStressStrainMethod::initialStrain, usingSecantMethod, mechMat);
             setPlasticMaterialCurve(plasticPlot, plasticityCurveMode, usingSecantMethod, mechMat);
        }
        //mechMat.mode = MechMaterialType::Elasticity;


        // шаги (на каждом шаге задаются вторые краевые условия)
        std::vector<MechGlobalStepParameters> *mechStep = new std::vector<MechGlobalStepParameters>;
        // нагружение - разгрузка
        MechGlobalStepParameters ms;  // шаг для МДТТ
        GlobalStep s;       // шаг
        std::vector<MechBoundaryCondition2Source_base *> *bc2Source;
        ms.slausolverParameters = slausolver_parameters;
        ms.timeMode = 0;    // 0 - квазистатическая задача
        ms.sigma0Mode = sigma0Mode;  // 1 - dSigma = integral(...), 2 - приращения
        ms.fixGrid = fixGrid;     // 1 - зафиксировать сетку
        ms.movingGridMode = 0;  // 0 - простые приращения, 1 - приращения Яумана
        ms.controlMode = 0; // 0 - игнорировать не достижение невязок(выдаются предупреждения в консоль) 1 - уменьшать шаг
        ms.iterLimit = nonlinearIterLimit;
        ms.slauResidualLimit = 1000;
        ms.epsResidualLimit = epsResidualLimit;
        ms.sigmaResidualLimit = 1000;
        // первый шаг (1 временной слой, заведомо упругое нагружение P1)
        s = (*step)[0];
        setBc2Sphere(bc2Source, 0, P1, s.t1, s.t2);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);
        // второй шаг (пластичное нагружение P2)
        s = (*step)[1];
        setBc2Sphere(bc2Source, P1, P2, s.t1, s.t2);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);
        // третий шаг (1 временной слой, начало разгрузки)
        s = (*step)[2];
        setBc2Sphere(bc2Source, P2, P3, s.t1, s.t2);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);
        // четвертый шаг (1 временной слой, заведомо упругая разгрузка до 0)
        s = (*step)[3];
        setBc2Sphere(bc2Source, P3, P4, s.t1, s.t2);
        ms.bc2Source = bc2Source;
        mechStep->push_back(ms);
        task.mechTask.mechStep = mechStep;
        // первые краевые условия
        std::vector<MechBoundaryCondition1Source> *mechBc1 = new std::vector<MechBoundaryCondition1Source>(5);
        (*mechBc1)[0].mode = {{ 0, -1, -1}};
        (*mechBc1)[0].u0 =   {{ 0, -1, -1}};
        (*mechBc1)[1].mode = {{-1,  0, -1}};
        (*mechBc1)[1].u0 =   {{-1,  0, -1}};
        (*mechBc1)[2].mode = {{-1, -1,  0}};
        (*mechBc1)[2].u0 =   {{-1, -1,  0}};
        (*mechBc1)[3].mode = {{-1, -1, -1}};
        (*mechBc1)[3].u0 =   {{-1, -1, -1}};
        (*mechBc1)[4].mode = {{-1, -1, -1}};
        (*mechBc1)[4].u0 =   {{-1, -1, -1}};
        //(*mechBc1)[4].mode = {{0, 0, 0}};
        //(*mechBc1)[4].u0 =   {{0, 0, 0}};
        task.mechTask.bc1Source = mechBc1;
        // инициализация начальных деформаций и напряжений (=0)
        task.mechTask.fe = new std::vector<MechFeData_base *>;
        task.mechTask.fe->resize(grid->fe.size());
        for (size_t i = 0; i < grid->fe.size(); i++)
        {
            MechFeData_LinearHexagon_homogeny_E_NLE_EP *feDataEl = new MechFeData_LinearHexagon_homogeny_E_NLE_EP;
            feDataEl->init(Integration::IntegrationType::Gauss3);
            (*task.mechTask.fe)[i] = feDataEl;
        }
        // инициализация начальных перемещений (=0)
        task.mechTask.vertex = new std::vector<MechVertexData>;
        for (size_t i = 0; i < grid->vertex.size(); i++)
        {
            MechVertexData v_el;
            for (int j = 0; j < 3; j++)
            {
                v_el.sum_du[j] = 0;
                v_el.du[j] = 0; // не имеет значения
            }
            task.mechTask.vertex->push_back(v_el);
        }

        // для криволинейных границ ###
        task.mechTask.vertexForCurvature = new std::vector<MechVertexData>;
        for (size_t i = 0; i < grid->vertexForCurvature.size(); i++)
        {
            MechVertexData v_el;
            for (int j = 0; j < 3; j++)
            {
                v_el.sum_du[j] = 0;
                v_el.du[j] = 0; // не имеет значения
            }
            task.mechTask.vertexForCurvature->push_back(v_el);
        }

        // инициализация начальных скоростей (=0) и ускорений (=0)
        task.mechTask.V0 = new Vector(grid->vertex.size()*3);
        task.mechTask.dV0 = new Vector(grid->vertex.size()*3);
        for (size_t i = 0; i < task.mechTask.V0->size(); i++)
        {
            (*task.mechTask.V0)[i] = 0;
            (*task.mechTask.dV0)[i] = 0;
        }

        // вторые краевые
        task.mechTask.bc2 = new std::vector<MechBoundaryCondition2>;
        MechBoundaryCondition2 bc2_el;
        bc2_el.FEsurfaceInd = 0;
        bc2_el.bc2SourceIndex = 0;
        (*task.mechTask.bc2).push_back(bc2_el);
        bc2_el.FEsurfaceInd = 1;
        bc2_el.bc2SourceIndex = 1;
        (*task.mechTask.bc2).push_back(bc2_el);

        // Контакт (отсутствует)
        // поверхность: жёсткий неподвижный цилиндр
        task.mechTask.rigidSurface = new std::vector<Surface_base *>;
        // контакты КЭ-поверхность - аналитическая поверхность
        task.mechTask.CC_FE_Rigid = new std::vector<ContactCondition_FE_Rigid>;
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
        // материал
        task.thermTask.material = new std::vector<ThermMaterialSource>(1);
        ThermMaterialSource &thermMat = (*task.thermTask.material)[0];
        thermMat.ro = 1;     // плотность
        thermMat.c = 5;      // удельная теплоёмкость
        thermMat.f = 0;      // мощность внутренних объёмных источников (стоков) тепла
        double L0 = 10;//10;
        thermMat.L =         // тензор теплопроводности
        {{
            {L0, 0, 0},
            {0, L0, 0},
            {0, 0, L0},
        }};
        int timeMode = 0;
        // Тест1
        if(thermTestIndex == 1)
            setBcSphereT_1(task.thermTask.bc1SourceT, task.thermTask.bc2SourceT,
                           1,   // температура на r = a
                           2    // температура на r = b
                           );
        // Тест2
        if(thermTestIndex == 2)
            setBcSphereT_2(task.thermTask.bc1SourceT, task.thermTask.bc2SourceT,
                           1,   // температура на r = a
                           2,   // температура окружающей среды r > b
                           1    // коэффициент конвективного теплообмена
                           );
        // Тест3
        if(thermTestIndex == 3)
            setBcSphereT_3(task.thermTask.bc1SourceT, task.thermTask.bc2SourceT,
                           1,   // температура окружающей среды r < a
                           10,  // коэффициент конвективного теплообмена
                           2,   // температура окружающей среды r > b
                           20   // коэффициент конвективного теплообмена
                           );
        // Тест4
        if(thermTestIndex == 4)
        {
            setBcSphereT_4(task.thermTask.bc1SourceT, task.thermTask.bc2SourceT, thermMat,
                           100, // плотность набегающего теплового потока на r = a
                           15,  // коэффициент теплопроводности
                           1    // температура на r = b
                           );
        }
        // Тест5
        if(thermTestIndex == 5)
        {
            setBcSphereT_5(task.thermTask.bc1SourceT, task.thermTask.bc2SourceT, thermMat,
                           10, // плотность набегающего теплового потока на r = a
                           1   // коэффициент теплопроводности
                           );
            timeMode = 1;
        }
        // шаги
        std::vector<ThermGlobalStep> *thermStep = new std::vector<ThermGlobalStep>;
        ThermGlobalStep ts;  // шаг для теплопроводности
        ts.slausolverParameters = slausolver_parameters;
        ts.timeMode = timeMode;    // статика
        thermStep->push_back(ts);
        thermStep->push_back(ts);
        thermStep->push_back(ts);
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
                double Q0 = (*task.thermTask.bc2SourceT)[0].q * 4*PI*SQR(a);
                T0 = T_a_ch + Q0/(4*PI*(*task.thermTask.material)[0].L.m[0][0])*(1./R - 1./a);
            }
            v_el.T = T0;
            v_el.newT = T0;
            task.thermTask.vertex->push_back(v_el);
        }
    }





    // сохренение информации о сетке в файл
    fprintf(f_inf, "r1 = %lf, r2 = %lf\n", gp.r1, gp.r2);
    fprintf(f_inf, "Grid: %dx%dx%d, q = %le\n", gp.N, gp.N, gp.Nparts, gp.q);
    if(gp.curvilinear == 0)
        fprintf(f_inf, "Linear\n");
    if(gp.curvilinear == 1)
        fprintf(f_inf, "Quadratic\n");
    if(gp.buildingMethod == 0)
        fprintf(f_inf, "buildingMethod = Divide arcs\n");
    if(gp.buildingMethod == 1)
        fprintf(f_inf, "Quadratic = Cube to sphere\n");
    fprintf(f_inf, "FENumber = %d\n", (int)grid->fe.size());
    fprintf(f_inf, "VertexesNumber = %d\n", (int)grid->vertex.size());
    fprintf(f_inf, "VertexesForCurvatureNumber = %d\n", (int)grid->vertexForCurvature.size());
    fprintf(f_inf, "bc1Number = %d\n", (int)grid->bc1.size());
    fprintf(f_inf, "FEsurface = %d\n", (int)grid->FESurface.size());
    for (size_t FEsurfaceInd = 0; FEsurfaceInd < grid->FESurface.size(); FEsurfaceInd++)
    {
        fprintf(f_inf, "FEsurface[%d]Size = %d\n", (int)FEsurfaceInd, (int)grid->FESurface[FEsurfaceInd].face.size());
    }
    //grid->buldOpenSCADModel("setka.scad");
    // сохренение информации о параметрах в файл
    fprintf(f_inf, "fixGrid = %d\n", fixGrid);
    fprintf(f_inf, "plasticPlot = %d\n", plasticPlot);
    fprintf(f_inf, "epsResidualLimit = %le\n", epsResidualLimit);
    fprintf(f_inf, "nonlinearIterLimit = %d\n", nonlinearIterLimit);
    fprintf(f_inf, "StepsNumber1 = %d\n", StepsNumber1);
    fprintf(f_inf, "StepsNumber2 = %d\n", StepsNumber2);
    fprintf(f_inf, "StepsNumber3 = %d\n", StepsNumber3);
    fprintf(f_inf, "StepsNumber4 = %d\n", StepsNumber4);
    fprintf(f_inf, "P = %lf\n", P);
    fprintf(f_inf, "P1 = %lf\n", P1);
    fprintf(f_inf, "P2 = %lf\n", P2);
    fprintf(f_inf, "P3 = %lf\n", P3);
    fprintf(f_inf, "P4 = %lf\n", P4);
    fprintf(f_inf, "sigma0Mode = %d\n", sigma0Mode);
    fprintf(f_inf, "HomogenyMode = %d\n", HomogenyMode);

    fclose(f_inf);
}
void Test_sphere::writeResults(const Task &task, const Solid::OutData &out)
{

    //task.grid->buldOpenSCADModel("setka.scad");
    using namespace Interpolation;
    //FILE * f = fopen("L.txt", "w");

    //bt.P = 1.e8;

    //SolverParameters sp = *task.sp;
    //Grid3D *grid = task.grid;
    //MassiveStack<Solid_Material> *materials = task.materials;
    //MassiveStack<Solid_BoundaryCondition1_source> *bc1 = task.bc1;
    //MassiveStack<Solid_BoundaryCondition2_source> *bc2 = task.bc2;
    //const std::vector<POINT3> *points = task.points;
        //BeamTask &bt = task.beamTask;

        //printf("Pmax = %le\n", bt.P);
    // инициализация
    // запуск решателя
    //Paint p;
    //p.buld_OpenSCAD_models("_OpenSCAD_modes.txt", g0, res0);
   //return;
    //return;
    //Paint p;
    //p.buld_OpenSCAD_models("_OpenSCAD_modes.txt", task.grid, &res);
    /*
    // Проверка ответа
    double Pnn = 0;
    // последовательное добавление поверхностных сил
    for (int nn = 0; nn <= sp.PN; nn++)
    {
        double Padd;
        if (nn == 0)
            Padd = bt.P * sp.proportion_of_P_elastic;
        else
            Padd = bt.P * (1 - sp.proportion_of_P_elastic) / sp.PN;
        Pnn += Padd;
        //printf("P = %le\n", Pnn);

        Solid_Material &m0 = (*task.materials)[0];

        // Проверка ответа
        double epsx, epsy, psixy;
        double epsx_numb, epsy_numb, psixy_numb;

        if (bt.type == BeamTask_x)
        {
            // одноосное растяжение
            epsx = (Pnn/3)*(1/(3*m0.K) * m0.k(Pnn/(9*m0.K)) +
                              (1/m0.G) * m0.g((2./9.) * Pnn*Pnn / m0.G/m0.G));
            epsy = (Pnn/3)*(1/(3*m0.K) * m0.k(Pnn/(9*m0.K)) -
                              (1/(2*m0.G)) * m0.g((2./9.) * Pnn*Pnn / m0.G/m0.G));
            // решение из книжки для 2/9, 8/9
            //epsx =
            //    (1. / 3.)*(1. / (3. * m0.K) +
            //    (1. + m0.g2 * (2. / 9.) * Pnn * Pnn / m0.G / m0.G) / m0.G)*Pnn;
            //epsy =
            //    (1. / 3.)*(1. / (3. * m0.K) -
            //    (1. / 2.)*
            //    (1. + m0.g2 * (2. / 9.) * Pnn * Pnn / m0.G / m0.G) / m0.G)*Pnn;
            // решение для произвольных коэффицентов интенсивностей
            //VECTOR6 sigma0 = res->fe[0].sigma0;
            //VECTOR6 sigma0 = { Pnn, 0, 0, 0, 0, 0 };
            //double eps0_eqv = m0.eps_sigma(m0.solve_sigma_eqv(sigma0));
            //epsx = Pnn / (9 * m0.K) + (2 * eps0_eqv) / (3 * sqrt(m0.c_M_eps));
            //epsy = Pnn / (9 * m0.K) - (eps0_eqv) / (3 * sqrt(m0.c_M_eps));
            //численное решение
            epsx_numb = out[0].pres[nn].eps.x[0];
            epsy_numb = out[0].pres[nn].eps.x[1];
        }
        if (bt.type == BeamTask_xyz)
        {
            // трехосное (гидростатическое) растяжение
            epsx = (Pnn/(3.*m0.K)) * m0.k(Pnn/(3*m0.K));
            // решение из книжки для 2/9, 8/9
            //epsx =
            //    (Pnn / (3. * m0.K)) * 1;
            //численное решение
            epsx_numb = out[0].pres[nn].eps.x[0];
        }
        if (bt.type == BeamTask_x_y)
        {
            // растяжение и сжатие по разным осям
            epsx = (Pnn/(2*m0.G)) * m0.g((2./3.) * Pnn*Pnn / m0.G/m0.G);
            epsy = -epsx;
            // решение из книжки для 2/9, 8/9
            //epsx =
            //    (Pnn / (2 * m0.G)) * ((1 + m0.g2*(2. / 3.)*Pnn*Pnn / m0.G / m0.G));
            //epsy = -epsx;
            //численное решение
            epsx_numb = out[0].pres[nn].eps.x[0];
            epsy_numb = out[0].pres[nn].eps.x[1];
        }
        if (bt.type == BeamTask_around_z)
        {
            // чистый сдвиг
            psixy = (Pnn/m0.G) * m0.g((2. / 3.) * Pnn*Pnn / m0.G/m0.G);
            // решение из книжки для 2/9, 8/9
            //psixy =
            //    (Pnn / m0.G) * ((1 + m0.g2*(2. / 3.)*Pnn*Pnn / m0.G / m0.G));
            //численное решение
            psixy_numb = out[0].pres[nn].eps.x[5];
        }

        // запись в файл
        if (bt.type == BeamTask_x || bt.type == BeamTask_xyz || bt.type == BeamTask_x_y)
        {
            fprintf(fepsx, "p = %le, %le vs %le, pogr = %le\n", Pnn, epsx_numb, epsx, ABS((epsx_numb - epsx) / epsx));
            fprintf(fepsx_graph, "%le %le\n", Pnn, epsx);
            fprintf(fepsx_pogr_graph, "%le %lf\n", Pnn, log10(ABS((epsx_numb - epsx) / epsx)));
        }
        if (bt.type == BeamTask_x || bt.type == BeamTask_x_y)
        {
            fprintf(fepsy, "p = %le, %le vs %le, pogr = %le\n", Pnn, epsy_numb, epsy, ABS((epsy_numb - epsy) / epsy));
            fprintf(fepsy_graph, "%le %le\n", Pnn, epsy);
            fprintf(fepsy_pogr_graph, "%le %lf\n", Pnn, log10(ABS((epsy_numb - epsy) / epsy)));

        }
        if (bt.type == BeamTask_around_z)
        {
            fprintf(fpsixy, "p = %le, %le vs %le, pogr = %le\n", Pnn, psixy_numb, psixy, ABS((psixy_numb - psixy) / psixy));
            fprintf(fpsixy_graph, "%le %le\n", Pnn, psixy);
            fprintf(fpsixy_pogr_graph, "%le %lf\n", Pnn, log10(ABS((psixy_numb - psixy) / psixy)));
        }
    }
    {
        res.save("out_u.txt", "out_eps.txt", "out_gamma.txt", "out_sigma.txt", "out_tau.txt");
        //Paint p;
        //p.buld_OpenSCAD_models("_OpenSCAD_modes.txt", &g, &res);
        //p.paint_pictures("_Paint.txt", &g, &res);
    }

    fclose(fepsx);
    fclose(fepsy);
    fclose(fpsixy);
    fclose(fepsx_graph);
    fclose(fepsy_graph);
    fclose(fpsixy_graph);
    fclose(fepsx_pogr_graph);
    fclose(fepsy_pogr_graph);
    fclose(fpsixy_pogr_graph);



    Solid_PointResult &r = out[0].pres[sp.PN];
    double L = r.u.x[0]*2;
    fprintf(f, "%d   %le\n", sp.PN, L);
    printf("%d   %le\n", sp.PN, L);
    */
    //}
    // Проверка ответа

    FILE *fT = fopen(fn_T.c_str(), "w");
    FILE *fT_ch = fopen(fn_T_ch.c_str(), "w");
    FILE *fT_pogr_abs = fopen(fn_T_pogr_abs.c_str(), "w");
    FILE *fT_pogr = fopen(fn_T_pogr.c_str(), "w");

    FILE *fu = fopen(fn_u.c_str(), "w");
    FILE *fu_ch = fopen(fn_u_ch.c_str(), "w");
    FILE *fu_pogr_abs = fopen(fn_u_pogr_abs.c_str(), "w");
    FILE *fu_pogr = fopen(fn_u_pogr.c_str(), "w");

    FILE *fsigma_r = fopen(fn_sigma_r.c_str(), "w");
    FILE *fsigma_r_ch = fopen(fn_sigma_r_ch.c_str(), "w");
    FILE *fsigma_fi = fopen(fn_sigma_fi.c_str(), "w");
    FILE *fsigma_fi_ch = fopen(fn_sigma_fi_ch.c_str(), "w");
    FILE *fsigma_r_0 = fopen(fn_sigma_r_0.c_str(), "w");
    FILE *fsigma_r_ch_0 = fopen(fn_sigma_r_ch_0.c_str(), "w");
    FILE *fsigma_fi_0 = fopen(fn_sigma_fi_0.c_str(), "w");
    FILE *fsigma_fi_ch_0 = fopen(fn_sigma_fi_ch_0.c_str(), "w");

    FILE *fsigma_r_pogr_abs = fopen(fn_sigma_r_pogr_abs.c_str(), "w");
    FILE *fsigma_r_pogr_0_abs = fopen(fn_sigma_r_pogr_0_abs.c_str(), "w");
    FILE *fsigma_fi_pogr_abs = fopen(fn_sigma_fi_pogr_abs.c_str(), "w");
    FILE *fsigma_fi_pogr_0_abs = fopen(fn_sigma_fi_pogr_0_abs.c_str(), "w");
    FILE *fsigma_r_pogr = fopen(fn_sigma_r_pogr.c_str(), "w");
    FILE *fsigma_r_pogr_0 = fopen(fn_sigma_r_pogr_0.c_str(), "w");
    FILE *fsigma_fi_pogr = fopen(fn_sigma_fi_pogr.c_str(), "w");
    FILE *fsigma_fi_pogr_0 = fopen(fn_sigma_fi_pogr_0.c_str(), "w");

    FILE *fres_maxpogr = fopen(fn_res_maxpogr.c_str(), "w");


    //FILE *fsigma_eqv_ch = fopen("__res_eqv_ch.txt", "w");
     //FILE *fs_t = fopen("__res_s_step0.txt", "w");
     //FILE *fs_dstep = fopen("__res_s_dstep.txt", "w");
    FILE *fs_iterationsNumber = fopen(fn_s_iterationsNumber.c_str(), "w");
    FILE *fs_slauResidual = fopen(fn_s_slauResidual.c_str(), "w");
    FILE *fs_stepSolvingTime = fopen(fn_s_stepSolvingTime.c_str(), "w");
    FILE *fs_nonlinearStateFENumber = fopen(fn_s_nonlinearStateFENumber.c_str(), "w");
    FILE *fs_slauNev = fopen(fn_s_slauNev.c_str(), "w");
    FILE *fs_epsNev = fopen(fn_s_epsNev.c_str(), "w");
    FILE *fs_sigmaNev = fopen(fn_s_sigmaNev.c_str(), "w");

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
        t1 = (*step)[globalStepNumber1].t2;
        t0 = (*step)[globalStepNumber0].t2;
    }
    if(step->size() == 2)
    {
        globalStepNumber1 = 0;
        globalStepNumber0 = 1;
        t1 = (*step)[globalStepNumber1].t2;
        t0 = (*step)[globalStepNumber0].t2;
    }
    if(step->size() == 1)
    {
        globalStepNumber1 = 0;
        globalStepNumber0 = -123;
        t1 = (*step)[globalStepNumber1].t2;
        t0 = -123;
    }

    // последовательное добавление поверхностных сил
    double P0 = 0;
    for (int globalStepNumber = 0; globalStepNumber < (int)task.step->size(); globalStepNumber++)
    {
        double T_pogr_max = 0;
        double T_pogr_max_abs = 0;
        double fi_pogr_max_abs = 0;
        double r_pogr_max_abs = 0;
        double fi_0_pogr_max_abs = 0;
        double r_0_pogr_max_abs = 0;
        double sigma_r_max = 0;
        double sigma_fi_max = 0;
        double sigma_r_0_max = 0;
        double sigma_fi_0_max = 0;
        double fi_pogr_max = 0;
        double r_pogr_max = 0;
        double fi_0_pogr_max = 0;
        double r_0_pogr_max = 0;

        // вывод картинок МДТТ
/*
        GridRectangleRegular3D gr;
        //double dz = task.sphere.r2;//2./task.sphere.N;
        gr.init({-task.sphere.r2/2,-task.sphere.r2/2,-task.sphere.r2/2},
                { task.sphere.r2/2, task.sphere.r2/2, task.sphere.r2/2},
                task.sphere.Nparts,
                task.sphere.Nparts,
                task.sphere.Nparts);
        Interpolator ip;

        {
            OutStepData &s = (*out)[nnn].step[(int)(*out)[nnn].step.size() - 1];
            if(s.t == t0 || s.t == t1)
            for(int resInd = 0; resInd < 2; resInd++)
            {
                OutData &o = (*out)[nnn];
                ip.init(&gr, InterpolatorType::Lagrange3D1, 0.01, 0, o.fe.size()*27);
                //ip.init(&gr, InterpolatorType::Hermite2D3, 0.01, 0, o.fe.size());
                for(size_t k = 0; k < o.fe.size(); k++)
                {
                    for (size_t pInd = 0; pInd < o.fe[k].pd.size(); pInd++)	// индекс точки Гаусса внутри конечного элемента
                    {
                        OutFePointData &outFePD = o.fe[k].pd[pInd];

                        POINT3 p = outFePD.p;
                        // главные напряжения
                        double max;
                        VECTOR3 ms;
                        outFePD.sumSigma.solveMain(ms);
                        int ind[3];
                        sortms(ms, ind);
                        max = 0;
                        double v0 = abs(ms[0]);
                        double v1 = abs(ms[1]);
                        double v2 = abs(ms[2]);
                        //if(v0 == v0 && v0 > max) max = v0;
                        //if(v1 == v1 && v1 > max) max = v1;
                        //if(v2 == v2 && v2 > max) max = v2;
                        if(!(v0 == v0) || !(v1 == v1) || !(v2 == v2))
                        {
                            print1("\ncomplex!\n");
                            max = 10e10;
                        }
                        // sigma_r
                        if(resInd == 0)
                        {
                            max = -ms[ind[0]];
                        }
                        // sigma_fi
                        if(resInd == 1)
                        {
                            max = -ms[ind[1]];//-(ms[ind[1]] + ms[ind[2]])/2;
                            //if(s.t == t0)// нагрузка

                            //if(s.t == t1)// разгрузка
                            //    max = ms[ind[1]];
                        }
                        //max = fe_el.sigma[0];
                        //if(abs(p[2]) < dz)
                        {
                            if(p[0] > gr.p1[0] && p[0] < gr.p2[0] &&
                               p[1] > gr.p1[1] && p[1] < gr.p2[1] &&
                               p[2] > gr.p1[2] && p[2] < gr.p2[2])
                                    ip.addPoint(p, max, UNKNOWN_FE_INDEX);
                            //ip.addPoint(p, log(max + 1), UNKNOWN_FE_INDEX);
                        }
                    }
                }
                ip.makeInterpolant();
                char fn[1000];
                if(resInd == 0)
                    sprintf(fn, "out%.2d_r.bmp", nnn);
                if(resInd == 1)
                    sprintf(fn, "out%.2d_fi.bmp", nnn);
                std::vector<Circle> circle(2);
                circle[0].x = 0;
                circle[0].y = 0;
                circle[0].r = 1;
                circle[0].mode = 0;
                circle[1].x = 0;
                circle[1].y = 0;
                circle[1].r = 4;
                circle[1].mode = 1;
                // ч/б, с градациями
                ip.saveInterpolant(0, fn, 321, 321, 1, 1, circle);
                if(resInd == 0)
                    sprintf(fn, "outc%.2d_r.bmp", nnn);
                if(resInd == 1)
                    sprintf(fn, "outc%.2d_fi.bmp", nnn);
                // цветной, с градациями
                ip.saveInterpolant(0, fn, 321, 321, 3, 1, circle);
                //OutStepData &s = o.step[nn];
                ip.release();
            }

        }
*/
        // вывод данных МДТТ
        if(task.mechTask.enabled)
        for (int numLocalStep = 0; numLocalStep < (int)(*out.mechOut)[globalStepNumber].step.size(); numLocalStep++)
        {
            MechOutStepData &mechOutStep_el = (*out.mechOut)[globalStepNumber].step[numLocalStep];
            double P;
            MechBoundaryCondition2Source_ScalarFunction *mbc2Source_sf = (MechBoundaryCondition2Source_ScalarFunction *)(*(*task.mechTask.mechStep)[globalStepNumber].bc2Source)[0];
            P = mbc2Source_sf->value(mechOutStep_el.t);
            if(mechOutStep_el.t == t1)
                P0 = P;
            //P0 += s.dt * (*(*step)[globalStepNumber].bc2Source)[0].dPdt[0];
            //fprintf(f, "p = %le\n", Pnn);
            //fprintf(f, "nn = %d, step0 = %le, P = %le\n", nn, s.step0, P0);
            //fprintf(fsigma_eqv_ch, "nn = %d, step0 = %le, P = %le\n", nn, s.step0, P0);
            // радиусы
            //double a = gp.r1;
            //double b = gp.r2;
            POINT3 v1 = grid->vertex[0];
            //POINT3 v2 = grid->vertex[task.sphere.Nparts];//16,32,64
            int N = gp.N/2+1;
            POINT3 v2 = grid->vertex[(3*N*N-3*N)*(gp.Nparts+1)];
            double a = v1.abs();
            double b = v2.abs();
            //OutFeData &r = (*out)[globalStepNumber].fe[j];
            //POINT3 C = r.p;
            //fprintf(fs_t, "%d %le\n", numLocalStep, mechOutStep_el.t);
            //fprintf(fs_dstep, "%le %le\n", ms.t, ms.dt);
            MechIterInf nlInf = mechOutStep_el.nlInf.iterInf.back();
            fprintf(fs_iterationsNumber, "%le %d\n", mechOutStep_el.t, nlInf.iterNumber);
            fprintf(fs_slauResidual, "%le %le\n", mechOutStep_el.t, nlInf.slau.slauResidual);
            //fprintf(fs_stepSolvingTime, "%le %le\n", mechOutStep_el.t, mechOutStep_el.nlInf.stepSolvingTime);
            fprintf(fs_nonlinearStateFENumber, "%le %d\n", mechOutStep_el.t, nlInf.plastic.nonlinearStateFENumber);
            fprintf(fs_slauNev, "%le %le\n", mechOutStep_el.t, nlInf.slau.slauResidualWithLastq);
            fprintf(fs_epsNev, "%le %le\n", mechOutStep_el.t, nlInf.plastic.maxEpsResidual);
            fprintf(fs_sigmaNev, "%le %le\n", mechOutStep_el.t, nlInf.plastic.maxSigmaResidual);


            if(mechOutStep_el.t == t0 || mechOutStep_el.t == t1)
            {
                MechMaterialSource &m0 = (*task.mechTask.material)[0];
                //P0 = bt.P;
                // расчет численного решения в узлах сетки (перемещения) и его погрешности
                if(mechOutStep_el.t == t1)
                {
                    for(size_t j = 0; j < (*out.mechOut)[globalStepNumber].vertex.size(); j++)
                    {
                        MechOutVertexData &v = (*out.mechOut)[globalStepNumber].vertex[j];
                        double R = grid->vertex[j].abs();
                        double u = m0.Talpha*30*R;
                        double u_ch = v.sum_du.abs();
                        double u_pogr_abs = u_ch - u;
                        double u_pogr = (u_ch - u)/fabs(u);
                        /*
                        if(ABS(T_pogr) > T_pogr_max)
                            T_pogr_max = ABS(T_pogr);
                        if(ABS(T_pogr_abs) > T_pogr_max_abs)
                            T_pogr_max_abs = ABS(T_pogr_abs);
                            */
                        fprintf(fu_ch, "%le %le\n", R, u_ch);
                        fprintf(fu_pogr, "%le %le\n", R, u_pogr*100);
                        fprintf(fu_pogr_abs, "%le %le\n", R, u_pogr_abs);
                    }
                    // расчет аналитического решения
                    int NumPoints = 1000;
                    for(int j = 0; j <= NumPoints; j++)
                    {
                        double R = a + (b-a)*j/NumPoints;
                        double u = m0.Talpha*30*R;
                        fprintf(fu, "%le %le\n", R, u);
                    }
                }


                // расчет численного решения в конечных элементах
                for(int j0 = 0; j0 < gp.Nparts; j0++)//int j0 = gp.Nparts - 1;
                {
                    double R = 0;
                    //double sigma_eqv = 0;
                    //double eps_eqv = 0;
                    double sigma_r_ch = 0;
                    double sigma_fi_ch = 0;
                    int resSize = (int)(*out.mechOut)[globalStepNumber].fe.size();
                    std::vector<ShperePoint> sp(resSize / gp.Nparts * 27);
                    int spInd = 0;
                    for(int j = j0; j < resSize; j += gp.Nparts)// выбераем КЭ на одном уровне по R
                    {
                        for (size_t pInd = 0; pInd < (*out.mechOut)[globalStepNumber].fe[j].pd.size(); pInd++)	// индекс точки Гаусса внутри конечного элемента
                        {
                            MechOutFePointData &r = (*out.mechOut)[globalStepNumber].fe[j].pd[pInd];

                            //OutFeData &r = (*out)[globalStepNumber].fe[j];
                            POINT3 C = r.p;
                            R = C.abs();
                            //sigma_eqv = m0.sigmaEqv(r.sigma);
                            //eps_eqv = m0.epsEqv(r.eps);

                            VECTOR3 ms;
                            if(sigmaSolvingType == 0)
                            {
                                // главные напряжения
                                r.mainStresses(ms);
                            }

                            if(sigmaSolvingType == 1)
                            {
                                // напряжения на площадке ###
                                VECTOR3 norm = -r.p;
                                norm = norm / norm.abs();  // направление к центру
                                MATR3x3 sigma3x3;
                                VECTOR3 vectSigma;
                                r.sumSigma.toMATR3x3(sigma3x3); // тензор напряжений
                                {
                                    vectSigma = sigma3x3*norm;// ms - вектор напряжений на площадке
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
                                /*VECTOR3 norm2 = norm*norm1;
                                {
                                    vectSigma = sigma3x3*norm2;// ms - вектор напряжений на площадке
                                    double absSigma = vectSigma*norm2;  // проекция на нормаль к площадке
                                    ms[2] = absSigma;
                                }*/
                            }



                            int ind[3];
                            sortms(ms, ind);

                            // ###########
                            //ms[ind[0]] = m0.epsEqv(r.sumEpsElastic);
                            //ms[ind[1]] = m0.epsEqv(r.sumEpsPlastic);



                            sigma_r_ch = ms[ind[0]];
                            sigma_fi_ch = ms[ind[1]];//(ms[ind[1]] + ms[ind[2]]) / 2.;


                            /*
                            double ttt;
                            POINT3 norm1 = r.p/r.p.abs();
                            solveSigmaFlat(r.sigma, norm1, sigma_r_ch, ttt);
                            POINT3 norm2 = {0,0,0};
                            if(abs(norm1[2]) > 1.e-5)
                                norm2 = {1, 1, -(norm1[0]+norm1[1]/norm1[2])};
                            else
                                if(abs(norm1[1]) > 1.e-5)
                                    norm2 = {1, -(norm1[0]+norm1[2]/norm1[1]), 1};
                                else
                                    norm2 = {-(norm1[1]+norm1[2]/norm1[0]), 1, 1};
                            norm2 /= norm2.abs();
                            solveSigmaFlat(r.sigma, norm2, sigma_fi_ch, ttt);
                            */
                            sp[spInd].R = R;
                            sp[spInd].sigma_r_ch = sigma_r_ch;
                            sp[spInd].sigma_fi_ch = sigma_fi_ch;
                            double c = solvec(a, b, P0, m0.elasticSigmaLimit);
                            double sigma_r;
                            double sigma_fi;
                            double sigma_r_0;
                            double sigma_fi_0;
                            solveAnalit(a,b,R,P0,m0.elasticSigmaLimit,c,sigma_r,sigma_fi,sigma_r_0,sigma_fi_0);

                            if(mechOutStep_el.t == t1)
                            {
                                sp[spInd].sigma_r = sigma_r;
                                sp[spInd].sigma_fi = sigma_fi;
                                sp[spInd].sigma_r_pogr_abs = -(sp[spInd].sigma_r - sp[spInd].sigma_r_ch);
                                sp[spInd].sigma_fi_pogr_abs = -(sp[spInd].sigma_fi - sp[spInd].sigma_fi_ch);
                                sp[spInd].sigma_r_pogr = -((sp[spInd].sigma_r - sp[spInd].sigma_r_ch)/fabs(sp[spInd].sigma_r));
                                sp[spInd].sigma_fi_pogr = -((sp[spInd].sigma_fi - sp[spInd].sigma_fi_ch)/fabs(sp[spInd].sigma_fi));

                                if(fabs(sp[spInd].sigma_fi_pogr_abs) > fi_pogr_max_abs)
                                    fi_pogr_max_abs = fabs(sp[spInd].sigma_fi_pogr_abs);
                                if(fabs(sp[spInd].sigma_r_pogr_abs) > r_pogr_max_abs)
                                    r_pogr_max_abs = fabs(sp[spInd].sigma_r_pogr_abs);
                                //if(R >= 1.2)  //## погрешность для sigma_fi считается для радиуса от 1.2 до 4
                                if(fabs(sp[spInd].sigma_fi_pogr) > fi_pogr_max)
                                    fi_pogr_max = fabs(sp[spInd].sigma_fi_pogr);
                                //if(R <= 3)  //## погрешность для sigma_r считается для радиуса от 1 до 3
                                if(fabs(sp[spInd].sigma_r_pogr) > r_pogr_max)
                                    r_pogr_max = fabs(sp[spInd].sigma_r_pogr);
                                fprintf(fsigma_r_ch, "%le %le\n", R, sp[spInd].sigma_r_ch);
                                fprintf(fsigma_fi_ch, "%le %le\n", R, sp[spInd].sigma_fi_ch);
                                fprintf(fsigma_r_pogr_abs, "%le %le\n", R, sp[spInd].sigma_r_pogr_abs*100);
                                fprintf(fsigma_fi_pogr_abs, "%le %le\n", R, sp[spInd].sigma_fi_pogr_abs*100);
                                fprintf(fsigma_r_pogr, "%le %le\n", R, fabs(sp[spInd].sigma_r_pogr)*100);
                                fprintf(fsigma_fi_pogr, "%le %le\n", R, fabs(sp[spInd].sigma_fi_pogr)*100);
                            }else
                            if(mechOutStep_el.t == t0)
                            {
                                sp[spInd].sigma_r = sigma_r_0;
                                sp[spInd].sigma_fi = sigma_fi_0;
                                sp[spInd].sigma_r_pogr_abs = -(sp[spInd].sigma_r - sp[spInd].sigma_r_ch);
                                sp[spInd].sigma_fi_pogr_abs = -(sp[spInd].sigma_fi - sp[spInd].sigma_fi_ch);
                                sp[spInd].sigma_r_pogr = -((sp[spInd].sigma_r - sp[spInd].sigma_r_ch)/fabs(sp[spInd].sigma_r));
                                sp[spInd].sigma_fi_pogr = -((sp[spInd].sigma_fi - sp[spInd].sigma_fi_ch)/fabs(sp[spInd].sigma_fi));
                                if(fabs(sp[spInd].sigma_fi_pogr_abs) > fi_0_pogr_max_abs)
                                    fi_0_pogr_max_abs = fabs(sp[spInd].sigma_fi_pogr_abs);
                                if(fabs(sp[spInd].sigma_r_pogr_abs) > r_0_pogr_max_abs)
                                    r_0_pogr_max_abs = fabs(sp[spInd].sigma_r_pogr_abs);
                                if(fabs(sp[spInd].sigma_fi_pogr) > fi_0_pogr_max)
                                    fi_0_pogr_max = fabs(sp[spInd].sigma_fi_pogr);
                                if(fabs(sp[spInd].sigma_r_pogr) > r_0_pogr_max)
                                    r_0_pogr_max = fabs(sp[spInd].sigma_r_pogr);
                                fprintf(fsigma_r_ch_0, "%le %le\n", R, sp[spInd].sigma_r_ch);
                                fprintf(fsigma_fi_ch_0, "%le %le\n", R, sp[spInd].sigma_fi_ch);
                                fprintf(fsigma_r_pogr_0_abs, "%le %le\n", R, sp[spInd].sigma_r_pogr_abs*100);
                                fprintf(fsigma_fi_pogr_0_abs, "%le %le\n", R, sp[spInd].sigma_fi_pogr_abs*100);
                                fprintf(fsigma_r_pogr_0, "%le %le\n", R, fabs(sp[spInd].sigma_r_pogr)*100);
                                fprintf(fsigma_fi_pogr_0, "%le %le\n", R, fabs(sp[spInd].sigma_fi_pogr)*100);
                            }
                            spInd++;
                        }
                    }
                    //double coef = sp.size();
                    // усредненные численные результаты
                    /*
                    R = 0;
                    sigma_r_ch = 0;
                    sigma_fi_ch = 0;
                    for(size_t spInd = 0; spInd < sp.size(); spInd++)
                    {
                        R += sp[spInd].R;
                        sigma_r_ch += sp[spInd].sigma_r_ch;
                        sigma_fi_ch += sp[spInd].sigma_fi_ch;
                    }
                    R /= coef;
                    sigma_r_ch /= coef;
                    sigma_fi_ch /= coef;
                    */
                    // максимальные sigma_r_ch
                    /*
                    int maxInd = 0;
                    for(size_t spInd = 1; spInd < sp.size(); spInd++)
                    {
                        if(sp[spInd].sigma_fi_ch < sp[maxInd].sigma_fi_ch)
                            maxInd = spInd;
                    }
                    R = sp[maxInd].R;
                    sigma_r_ch = sp[maxInd].sigma_r_ch;
                    sigma_fi_ch = sp[maxInd].sigma_fi_ch;

                    maxInd = 0;
                    for(size_t spInd = 1; spInd < sp.size(); spInd++)
                    {
                        if(sp[spInd].sigma_fi_ch > sp[maxInd].sigma_fi_ch)
                            maxInd = spInd;
                    }
                    R = sp[maxInd].R;
                    sigma_r_ch = sp[maxInd].sigma_r_ch;
                    sigma_fi_ch = sp[maxInd].sigma_fi_ch;
                    */
                    /*
                    //fprintf(f, "R = %le  SigmaEqv = %le EpsEqv = %le Sigma_r = %le   Sigma_fi = %le\n", R, sigma_eqv, eps_eqv, sigma_r_ch, sigma_fi_ch);
                    // аналитическое решение
                    double c = solvec(a, b, P0, m0.elasticSigmaLimit);
                    double sigma_r;
                    double sigma_fi;
                    double sigma_r_0;
                    double sigma_fi_0;
                    solveAnalit(a,b,R,P0,m0.elasticSigmaLimit,c,sigma_r,sigma_fi,sigma_r_0,sigma_fi_0);

                    if(s.t == t1)
                    {
                        fprintf(fsigma_r_ch, "%le %le\n", R, sigma_r_ch);
                        fprintf(fsigma_fi_ch, "%le %le\n", R, sigma_fi_ch);
                        fprintf(fsigma_r_pogr, "%le %le\n", R, abs((sigma_r_ch-sigma_r)));
                        fprintf(fsigma_fi_pogr, "%le %le\n", R, abs((sigma_fi_ch-sigma_fi)));

                    }else
                    if(s.t == t0)
                    {
                        fprintf(fsigma_r_ch_0, "%le %le\n", R, sigma_r_ch);
                        fprintf(fsigma_fi_ch_0, "%le %le\n", R, sigma_fi_ch);
                        fprintf(fsigma_r_pogr_0, "%le %le\n", R, abs((sigma_r_ch-sigma_r_0)));
                        fprintf(fsigma_fi_pogr_0, "%le %le\n", R, abs((sigma_fi_ch-sigma_fi_0)));
                    }
*/
                    //if(sigma_eqv >= elasticSigmaLimit*0.95)
                    //    fprintf(fsigma_eqv_ch, "## ");
                    //fprintf(fsigma_eqv_ch, "R = %le eps_eqv = %.15le vs %.15le sigma_eqv = %le sigma_fi-sigma_r = %le\n", R, eps_eqv, m0.eps_sigma(sigma_eqv, r.sigmaT), sigma_eqv, sigma_r_ch - sigma_fi_ch);
                }

                // расчет аналитического решения
                int NumPoints = 1000;
                double c = solvec(a, b, P0, m0.elasticSigmaLimit);
                for(int j = 0; j <= NumPoints; j++)
                {
                    double R = a + (b-a)*j/NumPoints;
                    double sigma_r;
                    double sigma_fi;
                    double sigma_r_0;
                    double sigma_fi_0;
                    solveAnalit(a,b,R,P0,m0.elasticSigmaLimit,c,sigma_r,sigma_fi,sigma_r_0,sigma_fi_0);
                    if(mechOutStep_el.t == t1)
                    {
                        fprintf(fsigma_r, "%le %le\n", R, sigma_r);
                        fprintf(fsigma_fi, "%le %le\n", R, sigma_fi);
                        if(fabs(sigma_r) > sigma_r_max)
                            sigma_r_max = fabs(sigma_r);
                        if(fabs(sigma_fi) > sigma_fi_max)
                            sigma_fi_max = fabs(sigma_fi);
                    }else
                    if(mechOutStep_el.t == t0)
                    {
                        fprintf(fsigma_r_0, "%le %le\n", R, sigma_r_0);
                        fprintf(fsigma_fi_0, "%le %le\n", R, sigma_fi_0);
                        if(fabs(sigma_r_0) > sigma_r_0_max)
                            sigma_r_0_max = fabs(sigma_r_0);
                        if(fabs(sigma_fi_0) > sigma_fi_0_max)
                            sigma_fi_0_max = fabs(sigma_fi_0);
                    }
                }
                //fprintf(fsigma_fi_0, "c = %le\n", c);
            }

        }

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
                ThermMaterialSource &m0 = (*task.thermTask.material)[0];
                // температура
                double T1 = (*task.thermTask.bc1SourceT)[3].T0;
                double T2 = (*task.thermTask.bc1SourceT)[4].T0;
                double Ta1 = (*task.thermTask.bc2SourceT)[0].Ta;
                double hi1 = (*task.thermTask.bc2SourceT)[0].hi;
                double h1 = hi1/m0.L.m[0][0];
                double Ta2 = (*task.thermTask.bc2SourceT)[1].Ta;
                double hi2 = (*task.thermTask.bc2SourceT)[1].hi;
                double h2 = hi2/m0.L.m[0][0];
                //double T_b_ch = (*out->thermOut)[globalStepNumber].vertex[task.sphere.Nparts].T;
                double T_b_ch = (*out.thermOut)[globalStepNumber].vertex[N_b].T;
                double T_a_ch = (*out.thermOut)[globalStepNumber].vertex[0].T;
                double Q0 = (*task.thermTask.bc2SourceT)[0].q * 4*PI*SQR(a);
                //int resSize = (int)(*out.thermOut)[globalStepNumber].vertex.size();
                // расчет численного решения и его погрешности
                for(size_t j = 0; j < (*out.thermOut)[globalStepNumber].vertex.size(); j++)
                {
                    ThermOutVertexData &v = (*out.thermOut)[globalStepNumber].vertex[j];
                    double R = grid->vertex[j].abs();
                    //if(R >= 1.05)
                    {
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
                    double T_pogr_abs = v.T-T;
                    double T_pogr = (v.T-T)/fabs(T);
                    if(fabs(T_pogr) > T_pogr_max)
                        T_pogr_max = fabs(T_pogr);
                    if(fabs(T_pogr_abs) > T_pogr_max_abs)
                        T_pogr_max_abs = fabs(T_pogr_abs);
                    fprintf(fT_ch, "%le %le\n", R, v.T);
                    fprintf(fT_pogr, "%le %le\n", R, T_pogr*100);
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
                    fprintf(fT, "%le %le\n", R, T);
                }
            }
        }



        // максимальные погрешности будут в конце файла
        if(globalStepNumber == globalStepNumber1)
        {
            fprintf(fres_maxpogr, "r_pogr_max_abs = %le\n", r_pogr_max_abs/sigma_r_max*100);
            fprintf(fres_maxpogr, "fi_pogr_max_abs = %le\n", fi_pogr_max_abs/sigma_fi_max*100);
            fprintf(fres_maxpogr, "r_pogr_max = %le\n", r_pogr_max*100);
            fprintf(fres_maxpogr, "fi_pogr_max = %le\n", fi_pogr_max*100);
            fprintf(fres_maxpogr, "T_pogr_max = %le\n", T_pogr_max*100);
            fprintf(fres_maxpogr, "T_pogr_max_abs = %le\n", T_pogr_max_abs);
        }
        if(globalStepNumber == globalStepNumber0)
        {
            fprintf(fres_maxpogr, "r_0_pogr_max_abs = %le\n", r_0_pogr_max_abs/sigma_r_0_max*100);
            fprintf(fres_maxpogr, "fi_0_pogr_max_abs = %le\n", fi_0_pogr_max_abs/sigma_fi_0_max*100);
            fprintf(fres_maxpogr, "r_0_pogr_max = %le\n", r_0_pogr_max*100);
            fprintf(fres_maxpogr, "fi_0_pogr_max = %le\n", fi_0_pogr_max*100);
        }

        //print("\nglobalStepNumber = %d outed\n", ARGS(globalStepNumber));
    }

    fclose(fT);
    fclose(fT_ch);
    fclose(fT_pogr);
    fclose(fT_pogr_abs);

    fclose(fu);
    fclose(fu_ch);
    fclose(fu_pogr_abs);
    fclose(fu_pogr);

    fclose(fsigma_r);
    fclose(fsigma_r_ch);
    fclose(fsigma_fi);
    fclose(fsigma_fi_ch);
    fclose(fsigma_r_0);
    fclose(fsigma_r_ch_0);
    fclose(fsigma_fi_0);
    fclose(fsigma_fi_ch_0);

    //fclose(fs_t);
    //fclose(fs_dstep);
    fclose(fs_iterationsNumber);
    fclose(fs_slauResidual);
    fclose(fs_stepSolvingTime);
    fclose(fs_nonlinearStateFENumber);
    fclose(fs_slauNev);
    fclose(fs_epsNev);
    fclose(fs_sigmaNev);

    fclose(fsigma_r_pogr_abs);
    fclose(fsigma_r_pogr_0_abs);
    fclose(fsigma_fi_pogr_abs);
    fclose(fsigma_fi_pogr_0_abs);
    fclose(fsigma_r_pogr);
    fclose(fsigma_r_pogr_0);
    fclose(fsigma_fi_pogr);
    fclose(fsigma_fi_pogr_0);
    fclose(fres_maxpogr);

}

#endif


std::vector<MechMaterialSource> *material;                      // материалы
std::vector<MechBoundaryCondition1Source> *bc1Source;           // 1-e краевые условия
std::vector<MechBoundaryCondition2Source_base *> *bc2Source;    // 2-е краевые
std::vector<Grid::MechBoundaryCondition2> *bc2;                 // 2-е краевые: индекс КЭ-поверхности + индекс MechBoundaryCondition2Source_base
ms.material = material;
ms.bc1Source = bc1Source;
//ms.bc2Source = bc2Source; - будет меняться
ms.bc2 = bc2;




// Сплайн
if(plasticPlot == 1)
{
    mechMat.type = MechMaterialType::Plasticity;
    mechMat.plasticityCurveMode = MechMaterialCurveType::Sigma_eps;
    mechMat.curveDependenceType = MechMaterialCurveDependenceType::None;
    int N = 10;
    // построение сплайнов
    double p1 = 0;
    double p2 = 0.01; // strToDouble(m0.vals[1]) - предел упругости
    std::vector<Elementary::VECTOR2> points;
    for(int i = 0; i < N; i++)
    {
        double eps = p1+(p2-p1)*i/N;
        double sigma = 3*mechMat.G*eps;
        if(sigma > mechMat.elasticSigmaLimit)
        {
            sigma = mechMat.elasticSigmaLimit + 3*mechMat.G*(eps - mechMat.elasticSigmaLimit/(3*mechMat.G))/1000;
        }
        VECTOR2 t;
        t.x[0] = eps;
        t.x[1] = sigma;
        points.push_back(t);
        if(eps > mechMat.elasticSigmaLimit/(3*mechMat.G)) break;
    }
    mechMat.epsFun.setSpline(points);
    // вывод кривой в файл
    FILE *f = fopen("___plastic.txt", "w");
    //p1 = 0;
    //p2 = 0.01; // strToDouble(m0.vals[1]) - предел упругости
    N = 1000;
    for(int i = 0; i < N; i++)
    {
        double eps = p1+(p2-p1)*i/N;
        double sigma = mechMat.eps(eps, 0, 0);
        fprintf(f, "%le %le\n", eps, sigma);
    }
    fclose(f);
}








// линейный спуск с горки
/*
if(epsEqv < epsEqv00)
    return 3*G*epsEqv;
else
    return 3*G*epsEqv00 - 3*G*(epsEqv - epsEqv00)*k00;
*/
// нелинейный подъём их краудерера
//return 3*G*epsEqv/(1-C*3*G*epsEqv);
// нелинейный спуск с горки
//return 3*G*epsEqv/(1+A*epsEqv*epsEqv);
// ползучесть
double f = pow(epsEqv/(c_A*t), 1./c_n);
double lin_f = 3*G*epsEqv;
if(lin_f < f)
    return lin_f;
else
    return f;






// пластичность
{
    double lin_f = 3*G*epsEqv;
    if(lin_f < elasticSigmaLimit)
    {
        return lin_f;
    }
    else
    {
        return elasticSigmaLimit;
    }
}
// ползучесть с линейным участком в начале
{
    double f = pow(epsEqv/(c_A*t), 1./c_n);
    double lin_f = 3*G*epsEqv;
    if(lin_f < f)
        return lin_f;
    else
        return f;
};



// линейный спуск с горки
/*
if(epsEqv < epsEqv00)
    return 1/(3*G);
else
    return 1/(-3*G*k00);
*/


// нелинейный подъём их краудерера
//double df = 3*G*(-3*C*G*epsEqv+1)-9*G*G*epsEqv*C;
// нелинейный спуск с горки
/*
double df = (3*G/(1+A*epsEqv*epsEqv) - 6*G*A*epsEqv*epsEqv/pow(1 + A*epsEqv*epsEqv, 2));
if(fabs(df) < df_min)
    df = (df_min) * SIGN(df);
return 1/df;
*/
// ползучесть
double df = 1 / (pow(epsEqv/(c_A*t), 1./c_n) / (c_n * epsEqv));
double dlin_f = 1 / (3*G);

double f = pow(epsEqv/(c_A*t), 1./c_n);
double lin_f = 3*G*epsEqv;
if(lin_f < f || t == 0 || epsEqv == 0)
    return dlin_f;
else
    return df;






FunParser::Function &fun = mechMat.epsFun;
char curveExpression[1000];

sprintf(curveExpression, "A*pow(sigma,n)*pow(t,m)\0");

fun.setExpression(std::string(curveExpression));
fun.args.clear();
fun.addArgument("sigma\0");  // 0
fun.addArgument("t\0");      // 1
fun.addArgument("A\0");      // 2
fun.addArgument("n\0");      // 3
fun.addArgument("m\0");      // 4
fun.setArgumentValue(2, A);
fun.setArgumentValue(3, n);
fun.setArgumentValue(4, m);
fun.parse();





// отладка
/*
fprintf(fcreep_debug, "iter = %d\n", nlInfCreep.iterNumber);
fprintf(fcreep_debug, "depsNonTerm = %le\n", m0.epsEqv(pd.depsNonTerm));
for(int i = 0; i < 6; i++)
{
    fprintf(fcreep_debug, "%.3le ", pd.depsNonTerm[i]);
}
fprintf(fcreep_debug, "\n");
fprintf(fcreep_debug, "creepEps = %le\n", m0.epsEqv(pd.itCreepIEpsData.depsCreep));
for(int i = 0; i < 6; i++)
{
    fprintf(fcreep_debug, "%.3le ", pd.itCreepIEpsData.depsCreep[i]);
}
fprintf(fcreep_debug, "\n");
*/







// производные перемещений
for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
        du[i][j] = set_du[i][j];
// расчет матрицы Яумена
/*
if(movingGridMode == 1)
{
    double dw01 = (du[1][0] - du[0][1])/2.;
    double dw12 = (du[2][1] - du[1][2])/2.;
    double dw20 = (du[0][2] - du[2][0])/2.;
    dw =
    {
        {
        {1,     dw01,   -dw20},
        {-dw01, 1,      dw12},
        {dw20,  -dw12,  1},
        }
    };
}
*/
// добавляются напряжения
// просто прибавляются
if(movingGridMode == 0)
{
    if(m0.plasticityMethodType == MechPlasticityMethodType::CreepInitialSigma)
        sumSigma += dsigma + itCreepISigmaData.dsigmaCreep;
    else
        sumSigma += dsigma;
}
// приращения Яумана##
if(movingGridMode == 1)
{
    /*MATR3x3 dsigma = VECTOR6toMATR3x3(dsigma) +
            VECTOR6toMATR3x3(sumSigma)*dw.transpose() +
            dw*VECTOR6toMATR3x3(sumSigma);
      sumSigma += MATR3x3toVECTOR6(dsigma);
            */
    sumSigma += dsigma;
    sumSigma = MATR3x3toVECTOR6(
                (dw*VECTOR6toMATR3x3(sumSigma))*dw.transpose());
}








// метод начальных напряжений
// расчёт ползучих деформаций
if(m0.plasticityMethodType == MechPlasticityMethodType::CreepInitialSigma)
{
    using namespace Elementary::Operations;
    //itCreepISigmaData.dsigmaCreep -> itCreepISigmaData.depsCreep
    MATR3x3x3x3 A;
    m0.solveInvertIsotropicD(A);
    MATR3x3 dsigmaCreep_tensor;
    itCreepISigmaData.dsigmaCreep.ToMATR3x3_sigma(dsigmaCreep_tensor);
    MATR3x3 depsCreep_tensor;
    depsCreep_tensor.clear();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    depsCreep_tensor.m[i][j] -= A.m[i][j][k][l]*dsigmaCreep_tensor.m[k][l];
    //itCreepISigmaData.depsCreep = dsigmaCreep_tensor;
    itCreepISigmaData.depsCreep[0] = depsCreep_tensor.m[0][0];
    itCreepISigmaData.depsCreep[1] = depsCreep_tensor.m[1][1];
    itCreepISigmaData.depsCreep[2] = depsCreep_tensor.m[2][2];
    itCreepISigmaData.depsCreep[3] = 2*depsCreep_tensor.m[1][2];
    itCreepISigmaData.depsCreep[4] = 2*depsCreep_tensor.m[2][0];
    itCreepISigmaData.depsCreep[5] = 2*depsCreep_tensor.m[0][1];
}














// метод начальных напряжений
// расчёт ползучих деформаций
if(m0.plasticityMethodType == MechPlasticityMethodType::CreepInitialSigma)
{
    using namespace Elementary::Operations;
    //itCreepISigmaData.dsigmaCreep -> itCreepISigmaData.depsCreep
    MATR3x3x3x3 A;
    m0.solveInvertIsotropicD(A);
    MATR3x3 dsigmaCreep_tensor;
    itCreepISigmaData.dsigmaCreep.ToMATR3x3_sigma(dsigmaCreep_tensor);
    MATR3x3 depsCreep_tensor;
    depsCreep_tensor.clear();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    depsCreep_tensor.m[i][j] -= A.m[i][j][k][l]*dsigmaCreep_tensor.m[k][l];
    //itCreepISigmaData.depsCreep = dsigmaCreep_tensor;
    itCreepISigmaData.depsCreep[0] = depsCreep_tensor.m[0][0];
    itCreepISigmaData.depsCreep[1] = depsCreep_tensor.m[1][1];
    itCreepISigmaData.depsCreep[2] = depsCreep_tensor.m[2][2];
    itCreepISigmaData.depsCreep[3] = 2*depsCreep_tensor.m[1][2];
    itCreepISigmaData.depsCreep[4] = 2*depsCreep_tensor.m[2][0];
    itCreepISigmaData.depsCreep[5] = 2*depsCreep_tensor.m[0][1];
}

// вычитаем температурные деформации
depsElastoplastic -= depsTerm;
if(m0.plasticityMethodType == MechPlasticityMethodType::CreepInitialEps)
    depsElastoplastic -= itCreepIEpsData.depsCreep;

// разделяем деформации depsNonTerm на пластичную и упругую составляющие
// расчёт приращения напряжений
using namespace Elementary::Operations;
switch (m0.plasticityMethodType)
{
// линейная упругость
case MechPlasticityMethodType::Elasticity:
// ползучесть методом начальных напряжений
case MechPlasticityMethodType::CreepInitialSigma:
// ползучесть методом начальных деформаций
case MechPlasticityMethodType::CreepInitialEps:
{
    // деформации
    // нет пластичнеских деформаций
    depsPlastic = {0,0,0,0,0,0};
    depsElastic = depsElastoplastic;
}break;
// пластичность
case MechPlasticityMethodType::Plasticity:
{
    // деформации
    // даны tan1, tan2Want
    if(itData.E_pl_is_inf)
    {
        // нет пластичнеских деформаций
        depsPlastic = {0,0,0,0,0,0};
        depsElastic = depsElastoplastic;
    }
    else
    {
        double sigma2Eqv;
        VECTOR6 sigma2;
        V1plusV26(sumSigma, dsigma, sigma2);
        // исходное положение - sigma1Eqv
        // положение после шага
        sigma2Eqv = m0.sigmaEqv(sigma2);
        // приращение эквивалентного напряжения
        double dSigmaEqv = sigma2Eqv - sigma1Eqv;
        // приращение эквивалентной пластической деформации
        double dEpsPlEqv = dSigmaEqv / itData.E_pl;
        VECTOR6 Z;
        MmulV6(m0.M_sigma, sumSigma, Z);  // M*sigma0 -> Z
        Vmulc6(Z, (2. / 3.) * sqrt(m0.Csigma/m0.Ceps) * dEpsPlEqv / sigma1Eqv, Z);    // M*sigma0*dEpsPlEqv/sigma1Eqv -> Z
        depsPlastic = Z;
        depsElastic = depsElastoplastic - depsPlastic;
    }
}break;
}
// напряжения
MmulV6(D_pl, depsElastoplastic, dsigma);    // D*depsNonTerm -> sigma





itData_new = itData;
using namespace Elementary::Operations;
double fx1, fy1, fx2, fy2;
double dfdx;
VECTOR6 eps2;
VECTOR6 sigma2;
V1plusV26(sumEpsElastoplastic, depsElastoplastic, eps2);
V1plusV26(sumSigma, dsigma, sigma2);
// исходное положение - eps1Eqv, sigma1Eqv
// реальное положение после шага
double eps2EqvReal_new = m0.epsEqv(eps2);
double sigma2EqvReal_new = m0.sigmaEqv(sigma2);
// желаемое положение после шага
double sigma2EqvWant_new;
double eps2EqvWant_new;
double dif1;
double dif2;
if(m0.PCType == MechPlasticityCurveType::Eps_sigma)
{
    sigma2EqvWant_new = m0.sigmaEqv(sigma2);
    eps2EqvWant_new = m0.eps(sigmaT, sigma2EqvWant_new, tl.t);
    dif1 = m0.difEps(sigmaT, sigma1Eqv, tl.t - tl.dt);
    dif2 = m0.difEps(sigmaT, sigma2EqvWant_new, tl.t);
}
if(m0.PCType == MechPlasticityCurveType::Sigma_eps)
{
    eps2EqvWant_new = m0.epsEqv(eps2);
    sigma2EqvWant_new = m0.sigma(epsT, eps2EqvWant_new, tl.t);
    dif1 = m0.difSigma(epsT, eps1Eqv, tl.t - tl.dt);
    dif2 = m0.difSigma(epsT, eps2EqvWant_new, tl.t);
}
// находимся на линейном участке?
if(dif1 == dif2)
{
    // в начале разгрузки difSigma(eps1Eqv, pd.epsT) даёт не ту производную!!
    //if(plasticIterNumber == 0)  // вдруг разгрузка? тогда на нулевой итерации можно сказать что находимся на линейном участке
    if(iterNumber == 0)  //###
        itData_new.needIterations = false;//true;//false;
    else
        itData_new.needIterations = true;
}
else
{
    itData_new.needIterations = true;
};






// подсчет точности
if(iterNumber >= 1)
{
    // погрешность - разница между реальным положением после шага и тем, которое хотели получить
    plasticResidual.eps = (eps2EqvReal_new - itData_new.eps2EqvWant)/((fabs(itData_new.eps2EqvWant) + fabs(eps1Eqv))/2);
    plasticResidual.sigma = (sigma2EqvReal_new - itData_new.sigma2EqvWant)/((fabs(itData_new.sigma2EqvWant) + fabs(sigma1Eqv))/2);
    if(m0.usingSecantMethod)
    {
        // метод секущих (приближенный метод Ньютона)
        if(m0.PCType == MechPlasticityCurveType::Eps_sigma)
        {
            fx2 = itData_new.sigma2EqvWant;
            fy2 = sigma2EqvWant_new;   // последняя пара значений
            if(iterNumber >= 2 /*&& !elastic*/)//&& deps > parameters.nonlinear_eps/10
            {
                fx1 = itData_new.fx; // предыдущая пара значений
                fy1 = itData_new.fy;
                dfdx = (fy2 - fy1) / (fx2 - fx1);   // ~производная
                double w = -1. / (dfdx - 1);
                if(w >= -1 && w <= 2)
                    sigma2EqvWant_new = (1-w)*itData_new.sigma2EqvWant + w*sigma2EqvWant_new;
            }
            itData_new.fx = fx2;
            itData_new.fy = fy2;
        }
        if(m0.PCType == MechPlasticityCurveType::Sigma_eps)
        {
            fx2 = itData_new.eps2EqvWant;
            fy2 = eps2EqvWant_new;   // последняя пара значений
            if(iterNumber >= 2 /*&& !elastic*/)
            {
                fx1 = itData_new.fx; // предыдущая пара значений
                fy1 = itData_new.fy;
                dfdx = (fy2 - fy1) / (fx2 - fx1);   // ~производная
                double w = -1. / (dfdx - 1);
                if(w >= -1 && w <= 2)
                    eps2EqvWant_new = (1-w)*itData_new.eps2EqvWant + w*eps2EqvWant_new;
            }
            itData_new.fx = fx2;
            itData_new.fy = fy2;
        }
        // направляем в соответствии с кривой деформирования
        if(m0.PCType == MechPlasticityCurveType::Eps_sigma)
        {
            eps2EqvWant_new = m0.eps(sigmaT, sigma2EqvWant_new, tl.t);
        }
        if(m0.PCType == MechPlasticityCurveType::Sigma_eps)
        {
            sigma2EqvWant_new = m0.sigma(epsT, eps2EqvWant_new, tl.t);
        }
    }
}
else
{
    plasticResidual.set_null();
}
itData_new.sigma2EqvWant = sigma2EqvWant_new;
itData_new.eps2EqvWant = eps2EqvWant_new;
itData_new.eps2EqvReal = eps2EqvReal_new;
itData_new.sigma2EqvReal = sigma2EqvReal_new;
//updateTan2WantReal(m0);

double delta_eps = itData_new.eps2EqvWant - eps1Eqv;
double delta_sigma = itData_new.sigma2EqvWant - sigma1Eqv;
if(delta_eps == 0 && delta_sigma == 0)
{
    // отсутствует движение по кривой - пассивное нагружение
    //itData_new.tan2Want = tan1;
    itData_new.E_pl_is_inf = true;
}
if(delta_eps != 0 && delta_sigma == 0)
{
    // максимальная пластичность
    //itData_new.tan2Want = бесконечность
    itData_new.tan2Want = 1.e100;
    itData_new.E_pl = 0;
    itData_new.E_pl_is_inf = false;
    //printf("E_pl = 0!\n");
}
if(delta_eps != 0 && delta_sigma != 0)
{
    // движение по кривой
    double dif1 = m0.difSigma(epsT, eps1Eqv, tl.t - tl.dt);
    double dif2 = m0.difSigma(epsT, eps2EqvWant_new, tl.t);
    if(dif1 == dif2 && dif1 == tan1)
    {
        // находимся на линейном участке
        //itData_new.tan2Want = tan1;
        itData_new.E_pl_is_inf = true;
    }
    else
    {
        itData_new.tan2Want = delta_eps/delta_sigma;
        //itData_new.E_pl = 1. / (delta_eps/delta_sigma - tan1);
        itData_new.E_pl = delta_sigma / (delta_eps - tan1*delta_sigma);
        //itData_new.E_pl = 1. / (delta_eps/delta_sigma - tan1);
        itData_new.E_pl_is_inf = false;
    }
}
/*
if(itData_new.eps2EqvWant - eps1Eqv == 0 || itData_new.sigma2EqvWant - sigma1Eqv == 0)    //##
    itData_new.tan2Want = tan1;  // В этом случае D_pl = D, т.к. тогда пластических деформаций нет
else
    itData_new.tan2Want = (itData_new.eps2EqvWant - eps1Eqv)/(itData_new.sigma2EqvWant - sigma1Eqv);

if(itData_new.eps2EqvReal - eps1Eqv == 0 || itData_new.sigma2EqvReal - sigma1Eqv == 0)
    itData_new.tan2Real = tan1;
else
    itData_new.tan2Real = (itData_new.eps2EqvReal - eps1Eqv)/(itData_new.sigma2EqvReal - sigma1Eqv);

itData_new.E_pl_is_inf = (itData_new.tan2Want == tan1);
if(!itData_new.E_pl_is_inf)
{
    itData_new.E_pl = 1. / (itData_new.tan2Want - tan1);
}
*/
//if(!itData_new.E_pl_is_inf && itData_new.E_pl < 0)
//    printf("E_pl = %le, delta_sigma = %le, delta_eps/delta_sigma - tan1 = %le\n", itData_new.E_pl, delta_sigma, delta_eps/delta_sigma - tan1);
    //printf("E_pl = %le, dsigma = %le, deps = %le\n", itData_new.E_pl, delta_sigma, delta_eps);













itCreepISigmaData_new = itCreepISigmaData;
using namespace Elementary::Operations;
// напряжение и деформации, полученные в результате упругого решения с начальными напряжениями
VECTOR6 sigma1;
VECTOR6 sigma2_elastic;
VECTOR6 eps1;
VECTOR6 eps2;
sigma1 = sumSigma;
V1plusV26(sigma1, dsigma, sigma2_elastic);
eps1 = sumEpsElastoplastic;
V1plusV26(eps1, depsElastoplastic, eps2);
// соответствующие эквивалентные напряжения и деформации
double sigma1Eqv = m0.sigmaEqv(sigma1);
double sigma2Eqv_elastic = m0.sigmaEqv(sigma2_elastic);
double eps1Eqv = m0.epsEqv(eps1);
double eps2Eqv = m0.epsEqv(eps2);
double sigma2Eqv_real;
double sigma2Eqv_want;
if(sigma1Eqv != 0)
{
    // полученные напряжения с учётом начальных напряжений
    VECTOR6 sigma2_real = sigma2_elastic + itCreepISigmaData.dsigmaCreep;
    sigma2Eqv_real = m0.sigmaEqv(sigma2_real);
    // направление девиатора напряжений
    //VECTOR6 sigma1_deviator_unit = sigma2_real;
    VECTOR6 sigma1_deviator_unit = sigma1;
    //VECTOR6 sigma1_deviator_unit = sigma2;
    //VECTOR6 sigma1_deviator_unit = (sigma1 + sigma2)/2;
    double sigma1_deviator_unitEqv = m0.sigmaEqv(sigma1_deviator_unit);
    double sigma1_0 = (sigma1_deviator_unit[0] + sigma1_deviator_unit[1] + sigma1_deviator_unit[2])/3;
    sigma1_deviator_unit[0] -= sigma1_0;
    sigma1_deviator_unit[1] -= sigma1_0;
    sigma1_deviator_unit[2] -= sigma1_0;
    sigma1_deviator_unit /= sigma1_deviator_unitEqv;//m0.sigmaEqv(sigma1_deviator_unit);
    // желаемое эквивалентное напряжение, которое соответствует(исходя из кривой ползучести) полученным деформациям
    sigma2Eqv_want = m0.sigma(0, eps2Eqv, tl.t);
    // приращение ползучего напряжения

    // без сглаживания
    double new_dsigmaCreepEqv = -m0.sigmaEqv(itCreepISigmaData.dsigmaCreep) + (sigma2Eqv_want - sigma2Eqv_real);

    // сглаживание
    /*
    //double w_add = 1;
    double w_add = 0.1;
    double new_dsigmaCreepEqv = -m0.sigmaEqv(itCreepISigmaData.dsigmaCreep) + (sigma2Eqv_want - sigma2Eqv_real)*w_add;
    */

    //double new_dsigmaCreepEqv = (sigma2Eqv_want - sigma2Eqv_real);
    // новое приближение
    VECTOR6 new_dsigmaCreep;
    Vmulc6(sigma1_deviator_unit, new_dsigmaCreepEqv, new_dsigmaCreep);
    itCreepISigmaData_new.dsigmaCreep = new_dsigmaCreep;

    // сглаживание
    /*
    //double w = 0.9;
    //double w = 0.5;
    double w = 0.9;
    itCreepISigmaData_new.dsigmaCreep = itCreepISigmaData.dsigmaCreep*w + itCreepISigmaData_new.dsigmaCreep*(1-w);
    */
}
else
{
    itCreepISigmaData_new.dsigmaCreep.clear();
}
itCreepISigmaData_new.eps2Eqv = eps2Eqv;

// расчет точности
if(iterNumber >= 1)
{
    // напряжения
    {
        // погрешность - изменение приращения эквивалентной деформации ползучести по сравнению с предыдущеё итерацией
        double dsigmaCreepEqv1 = m0.sigmaEqv(itCreepISigmaData.dsigmaCreep);
        double dsigmaCreepEqv2 = m0.sigmaEqv(itCreepISigmaData_new.dsigmaCreep);
        //double sum = (dsigmaCreepEqv1 + dsigmaCreepEqv2)/2;
        double sum = sigma1Eqv;
        //double sum = fabs(sigma2Eqv_elastic - sigma1Eqv); //+ (fabs(dsigmaCreepEqv1) + fabs(dsigmaCreepEqv2))/2;
        //double sum = (fabs(sigma1Eqv) + fabs(sigma2Eqv_elastic))/2 + (fabs(dsigmaCreepEqv1) + fabs(dsigmaCreepEqv2))/2;
        //double sum = (fabs(dsigmaCreepEqv1) + fabs(dsigmaCreepEqv2))/2;
        //double sum = (fabs(sigma1Eqv) + fabs(sigma2Eqv_elastic))/2 + (fabs(dsigmaCreepEqv1) + fabs(dsigmaCreepEqv2))/2;
        //double sum = fabs(sigma1Eqv) + (fabs(dsigmaCreepEqv1) + fabs(dsigmaCreepEqv2))/2;
        if(sum != 0)
        {
            //creepISigmaResidualLast.sigma = (dsigmaCreepEqv2 - dsigmaCreepEqv1) / sum;
            creepISigmaResidual.sigma = (sigma2Eqv_real - sigma2Eqv_want) / sum;
        }
        else
        {
            creepISigmaResidual.sigma = 0;
        }
    }
    // деформации
    {
        double deps1 = fabs(itCreepISigmaData.eps2Eqv - eps1Eqv);
        double deps2 = fabs(itCreepISigmaData_new.eps2Eqv - eps1Eqv);
        //double sum = (deps1 + deps2)/2;
        double sum = eps1Eqv;
        if(sum != 0)
        {
            //creepISigmaResidualLast.eps = (deps2 - deps1) / sum;
            creepISigmaResidual.eps = (itCreepISigmaData_new.eps2Eqv - itCreepISigmaData.eps2Eqv) / sum;
        }
        else
        {
            creepISigmaResidual.eps = 0;
        }
    }
}
else
{
    creepISigmaResidual.set_undefined();
}
itCreepISigmaData_new.needIterations = true;












itCreepIEpsData_new = itCreepIEpsData;
using namespace Elementary::Operations;
VECTOR6 sigma1;
VECTOR6 sigma2;
VECTOR6 sigma_srednee;
sigma1 = sumSigma;
V1plusV26(sigma1, dsigma, sigma2);
double w_sigma = 1;
//double w_sigma = 0;
//double w_sigma = 0.5;
//double w_sigma = 0.9;
sigma_srednee = sigma1*w_sigma + sigma2*(1-w_sigma);
//sigma_srednee = (sigma1 + sigma2)/2;
//sigma_srednee = sigma2;
//sigma_srednee = sigma1;
double sigma1Eqv = m0.sigmaEqv(sigma1);
double sigma2Eqv = m0.sigmaEqv(sigma2);
double sigma_srednee_Eqv = m0.sigmaEqv(sigma_srednee);
if(sigma_srednee_Eqv == 0)
//if(sigma1Eqv == 0)
//if(sigma2Eqv == 0)
{
    itCreepIEpsData_new.depsCreep.clear();   // 0 -> creepEps
    //Vmulc6(itCreepData_new.creepEps, 0, itCreepData_new.creepEps);  // M*0 / sigma2Eqv -> creepEps
}
else
{
    MmulV6(m0.M_sigma, sigma_srednee, itCreepIEpsData_new.depsCreep);                        // M*sigma_srednee -> creepEps
    Vmulc6(itCreepIEpsData_new.depsCreep, 1. / sigma_srednee_Eqv, itCreepIEpsData_new.depsCreep);  // M*sigma_srednee / sigma_srednee_Eqv -> creepEps
    //MmulV6(m0.M_sigma, sigma1, itCreepData_new.depsCreep);                        // M*sigma1 -> creepEps
    //Vmulc6(itCreepData_new.depsCreep, 1. / sigma1Eqv, itCreepData_new.depsCreep);  // M*sigma1 / sigma1Eqv -> creepEps
    //MmulV6(m0.M_sigma, sigma2, itCreepData_new.creepEps);                        // M*sigma2 -> creepEps
    //Vmulc6(itCreepData_new.creepEps, 1. / sigma2Eqv, itCreepData_new.creepEps);  // M*sigma2 / sigma2Eqv -> creepEps
}
double dF;
dF = m0.eps(0, sigma2Eqv, tl.t) -
     m0.eps(0, sigma1Eqv, tl.t - tl.dt);

// новое приближение
Vmulc6(itCreepIEpsData_new.depsCreep, dF, itCreepIEpsData_new.depsCreep);  // M*sigma_srednee / sigma_srednee_Eqv * dF -> creepEps



// сглаживание
double w_eps = 0;
//double w_eps = 0.5;
//double w_eps = 0.9;
VECTOR6 new_depsCreep = itCreepIEpsData.depsCreep*w_eps + itCreepIEpsData_new.depsCreep*(1-w_eps);
itCreepIEpsData_new.depsCreep = new_depsCreep;



// отладка
/*
fprintf(fcreep_debug, "sigma1Eqv = %le\n", sigma1Eqv);
for(int i = 0; i < 6; i++)
{
    fprintf(fcreep_debug, "%.3le ", sigma1[i]);
}
fprintf(fcreep_debug, "\n");
fprintf(fcreep_debug, "sigma2Eqv = %le\n", sigma2Eqv);
for(int i = 0; i < 6; i++)
{
    fprintf(fcreep_debug, "%.3le ", sigma2[i]);
}
fprintf(fcreep_debug, "\n");
fprintf(fcreep_debug, "dF = %le\n", dF);
fprintf(fcreep_debug, "creepEps_new = %le\n", m0.epsEqv(itCreepIEpsData_new.depsCreep));
for(int i = 0; i < 6; i++)
{
    fprintf(fcreep_debug, "%.3le ", itCreepIEpsData_new.depsCreep[i]);
}
fprintf(fcreep_debug, "\n\n");
*/

itCreepIEpsData_new.sigma2Eqv = sigma2Eqv;
// подсчет точности
if(iterNumber >= 1)
{
    // деформации
    {
        // погрешность - изменение приращения эквивалентной деформации ползучести по сравнению с предыдущеё итерацией
        double dEpsCreepEqv1 = m0.epsEqv(itCreepIEpsData.depsCreep);
        double dEpsCreepEqv2 = m0.epsEqv(itCreepIEpsData_new.depsCreep);
        //double sum = m0.epsEqv(sumEpsNonTerm);
        //double sum = m0.epsEqv(sumEpsNonTerm+depsNonTerm);
        //double sum = (fabs(dEpsCreepEqv1) + fabs(dEpsCreepEqv2))/2 +
        //        1*m0.epsEqv(sumEpsNonTerm);
        double sum = (fabs(dEpsCreepEqv1) + fabs(dEpsCreepEqv2))/2;
        if(sum != 0)
        {
            creepIEpsResidual.eps = (dEpsCreepEqv2 - dEpsCreepEqv1) / sum;
        }
        else
        {
            creepIEpsResidual.eps = 0;
        }
    }
    // напряжения
    {
        double sum = sigma1Eqv;
        if(sum != 0)
        {
            creepIEpsResidual.sigma = (itCreepIEpsData_new.sigma2Eqv - itCreepIEpsData.sigma2Eqv) / sum;
        }
        else
        {
            creepIEpsResidual.sigma = 0;
        }
    }

}
else
{
    creepIEpsResidual.set_undefined();
}
itCreepIEpsData_new.needIterations = true;

//Vmulc6(itCreepData_new.creepEps, 0, itCreepData_new.creepEps);
//itCreepData_new.needIterations = false;





















// 3D интерполянт (старая версия, новой версии пока нет)

// тип базисных функций при интерполировании
enum class InterpolatorType
{
    Lagrange3D1,  // билинейные Лагранжевы функции
    Hermite2D3,   // бикубические Эрмитовы базисные функции
};

// точка и скаляр
struct POINT3_VALUE
{
    POINT3 p;       // точка
    double value;   // значение скаляра в этой точке
    int fei;        // индекс конечного элемента, содержащего точку p
};
/*
struct POINT2_VALUE
{
    POINT2 p;       // точка
    double value;   // значение скаляра в этой точке
    int fei;        // индекс конечного элемента, содержащего точку p
};
struct POINT1_VALUE
{
    POINT1 p;       // точка
    double value;   // значение скаляра в этой точке
    int fei;        // индекс конечного элемента, содержащего точку p
};*/

// круг-заплатка
struct Circle
{
    double x, y;    // центр круга (z = 0)
    double r;       // радиус
    int mode;       // 0 - затирание вне окружности
                    // 1 - затирание внутри окружности
} ;
        /*
        // построитель интерполянта по заданному набору точек и сетке
        class Interpolator
        {
        public:
            //void init(const GridHexagone *set_grid, const InterpolationType set_type, const double set_alpha, const double set_betta, const int set_max_size);
            void init(const Grid::__GridRectangleRegular3D *set_grid, const InterpolatorType set_type, const double set_alpha, const double set_betta, const int set_max_size);
            void release();
            //void addPoint(const POINT1 &p, const double value, const int fe_index);
            //void addPoint(const POINT2 &p, const double value, const int fe_index);
            void addPoint(const POINT3 &p, const double value, const int fe_index);
            void makeInterpolant();
            void saveInterpolant(const double z0, const char *file_name, const int bmp_x, const int bmp_y, const int colorMode, const int pictureMode, const std::vector<Interpolation::Circle> circle);
            //double fun(const POINT1 &p)const;
            //double fun(const POINT2 &p)const;
            double fun(const POINT3 &p)const;
        private:
            void (Interpolator ::* make_ip_ptr)() = nullptr;
            //double (Interpolator ::* fun1_ptr)(const POINT1 &p)const = nullptr;
            //double (Interpolator ::* fun2_ptr)(const POINT2 &p)const = nullptr;
            double (Interpolator ::* fun3_ptr)(const POINT3 &p)const = nullptr;
            void (Interpolator ::* addPoint_ptr)(const POINT3 &p, const double value, const int fe_index) = nullptr;
            void addPoint_RectangleRegular3D_Lagrange1(const POINT3 &p, const double value, const int fe_index);
            void make_ip_RectangleRegular3D_Lagrange1();
            double fun_RectangleRegular3D_Lagrange1(const POINT3 &p)const;
            void addPoint_RectangleRegular3D_Hermite3(const POINT3 &p, const double value, const int fe_index);
            void make_ip_RectangleRegular3D_Hermite3();
            double fun_RectangleRegular3D_Hermite3(const POINT3 &p)const;
            const void *grid;   // указатель на сетку
            double alpha;       // параметр влияющий на точность интерполяции. Чем меньше, тем более точный, но тем больше точек необходимо для построения
            double betta;
            void *m = nullptr;      // точки и значения скаляров в каждой точке
            double *res;            // решение СЛАУ (после вызова make_ip)
            int size = 0;           // размер массива m (будет сортироваться и хранить только точки внутри области построения интерполянта)
            int max_size = 0;       // выделенная память для массива m
            int matrixSize = 0;     // размер матрицы СЛАУ
            double l[3];            // размерности прямоугольной регулярной сетки
            double h[3];            // размерности ячеек прямоугольной регулярной сетки
        };
        */
//#define RGB32(a, r, g, b)   ((b)+((g)<<8)+((r)<<16)+((a)<<24))			//возвращает индекс цвета (a, red,green,blue), каждая компонента <=2^8
//#define X_BMP(x_)	((int)(((x_) - grid3D->p1[0])/(grid3D->p2[0] - grid3D->p1[0])*bmp.x))
//#define Y_BMP(y_)	((int)(((y_) - grid3D->p1[1])/(grid3D->p2[1] - grid3D->p1[1])*bmp.y))
//#define X_REAL(x_)	(grid3D->p1[0]+(((double)(x_))/bmp.x)*(grid3D->p2[0] - grid3D->p1[0]))
#define Y_REAL(y_)	(grid3D->p1[1]+(((double)(y_))/bmp.y)*(grid3D->p2[1] - grid3D->p1[1]))
#define XI(i)	(2 * (((i) / 4) % 2) + ((i) % 2))
#define YI(i)	(2 * ((i) / 8) + (((i) / 2) % 2))


void Interpolator::init(const Grid::__GridRectangleRegular3D *set_grid, const InterpolatorType set_type, const double set_alpha, const double set_betta, const int set_max_size)
{
    release();
    grid = set_grid;
    alpha = set_alpha;
    betta = set_betta;
    max_size = set_max_size;
    switch (set_type)
    {
    case InterpolatorType::Lagrange3D1:
        matrixSize = (set_grid->N[0] + 1)*(set_grid->N[1] + 1)*(set_grid->N[2] + 1);
        m = new POINT3_VALUE[max_size];
        res = new double[matrixSize];
        addPoint_ptr = &Interpolator::addPoint_RectangleRegular3D_Lagrange1;
        make_ip_ptr = &Interpolator::make_ip_RectangleRegular3D_Lagrange1;
        fun3_ptr = &Interpolator::fun_RectangleRegular3D_Lagrange1;
        break;
    case InterpolatorType::Hermite2D3:
        matrixSize = (set_grid->N[0] + 1)*(set_grid->N[1] + 1)*4;
        m = new POINT3_VALUE[max_size];
        res = new double[matrixSize];
        addPoint_ptr = &Interpolator::addPoint_RectangleRegular3D_Hermite3;
        make_ip_ptr = &Interpolator::make_ip_RectangleRegular3D_Hermite3;
        fun3_ptr = &Interpolator::fun_RectangleRegular3D_Hermite3;
        break;
    default:
        break;
    }
    Grid::__GridRectangleRegular3D *g = (Grid::__GridRectangleRegular3D *)grid;
    l[0] = g->p2[0] - g->p1[0];
    l[1] = g->p2[1] - g->p1[1];
    l[2] = g->p2[2] - g->p1[2];
    h[0] = l[0] / g->N[0];
    h[1] = l[1] / g->N[1];
    h[2] = l[2] / g->N[2];
}
void Interpolator::release()
{
    if (m != nullptr)
    {
        delete []m;
        delete []res;
        m = nullptr;
        max_size = 0;
    }
    size = 0;
}
void Interpolator::addPoint(const POINT3 &p, const double value, const int fe_index)
{
    (this->*addPoint_ptr)(p, value, fe_index);
}
void Interpolator::makeInterpolant()
{
    (this->*make_ip_ptr)();
}
double Interpolator::fun(const POINT3 &p)const
{
    return (this->*fun3_ptr)(p);
}
void Interpolator::saveInterpolant(const double z0, const char *file_name, const int bmp_x, const int bmp_y, const int colorMode, const int pictureMode, const std::vector<Circle> circle)
{
    Grid::__GridRectangleRegular3D *grid3D = (Grid::__GridRectangleRegular3D *)grid;
    POINT3_VALUE *m3 = (POINT3_VALUE *)m;
    int r,g,b;
    double valueMin = +1.e+200;
    double valueMax = -1.e+200;
    for (int k = 0; k < size; k++)
    {
        double value = m3[k].value;
        if (value < valueMin) valueMin = value;
        if (value > valueMax) valueMax = value;
    }
    double valueNull = (valueMin + valueMax) / 2;
    BMP_INF bmp;		//
    // инициализация картинки
    bmp.x = bmp_x;
    bmp.y = bmp_y;
    bmp.c = new unsigned int[bmp.x * bmp.y];
    // рисование картинки
    for (int y = 0; y < bmp.y; y++)
        for (int x = 0; x < bmp.x; x++)
        {
            POINT3 p;
            p[0] = grid3D->p1[0] + l[0] * x / bmp.x;
            p[1] = grid3D->p1[1] + l[1] * y / bmp.y;
            p[2] = z0;
            //p[2] *= 2;
            solve_rgb(colorMode, fun(p), valueMin, valueNull, valueMax, r, g, b);
            PAINT(x, y, r, g, b);
        }
    // затирание областей внутри или вне набора оркужностей
    // цвет отверстия
    r = 255;
    g = 255;
    b = 255;
    for(size_t i = 0; i < circle.size(); i++)
    {
        for (int y = 0; y < bmp.y; y++)
            for (int x = 0; x < bmp.x; x++)
            {
                bool doit = false;
                double h = sqrt(SQR(circle[i].x - Y_REAL(x)) + SQR(circle[i].y - Y_REAL(y)));
                if(circle[i].mode == 0)
                    if(h < circle[i].r)
                        doit = true;
                if(circle[i].mode == 1)
                    if(h >= circle[i].r)
                        doit = true;
                if(doit)
                {
                    if(colorMode == 1)
                    {
                        int col = (x + y)%2;
                        r = col*255;
                        g = col*255;
                        b = col*255;
                        PAINT(x, y, r, g, b);
                    }
                    if(colorMode == 2)
                        PAINT(x, y, r, g, b);
                    if(colorMode == 3)
                        PAINT(x, y, r, g, b);

                }
            }
    }
    // пупырышки чтобы различать круг среди оттенков серого
    /*
    if(bmp_mode == 1)
    {
    // цвет отверстия
    r = 0;
    g = 0;
    b = 0;
    for(int i = 0; i < size_circle; i++)
    {
        double h, l;
        double x1, x2;
        double y1, y2;
        double x, y;
        y1 = Y_BMP(circle[i].y - circle[i].r);
        y2 = Y_BMP(circle[i].y + circle[i].r);
        for(y = y1; y <= y2; y+=2)
        {
            h = ABS(circle[i].y - Y_REAL(y));
            if(SQR(circle[i].r) - SQR(h) < 0) goto next2;
            l = sqrt(SQR(circle[i].r) - SQR(h));
            x1 = X_BMP(circle[i].x - l);
            x2 = X_BMP(circle[i].x + l);
            for(x = x1; x <= x2; x+=2)
            {
                //if(rand()%5 == 0)
                    PAINT((int)x, (int)y, r, g, b);
            }
            next2:;
        }
    }
    }*/
    // вывод картинки без градаций
    if(pictureMode == 0)
    {
        bmp_save(file_name, bmp);
        delete[] bmp.c;
        return;
    }
    // вывод картинки с градациями
    if(pictureMode == 1)
    {
        BMP_INF symbols;	// изображения символов
        BMP_INF bmp2;		// интерполянт + градации
        bmp_load("_symbols.bmp", symbols);
        // инициализация второй картинки
        bmp2.x = bmp_x+LEN*SIZE_X;
        bmp2.y = bmp_y;
        bmp2.c = new unsigned int[bmp2.x * bmp2.y];
        for(int y = 0; y < bmp.y; y++)
        for(int x = 0; x < bmp.x; x++)
        {
            bmp2.c[y*bmp2.x + x] = bmp.c[y*bmp.x + x];
        }
        for(int y = 0; y < bmp.y; y++)
        for(int x = bmp.x; x < bmp2.x; x++)
        {
            bmp2.c[y*bmp2.x + x] = RGB32(0, 255, 255, 255);;
        }

        for(int y = 0; y < SIZE_Y/2; y++)
        for(int x = bmp.x; x < bmp2.x; x++)
        {
            double value = valueMin;
            solve_rgb(colorMode, value, valueMin, valueNull, valueMax, r, g, b);
            bmp2.c[y*bmp2.x + x] = RGB32(0, r, g, b);;
        }
        for(int y = bmp.y-SIZE_Y/2; y < bmp.y; y++)
        for(int x = bmp.x; x < bmp2.x; x++)
        {
            double value = valueMax;
            solve_rgb(colorMode, value, valueMin, valueNull, valueMax, r, g, b);
            bmp2.c[y*bmp2.x + x] = RGB32(0, r, g, b);;
        }

        for(int y = SIZE_Y/2; y < bmp.y-SIZE_Y/2; y++)
        for(int x = bmp.x; x < bmp2.x; x++)
        {
            double value = valueMin + (valueMax - valueMin)*(y-SIZE_Y / 2)/(bmp.y - SIZE_Y);
            solve_rgb(colorMode, value, valueMin, valueNull, valueMax, r, g, b);
            bmp2.c[y*bmp2.x + x] = RGB32(0, r, g, b);;
        }
        double val = valueMin;
        double dval = (valueMax - valueMin)/K;
        for(int i = 0; i <= K; i++)
        {
            char str[100];
            str[0] = ' ';
            if(val >=0 )
                sprintf(str+1, "%.3le", val);
            else
                sprintf(str, "%.3le", val);
            print(str, bmp.x, SIZE_Y/2 + i*(bmp.y-SIZE_Y)/K, colorMode, symbols, bmp2);
            val += dval;
        }
        //print(0, 100, "-1.234e980", bmp2);
        // сохранение изображение
        bmp_save(file_name, bmp2);
        delete[] symbols.c;
        delete[] bmp2.c;
        delete[] bmp.c;
        return;
    }
}

int POINT3_VALUE_cmp(const void* x1, const void* x2)
{
    POINT3_VALUE *px1 = (POINT3_VALUE *)x1;
    POINT3_VALUE *px2 = (POINT3_VALUE *)x2;
    if (px1->fei < px2->fei) return -1;
    if (px1->fei > px2->fei) return  1;
    return 0;
}
void Interpolator::addPoint_RectangleRegular3D_Lagrange1(const POINT3 &p, const double value, const int fe_index)
{
    POINT3_VALUE *m3 = (POINT3_VALUE *)m;
    Grid::__GridRectangleRegular3D *g = (Grid::__GridRectangleRegular3D *)grid;
    if (fe_index != UNKNOWN_FE_INDEX)
    {
        m3[size] = { p, value, fe_index };
        size++;
    }
    else
    {
        int nk = (int)((p[2] - g->p1[2]) / l[2] * g->N[2]);
        int ni = (int)((p[1] - g->p1[1]) / l[1] * g->N[1]);
        int nj = (int)((p[0] - g->p1[0]) / l[0] * g->N[0]);
        int k = nj + ni * (g->N[0] + 1) + nk * (g->N[1] + 1) * (g->N[0] + 1);
        m3[size].p = p;
        m3[size].value = value;
        m3[size].fei = k;
        size++;
    }

}
void Interpolator::make_ip_RectangleRegular3D_Lagrange1()
{
    Grid::__GridRectangleRegular3D *g = (Grid::__GridRectangleRegular3D *)grid;
    POINT3_VALUE *m3 = (POINT3_VALUE *)m;
    // сортировка точек по параллелепипедам в которые они попали
    qsort(m, size, sizeof(POINT3_VALUE), POINT3_VALUE_cmp);
    double G[8][8];		// локальная матрица жесткости
    // построение локальной матрицы G
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
        {
            int xi = i % 2;
            int yi = (i / 2) % 2;
            int zi = (i / 4);
            int xj = j % 2;
            int yj = (j / 2) % 2;
            int zj = (j / 4);
            G[i][j] =
                (G_LINEAR_2D[xi][xj] / h[0])*(M_LINEAR_2D[yi][yj] * h[1])*(M_LINEAR_2D[zi][zj] * h[2]) +
                (M_LINEAR_2D[xi][xj] * h[0])*(G_LINEAR_2D[yi][yj] / h[1])*(M_LINEAR_2D[zi][zj] * h[2]) +
                (M_LINEAR_2D[xi][xj] * h[0])*(M_LINEAR_2D[yi][yj] * h[1])*(G_LINEAR_2D[zi][zj] / h[2]);
        }
    SlauSolving::SSCM matrix;                          // матрица
    matrix.init();
    SlauSolving::SSCMBulder mBulder;                    // сборщик
    Vector b;                                           // правая часть
    SlauSolving::SSCMPreconditioner mPreconditioner;    // предобусловливатель
    SlauSolving::SSCMSolver mSolver;                    // решатель СЛАУ
    Vector x;            // решение СЛАУ

    SolverParameters ssp;
    ssp.type = SlauMetod_LOS;
    ssp.preconditioning = SlauPreconditioning_none;//SlauPreconditioning_none;
    ssp.eps = 1.e-14;
    ssp.maxIter = 500;
    b.resize(matrixSize);
    x.resize(matrixSize);
    for (int i = 0; i < matrixSize; i++)
    {
        b[i] = 0;
        x[i] = 0;
    }
    mBulder.setMatrixSize(matrixSize);
    mBulder.start();
    // занесение элементов локальных матриц в глобальную
    int gcount = 0;
    int count;
    for (int k = 0; k < g->N[2]; k++)
        for (int i = 0; i < g->N[1]; i++)
            for (int j = 0; j < g->N[0]; j++)
            {					// выбран конечный элемент (k,i,j)
                int gind[8];	// глобальные номера локальных вершин конечного элемента (k,i,j)
                int v0 = j + i * (g->N[0] + 1) + k * (g->N[1] + 1) * (g->N[0] + 1);	// глобальный номер вершины 0 - он же идентификатор (key) конечного элемента
                POINT3 p0;
                p0[0] = g->p1[0] + j * l[0] / g->N[0];
                p0[1] = g->p1[1] + i * l[1] / g->N[1];
                p0[2] = g->p1[2] + k * l[2] / g->N[2];	// прямоугольный конечный элемент задается одной вершиной с локальным номером 0
                for (int li = 0; li < 8; li++)
                {
                    int xi = li % 2;
                    int yi = (li / 2) % 2;
                    int zi = (li / 4);
                    gind[li] = v0 + xi + yi * (g->N[0] + 1) + zi * (g->N[1] + 1) * (g->N[0] + 1);	// глобальный номер локальной вершины li
                }
                // подсчет суммы базисных функций от точек на этом конечном элементе
                // и сборка локальной матрицы конечного элемента (k,i,j)
                for (int li = 0; li < 8; li++)
                    for (int lj = 0; lj < 8; lj++)	// li,lj - локальные номера базисных функций (координаты в G)
                    {
                        // локальные номера вершин li, lj соответствуют глобальным номерам вершин gind[li], gind[lj]
                        double sum = 0;
                        count = gcount;
                        for (;;)
                        {
                            if (m3[count].fei > v0 || count >= size) break;
                            sum +=	Fem::rectangleScalarLinear3D(p0, m3[count].p, li, h[0], h[1], h[2])*
                                    Fem::rectangleScalarLinear3D(p0, m3[count].p, lj, h[0], h[1], h[2]);
                            count++;
                        }
                        if(gind[li] >= gind[lj])
                            mBulder.addElement_not_null(sum + alpha*G[li][lj], gind[li], gind[lj]);
                    }
                // сборка правой части
                for (int li = 0; li < 8; li++)	// li - локальные номера базисных функций
                {
                    // локальные номера вершин li соответствуют глобальным номерам вершин gind[li]
                    double sum = 0;
                    count = gcount;
                    for (;;)
                    {
                        if (m3[count].fei > v0 || count >= size) break;
                        sum += Fem::rectangleScalarLinear3D(p0, m3[count].p, li, h[0], h[1], h[2]) * m3[count].value;
                        count++;
                    }
                    b[gind[li]] += sum;
                }
                gcount = count;	// переходим к точкам следующего элемента
            }
    mBulder.complete(matrix);
    double residual;
    double relativeResidual;
    int iterations;
    double time;
    double timePreconditioner;
    mPreconditioner.init(matrix, ssp.preconditioning, timePreconditioner);
    mSolver.init(matrix.getMatrixSize());
    mSolver.solve(matrix, b,
                  mPreconditioner, x, ssp,
                  x, residual, relativeResidual, iterations, time);
    for (int i = 0; i < matrixSize; i++)
        res[i] = x[i];
    matrix.release();
    mBulder.release();
    b.clear();
    mPreconditioner.release();
    mSolver.release();
    x.clear();
}
double Interpolator::fun_RectangleRegular3D_Lagrange1(const POINT3 &p)const
{
    int k, i, j;
    int xi, yi, zi;
    int gind;
    double E = 0;
    Grid::__GridRectangleRegular3D *g = (Grid::__GridRectangleRegular3D *)grid;
    k = (int)((p[2] - g->p1[2]) / l[2] * g->N[2]);
    i = (int)((p[1] - g->p1[1]) / l[1] * g->N[1]);
    j = (int)((p[0] - g->p1[0]) / l[0] * g->N[0]);
    int v0 = j + i * (g->N[0] + 1) + k * (g->N[1] + 1) * (g->N[0] + 1);	// глобальный номер вершины 0 - он же идентификатор (key) конечного элемента
    POINT3 p0;
    p0[0] = g->p1[0] + j * l[0] / g->N[0];
    p0[1] = g->p1[1] + i * l[1] / g->N[1];
    p0[2] = g->p1[2] + k * l[2] / g->N[2];	// прямоугольный конечный элемент задается одной вершиной с локальным номером 0
    for (int li = 0; li < 8; li++)
    {
        xi = li % 2;
        yi = (li / 2) % 2;
        zi = (li / 4);
        gind = v0 + xi + yi * (g->N[0] + 1) + zi * (g->N[1] + 1) * (g->N[0] + 1);	// глобальный номер локальной вершины li
        E += res[gind] * Fem::rectangleScalarLinear3D(p0, p, li, h[0], h[1], h[2]);
    }
    return E;

}
void Interpolator::addPoint_RectangleRegular3D_Hermite3(const POINT3 &p, const double value, const int fe_index)
{
    POINT3_VALUE *m3 = (POINT3_VALUE *)m;
    Grid::__GridRectangleRegular3D *g = (Grid::__GridRectangleRegular3D *)grid;
    if (fe_index != UNKNOWN_FE_INDEX)
    {
        m3[size] = { p, value, fe_index };
        size++;
    }
    else
    {
        int ni = (int)((p[1] - g->p1[1]) / l[1] * g->N[1]);
        int nj = (int)((p[0] - g->p1[0]) / l[0] * g->N[0]);
        int k = nj * 4 + ni * 4 * (g->N[0] + 1);
        m3[size].p = p;
        m3[size].value = value;
        m3[size].fei = k;
        size++;
    }
}
// одномерная скалярная эрмитова базисная функция 3-го порядка в точке a
double hermite1D(const double a0, const double a1, const double a, const int i, const int dif)
{
    double h = (a1 - a0);
    double ksi = (a - a0) / h;
    if (dif == 0)
    {
        if (i == 0) return 1 - 3 * ksi*ksi + 2 * ksi*ksi*ksi;
        if (i == 1) return (ksi - 2 * ksi*ksi + ksi*ksi*ksi)*h;	// чтобы производная в левом узле = 1 (функции 1 и 3)
        if (i == 2) return 3 * ksi*ksi - 2 * ksi*ksi*ksi;
        if (i == 3) return (-ksi*ksi + ksi*ksi*ksi)*h;
    } /*else
    if (dif == 1)	// первая производная по a
    {
        double dksi = 1. / h;
        if (i == 0) return -6*ksi*dksi + 6 * ksi*ksi*dksi;
        if (i == 1) return (dksi - 4 * ksi*dksi + 3*ksi*ksi*dksi)*h;
        if (i == 2) return 6 * ksi*dksi - 6 * ksi*ksi*dksi;
        if (i == 3) return (-2*ksi*dksi + 3*ksi*ksi*dksi)*h;
    }*/
    return 0;
}
// двумерная скалярная эрмитова базисная функция 3-го порядка в точке p (i = 0..15)
double hermite2D(const POINT3 &p0, const double hx, const double hy, const POINT3 &p, const int i, const DIF_STATE3 &dif)
{
    int k1 = XI(i);
    int k2 = YI(i);
    return	hermite1D(p0[0], p0[0] + hx, p[0], k1, dif[0])*
            hermite1D(p0[1], p0[1] + hy, p[1], k2, dif[1]);
}
void Interpolator::make_ip_RectangleRegular3D_Hermite3()
{
    Grid::__GridRectangleRegular3D *g = (Grid::__GridRectangleRegular3D *)grid;
    POINT3_VALUE *m3 = (POINT3_VALUE *)m;
    // сортировка точек по параллелепипедам в которые они попали
    qsort(m, size, sizeof(POINT3_VALUE), POINT3_VALUE_cmp);
    double G1x[4][4];
    double M1x[4][4];
    double G1y[4][4];
    double M1y[4][4];
    double G[16][16];		// локальная матрица жесткости
    SlauSolving::SSCM matrix;                          // матрица
    matrix.init();
    SlauSolving::SSCMBulder mBulder;                    // сборщик
    Vector b;                                           // правая часть
    SlauSolving::SSCMPreconditioner mPreconditioner;    // предобусловливатель
    SlauSolving::SSCMSolver mSolver;                    // решатель СЛАУ
    Vector x;            // решение СЛАУ

    SolverParameters ssp;
    ssp.type = SolverType::Direct;
    ssp.preconditioning = Preconditioning::SlauPreconditioning_LLT;//SlauPreconditioning_none;
    ssp.eps = 1;//1.e-14;
    ssp.maxIter = 1;
    b.resize(matrixSize);
    x.resize(matrixSize);
    for (int i = 0; i < matrixSize; i++)
    {
        b[i] = 0;
        x[i] = 0;
    }
    mBulder.setMatrixSize(matrixSize);
    mBulder.start();
    // локальные матрицы для одномерных эрмитовых базисных функций 3-го порядка
     G1x[0][0] = 36;
     G1x[0][1] = 3 * h[0];
     G1x[0][2] = -36;
     G1x[0][3] = 3 * h[0];
      G1x[1][0] = 3 * h[0];
      G1x[1][1] = 4 * h[0]*h[0];
      G1x[1][2] = -3 * h[0];
      G1x[1][3] = -h[0]*h[0];
     G1x[2][0] = -36;
     G1x[2][1] = -3 * h[0];
     G1x[2][2] = 36;
     G1x[2][3] = -3 * h[0];
      G1x[3][0] = 3 * h[0];
      G1x[3][1] = -h[0]*h[0];
      G1x[3][2] = -3 * h[0];
      G1x[3][3] = 4 * h[0]*h[0];

     G1y[0][0] = 36;
     G1y[0][1] = 3 * h[1];
     G1y[0][2] = -36;
     G1y[0][3] = 3 * h[1];
      G1y[1][0] = 3 * h[1];
      G1y[1][1] = 4 * h[1]*h[1];
      G1y[1][2] = -3 * h[1];
      G1y[1][3] = -h[1]*h[1];
     G1y[2][0] = -36;
     G1y[2][1] = -3 * h[1];
     G1y[2][2] = 36;
     G1y[2][3] = -3 * h[1];
      G1y[3][0] = 3 * h[1];
      G1y[3][1] = -h[1]*h[1];
      G1y[3][2] = -3 * h[1];
      G1y[3][3] = 4 * h[1]*h[1];

    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
        G1x[i][j] *= 1.L / (30 * h[0]);
        G1y[i][j] *= 1.L / (30 * h[1]);
    }

     M1x[0][0] = 156;
     M1x[0][1] = 22 * h[0];
     M1x[0][2] = 54;
     M1x[0][3] = -13 * h[0];
      M1x[1][0] = 22 * h[0];
      M1x[1][1] = 4 * h[0]*h[0];
      M1x[1][2] = 13 * h[0];
      M1x[1][3] = -3*h[0]*h[0];
     M1x[2][0] = 54;
     M1x[2][1] = 13 * h[0];
     M1x[2][2] = 156;
     M1x[2][3] = -22 * h[0];
      M1x[3][0] = -13 * h[0];
      M1x[3][1] = -3*h[0]*h[0];
      M1x[3][2] = -22 * h[0];
      M1x[3][3] = 4 * h[0]*h[0];

     M1y[0][0] = 156;
     M1y[0][1] = 22 * h[1];
     M1y[0][2] = 54;
     M1y[0][3] = -13 * h[1];
      M1y[1][0] = 22 * h[1];
      M1y[1][1] = 4 * h[1]*h[1];
      M1y[1][2] = 13 * h[1];
      M1y[1][3] = -3*h[1]*h[1];
     M1y[2][0] = 54;
     M1y[2][1] = 13 * h[1];
     M1y[2][2] = 156;
     M1y[2][3] = -22 * h[1];
      M1y[3][0] = -13 * h[1];
      M1y[3][1] = -3*h[1]*h[1];
      M1y[3][2] = -22 * h[1];
      M1y[3][3] = 4 * h[1]*h[1];


    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
        M1x[i][j] *= h[0] / 420;
        M1y[i][j] *= h[1] / 420;
    }

    // построение локальной матрицы G для эрмитовых базисных функций 3-го порядка
    for(int i = 0; i < 16; i++)
    for(int j = 0; j < 16; j++)
    {
        int xi = XI(i);
        int yi = YI(i);
        int xj = XI(j);
        int yj = YI(j);
        G[i][j] =
            G1x[xi][xj]*M1y[yi][yj] +
            M1x[xi][xj]*G1y[yi][yj];
    }
    // занесение элементов локальных матриц в глобальную
    int gcount = 0;
    int count;
    for(int i = 0; i < g->N[1]; i++)
    for(int j = 0; j < g->N[0]; j++)
    {					// выбран конечный элемент (i,j)
        int gind[16];	// глобальные номера локальных функций конечного элемента (i,j)
        int v0 = j * 4 + i * 4 * (g->N[0] + 1);	// глобальный номер вершиной с локальным номером 0
        POINT3 p0;
        p0[0] = g->p1[0] + j * l[0] / g->N[0];
        p0[1] = g->p1[1] + i * l[0] / g->N[1];	// прямоугольный конечный элемент задается одной вершиной с локальным номером 0
        for(int li = 0; li < 16; li++)
        {
            if (li < 8)
                gind[li] = v0 + li;
            else
                gind[li] = v0 + 4  * (g->N[0] + 1) + (li - 8);
        }
        // подсчет суммы базисных функций от точек на этом конечном элементе
        // и сборка локальной матрицы конечного элемента (i,j)
        for(int li = 0; li < 16; li++)
        for(int lj = 0; lj < 16; lj++)	// li,lj - локальные номера базисных функций (координаты в G)
        {
            // локальные номера вершин li, lj соответствуют глобальным номерам вершин gind[li], gind[lj]
            double sum = 0;
            count = gcount;
            for (;;)
            {
                if (m3[count].fei > v0 || count >= size) break;
                sum +=  hermite2D(p0, h[0], h[1], m3[count].p, li, Fem::dif_NULL3) *
                        hermite2D(p0, h[0], h[1], m3[count].p, lj, Fem::dif_NULL3);
                count++;
            }
            if(gind[li] >= gind[lj])
                mBulder.addElement_not_null(sum + alpha*G[li][lj], gind[li], gind[lj]);
        }
        // сборка правой части
        for(int li = 0; li < 16; li++)	// li - локальные номера базисных функций
        {
            // локальные номера вершин li соответствуют глобальным номерам вершин gind[li]
            double sum = 0;
            count = gcount;
            for(;;)
            {
                if (m3[count].fei > v0 || count >= size) break;
                sum += hermite2D(p0, h[0], h[1], m3[count].p, li, dif_NULL3) * m3[count].value;
                count++;
            }
            b[gind[li]] += sum;
        }
        gcount = count;	// переходим к точкам следующего элемента
    }
    mBulder.complete(matrix);
    double residual;
    double relativeResidual;
    int iterations;
    double time;
    double timePreconditioner;
    mPreconditioner.init(matrix, ssp.preconditioning, timePreconditioner);
    mSolver.init(matrix.getMatrixSize());
    mSolver.solve(matrix, b,
                  mPreconditioner, x, ssp,
                  x, residual, relativeResidual, iterations, time);
    for (int i = 0; i < matrixSize; i++)
        res[i] = x[i];
    matrix.release();
    mBulder.release();
    b.clear();
    mPreconditioner.release();
    mSolver.release();
    x.clear();
}
double Interpolator::fun_RectangleRegular3D_Hermite3(const POINT3 &p)const
{
    Grid::__GridRectangleRegular3D *g = (Grid::__GridRectangleRegular3D *)grid;
    double E = 0;
    int i = (int)((p[1] - g->p1[1]) / l[1]*g->N[1]);
    int j = (int)((p[0] - g->p1[0]) / l[0]*g->N[0]);
    int v0 = j * 4 + i * 4 * (g->N[0] + 1);	// глобальный номер вершины 0
    POINT3 p0;
    p0[0] = g->p1[0] + j * l[0] / g->N[0];
    p0[1] = g->p1[1] + i * l[1] / g->N[1];		// прямоугольный конечный элемент задается одной вершиной с локальным номером 0
    for(int li = 0; li < 16; li++)
    {
        int gind;
        if (li < 8)
            gind = v0 + li;
        else
            gind = v0 + 4 * (g->N[0] + 1) + (li - 8);
        E += res[gind] * hermite2D(p0, h[0], h[1], p, li, dif_NULL3);
    }
    return E;
}

/*bool thereIsContact = (CC_FE_Rigid->size() != 0);
if(thereIsContact)
{
    // 1) запуск итераций контакта и пластичности по-очереди
    // начнём с контакта (пластичность по касательной)
    nlInf.contact.iterationMode = IterationMode::Free;
    nlInf.plastic.iterationMode = IterationMode::ReadOnly;
    // 2) запуск итераций контакта и пластичности параллельно
    //nlInf.contact.iterationMode = ContactIterationMode::Free;
    //nlInf.plastic.iterationMode = PlasticIterationMode::Free;

}
else
{
    // если контакта нет, то запускаем пластичность, переключений в этом случае не будет
    nlInf.plastic.iterationMode = IterationMode::Free;
    nlInf.contact.iterationMode = IterationMode::ReadOnly;
}*/
/*
if(nlInf.contact.iterationMode == IterationMode::Free &&
   nlInf.plastic.iterationMode == IterationMode::ReadOnly)
{
    if(nlInf.contact.accuracyAchieved)
    {
        // точность достигнута - переключаемся на пластичность
        nlInf.plastic.iterationMode = IterationMode::Free;
        nlInf.contact.iterationMode = IterationMode::ReadOnly;
    }
}
else
if(nlInf.plastic.iterationMode == IterationMode::Free &&
   nlInf.contact.iterationMode == IterationMode::ReadOnly &&
   thereIsContact)
{
    if(nlInf.plastic.accuracyAchieved)
    {
        // точность достигнута - переключаемся на контакт
        nlInf.contact.iterationMode = IterationMode::Free;
        nlInf.plastic.iterationMode = IterationMode::ReadOnly;
    }
}*/
// вычисление перемещений
/*
for (size_t vertexInd = 0; vertexInd < grid->vertex.size(); vertexInd++)
{
    (*vertex)[vertexInd].du[0] = q[3 * vertexInd + 0];
    (*vertex)[vertexInd].du[1] = q[3 * vertexInd + 1];
    (*vertex)[vertexInd].du[2] = q[3 * vertexInd + 2];
}*/

// регулярная прямоугольная сетка 3D
struct __GridRectangleRegular3D
{
    POINT3 p1;
    POINT3 p2;
    int N[3];
    void init(const POINT3 &set_p1, const POINT3 &set_p2, const int set_N1, const int set_N2, const int set_N3);
    int find_fe_index(const POINT3 &p)const;
};
void __GridRectangleRegular3D::init(const POINT3 &set_p1, const POINT3 &set_p2, const int set_N1, const int set_N2, const int set_N3)
{
    p1 = set_p1;
    p2 = set_p2;
    N[0] = set_N1;
    N[1] = set_N2;
    N[2] = set_N3;
}
int __GridRectangleRegular3D::find_fe_index(const POINT3 &p)const
{
    double l[3];
    l[0] = p2[0] - p1[0];
    l[1] = p2[1] - p1[1];
    l[2] = p2[2] - p1[2];
    int nk = (int)((p[2] - p1[2]) / l[2] * N[2]);
    int ni = (int)((p[1] - p1[1]) / l[1] * N[1]);
    int nj = (int)((p[0] - p1[0]) / l[0] * N[0]);
    return nj + ni * (N[0] + 1) + nk * (N[1] + 1) * (N[0] + 1);
}


// параметры двумерной регулярной сетки с тиражированием в 3D
/*struct Regular2DParameters
{
    // параметры для 2D сетки
    // координата 0 - x, координата 1 - y
    // vertexNumber[i] - количество вершин по координате i
    // vertex - массив вершин, размер=vertexNumber[0]*vertexNumber[1]
    // K[i] - массив количеств разбиений по координате i, размер=vertexNumber[i]-1
    // bc1SourceIndex[i] - массив индексов ресурсов 1-х краевых на отрезках, перпендикулярных координате i, размер=vertexNumber[i]*(vertexNumber[1-i]-1)
    // bc2SourceIndex[i] - массив индексов ресурсов 2-х краевых -//-
    int vertexNumber[2];
    std::vector<POINT2> vertex;
    std::vector<int> K[2];
    std::vector<int> bc1SourceIndex[2];
    std::vector<int> bc2SourceIndex[2];
    // параметры для тиражирования
    int vertexNumberZ;     // количество вершин по z
    double z1, z2;  // первая координата по z

    // поиск индекса вершины с координатами i, j
    //POINT2 getVertex(const int i, const int j)const
    //{
    //    return vertex[i*vertexNumber[0] + j];
    //}

};*/

/*
double d = 0.1;//+1.0;
int STEPS_PUSH = 4;  // количество шагов вдавливания и вынимания
int CORR = 0;        // количество шагов коррекции после каждого шага (1 или 0)


InterpolantSurface_base::update(time, grid, set_param);
int N = ceil(time/(1 + CORR) - 0.001);
VECTOR3 unitV = {0, -1, 0};
int N_tuda = STEPS_PUSH;//(1 + CORR)*STEPS_PUSH;
if(N <= N_tuda)
{
    C = C0 + unitV*d*N/N_tuda;
}
else
{
    if(N >= N_tuda*2)
        C = C0 + unitV*d*(2 - (double)N/N_tuda) - unitV*1;
    else
        C = C0 + unitV*d*(2 - (double)N/N_tuda);
}
*/
/*if(time <= STEPS_PUSH)
{
    C = C0 + (time)*V;
}
else
{
    if(time >= 99.9999)
        C = C0 + (100 - time - 50)*V;
    else
        C = C0 + (100 - time)*V;
}*/


// нумерация контактных вершин
for(size_t csIndex = 0; csIndex < contactFESurfaceIndex.size(); csIndex++)
{

    size_t FESurfaceInd = contactFESurfaceIndex[csIndex]; // индекс КЭ поверхности, которая будет учавствовать в контакте
    std::vector<Grid::FEFace> &face = FESurface[FESurfaceInd].face;    // набор граней поверхности
    for(size_t faceInd = 0; faceInd < face.size(); faceInd++)
    {
        // контактные узлы
        // только контактные узлы
        /*int vi_4[4];   // индексы вершин, принадлежащие данной грани
        int surfaceVertexesNumber = fe[face[faceInd].feIndex]->getFaceVertexIndexes(face[faceInd].faceIndex, vi_4);
        for(int i = 0; i < surfaceVertexesNumber; i++)
        {
            //vi_4[i] - индекс вершины грани
            if(node[vi_4[i]].level != -2)
            {
                node[vi_4[i]].level = -2;   // удаление
                newOrder_seporator.push_back(vi_4[i]);// сбрасывание
            }
        }*/
        // контактные узлы и узлы соседних КЭ
        /*
        int vi_8[8];
        fe[face[faceInd].feIndex]->getVertexIndexes(vi_8);
        for(int i = 0; i < 8; i++)
        {
            if(node[vi_8[i]].level != -2)
            {
                node[vi_8[i]].level = -2;   // удаление
                newOrder_seporator.push_back(vi_8[i]);// сбрасывание
            }
        }*/
    }
}

// копирование строки i из матрицы (p, e) в плотный вектор s
// элементы вектора, которые отсутствуют в строке матрицы, не меняются
/*void SSCMStrCopyToVector(const size_t i, const SSCMPortrait &p, const SSCMElements &e, std::vector<double> &s)
{
    size_t k = p.ind[i]; // k - индекс элемента i-й строки
    size_t k_max = p.ind[i + 1];
    while(k < k_max) // обход элементов строки
    {
        size_t j = p.ai[k]; // j - номер столбца элемента e.a[k]
                            // M2(i, j) = e.a[k]
        s[j] = e.a[k];
        k++;
    }
}*/
/*
// копирование плотного вектора s в строку i матрицы (p, e)
// элементы вектора, которые отсутствуют в строке матрицы, не копируются
void VectorCopyToSSCMstr(const size_t i, const std::vector<double> &s, const SSCMPortrait &p, SSCMElements &e)
{
    size_t k = p.ind[i]; // k - индекс элемента i-й строки
    size_t k_max = p.ind[i + 1];
    while(k < k_max) // обход элементов строки
    {
        size_t j = p.ai[k]; // j - номер столбца элемента e.a[k]
                            // M2(i, j) = e.a[k]
        e.a[k] = s[j];
        k++;
    }
}*/

openscad
// отрисовка вершин и/или номеров вершин
//if(mode.POINT3s == 1 || mode.text == 1)
for (size_t i = 0; i < vertex.size(); i++)
{
    POINT3 xyz = vertex[i];	// обычные координаты вершины без деформаций
    // вершина
    //if(mode.POINT3s == 1)
            fprintf(f, "color([0.0,0.0,0.0,1])translate([%.16le,%.16le,%.16le])cube(size=[%.16le,%.16le,%.16le],center=true);\n",
                xyz[0], xyz[1], xyz[2],
                0.01L, 0.01L, 0.01L);
        //fprintf(f, "color([0.0,0.0,0.0,1])translate([%lf,%lf,%lf])sphere(d = 0.03);", xyz[0], xyz[1], xyz[2]);
    // номер вершины
        /*
    //if(mode.text == 1)
    if (i%2 == 0)
    {
            double size = 0.1;
        fprintf(f, "\
translate([%.16le,%.16le,%.16le])rotate([%.16le,%.16le,%.16le])linear_extrude(height=%.16le)\
text(text=\"%d\",size=%.16le,font=\"Times New Roman:style=Bold\",valign=\"center\",halign=\"center\");",
                    xyz[0] + size, xyz[1] + size, xyz[2] + 0,
                    0., 0., 0.,
                    size,
                    i,
                    size
                   );

    }*/
}
//if(mode.projection == 1)	// проекция на плоскость z=0 (конец блока)
    //fprintf(f, "}");



/*
// ребро(ключ)
typedef
struct
{
    int vIndex[2];  // индексы вершин
}EDGE_KEY;
// ребро(значение)
typedef
struct
{
    int N;      // количество узлов на ребре
    double k;   // коэффициент сгущения
}EDGE_VALUE;
// грань
typedef
struct
{
    int vIndex[4];  // индексы вершин
}FACE_VALUE;*/

/*
// параметры сетки - вдавливание жёсткого шара в упругое полупространство
struct ContactSphereParameters
{
    SurfaceType surfaceType;                // способ представления поверхности
    DecompositionType decompositionType;    // способ декомпозиции поверхности

    double d;   // перемещение цилиндра по вертикали
    double y0;  // начальное положение по y
    double R0;  // радиус(константа)

    Interpolation::Interpolant1D_base *R_fun;
    Interpolation::Interpolant1D_base *x_fun;
    Interpolation::Interpolant1D_base *y_fun;
    Interpolation::Interpolant1D_base *z_fun;

    double E;              // модуль упругости
    double Nu;             // коэффициент Пуассона

    int STEPS_PUSH;

    std::vector<POINT3> v;// вершины
    std::vector<HexagonArea> vInd;// 6-гранники, заданные 8-ю вершинами


    // заданы:
    // массив вершин
    // массив 6-гранников (vIndex[8], N[3], k[3])
    // массив рёбер (vIndex[2], N, k) + хэш таблица для поиска индекса ребра по вершинам
    // массив граней (vIndex[4]) + хэш таблица для поиска индекса грани по вершинам
    //  (в этих хэш таблицах вершины упорядочены)

    // построение хэш таблицы узлов:
    // для заданного 6-гранника обходятся координаты узлов (i1, i2, i3) и составляется ID узла:
    //  для вершин тип=вершина, индекс1=индекс вершины
    //  для рёбер тип=ребро, индекс1=индекс ребра, индекс2=номер вершины ребра
    //  для граней тип=грань, индекс1=индекс грани, индекс2,индекс3=координаты вершины грани
    //  для 6-гранников тип=6-гранник, индекс1=индекс 6-гранника, индекс2,индекс3,индекс4=координаты вершины 6-гранника
    // если такого узла в хэш таблице узлов нет, то добавляется




    //POINT3 v[16]; // вершины
    //int vInd[5][8]; // 5 шестигранных областей, которые будут разбиваться на КЭ
    // узел может быть одним из следующих типов:
    // 1) вершина 6-гранника - индекс вершины из v[16]
    // 2) внетренность ребра 6-гранника
    // 3) внетренность грани 6-гранника
    // 4) внутренность 6-гранника
    //POINT3 a[5][8]; // 5 шестигранных областей, которые будут разбиваться на КЭ
    //VECTOR3_uint N[5];
};*/


/*
// подграф
struct Subgraph
{
    size_t i1;
    size_t i2;  // node[i1]..node[i2]-1 - подграф
};

// разделение подграфа g (из node) сечением
// заполняются индексы уровней level найденных вложенных сечений, возвращается индекс уровня разделителя
// если separatorLevel == 0, то весь подграф является разделителем (т.е. сортировка не нужна)
void selectSeparator(const Subgraph &g, std::vector<NDNode> &node,
                     size_t &separatorLevel)
{

}
// вырезание разделителя и сохранение порядка его узлов в массив
// при этом разделитель уменьшается, если это возможно
// вырезанные узлы помечаются признаком level = 0
void cutSeporator(const Subgraph &g, std::vector<NDNode> &node, const size_t separatorLevel,
                  std::vector<size_t> &newOrder)
{

}

// переупорядочение узлов подграфа в чледующем порядке:
// сначала нижние уровни, потом верхние уровни, затем разделитель
// связи между узлами обновляются, чтобы оставались корректными
// полученные подграфы отправляются в конец очереди subgraphsSet
void sortSeparatedSubgraph(const Subgraph &g, std::vector<NDNode> &node,
                           std::deque<Subgraph> &subgraphsSet)
{

}

void Grid3D::nestedDissection()
{
    const size_t size = vertex.size();
    // инициализация
    std::vector<NDNode> node;       // узлы со связями
    node.resize(size);

    // обход всех КЭ и составление массивов рёберных связей для каждого узла
    for (size_t feInd = 0; feInd < fe.size(); feInd++)   // индекс конечного элемента
    {
        int vi_8[8];
        fe[feInd]->getVertexIndexes(vi_8);
        // рёбра + диагонали
        for(int i = 0; i < 16; i++)	// i - индекс ребра edgeHexogen_extended[]
        {
            size_t vertexIndex1 = vi_8[edgeHexogen_extended[i][0]];
            size_t vertexIndex2 = vi_8[edgeHexogen_extended[i][1]]; // vertexIndex1 - vertexIndex2 - индексы вершин ребра
            node[vertexIndex1].addEdge(vertexIndex2);
            node[vertexIndex2].addEdge(vertexIndex1);
        }
    }

    std::vector<size_t> newOrder; // изменения индексов узлов (изначально пуст)
    newOrder.reserve(size);
    std::deque<Subgraph> subgraphsSet;  // множество несвязных подграфов
    Subgraph g;
    g.i1 = 0;
    g.i2 = size;
    subgraphsSet.push_back(g); // изначально 1 область - весь граф

    // нахождение сечений
    for(;;)
    {
        // 1) извлечение очередного подграфа
        if(subgraphsSet.empty())
            break;
        g = subgraphsSet.front();
        subgraphsSet.pop_front();
        // 2) разделение подграфа g (из node) сечением
        size_t separatorLevel;
        selectSeparator(g, node,
                        separatorLevel);
        // 3) вырезание разделителя и сохранение новой нумерации его узлов в массив
        cutSeporator(g, node, separatorLevel,
                     newOrder);
        // 4) сортировка разделённого подграфа g, чтобы его части лежали в 2-х последовательных массивах
        if(separatorLevel != 0)
            sortSeparatedSubgraph(g, node,
                                  subgraphsSet);
    }
}*/


struct GraphNode
{
    //size_t vertexIndex;       // индекс узла в исходной нумерации
    size_t newVertexIndex;      // индекс узла после сортировки
    bool isSorted;              // = true если узлу присвоен новый индекс
    std::vector<size_t> edge;   // рёбра: массив индексов соседних узлов (в исходной нумерации)
    GraphNode()
    {
        edge.reserve(32);
        isSorted = false;
    }
    void addEdge(const size_t newEdgeVertexIndex)
    {
        for(size_t edgeInd = 0; edgeInd < edge.size(); edgeInd++)
        {
            if(edge[edgeInd] == newEdgeVertexIndex)
                return;// такое ребро уже есть
        }
        edge.push_back(newEdgeVertexIndex);
    }
    void setNewVertexIndex(const size_t push_vertexIndex, std::deque<size_t> &activeVertexIndex, size_t &currentMaxVertexIndex)
    {
        if(!isSorted)
        {
            newVertexIndex = currentMaxVertexIndex;
            currentMaxVertexIndex--;
            isSorted = true;
            activeVertexIndex.push_back(push_vertexIndex);
        }
    }
};

void Grid3D::sortVertexes(const std::vector<int> &contactFESurfaceIndex)
{
    std::vector<GraphNode> node;
    node.resize(vertex.size());
    // обход всех КЭ и составление массивов рёберных связей для каждого узла
    for (size_t feInd = 0; feInd < fe.size(); feInd++)   // индекс конечного элемента
    {
        int vi_8[8];
        fe[feInd]->getVertexIndexes(vi_8);
        // просто рёбра
        /*for(int i = 0; i < 12; i++)	// i - индекс ребра edgeHexogen[]
        {
            size_t vertexIndex1 = vi_8[edgeHexogen[i][0]];
            size_t vertexIndex2 = vi_8[edgeHexogen[i][1]]; // vertexIndex1 - vertexIndex2 - индексы вершин ребра
            node[vertexIndex1].addEdge(vertexIndex2);
            node[vertexIndex2].addEdge(vertexIndex1);
        }*/
        // рёбра + диагонали
        for(int i = 0; i < 16; i++)	// i - индекс ребра edgeHexogen_extended[]
        {
            size_t vertexIndex1 = vi_8[edgeHexogen_extended[i][0]];
            size_t vertexIndex2 = vi_8[edgeHexogen_extended[i][1]]; // vertexIndex1 - vertexIndex2 - индексы вершин ребра
            node[vertexIndex1].addEdge(vertexIndex2);
            node[vertexIndex2].addEdge(vertexIndex1);
        }
    } // feInd

    // нумерация контактных вершин самыми большими индексами (начиная с vertex.size() - 1)
    // и инициализация массива активных узлов
    size_t currentMaxVertexIndex = vertex.size() - 1;
    std::deque<size_t> activeVertexIndex;       // массив активных узлов (индексы узлов в старой нумерации)
    for(size_t csIndex = 0; csIndex < contactFESurfaceIndex.size(); csIndex++)
    {
        size_t FESurfaceInd = contactFESurfaceIndex[csIndex]; // индекс КЭ поверхности, которая будет учавствовать в контакте
        std::vector<Grid::FEFace> &face = FESurface[FESurfaceInd].face;    // набор граней поверхности
        for(size_t faceInd = 0; faceInd < face.size(); faceInd++)
        {
            // контактные узлы и узлы соседних КЭ
            /*
            int vi_8[8];
            fe[face[faceInd].feIndex]->getVertexIndexes(vi_8);
            for(int i = 0; i < 8; i++)
            {
                node[vi_8[i]].setNewVertexIndex(vi_8[i], activeVertexIndex, currentMaxVertexIndex);
            }*/
            // только сами контактные узлы
            int vi_4[4];   // индексы вершин, принадлежащие данной грани
            int surfaceVertexesNumber = fe[face[faceInd].feIndex]->getFaceVertexIndexes(face[faceInd].faceIndex, vi_4);
            for(int i = 0; i < surfaceVertexesNumber; i++)
            {
                //vi_4[i] - индекс вершины грани
                node[vi_4[i]].setNewVertexIndex(vi_4[i], activeVertexIndex, currentMaxVertexIndex);
            }
        }
    }
    // если контактных узлов нет, то инициализируем узлом с индексом vertex.size() - 1
    if(activeVertexIndex.empty())
    {
        size_t firstVertexIndex = vertex.size() - 1;
        node[firstVertexIndex].setNewVertexIndex(firstVertexIndex, activeVertexIndex, currentMaxVertexIndex);
    }

    // обход в ширину
    for(;;)
    {
        if(activeVertexIndex.empty())
            break;
        // взятие из начала очереди очередной активной вершины
        size_t activeVertexIndex_el = activeVertexIndex.front();
        activeVertexIndex.pop_front();
        // добавление в конец очереди соседних вершин
        GraphNode &nodeEl = node[activeVertexIndex_el]; // активный узел
        for(size_t edgeInd = 0; edgeInd < nodeEl.edge.size(); edgeInd++)
        {
            size_t edgeEl = nodeEl.edge[edgeInd];
            node[edgeEl].setNewVertexIndex(edgeEl, activeVertexIndex, currentMaxVertexIndex);
        }
    }
    fprintf(stderr, "currentMaxVertexIndex = %d\n", (int)currentMaxVertexIndex);

    // перенумерация индексов вершин: vertexIndex -> node[vertexIndex].newVertexIndex
    // fe
    for (size_t feInd = 0; feInd < fe.size(); feInd++)   // индекс конечного элемента
    {
        int vi_8[8];
        fe[feInd]->getVertexIndexes(vi_8);
        int new_vi_8[8];
        for(int i = 0; i < 8; i++)
        {
            new_vi_8[i] = node[vi_8[i]].newVertexIndex;
        }
        fe[feInd]->updateVertexIndexes(new_vi_8);
    } // feInd
    // bc1
    for (size_t bc1Ind = 0; bc1Ind < bc1.size(); bc1Ind++)
    {
        bc1[bc1Ind].vertexIndex = node[bc1[bc1Ind].vertexIndex].newVertexIndex;;
    }

    // перестановка самих вершин
    std::vector<POINT3> vertex_copy;
    vertex_copy.resize(vertex.size());
    for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
    {
        vertex_copy[vertexIndex] = vertex[vertexIndex];
    }
    for(size_t vertexIndex = 0; vertexIndex < vertex.size(); vertexIndex++)
    {
        //vertex[vertexIndex] = vertex_copy[node[vertexIndex].newVertexIndex];
        vertex[node[vertexIndex].newVertexIndex] = vertex_copy[vertexIndex];
    }
}


// сечение (узел дерева разбиений)
/*struct SeparationNode
{
    std::vector<size_t> vi; // индексы вершин разделителя
    SeparationNode *parent; // предок
    std::vector<SeparationNode *> child; // потомки
    SeparationNode()
    {
        parent = nullptr;
        child.clear();
    }
};*/



/*
    if(condensation_coord == -1)
    {
        // сгущения нет
        POINT3 p01;
        POINT3 p23;
        POINT3 p01_23;
        POINT3 p45;
        POINT3 p67;
        POINT3 p45_67;
        findPointOnTheLine_3d(v[0], v[1], N[0], 1, i[0], p01);
        findPointOnTheLine_3d(v[2], v[3], N[0], 1, i[0], p23);
        findPointOnTheLine_3d(v[4], v[5], N[0], 1, i[0], p45);
        findPointOnTheLine_3d(v[6], v[7], N[0], 1, i[0], p67);
        findPointOnTheLine_3d(p01, p23, N[1], 1, i[1], p01_23);
        findPointOnTheLine_3d(p45, p67, N[1], 1, i[1], p45_67);
        findPointOnTheLine_3d(p01_23, p45_67, N[2], 1, i[2], x);
    }
    if(condensation_coord == 0)
    {
        // сгущение по x
        POINT3 p02;
        POINT3 p46;
        POINT3 p02_46;
        POINT3 p13;
        POINT3 p57;
        POINT3 p13_57;
        // y
        findPointOnTheLine_3d(v[0], v[2], N[1], 1, i[1], p02);
        findPointOnTheLine_3d(v[4], v[6], N[1], 1, i[1], p46);
        findPointOnTheLine_3d(v[1], v[3], N[1], 1, i[1], p13);
        findPointOnTheLine_3d(v[5], v[7], N[1], 1, i[1], p57);
        // z
        findPointOnTheLine_3d(p02, p46, N[2], 1, i[2], p02_46);
        findPointOnTheLine_3d(p13, p57, N[2], 1, i[2], p13_57);
        // x
        findPointOnTheLine_3d(p02_46, p13_57, N[0], condensation_q, i[0], x);
    }
    if(condensation_coord == 1)
    {
        // сгущение по y
        POINT3 p04;
        POINT3 p15;
        POINT3 p04_15;
        POINT3 p26;
        POINT3 p37;
        POINT3 p26_37;
        // z
        findPointOnTheLine_3d(v[0], v[4], N[2], 1, i[2], p04);
        findPointOnTheLine_3d(v[1], v[5], N[2], 1, i[2], p15);
        findPointOnTheLine_3d(v[2], v[6], N[2], 1, i[2], p26);
        findPointOnTheLine_3d(v[3], v[7], N[2], 1, i[2], p37);
        // x
        findPointOnTheLine_3d(p04, p15, N[0], 1, i[0], p04_15);
        findPointOnTheLine_3d(p26, p37, N[0], 1, i[0], p26_37);
        // y
        findPointOnTheLine_3d(p04_15, p26_37, N[1], condensation_q, i[1], x);
    }
    if(condensation_coord == 2)
    {
        // сгущение по z
        POINT3 p01;
        POINT3 p23;
        POINT3 p01_23;
        POINT3 p45;
        POINT3 p67;
        POINT3 p45_67;
        // x
        findPointOnTheLine_3d(v[0], v[1], N[0], 1, i[0], p01);
        findPointOnTheLine_3d(v[2], v[3], N[0], 1, i[0], p23);
        findPointOnTheLine_3d(v[4], v[5], N[0], 1, i[0], p45);
        findPointOnTheLine_3d(v[6], v[7], N[0], 1, i[0], p67);
        // y
        findPointOnTheLine_3d(p01, p23, N[1], 1, i[1], p01_23);
        findPointOnTheLine_3d(p45, p67, N[1], 1, i[1], p45_67);
        // z
        findPointOnTheLine_3d(p01_23, p45_67, N[2], condensation_q, i[2], x);
    }
*/

// E += a*b
/*inline static void mul_add(const MATR3x3 &a, const MATR3x3 &b,
                           MATR3x3 &E)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
                E.m[i][j] += a.m[i][k] * b.m[k][j];
        }
    }
}*/
// E -= a_T_transpose * b
/*inline static void leftMulTranspose_sub(const MATR3x3 &a_T, const MATR3x3 &b,
                                         MATR3x3 &E)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
                E.m[i][j] -= a_T.m[k][i] * b.m[k][j];
        }
    }
}*/
// res = a * this_transpose - умножение матрицы a справа на транспонированную матрицу b
/*inline static void rightMulTranspose_add(const MATR3x3 &a, const MATR3x3 &b,
                              MATR3x3 &res)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            res.m[i][j] = 0;
            for (int k = 0; k < 3; k++)
                res.m[i][j] += a.m[i][k] * b.m[j][k];
        }
    }
}*/

/*MATR3x3 transpose()const
{
    MATR3x3 t;
    t.m[0][0] = m[0][0];
    t.m[0][1] = m[1][0];
    t.m[0][2] = m[2][0];
    t.m[1][0] = m[0][1];
    t.m[1][1] = m[1][1];
    t.m[1][2] = m[2][1];
    t.m[2][0] = m[0][2];
    t.m[2][1] = m[1][2];
    t.m[2][2] = m[2][2];
    return t;

    //return MATR3x3(
    //    m[0][0], m[1][0], m[2][0],
    //    m[0][1], m[1][1], m[2][1],
    //    m[0][2], m[1][2], m[2][2]);
}*/

/*MATR3x3(const double m00, const double m01, const double m02,
        const double m10, const double m11, const double m12,
        const double m20, const double m21, const double m22)
{
    m[0][0] = m00;
    m[0][1] = m01;
    m[0][2] = m02;
    m[1][0] = m10;
    m[1][1] = m11;
    m[1][2] = m12;
    m[2][0] = m20;
    m[2][1] = m21;
    m[2][2] = m22;
}*/
/*
inline MATR3x3 sqrt()const
{
    return *this;//#####
}*/




// подсчёт среднего расстояния до пограничных контактных узлов
int local_contactPointsCount = 0;
double local_r_summ = 0;
int local_r_count = 0;
for(int i = 0; i < 4; i++)
{
    if(vertexData[vi[i]].contact)
        local_contactPointsCount++;
}
// граница?
if(local_contactPointsCount >= 1 && local_contactPointsCount <= 3)
{
    // учитываем пересечения рёбер со сферой
    local_r_summ = 0;
    local_r_count = 0;
    int p1p2[4][2] =
    {
        {0, 1},
        {0, 2},
        {1, 3},
        {2, 3},
    };
    for(int n = 0; n < 4; n++)
    {
        int i = p1p2[n][0];
        int j = p1p2[n][1];
        if(vertexData[vi[i]].contact != vertexData[vi[j]].contact)
        {
            POINT3 p1; // p1 - узел в состоянии контакта
            POINT3 p2; // p1 - узел не в состоянии контакта
            if(vertexData[vi[i]].contact)
            {
                p1 = task.grid->vertex[vi[i]];
                p2 = task.grid->vertex[vi[j]];
            }
            else
            {
                p1 = task.grid->vertex[vi[j]];
                p2 = task.grid->vertex[vi[i]];
            }
            // поиск пересечения p1-p2 со сферой
            //VECTOR3 a = p1;
            VECTOR3 d = p1 - ((Interpolation::AnaliticalSurface_Sphere*)(*task.mechTask.rigidSurface)[0])->C;
            VECTOR3 r = (p2 - p1) / (p2 - p1).abs();
            double rd = r * d;
            double rr = r * r;
            double dd = d * d;
            double koren = sqrt(rd*rd - rr*(dd-R*R));
            double k1 = (-rd + koren) / rr;
            double k2 = (-rd - koren) / rr; // = 0
            VECTOR3 p = p1 + k1*r;  // искомая точка пересечения

            local_r_summ += solve_r(p);
            local_r_count++;

            /*
            if(fabs(k1) > fabs(k2))
            {
                // пересечение находится на луче p1-p2
                VECTOR3 p = p1 + k1*r;  // искомая точка пересечения
                local_r_summ += solve_r(p);
                local_r_count++;
            }
            else
            {
                //local_r_summ += solve_r(p1);
                //local_r_count++;
            }
            */
        }
    }
    local_r_summ /= local_r_count;
    local_r_count = 1;
    // 0. учитываем самый контактный узел пограничной грани
    /*
    double local_r_max = -1;
    for(int i = 0; i < 4; i++)
    {
        if(vertexData[vi[i]].contact)
        {
            const POINT3 p = task.grid->vertex[vi[i]];
            double r = solve_r(p);
            if(r > local_r_max)
                local_r_max = r;
        }
    }
    local_r_summ = local_r_max;
    local_r_count = 1;
    */
    // 1. учитываем средний контактный узел пограничной грани
    /*
    for(int i = 0; i < 4; i++)
    {
        if(vertexData[vi[i]].contact)
        {
            const POINT3 p = task.grid->vertex[vi[i]];
            local_r_summ += solve_r(p);
        }
    }
    local_r_summ /= local_contactPointsCount;
    local_r_count = 1;
    */
    // 2. учитываем каждую граничную вершину
    // каждая вершина учитывается не более 1 раза
    /*
    for(int i = 0; i < 4; i++)
    {
        if(vertexData[vi[i]].contact && !vertexCounted[vi[i]])
        {
            const POINT3 p = task.grid->vertex[vi[i]];
            local_r_summ += solve_r(p);
            local_r_count++;
            vertexCounted[vi[i]] = true;
        }
    }*/
    // 3. учитываем каждую граничную вершину
    // вершина, которая принадлежит нескольким пограничным граням, прибавляется несколько раз
    /*
    for(int i = 0; i < 4; i++)
    {
        if(vertexData[vi[i]].contact)
        {
            const POINT3 p = task.grid->vertex[vi[i]];
            local_r_summ += solve_r(p);
            local_r_count++;
        }
    }
    */
    r_summ += local_r_summ;
    r_count += local_r_count;
}




// Расчёт радиуса(устарело)
/*
double local_r_summ = 0;
int local_r_count = 0;
// граница?
if(local_contactPointsCount >= 1 && local_contactPointsCount <= 3)
{
    // учитываем пересечения рёбер со сферой
    local_r_summ = 0;
    local_r_count = 0;
    int p1p2[4][2] =
    {
        {0, 1},
        {0, 2},
        {1, 3},
        {2, 3},
    };
    for(int n = 0; n < 4; n++)
    {
        int i = p1p2[n][0];
        int j = p1p2[n][1];
        if(vertexData[vi[i]].contact != vertexData[vi[j]].contact)
        {
            POINT3 p1; // p1 - узел в состоянии контакта
            POINT3 p2; // p1 - узел не в состоянии контакта
            if(vertexData[vi[i]].contact)
            {
                p1 = task.grid->vertex[vi[i]];
                p2 = task.grid->vertex[vi[j]];
            }
            else
            {
                p1 = task.grid->vertex[vi[j]];
                p2 = task.grid->vertex[vi[i]];
            }
            // поиск пересечения p1-p2 со сферой
            //VECTOR3 a = p1;
            VECTOR3 d = p1 - ((Interpolation::AnaliticalSurface_Sphere*)(*task.mechTask.rigidSurface)[0])->C;
            VECTOR3 r = (p2 - p1) / (p2 - p1).abs();
            double rd = r * d;
            double rr = r * r;
            double dd = d * d;
            double koren = sqrt(rd*rd - rr*(dd-R*R));
            double k1 = (-rd + koren) / rr;
            double k2 = (-rd - koren) / rr; // = 0
            VECTOR3 p = p1 + k1*r;  // искомая точка пересечения

            local_r_summ += solve_r(p);
            local_r_count++;
        }
    }
    local_r_summ /= local_r_count;
    local_r_count = 1;
    r_summ += local_r_summ;
    r_count += local_r_count;
}*/

#endif
#endif // _TRASH_CODE_H
