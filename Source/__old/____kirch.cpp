
void Solver::kirch_test(KirchTask &kt, KirschParameters &kp)
{
	FILE *fres = fopen("__res.txt", "w");
	FILE *fres_0 = fopen("__res_0.txt", "w");
	FILE *fres_45 = fopen("__res_45.txt", "w");
	FILE *fres_90 = fopen("__res_90.txt", "w");
	FILE *fres_delta_max1 = fopen("__res_delta_max1.txt", "w");
	FILE *fres_delta_max2 = fopen("__res_delta_max2.txt", "w");
	if (kt.material_index == 1)
	{
		kp.material.set_K_G(1.37*9.81*(1.e10), 0.46*9.81*(1.e10));
		kp.material.g2 = 0.18*1.e6;
	}
	if (kt.material_index == 2)
	{
		kp.material.set_K_G(1.324*9.81*(1.e10), 0.468*9.81*(1.e10));
		kp.material.g2 = 0.040*1.e6;
	}
	if (kt.material_index == 3)
	{
		kp.material.set_K_G(1.786*9.81*(1.e10), 0.853*9.81*(1.e10));
		kp.material.g2 = 0.085*1.e6;
	}
	kp.material.set_matrix_isotropic();
	

	op.init(kt.integration_mode_cube, kt.integration_foursquare, kt.integration_interval);
	{
		printf("\nq = %lf\n", kt.q);
		fprintf(fres, "\nq = %lf\n", kt.q);
		for (int N_r = kt.N_r1; N_r <= kt.N_r2; N_r += kt.N_r_inc, N_r *= kt.N_r_mul)
		{
			int N_fi = N_r * 2;
			// построение сетки для задачи Кирша
			kp.N_fi = N_fi;
			kp.N_r = N_r;
			kp.q = kt.q;
			set_kirch_grid(kp, g);
			// инициализация слау
			res.init(g.v.size, g.fe.size);
			slau.init(g.v.size * 3, g.fe.size * 24 * 24 * 2, "_slau_mode.txt");	// для каждой вершины 3 базисные функции (для u.x[0], u.x[1], u.x[2])
			
			double P = kt.P - kt.Pstart;	// общие силы за вычетом сил нулевого шага
			double Pnn = 0;
			// последовательное добавление поверхностных сил
			for (int nn = 0; nn <= kt.PN; nn++)
			{
				double max1 = 0, max2 = 0;
				double d_0, d_45, d_90;
				double Padd;
				if (nn == 0)
					Padd = kt.Pstart;
				else
					Padd = P / (double)kt.PN;
				Pnn += Padd;
				// изменение действующих сил
				Solid_BoundaryCondition2_source bc2_source0 = { Solid_BoundaryCondition2_source_scalar, { { 0, 0, 0 } } };
				Solid_BoundaryCondition2_source bc2_source1 = { Solid_BoundaryCondition2_source_scalar, { { -Padd, 0, 0 } } };
				kp.bc2_left = bc2_source1;
				kp.bc2_right = bc2_source1;
				if (kt.Pmode == 1)
				{
					kp.bc2_top = bc2_source0;
					kp.bc2_bottom = bc2_source0;
				}
				if (kt.Pmode == 2)
				{
					kp.bc2_top = bc2_source1;
					kp.bc2_bottom = bc2_source1;
				}
				set_kirch_sources(kp, m, bc1_source, bc2_source);

				printf("P = %le\n", Pnn);
				fprintf(fres, "\nP = %le\n", Pnn);
				// Итерации нелинейной задачи
				for (iter_pl = 0; iter_pl < kt.maxiter; iter_pl++)
				{
					int time1 = clock();
					printf("%3d %3d %3d genA..", nn, N_r, N_fi);
					slau.start();
					buld_A();
					printf("sbor..");
					slau.complete_A();
					buld_b();
					add_bc2();
					printf("bc1..");
					add_bc1();
					printf("N=%6d,size=%8d,", slau.N, slau.slau_count);
					printf("solv..");
					slau.solve();
					
					// проверка условия завершения итерационного процесса
					double nev = 0;
					if (iter_pl != 0)
					{
						slau.solve_nev(res.x, res.r);
						nev = sqrt(slau.scal_mul(res.r, res.r) / slau.scal_mul(slau.b, slau.b));
						if (nev < 1.e-9) iter_pl = 1000000;
					}
					for (int i = 0; i < slau.N; i++)
						res.x[i] = slau.x[i];

					// расчет результатов
					results_solve();

					// РУЧНОЙ ВЫВОД В ФАЙЛ
					// частный отладочный вывод
					max1 = 0;
					max2 = 0;
					double dfmin_0 = 100, dfmin_45 = 100, dfmin_90 = 100;
					for (int k = 0; k < g.fe.size; k++)
					{
						int vsize;
						bool granica;
						switch (g.fe.fe[k].type)
						{
						case FE_PRISM3:
							vsize = 6;
						break;
						case FE_HEXAGON:
							vsize = 8;
						break;
						}
						granica = false;
						for (int t = 0; t < vsize; t++)
							if ((g.fe.fe[k].vi[t] / 2) % (N_r + 1) == 0)
							{
								granica = true;
								break;
							}
						if (granica)
						{
							POINT3 p = VECTOR3_NULL;
							// центр конечного элемента
							for (int t = 0; t < vsize; t++)	// t - локальный номер вершины
							{
								p.x[0] += g.v.v[g.fe.fe[k].vi[t]].x[0];
								p.x[1] += g.v.v[g.fe.fe[k].vi[t]].x[1];
								p.x[2] += g.v.v[g.fe.fe[k].vi[t]].x[2];
							}
							p.x[0] /= vsize;
							p.x[1] /= vsize;
							p.x[2] /= vsize;
							if (p.x[0] > 0 && p.x[1] > 0)
							{
								// угол
								double f = atan(p.x[1] / p.x[0]);
								//printf("f(%lf,%lf)=%lf\n", p.x[0], p.x[1], f);
								VECTOR6 eps, sigma;
								for (int i = 0; i < 6; i++)
								{
									eps.x[i] = res.fe[k].eps0.x[i] + res.fe[k].eps.x[i];
									sigma.x[i] = res.fe[k].sigma0.x[i] + res.fe[k].sigma.x[i];
								}
								double O = (sigma.x[0] * SQR(sin(f)) + sigma.x[1] * SQR(cos(f)) - sigma.x[5] * sin(2 * f)) / (Pnn);
								double O_analit1;
								double tK;
								double tG;
								double tg2;
								double tL;
								double O_analit2;
								if (kt.Pmode == 1)
								{
									O_analit1 = (1. - 2.*cos(2 * f));
									tK = m.m[g.fe.fe[k].mi].K;
									tG = m.m[g.fe.fe[k].mi].G;
									tg2 = m.m[g.fe.fe[k].mi].g2;
									tL = tK*tg2 / (3 * tK + tG) / tG / tG;
									O_analit2 = 1 - 2 * cos(2 * f) + (4.388*cos(2 * f) - 3.066 - 2.107*cos(4 * f) + 0.775*cos(6 * f))*tL*Pnn*Pnn;
								}
								if (kt.Pmode == 2)
								{
									O_analit1 = 2;
									tK = m.m[g.fe.fe[k].mi].K;
									tG = m.m[g.fe.fe[k].mi].G;
									tg2 = m.m[g.fe.fe[k].mi].g2;
									tL = tK*tg2 / (3 * tK + tG) / tG / tG;
									O_analit2 = 2 * (1 - (1.500 + 10.605*tL*Pnn*Pnn)*tL*Pnn*Pnn);
								}

								if (ABS((O - O_analit1)) > max1) max1 = ABS(O - O_analit1);
								if (ABS((O - O_analit2)) > max2) max2 = ABS(O - O_analit2);
								fprintf(fres, "(%.3le: %.3le vs %.3le, %.3le), ", f, O, O_analit1, O_analit2);
								if (ABS(f - PI / 2) < dfmin_90)
								{
									d_90 = O;
									dfmin_90 = ABS(f - PI / 2);
								}
								if (ABS(f - PI / 4) < dfmin_45)
								{
									d_45 = O;
									dfmin_45 = ABS(f - PI / 4);
								}
								if (ABS(f - 0) < dfmin_0)
								{
									d_0 = O;
									dfmin_0 = ABS(f - 0);
								}
							}
						}
					}
					printf("%6d %.1le %.2le %.2le\n",
						slau.get_iter(), slau.get_nev(),
						max1, max2);
					fprintf(fres, "%4d	%4d	%6d	%.1le	%.1lf	%.5le	%.5le\n",
						N_r, N_fi,
						slau.get_iter(), slau.get_nev(),
						((double)clock() - time1) / CLOCKS_PER_SEC,
						max1, max2);
					//if (slau.iter >= 20000) break;
					printf("nev = %le\n", nev);
//break;
				}
				fprintf(fres_delta_max1, "%le	%le\n", Pnn, max1);
				fprintf(fres_delta_max2, "%le	%le\n", Pnn, max2);
				fprintf(fres_90, "%le	%le\n", Pnn, d_90);
				fprintf(fres_45, "%le	%le\n", Pnn, d_45);
				fprintf(fres_0, "%le	%le\n", Pnn, d_0);
				move_setka();
				m.release();
				bc2_source.release();
				bc1_source.release();
//break;
			}
			if (N_r == kt.N_r2)
			{
				res.save("out_u.txt", "out_eps.txt", "out_gamma.txt", "out_sigma.txt", "out_tau.txt");
				Paint p;
				p.buld_OpenSCAD_models("_OpenSCAD_modes.txt", &g, &res);
				p.paint_pictures("_Paint.txt", &g, &res);
			}
			slau.release();
			res.release();
			g.release();
		}
	}
	op.release();
	fclose(fres);
	fclose(fres_0);
	fclose(fres_45);
	fclose(fres_90);
	fclose(fres_delta_max1);
	fclose(fres_delta_max2);
}

void main()
{
	Solver s;
	KirchTask kt;
	KirschParameters kp;
	// неизменные параметры
	kp.fe_type = FE_HEXAGON;//FE_PRISM3;//FE_HEXAGON
	kp.a = 100;
	kp.b = 100;
	kp.c = 100;
	kp.r = 1;
	kp.material.c_M_sigma = tc_M_sigma;
	kp.material.c_Z = tc_Z;
	kp.material.M_sigma = tM_sigma;
	kp.material.F = { { 0, 0, 0 } };
	// изменяющиеся параметры
	kt.N_r1 = 64;
	kt.N_r2 = 64;
	kt.N_r_inc = 4;
	kt.N_r_mul = 1;
	kt.q = 1. + 1. / 8.;
	kt.Pstart = 1000;
	kt.P = 400.e5;
	kt.PN = 4;				// 0 - упругая задача с поверхностными силами kt.Pstart
	kt.maxiter = 8;
	kt.material_index = 1;		// 1,2,3
	kt.Pmode = 2;				// 1 - одноосное растяжение, 2 - двуосное растяжение
	kt.integration_mode_cube = INTEGRATION_MODE_GAUSS2;
	kt.integration_foursquare = INTEGRATION_MODE_GAUSS2;
	kt.integration_interval = INTEGRATION_MODE_GAUSS2;
	s.kirch_test(kt, kp);
	//#include "plast_material.h"
	//printf("%le %le %le\n", K, G, g2); 
	//return;
	/*
	Solid_Materials m;
	m.load("in_material.txt");
	for (double x = 1.e-10; x < 1; x *= 1.1)
	{
		//printf("sigma(%le) = %le\n", -x, m.m[0].sigma_eps(-x));
		printf("sigma_tan(%le) = %le\n", -x, m.m[0].sigma_eps_tan(-x));
	}
		return;
		*/
}
