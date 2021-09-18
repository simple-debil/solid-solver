#define _CRT_SECURE_NO_WARNINGS

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"

#include "paint.h"

#include "slausolving.h"
#include "integration.h"
#include "fem.h"
#include "interpolation.h"

using namespace SlauSolving;
using namespace Integration;
using namespace Fem;
using namespace Interpolation;

#define RGB32(a, r, g, b)   ((b)+((g)<<8)+((r)<<16)+((a)<<24))			//возвращает индекс цвета (a, red,green,blue), каждая компонента <=2^8
#define PAINT(x_, y_, r, g, b) {if(x_ >= 0 && x_ < bmp.x && y_ >= 0 && y_ < bmp.y) bmp.c[(y_)*bmp.x + (x_)] = RGB32(0, r, g, b);}
#define NORM(a)	((a) < 0 ? (0) : ((a) > 255 ? 255 : (a)))
#define X_BMP(x_)	((((x_) - p1.x[0])/(p2.x[0] - p1.x[0])*bmp.x))
#define Y_BMP(y_)	((((y_) - p1.x[0])/(p2.x[0] - p1.x[0])*bmp.x))
#define X_REAL(x_)	(p1.x[0]+(((double)(x_))/bmp.x)*(p2.x[0] - p1.x[0]))
#define Y_REAL(y_)	(p1.x[1]+(((double)(y_))/bmp.y)*(p2.x[1] - p1.x[1]))

#define MALENKOE_CHISLO (1.e-20L)

int rebroHexogen[12][2]=
{
    {0,1},{0,2},{0,4},
    {3,1},{3,2},{3,7},
    {5,1},{5,4},{5,7},
    {6,2},{6,4},{6,7},      // пары локальных номеров вершин 6-гренника - ребра
};

int rebroPrism3[9][2]=
{
    {0,1},{0,2},{1,2},
    {3,4},{3,5},{4,5},
    {0,3},{1,4},{2,5},      // пары локальных номеров вершин треугольной призмы - ребра
};

namespace Painting
{
// построение модели тела в формате OpenSCAD
void Paint::buld_OpenSCAD_model(const OpenSCAD_MOLDEL_MODE &mode)
{
	printf("%s: ", mode.file_name);
	FILE *f = fopen(mode.file_name, "w");
	if(mode.projection == 1)	// проекция на плоскость z=0
		fprintf(f, "projection(cut = false){");
	// отрисовка ребер
	if(mode.rebra == 1)
	{
    for(int k = 0; k < grid->fe.getCount(); k++)	// k - индекс конечного элемента
	{
        if(grid->fe[k].type == GridFEType::Hexagon)
        for(int i = 0; i < 12; i++)	// i - индекс ребра rebroHexogen[]
		{
			POINT3 p1, p2;
			VECTOR3 dp;
			double len, b, c;
            p1 = grid->v[grid->fe[k].vi[rebroHexogen[i][0]]];
            p2 = grid->v[grid->fe[k].vi[rebroHexogen[i][1]]];	// задано ребро p1-p2
			dp.x[0] = p2.x[0] - p1.x[0];
			dp.x[1] = p2.x[1] - p1.x[1];
			dp.x[2] = p2.x[2] - p1.x[2];	// вектор dp из точки p1 совпадет с ребром p1-p2
			len = sqrt(SQR(dp.x[0]) + SQR(dp.x[1]) + SQR(dp.x[2]));	// длина вектора
			b = acos(dp.x[2]/len)*180.L/PI;
			c = (ABS(dp.x[0]) < MALENKOE_CHISLO) ? (SIGN(dp.x[1])*90) : ((dp.x[0]>0) ? atan(dp.x[1]/dp.x[0])*180.L/PI : atan(dp.x[1]/dp.x[0])*180.L/PI+180); 
			fprintf(f, "color([0.0,0.0,0.0,1])translate([%.16le,%.16le,%.16le])rotate([0,%.16le,%.16le])translate([0,0,%.16le])cube(size=[%.16le,%.16le,%.16le],center=true);\n",
				p1.x[0], p1.x[1], p1.x[2],
				b, c, 
				len/2,
				mode.rebra_size, mode.rebra_size, len);
            //fprintf(f, "color([0.0,0.0,0.0,1])translate([%lf,%lf,%lf])rotate([0, %lf, %lf])cylinder(h=%lf, r=%lf);\n", p1.x[0], p1.x[1], p1.x[2], b, c, len, mode.getCount()_cylinder);
		}
        if(grid->fe[k].type == GridFEType::Prism)
        for(int i = 0; i < 9; i++)	// i - индекс ребра rebroHexogen[]
        {
            POINT3 p1, p2;
            VECTOR3 dp;
            double len, b, c;
            p1 = grid->v[grid->fe[k].vi[rebroPrism3[i][0]]];
            p2 = grid->v[grid->fe[k].vi[rebroPrism3[i][1]]];	// задано ребро p1-p2
            dp.x[0] = p2.x[0] - p1.x[0];
            dp.x[1] = p2.x[1] - p1.x[1];
            dp.x[2] = p2.x[2] - p1.x[2];	// вектор dp из точки p1 совпадет с ребром p1-p2
            len = sqrt(SQR(dp.x[0]) + SQR(dp.x[1]) + SQR(dp.x[2]));	// длина вектора
            b = acos(dp.x[2]/len)*180.L/PI;
            c = (ABS(dp.x[0]) < MALENKOE_CHISLO) ? (SIGN(dp.x[1])*90) : ((dp.x[0]>0) ? atan(dp.x[1]/dp.x[0])*180.L/PI : atan(dp.x[1]/dp.x[0])*180.L/PI+180);
            fprintf(f, "color([0.0,0.0,0.0,1])translate([%.16le,%.16le,%.16le])rotate([0,%.16le,%.16le])translate([0,0,%.16le])cube(size=[%.16le,%.16le,%.16le],center=true);\n",
                p1.x[0], p1.x[1], p1.x[2],
                b, c,
                len/2,
                mode.rebra_size, mode.rebra_size, len);
            //fprintf(f, "color([0.0,0.0,0.0,1])translate([%lf,%lf,%lf])rotate([0, %lf, %lf])cylinder(h=%lf, r=%lf);\n", p1.x[0], p1.x[1], p1.x[2], b, c, len, mode.getCount()_cylinder);
        }
    }
	}
	// отрисовка вершин и/или номеров вершин
	if(mode.POINT3s == 1 || mode.text == 1)
    for (int i = 0; i < grid->v.getCount(); i++)
	{
        POINT3 xyz = grid->v[i];	// обычные координаты вершины без деформаций
		// вершина
		if(mode.POINT3s == 1)
			fprintf(f, "color([0.0,0.0,0.0,1])translate([%.16le,%.16le,%.16le])cube(size=[%.16le,%.16le,%.16le],center=true);",
				xyz.x[0], xyz.x[1], xyz.x[2],
				mode.POINT3s_size, mode.POINT3s_size, mode.POINT3s_size);
			//fprintf(f, "color([0.0,0.0,0.0,1])translate([%lf,%lf,%lf])sphere(d = 0.03);", xyz.x[0], xyz.x[1], xyz.x[2]);
		// номер вершины
		if(mode.text == 1)
		{
			fprintf(f, "\
translate([%.16le,%.16le,%.16le])rotate([%.16le,%.16le,%.16le])linear_extrude(height=%.16le)\
text(text=\"%d\",size=%.16le,font=\"Times New Roman:style=Bold\",valign=\"center\",halign=\"center\");",
						xyz.x[0] + mode.text_dx, xyz.x[1] + mode.text_dy, xyz.x[2] + mode.text_dz,
						mode.text_a, mode.text_b, mode.text_c,
						mode.text_size_z,						
						i,
						mode.text_size_xy
					   );
		}		
	}
	// отрисовка граней
	if(mode.polygons == 1)
	{
	fprintf(f, "color([1,1,1,1])polyhedron(\nPOINT3s = [\n");
	// список точек-вершин (нумерация та же: 0,1,2..)
    for (int i = 0; i < grid->v.getCount(); i++)
	{
        POINT3 xyz = grid->v[i];	// обычные координаты вершины без деформаций
		fprintf(f,"[%.16le,%.16le,%.16le]", xyz.x[0], xyz.x[1], xyz.x[2]);
        if (i + 1 < grid->v.getCount()) // если это не последняя точка, то нужна запитая
			fprintf(f,",");
	}
	fprintf(f, "\n],\nfaces = [\n");
	// для каждого конечного элемента строятся 6*2 треугольников
    for (int k = 0; k < grid->fe.getCount(); k++)
	{
        FiniteElement el = grid->fe[k];
		fprintf(f,"[%d,%d,%d],", el.vi[0], el.vi[1], el.vi[2]);
		fprintf(f,"[%d,%d,%d],", el.vi[3], el.vi[1], el.vi[2]);

		fprintf(f,"[%d,%d,%d],", el.vi[0], el.vi[1], el.vi[4]);
		fprintf(f,"[%d,%d,%d],", el.vi[5], el.vi[1], el.vi[4]);

		fprintf(f,"[%d,%d,%d],", el.vi[1], el.vi[3], el.vi[5]);
		fprintf(f,"[%d,%d,%d],", el.vi[7], el.vi[3], el.vi[5]);

		fprintf(f,"[%d,%d,%d],", el.vi[2], el.vi[3], el.vi[6]);
		fprintf(f,"[%d,%d,%d],", el.vi[7], el.vi[3], el.vi[6]);

		fprintf(f,"[%d,%d,%d],", el.vi[0], el.vi[2], el.vi[4]);
		fprintf(f,"[%d,%d,%d],", el.vi[6], el.vi[2], el.vi[4]);
		
		fprintf(f,"[%d,%d,%d],", el.vi[4], el.vi[5], el.vi[6]);
		fprintf(f,"[%d,%d,%d]",  el.vi[7], el.vi[5], el.vi[6]);
        if (k + 1 < grid->fe.getCount()) // если это не последний конечный элемент, то нужна запитая
			fprintf(f,",");
	}
	fprintf(f, "\n]);\n");
	}
	if(mode.projection == 1)	// проекция на плоскость z=0 (конец блока)
		fprintf(f, "}");
	fclose(f);
	printf("\n");
}

// построение несколько моделей в OpenSCAD с параметрами из файла
void Paint::buld_OpenSCAD_models(const char *file_name, const Grid3D *grid_new, const Solid_Results *res_new)
{
	grid = grid_new;
	res = res_new;
	int N_models;
	OpenSCAD_MOLDEL_MODE mode;
	FILE *f = fopen(file_name, "r");
	fscanf(f, "%d", &N_models);
	for(int i = 0; i < N_models; i++)
	{
		fscanf(f, "%s", mode.file_name);
		fscanf(f, "%d", &mode.nagruzka);
		fscanf(f, "%d", &mode.projection);
		fscanf(f, "%d", &mode.rebra);
		fscanf(f, "%lf", &mode.rebra_size);
		fscanf(f, "%d", &mode.POINT3s);
		fscanf(f, "%lf", &mode.POINT3s_size);
		fscanf(f, "%d", &mode.polygons);
		fscanf(f, "%d", &mode.text);
		fscanf(f, "%lf", &mode.text_a);
		fscanf(f, "%lf", &mode.text_b);
		fscanf(f, "%lf", &mode.text_c);
		fscanf(f, "%lf", &mode.text_dx);
		fscanf(f, "%lf", &mode.text_dy);
		fscanf(f, "%lf", &mode.text_dz);
		fscanf(f, "%lf", &mode.text_size_xy);
		fscanf(f, "%lf", &mode.text_size_z);
		buld_OpenSCAD_model(mode);
	}
	fclose(f);
}

// сохраняет 32-битное изображение в формате .bmp
void Paint :: bmp_save(const char *file_name, const BMP_INF &inf)
{
	int size1 = 14, size2 = 12;
	BMP_INF1 inf1;
	BMP_INF2 inf2;
	FILE *f = fopen(file_name,"wb");
	inf1.tip = 0x4D42;
	inf1.sizef = size1 + size2 + 4 * inf.x * inf.y;
	inf1.rez1 = 0;
	inf1.rez2 = 0;
	inf1.gotobits = size1 + size2;
	inf2.lenstruct = size2;
	inf2.x = (unsigned short) inf.x;
	inf2.y = (unsigned short) inf.y;
	inf2.planes = 1;
	inf2.bits = 32;
	fseek(f, 0, SEEK_SET); fwrite(&inf1.tip, 1, 2, f);
	fseek(f, 2, SEEK_SET); fwrite(&inf1.sizef, 1, 4, f);
	fseek(f, 6, SEEK_SET); fwrite(&inf1.rez1, 1, 2, f);
	fseek(f, 8, SEEK_SET); fwrite(&inf1.rez2, 1, 2, f);
	fseek(f, 10, SEEK_SET); fwrite(&inf1.gotobits, 1, 4, f);
	fseek(f, 14, SEEK_SET); fwrite(&inf2.lenstruct, 1, 4, f);
	fseek(f, 18, SEEK_SET); fwrite(&inf2.x, 1, 2, f);
	fseek(f, 20, SEEK_SET); fwrite(&inf2.y, 1, 2, f);
	fseek(f, 22, SEEK_SET); fwrite(&inf2.planes, 1, 2, f);
	fseek(f, 24, SEEK_SET); fwrite(&inf2.bits, 1, 2, f);
	fseek(f, 26, SEEK_SET); fwrite(inf.c, 1, 4 * inf.x * inf.y, f);
	fclose(f);
};

// загружает 32-битное изображение в формате .bmp
void Paint :: bmp_load(const char *file_name, BMP_INF &inf)
{
	BMP_INF2 inf2;
	FILE *f = fopen(file_name,"rb");
	fseek(f, 18, SEEK_SET); fread(&inf2.x, 1, 2, f);
	fseek(f, 20, SEEK_SET); fread(&inf2.y, 1, 2, f);
	inf.x = inf2.x;
	inf.y = inf2.y;
	inf.c = new unsigned int[inf.x * inf.y];
	fseek(f, 26, SEEK_SET); fread(inf.c, 1, 4 * inf.x * inf.y, f);
	fclose(f);
};

// вычисление цвета по значению скаляра
void Paint :: solve_rgb(const double t, int &r, int &g, int &b)
{
	if(bmp_mode == 3)	// 3 цвета
	{
		if(t >= scalar_null)
		{
			r = NORM((int)((t - scalar_null) / (scalar_max - scalar_null) * (255)));
			b = 0;
			g = (255 - r)/2;
		}
		else
		{
			r = 0;
			b = NORM((int)((scalar_null - t) / (scalar_null - scalar_min) * (255)));
			g = (255 - b)/2;
		}
	}else
	if(bmp_mode == 2)	// 2 цвета
	{
		r = NORM((int)((t - scalar_min) / (scalar_max - scalar_min) * 255));
		g = 0;
		b = 255 - r;
	}else
	if(bmp_mode == 1)	// 1 цвет
	{
		double delta = ((t - scalar_min) / (scalar_max - scalar_min));
		r = NORM((int)(pow(delta, k_mono) * 255));
		r = 255 - r;
		g = r;
		b = r;
	}
}

// -1.234e-200
int SIZE_X = 12;
int SIZE_Y = 14;	// размер одного символа
int LEN = 11;		// длина числа в символах
int KK = 11;

// выводит символ номер s из картинки symbols в картинку out с координатами x, y
void Paint :: print(int x0, int y0, const char *s, BMP_INF out)
{
	unsigned int voidc = RGB32(0, 255, 255, 255);
	y0 -= SIZE_Y/2;
	int c, color_real;
	if(bmp_mode != 1)
	{
		color_real = RGB32(0, 0, 0, 0);
	}
	else
	{
		//color_real = ~out.c[y*out.x + x];
	}

	for(int i = 0; (c = s[i]) != 0 && i < LEN; i++)
	{
		if(c >= '0' && c <= '9')
			c -= '0';
		else
		if(c == '-')
			c = 10;
		else
		if(c == '.')
			c = 11;
		else
		if(c == 'e')
			c = 12;
		else
		if(c == '+')
			c = 13;
		else
		if(c == ' ')
			goto next;
		else
		{
			//printf("<%c>\n", c);
			return;
		}
		for(int y = y0; y < y0 + SIZE_Y; y++)
		for(int x = x0; x < x0 + SIZE_X; x++)
		{
			unsigned int color = symbols.c[(y-y0)*symbols.x + (x-x0) + c*SIZE_X];
			if(color != voidc)
			if(bmp_mode != 1)
				out.c[y*out.x + x] = color_real;
			else
			{
				int c = (out.c[y*out.x + x]%256 + 128)%256;
				//out.c[y*out.x + x] = out.c[y*out.x + x]^color_real;
				if(c >= 128)
					out.c[y*out.x + x] = RGB32(0, 255, 255, 255);// RGB32(0, c, c, c);
				else
					out.c[y*out.x + x] = RGB32(0, 0, 0, 0);
			}
		}
next:;
		x0 += SIZE_X;
	}
}

// построение скаляра на сечении z = 0 в виде картинки bmp
void Paint::paint_picture(void)
{
	printf("%s: ", out_file_name);
    const void *scalar_p;
	int j, k;
	int r, g, b;
	double lx, ly;
	lx = p2.x[0] - p1.x[0];
	ly = p2.x[1] - p1.x[1];
	// получение указателя на первый элемент скаляра
    if (scalar_id == 1) scalar_p = &res->fe[0].eps.x[0];
    if (scalar_id == 2) scalar_p = &res->fe[0].eps.x[1];
    if (scalar_id == 3) scalar_p = &res->fe[0].eps.x[2];
    if (scalar_id == 4) scalar_p = &res->fe[0].eps.x[3];
    if (scalar_id == 5) scalar_p = &res->fe[0].eps.x[4];
    if (scalar_id == 6) scalar_p = &res->fe[0].eps.x[5];
    if (scalar_id == 7) scalar_p = &res->fe[0].sigma.x[0];
    if (scalar_id == 8) scalar_p = &res->fe[0].sigma.x[1];
    if (scalar_id == 9) scalar_p = &res->fe[0].sigma.x[2];
    if (scalar_id == 10) scalar_p = &res->fe[0].sigma.x[3];
    if (scalar_id == 11) scalar_p = &res->fe[0].sigma.x[4];
    if (scalar_id == 12) scalar_p = &res->fe[0].sigma.x[5];
	GridRectangleRegular_3D gr;
	gr.init(p1, p2, Nx, Ny, Nz);
    Interpolator ip;
    ip.init(&gr, LinearInterpolation, alpha, 0, grid->fe.getCount());
	scalar_min = +1.e+200;
	scalar_max = -1.e+200;
    for (k = 0; k < grid->fe.getCount(); k++)
	{
        switch (grid->fe[k].type)
		{
        case GridFEType::Hexagon:
		{
			HEXAGON e;
			VECTOR3 v;
			for (j = 0; j < 8; j++)
                e.v[j] = grid->v[grid->fe[k].vi[j]]; // скопировали вершины 6-гранника
			v = VECTOR3_NULL;	// обнулили v
			for (j = 0; j < 8; j++)
			{
				v.x[0] += e.v[j].x[0];
				v.x[1] += e.v[j].x[1];
				v.x[2] += e.v[j].x[2];
			}
			v.x[0] /= 8;
			v.x[1] /= 8;
			v.x[2] /= 8;	// центр 6-гранника - среднее арифметическое его вершин
			double scalar = *(double*)((char*)scalar_p + k*sizeof(Solid_Result_Fe));
			// отбраковываем лишние
			if (
				v.x[0] >= p1.x[0] && v.x[0] <= p2.x[0] &&
				v.x[1] >= p1.x[1] && v.x[1] <= p2.x[1] &&
				v.x[2] >= p1.x[2] && v.x[2] <= p2.x[2])
			{
				ip.add_point(v, scalar, UNKNOWN_FE_INDEX);
				if (scalar < scalar_min) scalar_min = scalar;
				if (scalar > scalar_max) scalar_max = scalar;
			}
		}
		break;
        case GridFEType::Prism:
		{
			VECTOR3 v;
			double x1, x2, x3, y1, y2, y3, z1, z2;
            x1 = grid->v[grid->fe[k].vi[0]].x[0];
            x2 = grid->v[grid->fe[k].vi[1]].x[0];
            x3 = grid->v[grid->fe[k].vi[2]].x[0];
            y1 = grid->v[grid->fe[k].vi[0]].x[1];
            y2 = grid->v[grid->fe[k].vi[1]].x[1];
            y3 = grid->v[grid->fe[k].vi[2]].x[1];
            z1 = grid->v[grid->fe[k].vi[0]].x[2];
            z2 = grid->v[grid->fe[k].vi[3]].x[2];

			v.x[0] = (x1 + x2 + x3) / 3;
			v.x[1] = (y1 + y2 + y3) / 3;
			v.x[2] = (z1 + z2) / 2;
			double scalar = *(double*)((char*)scalar_p + k*sizeof(Solid_Result_Fe));
			if (k == 0)
			{
				scalar_min = scalar;
				scalar_max = scalar;
			}
			// отбраковываем лишние
			if (v.x[0] >= p1.x[0] && v.x[0] <= p2.x[0] &&
				v.x[1] >= p1.x[1] && v.x[1] <= p2.x[1] &&
				v.x[2] >= p1.x[2] && v.x[2] <= p2.x[2])
			{
				ip.add_point(v, scalar, UNKNOWN_FE_INDEX);
				if (scalar < scalar_min) scalar_min = scalar;
				if (scalar > scalar_max) scalar_max = scalar;
			}
		}
		break;
		}
	}
	scalar_null = (scalar_min + scalar_max) / 2;
	printf("solv..");
    ip.makeInterpolant();
	printf("paint..");
	// инициализация картинки
	bmp.x = bmp_x;
	bmp.y = bmp_y;
	bmp.c = new unsigned int[bmp.x * bmp.y];
	// рисование картинки
	for (int y = 0; y < bmp.y; y++)
		for (int x = 0; x < bmp.x; x++)
		{
			POINT3 p;
			p.x[0] = p1.x[0] + lx*x / bmp.x;
			p.x[1] = p1.x[1] + ly*y / bmp.y;
			p.x[2] = 0;
			solve_rgb(ip.fun(p), r, g, b);
			PAINT(x, y, r, g, b);
		}
	ip.release();
	// цвет отверстия
	r = 255;
	g = 255;
	b = 255;
	for(int i = 0; i < size_circle; i++)
	{
		double h, l;
		double x1, x2;
		double y1, y2;
		double x, y;
		y1 = Y_BMP(circle[i].y - circle[i].r);
		y2 = Y_BMP(circle[i].y + circle[i].r);
		for(y = y1; y <= y2; y++)
		{
			h = ABS(circle[i].y - Y_REAL(y));			
			if(SQR(circle[i].r) - SQR(h) < 0) goto next1;
			l = sqrt(SQR(circle[i].r) - SQR(h));
			x1 = X_BMP(circle[i].x - l);
			x2 = X_BMP(circle[i].x + l);
			for(x = x1; x <= x2; x++)
			{
				PAINT((int)x, (int)y, r, g, b);
			}
			next1:;
		}
	}
	// пупырышки чтобы различать круг среди оттенков серого
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
	}
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
		solve_rgb(scalar_min, r, g, b);
		bmp2.c[y*bmp2.x + x] = RGB32(0, r, g, b);;
	}
	for(int y = bmp.y-SIZE_Y/2; y < bmp.y; y++)
	for(int x = bmp.x; x < bmp2.x; x++)
	{
		solve_rgb(scalar_max, r, g, b);
		bmp2.c[y*bmp2.x + x] = RGB32(0, r, g, b);;
	}

	for(int y = SIZE_Y/2; y < bmp.y-SIZE_Y/2; y++)
	for(int x = bmp.x; x < bmp2.x; x++)
	{
		solve_rgb(scalar_min + (scalar_max - scalar_min)*(y-SIZE_Y / 2)/(bmp.y - SIZE_Y), r, g, b);
		bmp2.c[y*bmp2.x + x] = RGB32(0, r, g, b);;
	}
	double scal = scalar_min;
	double dscal = (scalar_max - scalar_min)/KK;
	for(int i = 0; i <= KK; i++)
	{
		char str[100];
		str[0] = ' ';
		if(scal >=0 )
			sprintf(str+1, "%.3le", scal);
		else
			sprintf(str, "%.3le", scal);
		print(bmp.x, SIZE_Y/2 + i*(bmp.y-SIZE_Y)/KK, str, bmp2);
		scal += dscal;
	}
	//print(0, 100, "-1.234e980", bmp2);
	// сохранение изображение
	bmp_save(out_file_name, bmp2);
	delete[] bmp2.c;
	delete[] bmp.c;
	printf("\n");
}

// рисование картинок
void Paint::paint_pictures(const char *file_name, const Grid3D *grid_new, const Solid_Results *res_new)
{
	grid = grid_new;
	res = res_new;
	bmp_load("_symbols.bmp", symbols);
	int N;
	FILE *f = fopen(file_name, "r");
	fscanf(f, "%lf%lf%lf", &p1.x[0], &p1.x[1], &p1.x[2]);
	fscanf(f, "%lf%lf%lf", &p2.x[0], &p2.x[1], &p2.x[2]);
	fscanf(f, "%d%d%d", &Nx, &Ny, &Nz);
	fscanf(f, "%lf", &alpha);
	fscanf(f, "%d%d", &bmp_x, &bmp_y);
	fscanf(f, "%d", &bmp_mode);
	//fscanf(f, "%lf", &k_mono);
		k_mono = 1;

	fscanf(f, "%d", &size_circle);
	circle = new CIRCLE[size_circle];
	for(int i = 0; i < size_circle; i++)
	{
		fscanf(f, "%lf%lf%lf", &circle[i].x, &circle[i].y, &circle[i].r);
	}

	fscanf(f, "%d", &N);
	for(int i = 0; i < N; i++)
	{
		fscanf(f, "%s", out_file_name);
		fscanf(f, "%d", &scalar_id);
		paint_picture();
	}
	fclose(f);
	delete[] symbols.c;
	delete[] circle;
}
}   // namespace Paint
