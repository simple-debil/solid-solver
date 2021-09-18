/* --------------------------------------------------------- */
// ВИЗУАЛИЗАТОР
/* --------------------------------------------------------- */

#ifndef _PAINT_H
#define _PAINT_H

#include "elementary.h"
#include "grid.h"
#include "solid.h"


namespace Painting
{
using namespace Elementary;
using namespace Grid;
using namespace Solid;

// параметры визуализации в OpenSCAD
typedef struct
{
	char file_name[100];	// имя выходного файла в формате OpenSCAD
	int nagruzka;			// учет деформаций
	int projection;			// проектирование
	int rebra;				// отрисовка ребер
	 double rebra_size;			// ширина ребра
	int POINT3s;				// отрисовка вершин
	 double POINT3s_size;		// диаметр узла
	int polygons;			// отрисовка граней
	int text;				// отрисовка номеров вершин
	 double text_a;
	 double text_b;
	 double text_c;			// углы поворота текста
	 double text_dx;
	 double text_dy;
	 double text_dz;		// смещение цифр от вершин
	 double text_size_xy;	// размер цифр на плоскости
	 double text_size_z;	// размер цифр в глубину
} OpenSCAD_MOLDEL_MODE;

// рисунак 32 бит x*y
typedef struct
{
	unsigned int *c;	// буффер
	int x;				// ширина
	int y;				// высота
} BMP_INF;

// вспомогательная структура 1
typedef struct
{
	unsigned short  tip;
	unsigned int sizef;
	unsigned short rez1;
	unsigned short rez2;
	unsigned int gotobits;
} BMP_INF1;

// вспомогательная структура 2
typedef struct
{
	unsigned int lenstruct;
	unsigned short x;
	unsigned short y;
	unsigned short planes;
	unsigned short bits;
} BMP_INF2;

// круг-заплатка
typedef struct
{
	double x, y;	// центр круга (z = 0)
	double r;		// радиус
} CIRCLE;

class Paint
{
public:
    const Grid3D *grid;
	const Solid_Results *res;
CIRCLE *circle;		// заплатки
int size_circle;	// размер массива circle
	char out_file_name[100];	// имя выходного файла в формате BMP
	int scalar_id;		// идентификатор скаляра, который нужно визуализировать
// RECT3 rect;	// прямоугольная область интерполянта (z=0), видимая на изображении
POINT3 p1;			// точка с наименьшими координатами
POINT3 p2;			// точка с наибольшими координатами
int Nx, Ny, Nz;		// размерность сетки
	double alpha;
	int bmp_x, bmp_y;	// требуемая размерность картинки
	int bmp_mode;		// режим отображения картинки: 1-черно-белая, 2-красный и синий, 3-красный, зеленый, синий
	double k_mono;		// коэффициент для искажения черно-белой картинки (1- без искажения)
	double scalar_min, scalar_max, scalar_null;	// min, max, среднее между min и max значения скаляра
	BMP_INF bmp;		//
	BMP_INF bmp2;		//
	BMP_INF symbols;	// изображения
	// визуализация расслабленного и нагруженного тела в OpenScad
	void buld_OpenSCAD_model(const OpenSCAD_MOLDEL_MODE &mode);
    void buld_OpenSCAD_models(const char *file_name, const Grid3D *grid_new, const Solid_Results *res_new);
	// визуализация скаляра на сечении Z = 0
	void bmp_save(const char *file_name, const BMP_INF &inf);
	void bmp_load(const char *file_name, BMP_INF &inf);
	void solve_rgb(const double t, int &r, int &g, int &b);
	void paint_picture(void);
    void paint_pictures(const char *file_name, const Grid3D *grid_new, const Solid_Results *res_new);
	void print(int x0, int y0, const char *s, BMP_INF out);
};
}   // namespace Painting
#endif  // _PAINT_H
