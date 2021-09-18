#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "glut.h"
#include "glaux.h"



#include <cstdlib>
#include <cmath>
#include <limits>


struct Point
{
	double x;
	double y;
	double z;
};

Point pts[1024];

////////////////////////////////////////////////////////////////////////////////////////////////
class cubic_spline
{
private:
	// Структура, описывающая сплайн на каждом сегменте сетки
	struct spline_tuple
	{
		double a, b, c, d, x;
	};
 
	spline_tuple *splines; // Сплайн
	std::size_t n; // Количество узлов сетки
 
	void free_mem(); // Освобождение памяти
 
public:
	cubic_spline(); //конструктор
	~cubic_spline(); //деструктор
 
	// Построение сплайна
	// x - узлы сетки, должны быть упорядочены по возрастанию, кратные узлы запрещены
	// y - значения функции в узлах сетки
	// n - количество узлов сетки
	void build_spline(const Point *p, std::size_t n);
 
	// Вычисление значения интерполированной функции в произвольной точке
	double f(double x) const;
};
 
cubic_spline::cubic_spline() : splines(NULL)
{
 
}
 
cubic_spline::~cubic_spline()
{
	free_mem();
}
 
void cubic_spline::build_spline(const Point *p, std::size_t n)
{
	free_mem();
 
	this->n = n;
 
	// Инициализация массива сплайнов
	splines = new spline_tuple[n];
	for (std::size_t i = 0; i < n; ++i)
	{
		splines[i].x = p[i].x;
		splines[i].a = p[i].y;
	}
	splines[0].c = 0.;
 
	// Решение СЛАУ относительно коэффициентов сплайнов c[i] методом прогонки для трехдиагональных матриц
	// Вычисление прогоночных коэффициентов - прямой ход метода прогонки
	double *alpha = new double[n - 1];
	double *beta = new double[n - 1];
	double A, B, C, F, h_i, h_i1, z;
	alpha[0] = beta[0] = 0.;
	for (std::size_t i = 1; i < n - 1; ++i)
	{
		h_i = p[i].x - p[i - 1].x, h_i1 = p[i + 1].x - p[i].x;
		A = h_i;
		C = 2. * (h_i + h_i1);
		B = h_i1;
		F = 6. * ((p[i + 1].y - p[i].y) / h_i1 - (p[i].y - p[i - 1].y) / h_i);
		z = (A * alpha[i - 1] + C);
		alpha[i] = -B / z;
		beta[i] = (F - A * beta[i - 1]) / z;
	}
 
	splines[n - 1].c = (F - A * beta[n - 2]) / (C + A * alpha[n - 2]);
 
	// Нахождение решения - обратный ход метода прогонки
	for (std::size_t i = n - 2; i > 0; --i)
		splines[i].c = alpha[i] * splines[i + 1].c + beta[i];
 
	// Освобождение памяти, занимаемой прогоночными коэффициентами
	delete[] beta;
	delete[] alpha;
 
	// По известным коэффициентам c[i] находим значения b[i] и d[i]
	for (std::size_t i = n - 1; i > 0; --i)
	{
		double h_i = p[i].x - p[i - 1].x;
		splines[i].d = (splines[i].c - splines[i - 1].c) / h_i;
		splines[i].b = h_i * (2. * splines[i].c + splines[i - 1].c) / 6. + (p[i].y - p[i - 1].y) / h_i;
	}
}
 
double cubic_spline::f(double x) const
{
	if (!splines)
		return std::numeric_limits<double>::quiet_NaN(); // Если сплайны ещё не построены - возвращаем NaN
 
	spline_tuple *s;
	if (x <= splines[0].x) // Если x меньше точки сетки x[0] - пользуемся первым эл-тов массива
		s = splines + 1;
	else if (x >= splines[n - 1].x) // Если x больше точки сетки x[n - 1] - пользуемся последним эл-том массива
		s = splines + n - 1;
	else // Иначе x лежит между граничными точками сетки - производим бинарный поиск нужного эл-та массива
	{
		std::size_t i = 0, j = n - 1;
		while (i + 1 < j)
		{
			std::size_t k = i + (j - i) / 2;
			if (x <= splines[k].x)
				j = k;
			else
				i = k;
		}
		s = splines + j;
	}
 
	double dx = (x - s->x);
	return s->a + (s->b + (s->c / 2. + s->d * dx / 6.) * dx) * dx; // Вычисляем значение сплайна в заданной точке.
}
 
void cubic_spline::free_mem()
{
	delete[] splines;
	splines = NULL;
}
////////////////////////////////////////////////////////////////////////////////////////////////







#define SIZE	32		//размер ячейки 
int N = 0;				//текущее число точек
int Width = 512, Height = 512, Deg_flag = 1, Mouse_mov = 0, pos_x = 0, pos_y = 0;
double MyScale = SIZE;
double OX = 0, OY = 0;
double delta;

int cmp(const void *p1, const void *p2)
{
	const Point *a = (const Point *)p1, *b = (const Point *)p2;
	if(a->x < b->x)return -1;else
	if(a->x > b->x)return 1;else
	return 0;
}

void display(void)
{
	double w,h;
	char buff[16];
	w = Width/2./MyScale;
	h = Height/2./MyScale;
	delta = SIZE/MyScale;
    glClear(GL_COLOR_BUFFER_BIT);//заполняет соответствующие буферы значениями, выбранными функциями glClearColor, glClearIndex, glClearDepth, glClearStencil и glClearAccum
	glLoadIdentity();//Функция glLoadIdentity заменяет текущую матрицу на единичную
	glScaled(MyScale,MyScale,MyScale);
	glTranslated(OX,OY,0);//Функция glTranslate выполняет сдвиг текущей матрицы на вектор (x, y, z)
	glColor3f(0.0, 0.0, 0.0);
    glBegin(GL_LINES);
	for (double y = 0; y <= 1024; y+=delta)
	{
		glVertex3f(-w-OX,y-OY,0);
		glVertex3f( w-OX,y-OY,0);
		glVertex3f(-w-OX,-y-OY,0);
		glVertex3f( w-OX,-y-OY,0);
	}
	for (double x = 0; x <= 1024; x+=delta)
	{
		glVertex3f( x-OX,-h-OY,0);
		glVertex3f( x-OX,h-OY,0);
		glVertex3f(-x-OX,-h-OY,0);
		glVertex3f(-x-OX,h-OY,0);
	}
	glEnd();
	if(Deg_flag)
	for (double y = 0; y <= 1024; y+=delta)
	{
		if(y+1/MyScale-OY >-h-OY && y+1/MyScale-OY < h-OY)
		{
		memset(buff,0,5);
		sprintf(buff,"%.3lf\0",y-OY);
		glRasterPos3f(1/MyScale-OX,y+1/MyScale-OY,0);
		for (int k = 0; k < 6; k++)
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10,*(buff+k));
		if(y!=0)
		{
		memset(buff,0,5);
		sprintf(buff,"%.3lf\0",-y-OY);
		glRasterPos3f(1/MyScale-OX,-y+1/MyScale-OY,0);
		for (int k = 0; k < 6; k++)
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10,*(buff+k));
		}
		}
	}
	if(Deg_flag)
	for (double x = 0; x <= 1024; x+=delta)
	{
		if(x+1/MyScale-OX >-w-OX && x+1/MyScale-OX < w-OX)
		{
		memset(buff,0,5);
		sprintf(buff,"%.3lf\0",x-OX);
		glRasterPos3f(x+1/MyScale-OX,-8/MyScale-OY,0);
		for (int k = 0; k < 6; k++)
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10,*(buff+k));			
		if(x!=0)
		{
		memset(buff,0,5);
		sprintf(buff,"%.3lf\0",-x-OX);
		glRasterPos3f(-x+1/MyScale-OX,-8/MyScale-OY,0);
		for (int k = 0; k < 6; k++)
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10,*(buff+k));
		}
		}
	}
	if(N >= 3)
    {
		glColor3f(0.0, 0.0, 1.0);
			cubic_spline spline;	//класс для построения сплайна
			spline.build_spline(pts, N);
		glPointSize(5.0);
		glColor3f(1.0, 0.0, 1.0);
		glBegin(GL_POINTS);
        for (int i = 0; i < N; i++)
            glVertex3dv((double*)&pts[i]);
		glEnd();
		glColor3f(1.0, 0.0, 0.0);
		glBegin(GL_LINE_STRIP);
        for (int i = 0; i < N-1; i++)
		{
			for(double t = 0; t <= 1.01; t += 0.01)
			{
				Point p;
				p.x = pts[i].x + t*(pts[i+1].x - pts[i].x);
				p.y = spline.f(p.x);
				p.z = 0;
				glVertex3dv((double*)&p);
			}
		}
		glEnd();
	}
	glFinish();//Функция glFinish блокирует выполнение программы до тех пор, пока не выполнятся все вызванные до этого команды OpenGL
	//glutPostRedisplay();
}

void myReshape(GLsizei w, GLsizei h)
{
	Width=w, Height=h;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
        glOrtho(-w/2, w/2, -h/2, h/2, -5.0, 5.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void mouse(int button, int state, int x, int y)
{
	Point p;
	p.x=((double)x-Width/2)/MyScale-OX;
	p.y=((double)Height/2-y)/MyScale-OY;
	p.z=0.0;
	if (state==GLUT_UP&&button==GLUT_LEFT_BUTTON)
	{
		pts[N] = p;
		N++;
		qsort(pts, N, sizeof(Point), cmp);
		Mouse_mov = 0;
	}
	if (state==GLUT_DOWN&&button==GLUT_RIGHT_BUTTON)
		{

			Mouse_mov = 1;
			pos_x=x;
			pos_y=y;
		}
	if (state==GLUT_UP&&button==GLUT_RIGHT_BUTTON)
	Mouse_mov = 0;
	glutPostRedisplay();
}

void motion(int x, int y)
{
	if (Mouse_mov)
	{	
		OX += (x - pos_x)/MyScale;
		OY += (pos_y - y)/MyScale;
		pos_x = x;
		pos_y = y;
	}
	glutPostRedisplay();
}

void Keyboard( unsigned char key, int x, int y )
{
	switch(key)
	{
	case 27: exit(0);
	case 'q': Deg_flag = !Deg_flag;			glutPostRedisplay(); break;
	case 'z': if(MyScale<1024)MyScale *= 1.1;						glutPostRedisplay(); break;
	case 'x': if(MyScale>1./1024)MyScale /= 1.1;					glutPostRedisplay(); break;
	case ' ': N=0; OX=0; OY=0; MyScale = SIZE;
											glutPostRedisplay(); break;
	case 'w': OY+=SIZE/MyScale;				glutPostRedisplay(); break;
	case 'a': OX-=SIZE/MyScale;				glutPostRedisplay(); break;
	case 's': OY-=SIZE/MyScale;				glutPostRedisplay(); break;
	case 'd': OX+=SIZE/MyScale;				glutPostRedisplay(); break;
	}
}


int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);
	glutInitWindowSize(Width, Height);
	glutCreateWindow("rgz");
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glEnable(GL_MAP1_VERTEX_3);//Вызов команды glEnable() активизирует одномерный вычислитель для трехмерных вершин
	glShadeModel(GL_FLAT);//Эта процедура указывает, сглаживать или нет углы между смежными полигонами
	glutDisplayFunc(display);
	glutReshapeFunc(myReshape);
	glutKeyboardFunc(Keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutMainLoop();
}