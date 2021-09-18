#include <QPainter>

#include "elementary.h"
#include "interpolation.h"
#include "visualization2d.h"


using namespace Elementary;
using namespace FunParser;
using namespace Grid;
using namespace Solid;
using namespace Tests;

namespace Visualization2d
{
POINT3 pInMin = {1.e100, 1.e100, 1.e100};
POINT3 pInMax = {-1.e100, -1.e100, -1.e100};
POINT3 pOutMin;
POINT3 pOutMax;     // область сетки
POINT3 sizes1;
POINT3 sizes2;
double textSize = 16;       // ###
int x_add_left = textSize*3;    // ### место слева
int y_add_bottom = textSize*2;  // ### место снизу
int x_add_right = textSize*7;   // ### место справа
int y_add_top = textSize*1;     // ### место сверху

#define X1_TO_X2(x) (x_add_left + (x - pInMin[0]) / (pInMax[0] - pInMin[0]) * (pOutMax[0] - pOutMin[0]))
#define Y1_TO_Y2(y) (y_add_top + sizes2[1] - (y - pInMin[1]) / (pInMax[1] - pInMin[1]) * (pOutMax[1] - pOutMin[1]))

double getValue(Tests::Test_base *testBuilder, Task &task, const Solid::OutData &out, const int globalStepIndex, const size_t feInd, const size_t subHexInd, const size_t valueIndex)
{
    MechOutFePointData &mechRes = (*out.mechOut)[globalStepIndex].fe[feInd].pd[subHexInd];
    //Heat::ThermOutFePointData &thermRes = (*out.thermOut)[globalStepIndex].fe[feInd].pd[0];
    //MechMaterialSource &mechMaterial = (*task.mechTask.material)[0];
    MechMaterialSource &mechMaterial = (*(*task.mechTask.mechStep)[globalStepIndex].material)[0];
    // сетка
    if(valueIndex == 0)
    {
        return 0;
    }
    // пластичность(выход за пределы упругости)
    if(valueIndex == 1)
    {
        //double sigmaEqv = mechMaterial.sigmaEqv(mechRes.sumSigma);
        if(mechRes.isYelding == true)
            return 1;
        else
            return 0;
        /*
        if(mechMaterial.sigmaEqv(mechRes.sumSigma) >= mechMaterial.elasticSigmaLimit)
            return 1;
        else
            return 0;
        */
        /*
        if(mechMaterial.epsEqv(mechRes.sumEpsElastoplastic) >= mechMaterial.elasticEpsLimit)
            return 1;
        else
            return 0;
        */
        /*
        if(mechRes.sigmaYield >= mechMaterial.elasticSigmaLimit)
            return 1;
        else
            return 0;
        */
        /*if(mechRes.epsYield>= mechMaterial.elasticEpsLimit)
            return 1;
        else
            return 0;*/
    }
    // эквивалентное напряжение
    if(valueIndex == 2)
    {
        double value = mechRes.sumSigma.eqv_SIGMA();
        double max = testBuilder->max_value(valueIndex);
        double min = testBuilder->min_value(valueIndex);
        if(value > max)
        {
            value = max;
        }
        if(value < min)
        {
            value = min;
        }
        return value;
    }
    // среднее напряжение
    if(valueIndex == 3)
    {
        return mechRes.sumSigma.hydrostatic();
    }
    // главные напряжения
    if(valueIndex >= 4 && valueIndex <= 6)
    {
        VECTOR3 ms;
        mechRes.mainStresses(ms);
        return ms[valueIndex - 4];
    }
    // напряжения
    if(valueIndex >= 7 && valueIndex <= 12)// разница 5 (6 штук)
    {
        if(valueIndex - 7 == 0) return mechRes.sumSigma.m[0][0];
        if(valueIndex - 7 == 1) return mechRes.sumSigma.m[1][1];
        if(valueIndex - 7 == 2) return mechRes.sumSigma.m[2][2];
        if(valueIndex - 7 == 3) return mechRes.sumSigma.m[1][2];
        if(valueIndex - 7 == 4) return mechRes.sumSigma.m[2][0];
        if(valueIndex - 7 == 5) return mechRes.sumSigma.m[0][1];
    }

    // эквивалентная упугопластическая деформация
    if(valueIndex == 13)
    {
        return (mechRes.sumEpsElastic + mechRes.sumEpsPlastic).eqv_EPS();
    }
    // эквивалентная упругая деформация
    if(valueIndex == 14)
    {
        return mechRes.sumEpsElastic.eqv_EPS();
    }
    // эквивалентная пластическая деформация
    if(valueIndex == 15)
    {
        return mechRes.sumEpsPlastic.eqv_EPS();
    }
    // средняя упугопластическая деформация
    if(valueIndex == 16)
    {
        return (mechRes.sumEpsElastic + mechRes.sumEpsPlastic).hydrostatic();
    }
    // средняя упругая деформация
    if(valueIndex == 17)
    {
        return mechRes.sumEpsElastic.hydrostatic();
    }
    // средняя пластическая деформация
    if(valueIndex == 18)
    {
        return mechRes.sumEpsPlastic.hydrostatic();
    }
    // деформации
    if(valueIndex >= 19 && valueIndex <= 24)// разница 5 (6 штук)
    {
        TENSOR2 sumEps_ep = mechRes.sumEpsElastic + mechRes.sumEpsPlastic;
        if(valueIndex - 19 == 0) return sumEps_ep.m[0][0];
        if(valueIndex - 19 == 1) return sumEps_ep.m[1][1];
        if(valueIndex - 19 == 2) return sumEps_ep.m[2][2];
        if(valueIndex - 19 == 3) return sumEps_ep.m[1][2];
        if(valueIndex - 19 == 4) return sumEps_ep.m[2][0];
        if(valueIndex - 19 == 5) return sumEps_ep.m[0][1];
    }

    // температура
    if(valueIndex == 25)
    {
        //return thermRes.T;
        if(task.thermTask.enabled)
            return (*out.thermOut)[globalStepIndex].fe[feInd].pd[0].T;
        else
            return 0;
    }

    // невязка напряжений
    if(valueIndex == 26)
    {
        return abs(mechRes.residual.sigma);
    }
    // невязка деформаций
    if(valueIndex == 27)
    {
        return abs(mechRes.residual.eps);
    }
    // коэффициент Одквиста
    if(valueIndex == 28)
    {
        return mechRes.q;
    }

    return 0;
}
QColor genColor(const double value, const int colorModeIndex)
{
    if(colorModeIndex == 0)
    {
        int r, g, b;
        Interpolation::calc_rgb(2, value, 0, 0.5, 1, r, g, b);
        return QColor(r, g, b, 255);
    }
    else
    if(colorModeIndex == 1)
    {
        int c = 255 - (int)(value*255);
        return QColor(c, c, c, 255);
    }
    else
    {
        return QColor(255, 255, 255, 255);
    }
}

void drawRegionIndexation1(QPainter &p, RegionIndexation &ri, RegionDecompositionParameters &param, RegionIndexationVertexDataMap *m)
{
    if(ri.down == nullptr) // область листовая?
    {
        // область листовая
        CUBE q0;      // куб
        ri.getCube(param,
                   q0);

        // рисуем куб
        QPointF points1[4];

        points1[0].setY(Y1_TO_Y2(q0.i[1][0]));
        points1[0].setX(X1_TO_X2(q0.i[0][0]));
        points1[1].setY(Y1_TO_Y2(q0.i[1][0]));
        points1[1].setX(X1_TO_X2(q0.i[0][1]));
        points1[2].setY(Y1_TO_Y2(q0.i[1][1]));
        points1[2].setX(X1_TO_X2(q0.i[0][1]));
        points1[3].setY(Y1_TO_Y2(q0.i[1][1]));
        points1[3].setX(X1_TO_X2(q0.i[0][0]));
        if(ri.defined && ri.state.side == 1 && ri.state.onBorder == false)
        {
            p.setBrush(QBrush(QColor(0, 0, 0, 128), Qt::Dense6Pattern));
            p.setPen(QColor(255, 0, 0, 255));
            p.drawPolygon(points1, 4);
        }
        else
        {
            p.setBrush(QBrush(QColor(0, 0, 0, 0), Qt::SolidPattern));
            p.setPen(QColor(255, 0, 0, 255));
            p.drawPolygon(points1, 4);
        }
    }
    else
    {
        // область не листовая
        // рисуем подобласти
        for(int downInd = 0; downInd < 8; downInd++)
        {
            drawRegionIndexation1(p, ri.down[downInd], param, m);
        }
    }
}
void drawRegion(QPainter &p, RegionWithSpheres &r)
{
    const double z0 = -0.75;
    // в области не содержатся тела
    if(r.bodysNomber == 0)
        return;
    // количество тел != 0
    if(r.down == nullptr) // область листовая?
    {
        /*
        // область листовая
        for(int sbInd = 0; sbInd < r.bodysNomber; sbInd++)
        {
            r.sb[sbInd].s
        }
        */
    }
    else
    {
        // область не листовая
        for(int downInd = 0; downInd < 8; downInd++)
        {
            drawRegion(p, r.down[downInd]);
        }
    }
    QPointF points1[4];

    points1[0].setY(Y1_TO_Y2(r.q0.i[1][0]));
    points1[0].setX(X1_TO_X2(r.q0.i[0][0]));
    points1[1].setY(Y1_TO_Y2(r.q0.i[1][0]));
    points1[1].setX(X1_TO_X2(r.q0.i[0][1]));
    points1[2].setY(Y1_TO_Y2(r.q0.i[1][1]));
    points1[2].setX(X1_TO_X2(r.q0.i[0][1]));
    points1[3].setY(Y1_TO_Y2(r.q0.i[1][1]));
    points1[3].setX(X1_TO_X2(r.q0.i[0][0]));

    p.setPen(QColor(255, 0, 0, 255));
    p.drawPolygon(points1, 4);


    if(r.s0.O[2] - r.s0.R < z0 && r.s0.O[2] + r.s0.R > z0)
    {
        int x1 = X1_TO_X2(r.s0.O[0] - r.s0.R);
        int x2 = X1_TO_X2(r.s0.O[0] + r.s0.R);
        int y1 = Y1_TO_Y2(r.s0.O[1] + r.s0.R);
        int y2 = Y1_TO_Y2(r.s0.O[1] - r.s0.R);

        if(r.down == nullptr)
        {
            //p.setPen(QColor(255, 255, 0, 255));
            QPen pen(QColor(0, 0, 0, 255), 2);
            pen.setStyle(Qt::DotLine);
            p.setPen(pen);
            p.drawEllipse(x1, y1, x2-x1+1, y2-y1+1);
        }
        else
        {
            /*
            QPen pen(QColor(0, 255, 0, 255), 2);
            pen.setStyle(Qt::DotLine);
            p.setPen(pen);
            p.drawEllipse(x1, y1, x2-x1+1, y2-y1+1);
            */
        }
    }
}
void drawRegionIndexation(QPainter &p, Surface_base *surface)
{
    if(surface->surfaceType == SurfaceType::AnaliticalSurface_Cylinder)
    {
        Interpolation::AnaliticalSurface_Cylinder *s = (Interpolation::AnaliticalSurface_Cylinder *)surface;
        if(s->decompositionType == Grid::RegionDecompositionType::cubeVertexIndexation && s->regionIndexationInited)
        {
            Interpolation::InterpolantSurfaceDecomposition_cubeVertexIndexation *decomposition = (Interpolation::InterpolantSurfaceDecomposition_cubeVertexIndexation*) s->decomposition;
            drawRegionIndexation1(p, decomposition->ri, decomposition->param, &decomposition->map);
        }
    }
    if(surface->surfaceType == SurfaceType::InterpolantSurface_Hermite3)
    {
        Interpolation::InterpolantSurface_Hermite3 *s = (Interpolation::InterpolantSurface_Hermite3 *)surface;
        if(s->decompositionType == Grid::RegionDecompositionType::cubeVertexIndexation && s->regionIndexationInited)
        {
            Interpolation::InterpolantSurfaceDecomposition_cubeVertexIndexation *decomposition = (Interpolation::InterpolantSurfaceDecomposition_cubeVertexIndexation*) s->decomposition;
            drawRegionIndexation1(p, decomposition->ri, decomposition->param, &decomposition->map);
        }
        if(s->decompositionType == Grid::RegionDecompositionType::surfaceAsSpheres && s->regionIndexationInited)
        {
            Interpolation::InterpolantSurfaceDecomposition_surfaceAsSpheres *decomposition = (Interpolation::InterpolantSurfaceDecomposition_surfaceAsSpheres*) s->decomposition;
            drawRegion(p, decomposition->ri);
        }
    }
    if(surface->surfaceType == SurfaceType::InterpolantSurface_Lagrange3)
    {
        Interpolation::InterpolantSurface_Lagrange3 *s = (Interpolation::InterpolantSurface_Lagrange3 *)surface;
        if(s->decompositionType == Grid::RegionDecompositionType::cubeVertexIndexation && s->regionIndexationInited)
        {
            Interpolation::InterpolantSurfaceDecomposition_cubeVertexIndexation *decomposition = (Interpolation::InterpolantSurfaceDecomposition_cubeVertexIndexation*) s->decomposition;
            drawRegionIndexation1(p, decomposition->ri, decomposition->param, &decomposition->map);
        }
        if(s->decompositionType == Grid::RegionDecompositionType::surfaceAsSpheres && s->regionIndexationInited)
        {
            Interpolation::InterpolantSurfaceDecomposition_surfaceAsSpheres *decomposition = (Interpolation::InterpolantSurfaceDecomposition_surfaceAsSpheres*) s->decomposition;
            drawRegion(p, decomposition->ri);
        }
    }
}
//void drawIntersections(QPainter &p, Region &r)
void drawProjections(QPainter &painter, Task &task, int fesInd1, int fesInd2)
{
    //FILE *fff = fopen("111111", "w");
    FiniteElementSurface &fes1 = task.grid->FESurface[fesInd1]; // подвижная поверхность
    FiniteElementSurface &fes2 = task.grid->FESurface[fesInd2]; // жёсткая поверхность
    // узлы первой поверхности (подвижной)
    // проецируются на грани второй поверхности (жёсткой)
    RegionWithSpheres &r = *fes2.r0;

        // drawRegion(painter, r);

    std::vector<Grid::FEFace> &face = fes1.face;    // набор граней, из которых нам понадобится набор узлов
    std::set<int> pushedVertexIndexes;
    pushedVertexIndexes.clear();
    for(size_t faceInd = 0; faceInd < face.size(); faceInd++)  // грани подвижной поверхности
    {
        int vi[27];   // индексы вершин грани
        int surfaceVertexesNumber = (*task.grid).fe[face[faceInd].feIndex]->getFaceVertexIndexes(face[faceInd].faceIndex, vi);
        for(int i = 0; i < surfaceVertexesNumber; i++)  // локальные номера вершин некоторой грани
        {
            if(pushedVertexIndexes.count(vi[i]) == 0)   // дважды одну и ту же вершину добавлять не будем
            {
                //ContactSurfaceVertexData_FE_Rigid CVD_FE_R;
                // CVD_FE_R.init(vi[i], CCSource_FE_Rigid_ind);
                // contactVertex_FE_Rigid.push_back(CVD_FE_R);
                pushedVertexIndexes.insert(vi[i]);
                // ищем, на какие грани жёсткой поверхнеости fes2 можно спроецировать вершину с индексом v[i] (которая принадлежит некоторой грани поверхности fes1)
                POINT3 p = task.grid->vertex[vi[i]];    // вершина подвижной поверхности fes1
                Sphere s;
                s.init0(p);     // сфера 0-го радиуса, соответствующая вершине подвижной поверхности
                                // (##а должен быть отрезок v -> (v + u))
                SpheredBodyArray body;  // индексы граней жёсткой поверхности fes2, с которыми возможно пересечение подвижного узла (сферы s)
                body.clear();
                r.findIntersectingBodies(s, body);
                for(size_t bodyInd = 0; bodyInd < body.size(); bodyInd++)
                {
                    size_t fes2Ind = body[bodyInd];    // индекс грани жёсткой поверхности, с которой возможен контакт подвижного узла
                    // p - координаты подвижного узла
                    int feInd = fes2.face[fes2Ind].feIndex; // индекс конечного элемента
                    int faceLocalInd = fes2.face[fes2Ind].faceIndex; // индекс грани



                    /*
                    // узел -> середина грани
                    {
                        POINT3 gran[4];
                        task.grid->fe[feInd]->getFaceVertexes(task.grid->vertex, faceLocalInd, gran);
                        POINT3 p1 = p;
                        POINT3 p2 = (gran[0] + gran[1] + gran[2] + gran[3]) / 4;
                        QLineF line(X1_TO_X2(p1[0]), Y1_TO_Y2(p1[1]),
                                    X1_TO_X2(p2[0]), Y1_TO_Y2(p2[1]));
                        QPen pen(QColor(0, 255, 0, 255), 2);
                        pen.setStyle(Qt::DotLine);
                        painter.setPen(pen);
                        painter.drawLine(line);
                    }
                    */

                    //faceLocalInd = 2;   // X = -1 ok
                    //faceLocalInd = 3;   // X = +1 ok
                    //faceLocalInd = 4;   // Y = -1
                    //faceLocalInd = 5;   // Y = +1


                    // // узел -> узел + нормаль
                    // нормали из середин граней

                    /*
                    {
                        POINT3 gran[4];
                        task.grid->fe[feInd]->getFaceVertexes(task.grid->vertex, faceLocalInd,
                                                              gran);
                        POINT3 v[27];
                        task.grid->fe[feInd]->getVertexes(task.grid->vertex, task.grid->vertexForCurvature,
                                                          v);
                        VECTOR3 normal;
                        task.grid->fe[feInd]->getFaceNormal(v, faceLocalInd, {0., 0.},
                                                            normal);

                        for(int i = 0; i < 8; i++)
                            fprintf(fff, "(%le, %le, %le) ", v[i][0], v[i][1], v[i][2]);
                        fprintf(fff, "\n");

                        //POINT3 nearestPoint = p + normal;
                        //QLineF line(X1_TO_X2(p[0]), Y1_TO_Y2(p[1]),
                        //            X1_TO_X2(nearestPoint[0]), Y1_TO_Y2(nearestPoint[1]));
                        POINT3 p1 = (gran[0] + gran[1] + gran[2] + gran[3]) / 4;
                        POINT3 p2 = p1 + normal / 10;
                        QLineF line(X1_TO_X2(p1[0]), Y1_TO_Y2(p1[1]),
                                    X1_TO_X2(p2[0]), Y1_TO_Y2(p2[1]));
                        QPen pen(QColor(255, 0, 0, 255), 2);
                        painter.setPen(pen);
                        painter.drawLine(line);
                    }
                    */




                    {
                        POINT3 nearestPoint;
                        VECTOR3 normal;
                        double h;
                        bool res = task.grid->fe[feInd]->project(task.grid->vertex, task.grid->vertexForCurvature, faceLocalInd, p,
                                                                 nearestPoint, normal, h);
                        if(res)
                        {
                            // если проекция есть на грань, то рисуем её (0 -> nearestPoint)
                            //nearestPoint = p + normal;
                            {
                                POINT3 p1 = p;
                                POINT3 p2 = nearestPoint;
                                QLineF line(X1_TO_X2(p1[0]), Y1_TO_Y2(p1[1]),
                                            X1_TO_X2(p2[0]), Y1_TO_Y2(p2[1]));
                                QPen pen(QColor(255, 0, 0, 255), 1);
                                painter.setPen(pen);
                                painter.drawLine(line);
                            }
                            {
                                POINT3 p1 = nearestPoint;
                                POINT3 p2 = nearestPoint + normal/10;
                                QLineF line(X1_TO_X2(p1[0]), Y1_TO_Y2(p1[1]),
                                            X1_TO_X2(p2[0]), Y1_TO_Y2(p2[1]));
                                QPen pen(QColor(0, 255, 0, 255), 1);
                                pen.setStyle(Qt::DotLine);
                                painter.setPen(pen);
                                painter.drawLine(line);
                            }
                            /*{
                                POINT3 gran[4];
                                task.grid->fe[feInd]->getFaceVertexes(task.grid->vertex, faceLocalInd,
                                                                      gran);
                                POINT3 p1 = (gran[0] + gran[1] + gran[2] + gran[3]) / 4;
                                POINT3 p2 = p1 + normal/10;
                                QLineF line(X1_TO_X2(p1[0]), Y1_TO_Y2(p1[1]),
                                            X1_TO_X2(p2[0]), Y1_TO_Y2(p2[1]));
                                QPen pen(QColor(255, 0, 0, 255), 1);
                                painter.setPen(pen);
                                painter.drawLine(line);
                            }*/
                        }
                    }


                }
            }
        }
    }


}

void myWidget_out::paintEvent(QPaintEvent *)
{
    QPainter p(this);
    // заливка фона
    p.setBrush(QBrush(QColor(255, 255, 255, 255), Qt::SolidPattern));
    p.setPen(QColor(0, 0, 0, 255));
    p.drawRect(QRect(0, 0, width() - 1, height() - 1));
    if(out == nullptr)
        return;
    if(!testBuilder->possibleToShow2d())
        return;
    if(!task->mechTask.enabled)
        return;

    double circleR;
    POINT2 circleCenter;
    bool thereIsContact = testBuilder->getContactSurfaceCircle(*task, vis2d_parameters->globalStepIndex, circleCenter, circleR);


    // поиск области ввода (координаты сетки)
    pInMin = {1.e100, 1.e100, 1.e100};
    pInMax = {-1.e100, -1.e100, -1.e100};
    for(size_t vertexInd = 0; vertexInd < task->grid->vertex.size(); vertexInd++)
    {
        POINT3 p0 = (*out->mechOut)[vis2d_parameters->globalStepIndex].vertex[vertexInd].p;
        // минимальные координаты узлов сетки
        if(p0[0] < pInMin[0]) pInMin[0] = p0[0];
        if(p0[1] < pInMin[1]) pInMin[1] = p0[1];
        if(p0[2] < pInMin[2]) pInMin[2] = p0[2];
        // максимальные координаты узлов сетки
        if(p0[0] > pInMax[0]) pInMax[0] = p0[0];
        if(p0[1] > pInMax[1]) pInMax[1] = p0[1];
        if(p0[2] > pInMax[2]) pInMax[2] = p0[2];
    }
    for(int i = 0; i <= 1; i++)
    {
        if(pInMin[i] > task->grid->regionIndexationParameters.q0.i[i][0])
            pInMin[i] = task->grid->regionIndexationParameters.q0.i[i][0];
        if(pInMax[i] < task->grid->regionIndexationParameters.q0.i[i][1])
            pInMax[i] = task->grid->regionIndexationParameters.q0.i[i][1];
    }

    /*
    if(pInMin[1] > task.grid->regionIndexationParameters.q0.i[1][0])//-16
        pInMin[1] = task.grid->regionIndexationParameters.q0.i[1][0];//-16;
    if(pInMax[1] < task.grid->regionIndexationParameters.q0.i[1][1])//4
        pInMax[1] = task.grid->regionIndexationParameters.q0.i[1][1];//4;
    if(pInMin[0] > task.grid->regionIndexationParameters.q0.i[0][0])
        pInMin[0] = task.grid->regionIndexationParameters.q0.i[0][0];
    if(pInMax[0] < task.grid->regionIndexationParameters.q0.i[0][1])
        pInMax[0] = task.grid->regionIndexationParameters.q0.i[0][1];
    */
    // настройка области вывода (координаты виджета)
    int w = vis2d_parameters->pictureSize;
        // выравниваем исходную область по целым (для осей)
    double OXY_start_x;
    double OXY_start_y;
    double dx = 1;          // ###
    double dy = 1;          // ###
    {
        int minInt_x = floor(pInMin[0]);
        int minInt_y = floor(pInMin[1]);
        if(fabs(minInt_x - pInMin[0]) < 0.3) minInt_x--;   // ###
        if(fabs(minInt_y - pInMin[1]) < 0.3) minInt_y--;   // ###
        OXY_start_x = minInt_x;        // ###
        OXY_start_y = minInt_y;        // ###
        pInMin[0] = OXY_start_x;
        pInMin[1] = OXY_start_y;
    }
    {
        int maxInt_x = ceil(pInMax[0]);
        int maxInt_y = ceil(pInMax[1]);
        if(fabs(maxInt_x - pInMax[0]) < 0.3) maxInt_x++;   // ###
        if(fabs(maxInt_y - pInMax[1]) < 0.3) maxInt_y++;   // ###
        pInMax[0] = maxInt_x;
        pInMax[1] = maxInt_y;
    }




    sizes1 = pInMax - pInMin;
    sizes2[0] = (double)w;
    sizes2[1] = (double)w*sizes1[1]/sizes1[0];
    pOutMin = {0, 0, 0};                             //
    pOutMax = {sizes2[0] - 1, sizes2[1] - 1, 0};     // область сетки

    resize(sizes2[0] + x_add_left + x_add_right, sizes2[1] + y_add_bottom + y_add_top);
    //p.translate(sizes2[0] * (0 - pMin1[0])/sizes1[0] + w_add, sizes2[1] * (0 - pMin1[1])/sizes1[1] + h_add);
    // поиск минимального и максимального значений отображаемого поля
    double valueMin = 1.e100;
    double valueMax = -1.e100;
    for(size_t feInd = 0; feInd < task->grid->fe.size(); feInd++)
    {
        double value = 0;
        if(task->grid->fe[feInd]->get_FEType() == Grid::FEType::LinearHexagon_XFEM)
        {
            size_t pd_size = (*out->mechOut)[vis2d_parameters->globalStepIndex].fe[feInd].pd.size();
            //Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<Grid::FE_LinearHexagon_XFEM *>(task->grid->fe[feInd]);
            //for(size_t subHexInd = 0; subHexInd < feEl_XFEM->subHexagonArea.size(); subHexInd++) // индекс 6-гранной подобласти КЭ
            for(size_t subHexInd = 0; subHexInd < pd_size; subHexInd++) // индекс 6-гранной подобласти КЭ
            {
                value = getValue(testBuilder, *task, *out, vis2d_parameters->globalStepIndex, feInd, subHexInd, vis2d_parameters->valueIndex);
                //value = 1.0*1.e4; //#####
                //value = 3.4*1.e4;
                // минимальное значение
                if(value < valueMin) valueMin = value;
                // максимальное значение
                if(value > valueMax) valueMax = value;
            }
        }
        else
        {
            value = getValue(testBuilder, *task, *out, vis2d_parameters->globalStepIndex, feInd, 0, vis2d_parameters->valueIndex);
            // минимальное значение
            if(value < valueMin) valueMin = value;
            // максимальное значение
            if(value > valueMax) valueMax = value;
        }
    }

    // подкраска области, в которую заключена сетка
    p.setBrush(QBrush(QColor(240, 240, 240, 255), Qt::SolidPattern));
    p.setPen(QColor(0, 0, 0, 0));
    p.drawRect(X1_TO_X2(pInMin[0]), Y1_TO_Y2(pInMin[1]), X1_TO_X2(pInMax[0]) - X1_TO_X2(pInMin[0]), Y1_TO_Y2(pInMax[1]) - Y1_TO_Y2(pInMin[1]));
    p.setPen(QColor(0, 0, 0, 255));

    int cn1, cn2;
    testBuilder->notFixedCoordinates(cn1, cn2);
    // вывод четырёхугольников
    //if(false)
    for(size_t feInd = 0; feInd < task->grid->fe.size(); feInd++)
    {
        std::array<bool, 6> faceState;
        if(task->grid->fe[feInd]->get_FEType() == Grid::FEType::LinearHexagon_XFEM)
        {
            Grid::FE_LinearHexagon_XFEM *feEl_XFEM = dynamic_cast<Grid::FE_LinearHexagon_XFEM *>(task->grid->fe[feInd]);
            for(size_t subHexInd = 0; subHexInd < feEl_XFEM->subHex.size(); subHexInd++) // индекс 6-гранной подобласти КЭ
            {
                // отрисовка подобласти (XFEM)
                testBuilder->needToDrawSubFe(*task, *out, feInd, subHexInd, vis2d_parameters->z_index, faceState);
                for(int faceIndex = 0; faceIndex < 6; faceIndex++)
                {
                    if(faceState[faceIndex])
                    {
                        POINT3 v[4];
                        feEl_XFEM->gn_getSubFaceVertexes(subHexInd, faceIndex, v);
                        std::swap(v[2], v[3]);
                        QPointF points[4];
                        for(int i = 0; i < 4; i++)
                        {
                            //int vertexInd = vi[i];
                            //POINT3 p0 = (*out->mechOut)[vis2d_parameters->globalStepIndex].vertex[vertexInd].p;
                            POINT3 p0 = v[i];
                            points[i].setX(X1_TO_X2(p0[cn1]));
                            points[i].setY(Y1_TO_Y2(p0[cn2]));
                        }
                        double value = getValue(testBuilder, *task, *out, vis2d_parameters->globalStepIndex, feInd, subHexInd, vis2d_parameters->valueIndex);
                        if(valueMax - valueMin == 0)
                            p.setBrush(QBrush(genColor(0, vis2d_parameters->colorModeIndex), Qt::SolidPattern));
                        else
                            p.setBrush(QBrush(genColor((value - valueMin) / (valueMax - valueMin), vis2d_parameters->colorModeIndex), Qt::SolidPattern));
                        //p.setBrush(QBrush(QColor(0, 0, 0, 0), Qt::SolidPattern));
                        p.setPen(QColor(0, 0, 0, 255));
                        p.drawPolygon(points, 4);

                    }
                }
            }
        }
        else
        {
            // отрисовка КЭ
            testBuilder->needToDrawFe(*task, *out, feInd, vis2d_parameters->z_index, faceState);
            for(int faceIndex = 0; faceIndex < 6; faceIndex++)
            {
                if(faceState[faceIndex])
                {
                    int vi[4];
                    task->grid->fe[feInd]->getFaceVertexIndexes(faceIndex, vi);
                    std::swap(vi[2], vi[3]);
                    QPointF points[4];
                    for(int i = 0; i < 4; i++)
                    {
                        int vertexInd = vi[i];
                        POINT3 p0 = (*out->mechOut)[vis2d_parameters->globalStepIndex].vertex[vertexInd].p;
                        points[i].setX(X1_TO_X2(p0[cn1]));
                        points[i].setY(Y1_TO_Y2(p0[cn2]));
                    }
                    double value = getValue(testBuilder, *task, *out, vis2d_parameters->globalStepIndex, feInd, 0, vis2d_parameters->valueIndex);
                    if(valueMax - valueMin == 0)
                        p.setBrush(QBrush(genColor(0, vis2d_parameters->colorModeIndex), Qt::SolidPattern));
                    else
                        p.setBrush(QBrush(genColor((value - valueMin) / (valueMax - valueMin), vis2d_parameters->colorModeIndex), Qt::SolidPattern));
                    //p.setBrush(QBrush(QColor(0, 0, 0, 0), Qt::SolidPattern));
                    p.setPen(QColor(0, 0, 0, 255));
                    p.drawPolygon(points, 4);

                }
            }
        }
        /*if(paintGrid)
        {
            int reindex[8] = {0, 1, 3, 2, 4,5,7,6};
            QPointF points1[4];
            QPointF points2[4];
            //FiniteElement fe_el = task.grid->fe[feInd];
            for(int i = 0; i < 4; i++)
            {
                int vertexInd = task->grid->fe[feInd]->globalVertexIndex(reindex[i]);
                //POINT3 dp = {0,0,0};
                POINT3 dp = (*out->mechOut)[vis2d_parameters->globalStepIndex].vertex[vertexInd].sum_du -
                            (*out->mechOut)[(*out->mechOut).size() - 2].vertex[vertexInd].sum_du;
                POINT3 p = task->grid->vertex[vertexInd] + dp;
                points1[i].setX(X1_TO_X2(p[0]));
                points1[i].setY(Y1_TO_Y2(p[1]));
            }
            double value = getValue(*task, *out, vis2d_parameters->globalStepIndex, feInd, vis2d_parameters->valueIndex);
            if(valueMax - valueMin == 0)
                p.setBrush(QBrush(genColor(0, vis2d_parameters->colorModeIndex), Qt::SolidPattern));
            else
                p.setBrush(QBrush(genColor((value - valueMin) / (valueMax - valueMin), vis2d_parameters->colorModeIndex), Qt::SolidPattern));
            p.drawPolygon(points1, 4);
            for(int i = 4; i < 8; i++)
            {
                int vertexInd = task->grid->fe[feInd]->globalVertexIndex(reindex[i]);
                //POINT3 dp = {0,0,0};
                POINT3 dp = (*out->mechOut)[vis2d_parameters->globalStepIndex].vertex[vertexInd].sum_du -
                            (*out->mechOut)[(*out->mechOut).size() - 2].vertex[vertexInd].sum_du;
                POINT3 p = task->grid->vertex[vertexInd] + dp;
                points2[i-4].setX(X1_TO_X2(p[0]));
                points2[i-4].setY(Y1_TO_Y2(p[1]));
            }
            p.setBrush(QBrush(QColor(0, 0, 0, 0), Qt::SolidPattern));
            p.setPen(QColor(0, 0, 0, 255));
            p.drawPolygon(points2, 4);
        }*/
    }


    // вывод стрелочек из узлов, на которые давит опора
    if(false)
    if(thereIsContact)
    for (size_t vertexInd = 0; vertexInd < task->grid->vertex.size(); vertexInd++)
    {
        bool contacted = (*out->mechOut)[vis2d_parameters->globalStepIndex].vertex[vertexInd].contact;
        if(contacted)
        {
            POINT3 p0 = (*out->mechOut)[vis2d_parameters->globalStepIndex].vertex[vertexInd].p;
            VECTOR3 F = (*out->mechOut)[vis2d_parameters->globalStepIndex].vertex[vertexInd].F_sum;
            VECTOR3 direction = F/F.abs();
            QPointF p1;
            QPointF p2;
            p1 = {X1_TO_X2(p0[0]), Y1_TO_Y2(p0[1])};
            p2 = {X1_TO_X2(p0[0] + direction[0]), Y1_TO_Y2(p0[1] + direction[1])};
            p.setBrush(QBrush(QColor(0, 0, 255, 255), Qt::SolidPattern));
            p.setPen(QColor(0, 0, 255, 255));
            //if(sumF != 0)
            //int size = textSize/4;
            //int x = p1.x() - size/2;
            //int y = p1.y() - size/2;
            //p.drawRect(x, y, size, size);
                p.drawLine(p1, p2);
        }
    }


    /*
    // первые краевые
    for (size_t bc1Index = 0; bc1Index < task.grid->bc1.size(); bc1Index++)
    {
        Grid::BoundaryCondition1 bc_el = task.grid->bc1[bc1Index];
        int vertexIndex = bc_el.vi;
        if(bc_el.si == 2)
        {
            POINT3 dp = (*out.mechOut)[globalStepIndex].vertex[vertexIndex].sum_du -
                        (*out.mechOut)[(*out.mechOut).size() - 2].vertex[vertexIndex].sum_du;
            POINT3 p0 = task.grid->vertex[vertexIndex] + dp;
            QPointF p1;
            p1 = {X1_TO_X2(p0[0]), Y1_TO_Y2(p0[1])};
            p.setBrush(QBrush(QColor(0, 0, 255, 255), Qt::SolidPattern));
            p.setPen(QColor(0, 0, 255, 255));
            //if(sumF != 0)
            int size = textSize/4;
            int x = p1.x() - size/2;
            int y = p1.y() - size/2;
            p.drawRect(x, y, size, size);
        }
    }
    */
    // отрисовка контактной КЭ-поверхности
    /*
    if(task.mechTask.contactElasticType == 1)
    for (size_t faceInd = 0; faceInd < (*task.mechTask.CFESurface)[0].face.size(); faceInd++)
    {
        break;
        int vi[4];
        POINT3 v[4];        // вершины грани
        int feInd = (*task.mechTask.CFESurface)[0].face[faceInd].feIndex;
        int faceLocalInd = (*task.mechTask.CFESurface)[0].face[faceInd].faceIndex;
        task.grid->fe[feInd]->getFaceVertexIndexes(faceLocalInd, vi);
        //task.grid->fe[feInd]->getFaceVertexes(task.grid->vertex, faceLocalInd, v);
        for(int i = 0; i < 4; i++)
        {
            int vertexInd = vi[i];
            POINT3 dp = (*out.mechOut)[globalStepIndex].vertex[vertexInd].sum_du -
                        (*out.mechOut)[(*out.mechOut).size() - 2].vertex[vertexInd].sum_du;
            v[i] = task.grid->vertex[vertexInd] + dp;
            POINT3 p0 = v[i];
            QPointF p1;
            p1 = {X1_TO_X2(p0[0]), Y1_TO_Y2(p0[1])};
            p.setBrush(QBrush(QColor(0, 0, 255, 255), Qt::SolidPattern));
            p.setPen(QColor(255, 0, 255, 255));
            int size = textSize/8;
            int x = p1.x() - size/2;
            int y = p1.y() - size/2;
            p.drawRect(x, y, size, size);
        }
    }*/

    // вывод осей и шкалы
    QPointF p1;
    QPointF p2;
    /*
    double x1 = -11;        // ###
    double y1 = -4;         // ###
    double dx = 1;          // ###
    double dy = 1;          // ###
    */
    double xMarkSize = textSize/2;
    double yMarkSize = textSize/2;
    QPointF dpTextX = {-textSize*1.5, textSize*2};
    QPointF dpTextY = {-textSize*3.5, textSize*0.4};
    // Ось ОХ
    double x = OXY_start_x;
    for(;;)
    {
        if(x >= pInMax[0])
            break;
        p1 = {X1_TO_X2(x), Y1_TO_Y2(OXY_start_y) - yMarkSize};
        p2 = {X1_TO_X2(x), Y1_TO_Y2(OXY_start_y) + yMarkSize};
        p.drawLine(p1, p2);
        char strt[1000];
        sprintf(strt, "%6.1lf", x);
        QFont font;
        font.setPixelSize(textSize);
        p.setFont(font);
        QPointF pText = p1;
        pText += dpTextX;
        p.setBrush(QBrush(QColor(0, 0, 0, 0), Qt::SolidPattern));
        p.setPen(QColor(0, 0, 0, 255));
        p.drawText(pText, strt);
        x += dx;
    }
    p1 = {X1_TO_X2(OXY_start_x), Y1_TO_Y2(OXY_start_y)};
    p2 = {X1_TO_X2(pInMax[0]), Y1_TO_Y2(OXY_start_y)};
    p.setBrush(QBrush(QColor(0, 0, 0, 0), Qt::SolidPattern));
    p.setPen(QColor(0, 0, 0, 255));
    p.drawLine(p1, p2);
    // Ось ОХ
    double y = OXY_start_y;
    for(;;)
    {
        if(y >= pInMax[1])
            break;
        p1 = {X1_TO_X2(OXY_start_x) + xMarkSize, Y1_TO_Y2(y)};
        p2 = {X1_TO_X2(OXY_start_x) - xMarkSize, Y1_TO_Y2(y)};
        p.drawLine(p1, p2);
        char strt[1000];
        sprintf(strt, "%6.1lf", y);
        QFont font;
        font.setPixelSize(textSize);
        p.setFont(font);
        QPointF pText = p1;
        pText += dpTextY;
        p.drawText(pText, strt);
        y += dy;
    }
    p1 = {X1_TO_X2(OXY_start_x), Y1_TO_Y2(OXY_start_y)};
    p2 = {X1_TO_X2(OXY_start_x), Y1_TO_Y2(pInMax[1])};
    p.drawLine(p1, p2);
    // шкала
    int N = 6;
    y = OXY_start_y;
    for(int i = 0; i <= N; i++)
    {
        double volume = (double)i / N;
        double value = valueMin + (valueMax - valueMin) * volume;
        double y1 = pInMin[1] + (pInMax[1] - pInMin[1]) * ((double)i - 0.5) / N;
        double y2 = pInMin[1] + (pInMax[1] - pInMin[1]) * ((double)i + 0.5) / N;
        int rect_x0 = X1_TO_X2(pInMax[0]) + textSize*1.5;
        int dx = textSize;
        QRect rect;
        rect.setCoords(rect_x0 - dx/2, Y1_TO_Y2(y1),
                       rect_x0 + dx/2, Y1_TO_Y2(y2));
        //p.setPen(QColor(0, 0, 0, 0));
        p.setPen(QColor(0, 0, 0, 255));
        p.setBrush(QBrush(genColor(volume, vis2d_parameters->colorModeIndex), Qt::SolidPattern));
        p.drawRect(rect);
        //p.setPen(QColor(0, 0, 0, 255));

        char strt[1000];
        sprintf(strt, "%6.1le", value);
        QFont font;
        font.setPixelSize(textSize);
        p.setFont(font);
        QPointF pText = {rect_x0 + textSize, Y1_TO_Y2((y1 + y2) / 2) + 1*textSize/2};
        p.drawText(pText, strt);
    }

    // окружность
    if(thereIsContact)
    {
        //p.setBrush(QBrush(QColor(255, 255, 255, 255), Qt::SolidPattern));
        p.setBrush(Qt::NoBrush);
        //p.setPen(QColor(0, 0, 0, 255));
        //p.setPen(Qt::DashLine);
        QPen pen(QColor(0, 0, 0, 255), 2);
        pen.setStyle(Qt::DotLine);
        p.setPen(pen);

        int x1 = X1_TO_X2(circleCenter[0] - circleR);
        int x2 = X1_TO_X2(circleCenter[0] + circleR);
        int y1 = Y1_TO_Y2(circleCenter[1] + circleR);
        int y2 = Y1_TO_Y2(circleCenter[1] - circleR);
        p.drawEllipse(x1, y1, x2-x1+1, y2-y1+1);
        //int R_2 = Y1_TO_Y2(-1 - R_1) - Y1_TO_Y2(-1);
        //p.drawEllipse(X1_TO_X2(0), Y1_TO_Y2(-1-R_1), R_2, R_2);
    }

    // отрисовка декомпозиции пространства
    bool paintFiniteElementSurface = testBuilder->needToPaintFiniteElementSurface();
    if(paintFiniteElementSurface)
    {
        //p.setBrush(QBrush(QColor(0, 0, 0, 0), Qt::SolidPattern));
        p.setPen(QColor(255, 0, 0, 255));
        p.setBrush(Qt::NoBrush);
        //FiniteElementSurface &fes = task.grid->FEsurface[3];
        //fes.buildRegions(task.grid);
        //Region &r = *fes.r0;
        //drawRegion(p, r);
        task->grid->FESurface[3].init(true, task->grid);
        task->grid->FESurface[3].update(0, nullptr);
        drawProjections(p, *task, 2, 3);
    }

    if(thereIsContact)
    {
        Surface_base *s = (*task->mechTask.rigidSurface)[(*task->mechTask.CC_FE_Rigid)[0].RigidSurfaceInd];
        drawRegionIndexation(p, s);
    }

}
} //namespace Visualization2d
