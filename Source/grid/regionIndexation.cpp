#define _CRT_SECURE_NO_WARNINGS

#include "regionIndexation.h"

const int SpheredBodyMaxSize = 1;

namespace Grid
{
void RegionWithSpheres::init(const CUBE &set_q0)
{
    q0 = set_q0;                            // инициализация куба
    down = nullptr;                         // изначально область листовая
    bodysNomber = 0;                        // не содержит тел
    // значение s0 не важно, т.к. область не содержит тел
    sb = nullptr;                           // тел нет
}
void RegionWithSpheres::release()
{
    // если в области не содержатся тела, то память освобождать не нужно
    if(bodysNomber == 0)
        return;
    if(down == nullptr) // область листовая?
    {
        // область листовая
        // нужно освободить sb, т.к. в области содержится хотя бы 1 тело
        delete []sb;
    }
    else
    {
        // область не листовая
        // освобождение памяти всех подобластей
        for(int downInd = 0; downInd < 8; downInd++)
            down[downInd].release();
        delete []down;
    }
}
void RegionWithSpheres::insert(const SpheredBody &sbEl)
{
    // область пустая?
    if(bodysNomber == 0)
    {
        // область пустая
        // инициализируем s0
        s0 = sbEl.s;
        // инициализируем массив тел
        sb = new SpheredBody[SpheredBodyMaxSize];
    }
    else
    {
        // область не пустая - расширяем s0
        s0.expand(sbEl.s);
    }
    if(down == nullptr) // область листовая?
    {
        // область листовая
        // не слишком ли много в области тел?
        if(bodysNomber == SpheredBodyMaxSize)
        {
            // слишком много тел
            // разбиваем область на подобласти
            POINT3 O;           // центр q0
            q0.getCenter(O);
            double x[3][3];             // x[k][i] - координаты номер k (0-x, 1-y, 2-z)
            for(int k = 0; k < 3; k++)
            {
                x[k][0] = q0.i[k][0];
                x[k][1] = O[k];
                x[k][2] = q0.i[k][1];
            }
            down = new RegionWithSpheres[8];
            int regionInd = 0;
            for(int i0 = 0; i0 <= 1; i0++)
                for(int i1 = 0; i1 <= 1; i1++)
                    for(int i2 = 0; i2 <= 1; i2++)
                    {
                        CUBE down_q0;
                        down_q0.i[0][0] = x[0][i0];
                        down_q0.i[0][1] = x[0][i0 + 1];
                        down_q0.i[1][0] = x[1][i1];
                        down_q0.i[1][1] = x[1][i1 + 1];
                        down_q0.i[2][0] = x[2][i2];
                        down_q0.i[2][1] = x[2][i2 + 1];
                        down[regionInd].init(down_q0);
                        regionInd++;
                    }
            // подобласти построены и инициализированы
            // распределяем тела исходной области по подобластям
            for(int sbInd = 0; sbInd < bodysNomber; sbInd++)
            {
                for(int downInd = 0; downInd < 8; downInd++)
                {
                    if(down[downInd].q0.intersect(sb[sbInd].s.O))
                    {
                        // центр тела внутри кубика - подобласть найдена
                        down[downInd].insert(sb[sbInd]);
                        break;
                    }
                }
                //addElement(sb[sbInd]);
            }
            for(int downInd = 0; downInd < 8; downInd++)
            {
                if(down[downInd].q0.intersect(sbEl.s.O))
                {
                    // центр тела внутри кубика - подобласть найдена
                    down[downInd].insert(sbEl);
                    break;
                }
            }
            //addElement(sbEl);
            // удаляем тела
            delete[] sb;
            sb = nullptr;
            // теперь область стала не листовая
        }
        else
        {
            // не слишком много тел - просто добавляем новое тело в массив тел
            sb[bodysNomber] = sbEl;
        }
    }
    else
    {
        // область не листовая
        // добавляем тело в одну из подобластей
        for(int downInd = 0; downInd < 8; downInd++)
        {
            if(down[downInd].q0.intersect(sbEl.s.O))
            {
                // центр тела внутри кубика - подобласть найдена
                down[downInd].insert(sbEl);
                break;
            }
        }
    }
    bodysNomber++;
}
void RegionWithSpheres::findIntersectingBodies(const Sphere &s, Grid::SpheredBodyArray &body) const
{
    // если в области не содержатся тела, то пересечений нет
    if(bodysNomber == 0)
        return;
    // если s0 и s не пересекаются, то пересечений нет
    if(!s0.intersect(s))
        return;
    // количество тел != 0 и s пересекается с s0
    if(down == nullptr) // область листовая?
    {
        // область листовая
        // поиск пересечений с каждым из шаров, которые содержатся в области
        for(int sbInd = 0; sbInd < bodysNomber; sbInd++)
        {
            if(sb[sbInd].s.intersect(s))
            {
                // нашли тело, с которым возможно есть пересечение
                // добавляем в массив индекс найденного тела
                body.push_back(sb[sbInd].bodyInd);
            }
        }
    }
    else
    {
        // область не листовая
        // поиск пересечений в подобластях
        for(int downInd = 0; downInd < 8; downInd++)
            down[downInd].findIntersectingBodies(s, body);
    }
}
void RegionWithSpheres::findIntersectingBodies(const POINT3 &p1, const POINT3 &p2, Grid::SpheredBodyArray &body) const
{
    // если в области не содержатся тела, то пересечений нет
    if(bodysNomber == 0)
        return;
    // если s0 и p1-p2 не пересекаются, то пересечений нет
    if(!s0.intersect(p1, p2))
        return;
    // количество тел != 0 и s пересекается с s0
    if(down == nullptr) // область листовая?
    {
        // область листовая
        // поиск пересечений с каждым из шаров, которые содержатся в области
        for(int sbInd = 0; sbInd < bodysNomber; sbInd++)
        {
            if(sb[sbInd].s.intersect(p1, p2))
            {
                // нашли тело, с которым возможно есть пересечение
                // добавляем в массив индекс найденного тела
                body.push_back(sb[sbInd].bodyInd);
            }
        }
    }
    else
    {
        // область не листовая
        // поиск пересечений в подобластях
        for(int downInd = 0; downInd < 8; downInd++)
            down[downInd].findIntersectingBodies(p1, p2, body);
    }
}

bool RegionIndexationVertexIndex_equal::operator()(const RegionIndexation_vertexIndex &vIndex1, const RegionIndexation_vertexIndex &vIndex2) const
{
    return vIndex1[0] == vIndex2[0] && vIndex1[1] == vIndex2[1] && vIndex1[2] == vIndex2[2];
    /*
    if(s1.x[2] < s2.x[2])
        return true;
    if(s1.x[2] > s2.x[2])
        return false;
    if(s1.x[1] < s2.x[1])
        return true;
    if(s1.x[1] > s2.x[1])
        return false;
    if(s1.x[0] < s2.x[0])
        return true;
    if(s1.x[0] > s2.x[0])
        return false;
    return false;   // равенство
    */
}
size_t RegionIndexation_vertexIndex_hash::operator()(const RegionIndexation_vertexIndex &vIndex) const
{
    size_t h0 = std::hash<unsigned int>()(vIndex[0]);
    size_t h1 = std::hash<unsigned int>()(vIndex[1]);
    size_t h2 = std::hash<unsigned int>()(vIndex[2]);
    return h0 ^ ((h1 ^ (h2 << 1)) << 1);
    //return h1 ^ (h2 << 1);
}

void RegionIndexation::init(const RegionDecompositionParameters &param, const RegionIndexation_vertexIndex &set_v0Index, const unsigned int set_h, const Grid::Surface_base *s, const Grid::SurfacePositionData &set_prevSolutionData,
                            Grid::RegionIndexationVertexDataMap *m)
{
    v0Index = set_v0Index;
    h = set_h;
    unsigned int deltaIndex = (1U << (param.h_max - h)); // "длина" стороны если мерить индексами
    double diag = getDiag(param, h);                     // длина диагонали куба

    RegionIndexation_State state0; // состояние нулевой вершины
    bool homogony = true;       // true если состояния всех вершин одинаковые
    bool in = false;            // true если хотя бы 1 шар содержит внутри куб

    // определение данных в вершинах
    // индексы вершин куба: v0[0] + deltaIndex*x0, v0[0] + deltaIndex*x1, v0[0] + deltaIndex*x2
    int localVertexIndex = 0;
    for(unsigned int x2 = 0; x2 <= 1; x2++)
    for(unsigned int x1 = 0; x1 <= 1; x1++)
    for(unsigned int x0 = 0; x0 <= 1; x0++)
    {
        RegionIndexation_vertexIndex vIndex;  // индекс вершины
        RegionIndexation_vertexData vd;       // данные в вершине
        vIndex[2] = v0Index[2] + deltaIndex*x2;
        vIndex[1] = v0Index[1] + deltaIndex*x1;
        vIndex[0] = v0Index[0] + deltaIndex*x0;
        // поиск готовых данных в вершине vIndex или их расчёт и сохранение
        initOrGetvertexData(param, s, vIndex, set_prevSolutionData, m,
                            vd);
        /*RegionIndexationVertexDataMap::iterator it = m->find(vIndex);
        if(it == m->end())
        {
            // информация о вершине с индексом vIndex не найдена
            // расчёт координат вершины
            POINT3 v;
            getVertex(param, vIndex,
                      v);
            // определение информации о ближайшей к v точке поверхности, используя начальное приближение set_ssd
            POINT3 v_nearestPoint;
            POINT3 v_normal;
            int v_side;
            bool v_onBorder;
            Grid::SurfaceSolutionData v_ssd;
            s->findNearestPoint(v, VECTOR3_NULL, set_ssd,
                                v_ssd, v_nearestPoint, v_normal, v_side, v_onBorder);
            // добавление информации о вершине в контейнер
            vd.R = (v - v_nearestPoint).abs();
            vd.state.side = v_side;
            vd.state.onBorder = v_onBorder;
            vd.ssd = v_ssd;
            std::pair<const RegionIndexation_vertexIndex, RegionIndexation_vertexData> m_value = {vIndex, vd};
            m->insert(m_value);
        }
        else
        {
            // информация о вершине с индексом vIndex найдена
            vd = (*it).second;
        }*/
        // итак, теперь есть информация vd о вершине с индексом vIndex
        if(localVertexIndex == 0)
        {
            state0 = vd.state;
        }
        if(homogony)
        {
            if(vd.state == state0)
            {
                if(vd.R >= diag)
                {
                    // куб находится внутри шара
                    in = true;
                }
            }
            else
            {
                // не во всех вершинах куба состояния одинаковые
                homogony = false;
            }
        }
        localVertexIndex++;
    }
    // определение общего состояния куба
    if(homogony && in)
    {
        // состояние любой точки внутри куба определено
        defined = true;
        state = state0;
    }
    else
    {
        // состояние точек внутри куба могут отличаться
        defined = false;
    }

    if((defined == false && h < param.h_max) || (h < param.h_min))
    {
        // куб делится на 8 частей
        RegionIndexation_vertexIndex vCenterIndex;  // индекс вершины в центре куба
        RegionIndexation_vertexData vCenterData;    // данные в вершине
        vCenterIndex[0] = v0Index[0] + deltaIndex/2;
        vCenterIndex[1] = v0Index[1] + deltaIndex/2;
        vCenterIndex[2] = v0Index[2] + deltaIndex/2;
        initOrGetvertexData(param, s, vCenterIndex, set_prevSolutionData, m,
                            vCenterData);
        if(h < param.h_min)
            vCenterData.ssd.index = -1;     // начальные приближения начинают использоваться только начиная с минимальной глубины
        down = new RegionIndexation[8];
        int localVertexIndex = 0;
        for(unsigned int x2 = 0; x2 <= 1; x2++)
        for(unsigned int x1 = 0; x1 <= 1; x1++)
        for(unsigned int x0 = 0; x0 <= 1; x0++)
        {
            RegionIndexation_vertexIndex v0Index_new;  // индекс вершины, который определяет 8-ю часть куба
            v0Index_new[2] = v0Index[2] + deltaIndex/2*x2;
            v0Index_new[1] = v0Index[1] + deltaIndex/2*x1;
            v0Index_new[0] = v0Index[0] + deltaIndex/2*x0;
            down[localVertexIndex].init(param, v0Index_new, h + 1, s, vCenterData.ssd,
                                        m);
            localVertexIndex++;
        }
    }
    else
    {
        // листовой узел
        down = nullptr;
    }
}
void RegionIndexation::release()
{
    if(down != nullptr)
    {
        for(int localVertexIndex = 0; localVertexIndex < 8; localVertexIndex++)
        {
            down[localVertexIndex].release();
        }
    }
}
void RegionIndexation::getInformation(const RegionDecompositionParameters &param, const Grid::RegionIndexationVertexDataMap *m, const POINT3 &p,
                                      bool &p_defined, int &p_side, bool &p_onBorder, Grid::SurfacePositionData &approximateSolutionData)const
{
    if(defined)
    {
        // состояние куба определено, возвращаем состояние
        p_defined = true;
        p_side = state.side;
        p_onBorder = state.onBorder;
    }
    else
    {
        // состояние куба не определено
        if(down == nullptr) // область листовая?
        {
            // область листовая
            // состояние определить не удалось, возвращаем первое приближение
            p_defined = false;
            RegionIndexationVertexDataMap::const_iterator it = m->find(v0Index);
            approximateSolutionData = (*it).second.ssd;
        }
        else
        {
            // область не листовая
            // продолжаем поиск состояния точки p в одной из подобластей
            for(int downInd = 0; downInd < 8; downInd++)
            {
                CUBE cube;      // куб - подобласть
                down[downInd].getCube(param,
                        cube);
                if(cube.intersect(p))
                {
                    // точка внутри подобласти - подобласть найдена
                    down[downInd].getInformation(param, m, p,
                                                 p_defined, p_side, p_onBorder, approximateSolutionData);
                    break;
                }
            }
        }
    }
}
void RegionIndexation::getVertex(const RegionDecompositionParameters &param, const RegionIndexation_vertexIndex &vIndex,
                                 POINT3 &v)const
{
    v[0] = param.q0.i[0][0] + (param.q0.i[0][1] - param.q0.i[0][0])/(1U << param.h_max)*vIndex[0];
    v[1] = param.q0.i[1][0] + (param.q0.i[1][1] - param.q0.i[1][0])/(1U << param.h_max)*vIndex[1];
    v[2] = param.q0.i[2][0] + (param.q0.i[2][1] - param.q0.i[2][0])/(1U << param.h_max)*vIndex[2];
}
void RegionIndexation::getCube(const RegionDecompositionParameters &param,
                               CUBE &cube)const
{
    unsigned int deltaIndex = (1U << (param.h_max - h));
    RegionIndexation_vertexIndex vIndex;  // индекс вершины куба
    POINT3 v;                             // координаты вершины куба
    // v0 - индекс вершины с наименьшими координатами
    getVertex(param, v0Index,
              v);
    cube.i[0][0] = v[0];
    cube.i[1][0] = v[1];
    cube.i[2][0] = v[2];
    // vIndex - индекс вершины с наибольшими координатами
    vIndex[0] = v0Index[0] + deltaIndex;
    vIndex[1] = v0Index[1] + deltaIndex;
    vIndex[2] = v0Index[2] + deltaIndex;
    getVertex(param, vIndex,
              v);
    cube.i[0][1] = v[0];
    cube.i[1][1] = v[1];
    cube.i[2][1] = v[2];
}
double RegionIndexation::getDiag(const RegionDecompositionParameters &param, const unsigned int h)const
{
    double a = (param.q0.i[0][1] - param.q0.i[0][0]) / (1U << h);
    double b = (param.q0.i[1][1] - param.q0.i[1][0]) / (1U << h);
    double c = (param.q0.i[2][1] - param.q0.i[2][0]) / (1U << h);
    return sqrt(a*a + b*b + c*c);
}
void RegionIndexation::initOrGetvertexData(const RegionDecompositionParameters &param, const Grid::Surface_base *s, const RegionIndexation_vertexIndex &vIndex, const SurfacePositionData &prevSolutionData,
                                           Grid::RegionIndexationVertexDataMap *m, RegionIndexation_vertexData &vd)const
{
    RegionIndexationVertexDataMap::iterator it = m->find(vIndex);
    if(it == m->end())
    {
        // информация о вершине с индексом vIndex не найдена
        // расчёт координат вершины
        POINT3 v;
        getVertex(param, vIndex,
                  v);
        // определение информации о ближайшей к v точке поверхности, используя начальное приближение prevSolutionData
        POINT3 v_nearestPoint;
        POINT3 v_normal;
        int v_side;
        bool v_onBorder;
        Grid::SurfacePositionData v_ssd;
        s->findNearestPoint(v, VECTOR3_NULL, prevSolutionData,
                            v_ssd, v_nearestPoint, v_normal, v_side, v_onBorder);
        // добавление информации о вершине в контейнер
        vd.R = (v - v_nearestPoint).abs();
        vd.state.side = v_side;
        vd.state.onBorder = v_onBorder;
        vd.ssd = v_ssd;
        std::pair<const RegionIndexation_vertexIndex, RegionIndexation_vertexData> m_value = {vIndex, vd};
        m->insert(m_value);
    }
    else
    {
        // информация о вершине с индексом vIndex найдена
        vd = (*it).second;
    }
}
}   // namespace Grid
