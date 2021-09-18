/* --------------------------------------------------------- */
// ОБРАБОТКА РЕЗУЛЬТАТОВ
/* --------------------------------------------------------- */

#ifndef POST_H
#define POST_H

#include "elementary.h"
#include "solid_base.h"

namespace Post
{
using namespace Elementary;

// дан массив узловых сил (узлы должны принадлежать поверхности)
// возвращает коэффициенты разложения давления по базису (ненулевые коэффициенты только перед базисными функциями, которые на данной поверхности могут становиться ненулевыми)
void calcPressureByNodeForces(const Grid::FiniteElementSurface &contactSurface, const Grid::Grid3D &grid, const std::vector<Solid::MechOutVertexData> &vertexData,
                               std::vector<VECTOR3> &vertexP);
void calcContactArea(const Grid::FiniteElementSurface &contactSurface, const Grid::Grid3D &grid, const std::vector<Solid::MechOutVertexData> &vertexData,
                      std::vector<VECTOR3> &vertexP, double &S_true, double &S_geom, double &S_shamanstvo);

}   // namespace Post

#endif // POST_H
