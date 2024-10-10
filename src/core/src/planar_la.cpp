#include "planar_la.h"


PlanarLA::PlanarLA(const double _alpha, const double _lambda)
  : alpha(_alpha)
{
  c0 = _alpha * _lambda;
  c1 = 1;
  c2 = -1.0 / ((_alpha - 1.0) * _lambda);
}
