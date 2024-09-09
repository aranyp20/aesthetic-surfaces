#include "logaesthetic_spline.h"
#include <cmath>


namespace core {
  Eigen::Vector3d LogAestheticSpline::f(double u) const
  {
    u *= 2 * M_PI;
    return {0, std::sinf(u), std::cosf(u)};
  }
}


