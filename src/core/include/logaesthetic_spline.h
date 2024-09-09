#pragma once

#include "parametric_spline_3d.h"

namespace core {

  class LogAestheticSpline : public ParametricSpline3d
  {
  protected:
    Eigen::Vector3d f(double u) const override;
  };

  
}
