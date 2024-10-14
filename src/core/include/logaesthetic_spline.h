#pragma once

#include "mesh.h"
#include "parametric_spline_3d.h"

namespace core {

  class LogAestheticSpline : public ParametricSpline3d
  {
  protected:
    Eigen::Vector3d f(double u) const override;
  public:
    common::BaseMesh execute() const;
    common::BaseMesh execute2() const;
    common::BaseMesh execute3() const;
  };

  
}
