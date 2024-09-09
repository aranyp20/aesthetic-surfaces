#pragma once

#include <eigen3/Eigen/Dense>
#include "mesh.h"


namespace core {

  class ParametricSpline3d
  {
  protected:
    virtual Eigen::Vector3d f(double u) const = 0;
  public:
    common::BaseMesh tessellate(size_t resolution) const;
  };

}
