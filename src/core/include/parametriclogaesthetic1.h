#pragma once

#include "parametricsurface.h"

namespace core {

  class ParametricLogAesthetic1 : public ParametricSurface
  {
  protected:
    Eigen::Vector3d f(double u, double v) const override;
  };

  
}
