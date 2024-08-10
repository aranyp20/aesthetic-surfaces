#pragma once
#include <eigen3/Eigen/Dense>
#include "mesh.h"

namespace common {
  namespace converter {

    inline Eigen::Vector3d meshPointToEigen(const common::MyMesh::Point& vp)
    {
      return {vp[0], vp[1], vp[2]};
    }
    
  }
}
