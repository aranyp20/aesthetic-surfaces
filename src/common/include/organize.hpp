#pragma once

#include <eigen3/Eigen/Dense>
#include <numeric>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <iostream>

namespace common {
namespace color {



  

    inline Eigen::Vector3d hsvToRgb(const Eigen::Vector3d& hsv)
    {
        const auto& h = hsv[0];
        const auto& s = hsv[1];
        const auto& v = hsv[2];
        
        int h_i = int(h * 6);
        float f = h * 6 - h_i;
        float p = v * (1 - s);
        float q = v * (1 - f * s);
        float t = v * (1 - (1 - f) * s);

        switch(h_i) {
            case 0:
                return {v, t, p};
            case 1:
                return {q, v, p};
            case 2:
                return {p, v, t};
            case 3:
                return {p, q, v};
            case 4:
                return {t, p, v};
            default:
                return {v, p, q};
        }
    }



}

    inline Eigen::Vector3d transformAsPoint(const Eigen::Matrix4d& trafo, const Eigen::Vector3d& point)
    {
        const Eigen::Vector4d result = trafo * Eigen::Vector4d(point[0], point[1], point[2], 1);
        return Eigen::Vector3d(result[0], result[1], result[2]);
    }


// average //

  inline Eigen::Vector3d average(const std::vector<Eigen::Vector3d>& v)
  {
    Eigen::Vector3d sum(0,0,0);
    for(const auto u : v) {
      sum+=u;
    }
    return sum/v.size();
  }
////////////

    inline Eigen::Vector3d toEigenV3(const OpenMesh::DefaultTraits::Point& from)
    {
        return Eigen::Vector3d(from[0], from[1], from[2]);
    }

    inline OpenMesh::DefaultTraits::Point toOpenMeshV3(const Eigen::Vector3d& from)
    {
        return OpenMesh::DefaultTraits::Point(from[0], from[1], from[2]);
    }

}
