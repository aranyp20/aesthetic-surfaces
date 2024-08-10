#pragma once
#include <eigen3/Eigen/Dense>


namespace common {
  namespace math {

    //todo too big for inline, put in cpp
    inline double heron(const Eigen::Vector3d& p1,
			const Eigen::Vector3d& p2,
			const Eigen::Vector3d& p3)
    {
      const double a = (p1 - p2).norm();
      const double b = (p1 - p3).norm();
      const double c = (p2 - p3).norm();
      const double s = (a + b + c) / 2;

      return std::sqrt(s * (s-a) * (s-b) * (s-c));
    }


    
  }


  
}
