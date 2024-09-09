#include "parametriclogaesthetic1.h"
#include <cmath>


namespace core {

  Eigen::Vector3d ParametricLogAesthetic1::f(double u, double v) const
  {
    constexpr double R = 5;
    constexpr double r = 0.2;
    
    const auto rho = u * 2 * M_PI * 1; 
    const auto theta = v * 2 * M_PI * 0.3;

    const double bp = 0.2;
    const double bg = 0.2;

    
    Eigen::Matrix3d Rz;
    Rz << std::cosf(theta), std::sinf(-theta), 0,
      std::sinf(theta), std::cosf(theta), 0,
      0, 0, 1;

    Eigen::Matrix3d Sc = Eigen::Matrix3d::Identity() * std::powf(M_E, bg * theta);

    Eigen::Matrix<double, 3, 1> M;
    M << R + std::powf(M_E, bp * rho) * std::cosf(rho), 0, std::powf(M_E, bp * rho) * std::sinf(rho);

    return Rz * Sc * M;
  }

  
}
