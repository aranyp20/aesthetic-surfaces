#include "logaesthetic_spline.h"
#include "mesh.h"
#include "settings.h"
#include "type_converter.h"
#include <cmath>
#include <cstddef>
#include <vector>
#include <complex>


namespace core {

  Eigen::Vector3d LogAestheticSpline::f(double u) const
  {
    u *= 2 * M_PI;
    return {0, std::sinf(u), std::cosf(u)};
  }


  namespace {
    double d_theta(const double s)
    {
      const double alpha = common::settings::log_aesthetic_alpha;
      constexpr double c0 = 1;
      constexpr double c1 = 0;

      if(s==0) {
	return 0;
      }
      
      return (alpha * std::pow(c0 * s + c1, 1.0 - 1.0 / alpha)) / ((alpha - 1) * c0);
    }

    double d_theta2(const double s, const double c0, const double c1, const double c2)
    {
      const double alpha = common::settings::log_aesthetic_alpha;

      if(s==0) {
	return 0;
      }
      
      return (alpha * std::pow(c0 * s + c1, 1.0 - 1.0 / alpha)) / ((alpha - 1) * c0) + c2;
    }
    
  }

  
  common::BaseMesh LogAestheticSpline::execute() const
  {
    std::vector<Eigen::Vector2d> spline;

    std::vector<common::BaseMesh::VertexHandle> vhs;

    std::complex<double> P(0,0); //lehet nem is kene complex
    double theta = 0;

    double s = 0;
    const double s_max = 2;
    const double ds = 1e-5;
    const size_t step_count = s_max / ds;


    common::BaseMesh retval;

    const auto dummy_vertex = retval.add_vertex(common::BaseMesh::Point(0,-10,-10));

    vhs.push_back(retval.add_vertex(common::BaseMesh::Point(0, P.real(), P.imag())));

    for(size_t i = 0 ; i < step_count; i++){
      s = i * ds;
      
      P.real(P.real() + std::cosf(theta) * ds);
      P.imag(P.imag() + std::sinf(theta) * ds);

      theta = d_theta(s);
      
      spline.emplace_back(P.real(), P.imag());

      vhs.push_back(retval.add_vertex(common::BaseMesh::Point(0, spline[i][0], spline[i][1])));
    }
    


    
    
    for(size_t i = 1; i < spline.size(); i++) {
      std::vector<common::BaseMesh::VertexHandle> vhs2;
      vhs2.push_back(dummy_vertex);
      vhs2.push_back(vhs[i-1]);
      vhs2.push_back(vhs[i]);
      retval.add_face(vhs2);
    }

    return retval;
  }


  common::BaseMesh LogAestheticSpline::execute2() const
  {
    common::BaseMesh retval;
    const auto dummy_vertex = retval.add_vertex(common::BaseMesh::Point(0,-10,-10));

    
    std::vector<Eigen::Vector2d> spline;

    std::vector<common::BaseMesh::VertexHandle> vhs;

    std::complex<double> P(0,0); //lehet nem is kene complex
    vhs.push_back(retval.add_vertex(common::BaseMesh::Point(0, P.real(), P.imag())));



    const auto complex_to_vec2 = [](const std::complex<double>& c)
    {
      return Eigen::Vector2d(c.real(), c.imag());
    };

    double s = 0;
    const double ds = 1e-5;


    double theta = 0;
    double lambda = 1;

    std::complex<double> _P0(0,0);
    std::complex<double> _P2(3,3);
    std::complex<double> _P1(1,0);

    const double theta_d = std::acos((complex_to_vec2(_P2) - complex_to_vec2(_P1)).normalized().dot(complex_to_vec2(_P1).normalized()));

    const double s0 = 2; //todo

    
    double c0, c1, c2;
    
    for(size_t i = 0; i < 10; i++) {
      c0 = common::settings::log_aesthetic_alpha * lambda;
      c1 = 1;
      c2 = 1.0 / ((common::settings::log_aesthetic_alpha - 1) * lambda);

      
    }
    


    const size_t step_count = s0 / ds;
    
    for(size_t i = 0 ; i < step_count; i++){
      s = i * ds;
      
      
      P.real(P.real() + std::cosf(theta) * ds);
      P.imag(P.imag() + std::sinf(theta) * ds);

      theta = d_theta2(s, c0, c1, c2);
      
      spline.emplace_back(P.real(), P.imag());

      vhs.push_back(retval.add_vertex(common::BaseMesh::Point(0, spline[i][0], spline[i][1])));
    }
    


    
    
    for(size_t i = 1; i < spline.size(); i++) {
      std::vector<common::BaseMesh::VertexHandle> vhs2;
      vhs2.push_back(dummy_vertex);
      vhs2.push_back(vhs[i-1]);
      vhs2.push_back(vhs[i]);
      retval.add_face(vhs2);
    }

    return retval;
  }
}


