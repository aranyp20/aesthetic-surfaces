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

    struct LAParams
    {
      double c0, c1, c2, lambda;
      double s0;
      double alpha = common::settings::log_aesthetic_alpha;
    };

    double calcTheta(const LAParams& la, const double s)
    {
      return la.alpha * std::powf(la.c0 * s + la.c1, 1.0 - 1.0/la.alpha) / ((la.alpha - 1.0) * la.c0) + la.c2;
    }
    
    Eigen::Vector2d logAestheticIncrementer(const LAParams& la_params, const double s, const double ds)
    {
      const auto theta = calcTheta(la_params, s);
      return ds * Eigen::Vector2d(std::cosf(theta), std::sinf(theta));
    }
    
    double bisection_computeThetaEDash(const double lambda, const double s0)
    {
      const double ds = 0.01;
      const size_t step_count = s0 / ds;

      Eigen::Vector2d P0(0,0);
      for(size_t i = 0; i < step_count; i++) {
	
      }
    }

    double bisection_computeThetaFDash(const double lambda)
    {
      
    }

    double bisection_computeThetaD(const Eigen::Vector2d& p0,
				   const Eigen::Vector2d& p1,
				   const Eigen::Vector2d& p2)
    {
      //alpha>1 and in form I
      return std::acosf((p0-p1).normalized().dot(Eigen::Vector2d(-1,0)));
    }

    double bisection_computeThetaE(const Eigen::Vector2d& p0,
				   const Eigen::Vector2d& p1,
				   const Eigen::Vector2d& p2)
    {
      //alpha>1 and in form I
      return std::acosf((p2-p0).normalized().dot((p1-p0).normalized()));
    }

    double bisection_computeThetaF(const Eigen::Vector2d& p0,
				   const Eigen::Vector2d& p1,
				   const Eigen::Vector2d& p2)
    {
      //alpha>1 and in form I
      return std::acosf((p0-p2).normalized().dot((p1-p2).normalized()));
    }
    
    double bisection(const double alpha, const int maxIteration, const double thetaD, const double thetaE, const double thetaF, const double s0)
    {
      constexpr double EPS = 1e-5;
      
      double lmin = 0.0,lmax = 1.0, f;
      int i = 0, enlarge = 0;
      if ( alpha == 1.0 ) enlarge = 1;
      else if ( alpha < 1.0 )lmax = 1 / (thetaD * (1 - alpha));
      else if ( alpha > 1.0 )lmax = 1 / (-thetaD * (1 - alpha));
      double lambda = (lmin + lmax) * 0.5;
      do {
	if (alpha <= 1)
	  f = thetaE - bisection_computeThetaEDash(lambda, s0);
	else f = thetaF - bisection_computeThetaFDash(lambda);
	if (std::fabs(f) < EPS)
	  return lambda; // found
	double pLambda = lambda;
	if ((f < 0.0 && alpha <= 1.0) ||(f > 0.0 && alpha > 1.0)) {
	  if (enlarge)
	    lmax *= 10.0;
	  lambda += (lmax - lambda) * 0.5;
	  lmin = pLambda;
	}
	else {
	  enlarge = 0;
	  lambda -= (lambda - lmin) * 0.5;
	  lmax = pLambda;
	}
	i++;
      } while ( i < maxIteration );

      return -1; // not found
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

    const auto& alpha = common::settings::log_aesthetic_alpha;
    
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

    std::complex<double> _P0(0,0);
    std::complex<double> _P2(3,3);
    std::complex<double> _P1(1,0);

    const auto P0 = complex_to_vec2(_P0);
    const auto P1 = complex_to_vec2(_P1);
    const auto P2 = complex_to_vec2(_P2);

    const auto thetaD = bisection_computeThetaD(P0, P1, P2);
    const auto thetaE = bisection_computeThetaE(P0, P1, P2);
    const auto thetaF = bisection_computeThetaF(P0, P1, P2);

    
    double c0, c1, c2;

    double lambda = 1; //todo
    
    c0 = alpha * lambda;
    c1 = 1;
    c2 = 1.0 / ((alpha - 1) * lambda);

      
    const double s0 = (std::powf(((thetaD - c2) * (alpha - 1) * c0 / alpha), 1.0 / (1.0 - 1.0 / alpha)) - c1) / c0;

    lambda = bisection(alpha, 10, thetaD, thetaE, thetaF, s0);


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


