#include "logaesthetic_spline.h"
#include "mesh.h"
#include "settings.h"
#include "type_converter.h"
#include <cmath>
#include <cstddef>
#include <vector>
#include <complex>
#include "planar_la.h"

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


    double calcTheta(const PlanarLA& la, const double s)
    {
      return la.alpha * std::powf(la.c0 * s + la.c1, 1.0 - 1.0/la.alpha) / ((la.alpha - 1.0) * la.c0) + la.c2;
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
    
    Eigen::Vector2d logAestheticIncrementer(const PlanarLA& la_params, const double s, const double ds)
    {
      const auto theta = calcTheta(la_params, s);
      return ds * Eigen::Vector2d(std::cosf(theta), std::sinf(theta));
    }

    double computeS0(const double lambda, const double thetaD)
    {
      const PlanarLA la(common::settings::log_aesthetic_alpha, lambda);
      const double exponent = -1. + 1./la.alpha;
      const double base = (thetaD - la.c2) * (la.alpha - 1) * la.c0 / la.alpha;
      return (std::powf(base, exponent) - la.c1) / la.c0;
	
    }

    Eigen::Vector2d computeP1(const Eigen::Vector2d P2, const double thetaD)
    {
      return Eigen::Vector2d((P2 + Eigen::Vector2d(std::cosf(thetaD), std::sinf(thetaD)) * (-P2[1] / std::sinf(thetaD)))[0], 0);
    }

    std::array<Eigen::Vector2d, 3> computeTriangle(const double lambda, const double thetaD)
    {
      const double s0 = computeS0(lambda, thetaD);

      const PlanarLA la(common::settings::log_aesthetic_alpha, lambda);

      const double ds = 0.001;
      const size_t step_count = s0 / ds;


      Eigen::Vector2d P2(0,0);

      double s = 0;
      for(size_t i = 0; i < step_count; i++) {
	s = i * ds;
	P2 += logAestheticIncrementer(la, s, ds);
      }

      const auto P1 = computeP1(P2, thetaD);

      const Eigen::Vector2d P0(0,0);

      return {P0, P1, P2};
    }
    
    double bisection_computeThetaEDash(const double lambda, const double thetaD)
    {
      const auto triangle = computeTriangle(lambda, thetaD);

      return bisection_computeThetaE(triangle[0], triangle[1], triangle[2]);
    }

    double bisection_computeThetaFDash(const double lambda, const double thetaD)
    {
      const auto triangle = computeTriangle(lambda, thetaD);

      return bisection_computeThetaF(triangle[0], triangle[1], triangle[2]);
    }
    
    double bisection(const double alpha, const int maxIteration, const double thetaD, const double thetaE, const double thetaF)
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
	  f = thetaE - bisection_computeThetaEDash(lambda, thetaD);
	else f = thetaF - bisection_computeThetaFDash(lambda, thetaD);
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

    vhs.push_back(retval.add_vertex(common::BaseMesh::Point(0, 0, 0)));



    const auto complex_to_vec2 = [](const std::complex<double>& c)
    {
      return Eigen::Vector2d(c.real(), c.imag());
    };

    double s = 0;
    const double ds = 1e-2;


    double theta = 0;

    std::complex<double> _P0(-3,3);
    std::complex<double> _P2(0,0);
    std::complex<double> _P1(-1,0);

    const auto P0 = complex_to_vec2(_P0);
    const auto P1 = complex_to_vec2(_P1);
    const auto P2 = complex_to_vec2(_P2);

    const auto thetaD = bisection_computeThetaD(P0, P1, P2);
    const auto thetaE = bisection_computeThetaE(P0, P1, P2);
    const auto thetaF = bisection_computeThetaF(P0, P1, P2);

    

    

      

    const auto lambda = 5.0;// bisection(alpha, 10, thetaD, thetaE, thetaF);
    const PlanarLA la(alpha, lambda);
    const auto s0 = 1.0;//computeS0(lambda, thetaD);
    const size_t step_count = s0 / ds;
    std::cout << thetaD<<" "<<computeS0(lambda, thetaD)<<" "<<step_count << std::endl;
    Eigen::Vector2d eval_P;
    
    for(size_t i = 0 ; i < step_count; i++){
      s = i * ds;
      

      eval_P += logAestheticIncrementer(la, s, ds);

      spline.push_back(eval_P);

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


  common::BaseMesh uvPlaneVhsToMesh(const std::vector<std::vector<common::BaseMesh::VertexHandle>>& vhs)
  {
    common::BaseMesh retval;

    return retval;
  }
  
  common::BaseMesh integralPointsToMesh(const std::vector<std::vector<Eigen::Vector3d>>& ip, const size_t width, const size_t height)
  {
    common::BaseMesh retval;

    //width<height>

    constexpr size_t resoltuion = 20;
    const double width_step = (double)width / resoltuion;
    const double height_step = (double)height / resoltuion; 

    std::vector<std::vector<common::BaseMesh::VertexHandle>>
      vhs(resoltuion, std::vector<common::BaseMesh::VertexHandle>(resoltuion));

    for(size_t i = 0; i < resoltuion; i++) {
      for(size_t j = 0; j < resoltuion; j++) {
	const Eigen::Vector3d& c_ip = ip[(size_t)(i * width_step)][(size_t)(j *height_step)];
	vhs[i][j] = retval.add_vertex(common::BaseMesh::Point(c_ip[0], c_ip[1], c_ip[2]));
      }
    }

    

    return retval;
  }
  
  common::BaseMesh LogAestheticSpline::execute3() const
  {
    common::BaseMesh retval;
    const auto& alpha = common::settings::log_aesthetic_alpha;
    


    const auto complex_to_vec2 = [](const std::complex<double>& c)
    {
      return Eigen::Vector2d(c.real(), c.imag());
    };

    double s_a = 0;
    double s_b = 0;
    const double ds = 1e-2;



    std::complex<double> _P0(-3,3);
    std::complex<double> _P2(0,0);
    std::complex<double> _P1(-1,0);

    const auto P0 = complex_to_vec2(_P0);
    const auto P1 = complex_to_vec2(_P1);
    const auto P2 = complex_to_vec2(_P2);

    const auto thetaD = bisection_computeThetaD(P0, P1, P2);
    const auto thetaE = bisection_computeThetaE(P0, P1, P2);
    const auto thetaF = bisection_computeThetaF(P0, P1, P2);

    

    

      

    const auto lambda = 5.0;// bisection(alpha, 10, thetaD, thetaE, thetaF);
    const PlanarLA la(alpha, lambda);
    const auto s0 = 1.0;//computeS0(lambda, thetaD);
    const size_t step_count = s0 / ds;

    Eigen::Vector2d eval_P_a(0,0);
    Eigen::Vector2d eval_P_b(0,0);


    const Eigen::Vector3d u(1,0,0);
    const Eigen::Vector3d v(0,1,0);
    const Eigen::Vector3d uxv(0,0,1);

    const auto to3d = [](const Eigen::Vector3d& i_hat, const Eigen::Vector3d& j_hat, const Eigen::Vector2d& p)
    {
      return i_hat * p[0] + j_hat * p[1];
    };

    std::vector<std::vector<Eigen::Vector3d>> all_points(step_count, std::vector<Eigen::Vector3d>(step_count));

    
    for(size_t i = 0 ; i < step_count; i++){
      s_a = i * ds;
      
      eval_P_a += logAestheticIncrementer(la, s_a, ds);

      eval_P_b = Eigen::Vector2d(0,0);
      
      for(size_t j = 0; j < step_count; j++) {
	s_b = j * ds;

	eval_P_b += logAestheticIncrementer(la, s_b, ds);

	all_points[i][j] = to3d(u, v, eval_P_a) + to3d(u, uxv, eval_P_b);
	
      }

    }
    


    
    return retval;
  }
}


