#pragma once

#include <map>
#include <string>

namespace common {
  namespace settings{

    enum Algorithm : unsigned int
      {
	BASIC = 0u,
	BEZIER = 1u,
	LOG_AESTHETIC = 2u
      };

    inline const std::map<Algorithm, std::string> alg_name {
      {BASIC, "Basic"},
      {BEZIER, "Bezier"},
      {LOG_AESTHETIC, "Log-aesthetic"}
    };

    enum CurvatureType
      {
	MEAN = 0u,
	GAUSSIAN = 1u
      };
    inline const std::map<CurvatureType, std::string> curvature_name {
      {MEAN, "Mean"},
      {GAUSSIAN, "Gaussian"}
    };
    
    inline double log_aesthetic_alpha = 2.1;
    inline Algorithm selected_alg = Algorithm::LOG_AESTHETIC;
    inline CurvatureType selected_curvature = CurvatureType::MEAN;
    inline bool sync = false;
    inline bool show_vertex_ids = false;

    inline std::string default_model = "deformed_ico.obj";
  }
}
