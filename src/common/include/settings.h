#pragma once

#include <map>

namespace common {
  namespace settings{

    enum Algorithm
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
    
    inline double log_aesthetic_alpha = -1.0;
    inline Algorithm selected_alg = Algorithm::BASIC;
    inline bool sync = false;
    inline bool show_vertex_ids = true;
  }
}
