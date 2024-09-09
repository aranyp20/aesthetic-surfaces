#include "parametricsurface.h"
#include "mesh.h"



namespace core {

  common::BaseMesh ParametricSurface::tessellate(size_t resolution) const
  {
    //resolution = 3;
    
    common::BaseMesh retval;

	  
    constexpr double u_min = 0;
    constexpr double u_max = 1;
    constexpr double v_min = 0;
    constexpr double v_max = 1;

    const double u_step_distance = (u_max - u_min) / resolution;
    const double v_step_distance = (v_max - v_min) / resolution;

    std::vector<common::BaseMesh::VertexHandle> vhs;
    for (size_t i = 0; i < resolution + 1; i++) {
      const double u = u_min + i * u_step_distance;
      for (size_t j = 0; j < resolution + 1; j++) {
	const double v = v_min + j * v_step_distance;

	const auto pos = f(u, v);
	common::BaseMesh::Point om_pos(pos[0], pos[1], pos[2]);
	
	vhs.push_back(retval.add_vertex(om_pos));
      }
    }

    // TODO refactor
    
    for (size_t i = 0; i < (resolution) * (resolution); i++) {
      const size_t left_right_index = std::floor((double)i / resolution) + i;

      //todo order might not be good so normals might be flipped
      std::vector<common::BaseMesh::VertexHandle> f_vhs;
      f_vhs.push_back(vhs[left_right_index]);
      f_vhs.push_back(vhs[left_right_index + 1]);
      f_vhs.push_back(vhs[left_right_index + resolution + 2]);
      f_vhs.push_back(vhs[left_right_index + resolution + 1]);
      retval.add_face(f_vhs);
    }

    return retval;
  }
  
}
