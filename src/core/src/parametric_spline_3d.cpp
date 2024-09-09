#include "parametric_spline_3d.h"
#include "mesh.h"

namespace core {

  common::BaseMesh ParametricSpline3d::tessellate(size_t resolution) const
  {
    common::BaseMesh retval;

    const auto dummy_vertex = retval.add_vertex(common::BaseMesh::Point(0,-10,-10));
	  
    constexpr double u_min = 0;
    constexpr double u_max = 1;

    const double u_step_distance = (u_max - u_min) / resolution;

    std::vector<common::BaseMesh::VertexHandle> vhs;
    for (size_t i = 0; i < resolution + 1; i++) {
      const double u = u_min + i * u_step_distance;

	const auto pos = f(u);
	common::BaseMesh::Point om_pos(pos[0], pos[1], pos[2]);
	
	vhs.push_back(retval.add_vertex(om_pos));

	if(i > 0) {
	  std::vector<common::BaseMesh::VertexHandle> vhs2;
	  vhs2.push_back(dummy_vertex);
	  vhs2.push_back(vhs[i-1]);
	  vhs2.push_back(vhs[i]);
	  retval.add_face(vhs2);
	}
    }
    
    return retval;
  }


  
}
