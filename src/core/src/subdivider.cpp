#include "subdivider.h"
#include "common_defines.h"
#include <stdexcept>
#include <sys/types.h>




namespace core {

  common::MyMesh::EdgeHandle Subdivider::findEdgeConnectingVertices(common::MyMesh& mesh, common::MyMesh::VertexHandle v1, common::MyMesh::VertexHandle v2) const {
    for (common::MyMesh::ConstEdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
      common::MyMesh::HalfedgeHandle heh0 = mesh.halfedge_handle(*e_it, 0);
      common::MyMesh::HalfedgeHandle heh1 = mesh.halfedge_handle(*e_it, 1);
      common::MyMesh::VertexHandle v0 = mesh.to_vertex_handle(heh0);
      common::MyMesh::VertexHandle v3 = mesh.to_vertex_handle(heh1);
      if ((v1 == v0 && v2 == v3) || (v1 == v3 && v2 == v0)) {
	// Found the edge connecting the two vertices
	return *e_it;
      }
    }
    
    // If no such edge found, return an invalid handle
    return common::MyMesh::EdgeHandle();
  }


  
  void Subdivider::processQuadFace(common::MyMesh::FaceHandle& face, common::MyMesh& mesh, ChildrenParents& children_parents_map, const std::map<common::MyMesh::EdgeHandle, common::MyMesh::VertexHandle>& halfpoint_map) const
  {
    std::vector<common::MyMesh::VertexHandle> border_vertices; /* layout: original-halfpoint-original...halfpoint */
    for(const common::MyMesh::HalfedgeHandle f_hedge : mesh.fh_range(face)) {
      border_vertices.push_back(mesh.from_vertex_handle(f_hedge));
      border_vertices.push_back(halfpoint_map.at(mesh.edge_handle(f_hedge)));
    }

    common::MyMesh::Point middle_point = (mesh.point(border_vertices[1]) + mesh.point(border_vertices[5])) / 2.0;
    common::MyMesh::VertexHandle middle_vertex = mesh.add_vertex(middle_point);

    children_parents_map[middle_vertex] = {border_vertices[1], border_vertices[5]};

    mesh.delete_face(face, false);

    mesh.add_face({border_vertices[7], border_vertices[0], border_vertices[1], middle_vertex});
    mesh.add_face({border_vertices[1], border_vertices[2], border_vertices[3], middle_vertex});
    mesh.add_face({border_vertices[3], border_vertices[4], border_vertices[5], middle_vertex});
    mesh.add_face({border_vertices[5], border_vertices[6], border_vertices[7], middle_vertex});
  }


  void Subdivider::processTriFace(common::MyMesh::FaceHandle& face, common::MyMesh& mesh, const std::map<common::MyMesh::EdgeHandle, common::MyMesh::VertexHandle>& halfpoint_map) const
  {
    std::vector<common::MyMesh::VertexHandle> border_vertices; /* layout: original-halfpoint-original...halfpoint */
    for(const common::MyMesh::HalfedgeHandle f_hedge : mesh.fh_range(face)) {
      border_vertices.push_back(mesh.from_vertex_handle(f_hedge));
      border_vertices.push_back(halfpoint_map.at(mesh.edge_handle(f_hedge)));
    }

    mesh.delete_face(face, false);

    mesh.add_face({border_vertices[0], border_vertices[1], border_vertices[5]});
    mesh.add_face({border_vertices[1], border_vertices[2], border_vertices[3]});
    mesh.add_face({border_vertices[3], border_vertices[4], border_vertices[5]});
    mesh.add_face({border_vertices[1], border_vertices[3], border_vertices[5]});
  }

  
  // Vertexhandle not present as child in this map must be original vertex
  //TODO refactor
  void Subdivider::execute(common::MyMesh& mesh, ChildrenParents& children_parents_map) const
  {
    using namespace common;

    std::vector<common::MyMesh::EdgeHandle> original_edges;
    for (common::MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
      original_edges.push_back(*e_it);
    }

    std::map<common::MyMesh::EdgeHandle, common::MyMesh::VertexHandle> halfpoint_map;

    for (const auto& original_edge : original_edges) {
      
      MyMesh::HalfedgeHandle heh = mesh.halfedge_handle(original_edge, 0);
      MyMesh::HalfedgeHandle oheh = mesh.opposite_halfedge_handle(heh);
      // Get the vertices
      MyMesh::VertexHandle v1 = mesh.to_vertex_handle(heh);
      MyMesh::VertexHandle v2 = mesh.to_vertex_handle(oheh);

      // Create the new vertex at the midpoint
      MyMesh::Point new_point = (mesh.point(v1) + mesh.point(v2)) / 2.0;
      MyMesh::VertexHandle vnew = mesh.add_vertex(new_point);
      halfpoint_map[original_edge] = vnew;
      children_parents_map[vnew] = {v1, v2};
    }

    std::vector<common::MyMesh::FaceHandle> original_faces;
    for (const auto& face : mesh.faces()) {
      original_faces.push_back(face);
    }

    for (auto& original_face : original_faces) {
      switch (mesh.valence(original_face)) {
      case 3: {
	processTriFace(original_face, mesh, halfpoint_map);
	break;
      }
      case 4: {
	processQuadFace(original_face, mesh, children_parents_map, halfpoint_map);
	break;
      }
      default: {
	throw std::runtime_error("Subdivision undefined for a face in the processed mesh.");
      }
      }
    }

    for (auto& original_edge : original_edges) {
      //mesh.delete_edge(original_edge, false);
    }
    mesh.garbage_collection();
    
    std::cout << mesh.n_faces() << std::endl;
    std::cout << mesh.n_edges() << std::endl;
    std::cout << mesh.n_vertices() << std::endl;
    
    return;
    //////////////////////////////////////////
    
    for (auto& original_edge : original_edges) {
      if(!(mesh.face_handle(mesh.halfedge_handle(original_edge,0)).is_valid() && mesh.face_handle(mesh.halfedge_handle(original_edge,1)).is_valid())) {
	continue;
      }

      
      MyMesh::HalfedgeHandle heh = mesh.halfedge_handle(original_edge, 0);
      MyMesh::HalfedgeHandle oheh = mesh.opposite_halfedge_handle(heh);

      // Get the vertices
      MyMesh::VertexHandle v1 = mesh.to_vertex_handle(heh);
      MyMesh::VertexHandle v2 = mesh.to_vertex_handle(oheh);

      // Create the new vertex at the midpoint
      MyMesh::Point new_point = (mesh.point(v1) + mesh.point(v2)) / 2.0;
      MyMesh::VertexHandle vnew = mesh.add_vertex(new_point);

      // Get the faces
      MyMesh::FaceHandle fh1 = mesh.face_handle(heh);
      MyMesh::FaceHandle fh2 = mesh.face_handle(oheh);

      // Store the vertices of the faces
      std::vector<MyMesh::VertexHandle> quad_vertices;
      std::vector<MyMesh::VertexHandle> tri_vertices;
 
      // Store the vertices of the quad
      if (mesh.valence(fh1) == 3) {
        for (MyMesh::VertexHandle vh : mesh.fv_range(fh1)) {
	  if (vh != v1 && vh != v2) {
	    tri_vertices.push_back(vh);
	  }
        }
	tri_vertices.push_back(v1);
	tri_vertices.push_back(v2);
      }

      // Store the vertices of the triangle
      if (mesh.valence(fh2) == 4) {
	quad_vertices.push_back(v1);
	quad_vertices.push_back(v2);


	for (MyMesh::VertexHandle vh : mesh.fv_range(fh2)) {
	  if (vh != v1 && vh != v2) {
	    quad_vertices.push_back(vh);
	  }
        }
      }



      // Delete the old faces
      mesh.delete_face(fh1, false);
      mesh.delete_face(fh2, false);

      MyMesh::Point new_point_opposite = (mesh.point(quad_vertices[2]) + mesh.point(quad_vertices[3])) / 2.0;
      MyMesh::VertexHandle vnew_opposite = mesh.add_vertex(new_point_opposite);
      
      // Add the new faces
      std::vector<MyMesh::VertexHandle> new_face1 = {quad_vertices[2], vnew_opposite, vnew, quad_vertices[1]};
      std::vector<MyMesh::VertexHandle> new_face2 = {vnew_opposite, quad_vertices[3], quad_vertices[0], vnew};
      std::vector<MyMesh::VertexHandle> new_face3 = {tri_vertices[2], vnew, tri_vertices[0]};
      std::vector<MyMesh::VertexHandle> new_face4 = {vnew, tri_vertices[1], tri_vertices[0]};
      mesh.add_face(new_face4);
      mesh.add_face(new_face3);
      mesh.add_face(new_face2);
      mesh.add_face(new_face1);

      mesh.garbage_collection();
      
      std::cout << mesh.n_faces() << std::endl;
      std::cout << mesh.n_edges() << std::endl;
      std::cout << mesh.n_vertices() << std::endl;
    }

    /////////////////////////////////

    /*
    std::vector<common::MyMesh::EdgeHandle> original_edges;
    for (common::MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
      original_edges.push_back(*e_it);
    }

    std::vector<common::MyMesh::VertexHandle> new_vertices;
    std::vector<common::MyMesh::EdgeHandle> flippable_edges;

    for(auto& original_edge_h : original_edges) {
      auto heh1 = mesh.halfedge_handle(original_edge_h, 0);
      auto heh2 = mesh.halfedge_handle(original_edge_h, 1);


      auto face1 = mesh.face_handle(heh1);
      auto face2 = mesh.face_handle(heh2);

	

      auto v1 = heh1.is_valid() ? mesh.from_vertex_handle(heh1) : mesh.from_vertex_handle(heh2);
      auto v2 = heh1.is_valid() ? mesh.to_vertex_handle(heh1) : mesh.to_vertex_handle(heh2);

      common::MyMesh::VertexHandle opposite_vertex1;
      common::MyMesh::VertexHandle opposite_vertex2;


      if (face1.is_valid()) {
	for (common::MyMesh::FaceVertexIter fv_it = mesh.fv_iter(face1); fv_it.is_valid(); ++fv_it) {
	  common::MyMesh::VertexHandle v_face = *fv_it;

	  if (v_face != v1 && v_face != v2) {
	    opposite_vertex1 = v_face;
	    break;
	  }
	}
      }

      if (face2.is_valid()) {
	for (common::MyMesh::FaceVertexIter fv_it = mesh.fv_iter(face2); fv_it.is_valid(); ++fv_it) {
	  common::MyMesh::VertexHandle v_face = *fv_it;

	  if (v_face != v1 && v_face != v2) {
	    opposite_vertex2 = v_face;
	    break;
	  }
	}
      }


      auto new_vertex = mesh.split(original_edge_h, mesh.calc_edge_midpoint(original_edge_h));

      new_vertices.push_back(new_vertex);

      children_parents_map.insert({new_vertex, {v1, v2}});
	

      if(face1.is_valid()){
	if(std::find(new_vertices.begin(), new_vertices.end(), opposite_vertex1) == new_vertices.end()){
	  flippable_edges.push_back(findEdgeConnectingVertices(mesh, new_vertex, opposite_vertex1));     

	}
      }

      if(face2.is_valid()){
	if(std::find(new_vertices.begin(), new_vertices.end(), opposite_vertex2) == new_vertices.end()){
	  flippable_edges.push_back(findEdgeConnectingVertices(mesh, new_vertex, opposite_vertex2));     

	}
      }
    }

    for (auto& flippable_edge : flippable_edges) {


      mesh.flip(flippable_edge);

    }
    */
  }


  Subdivider::ChildrenParents Subdivider::execute(common::MyMesh& mesh, const size_t iteration_count) const
  {
    ChildrenParents retval;
    for (size_t i = 0; i < iteration_count; i++) {
      execute(mesh, retval);
    }
    return retval;
  }


  
} //namespace core
