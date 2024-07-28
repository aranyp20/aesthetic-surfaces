#include "subdivider.h"
#include "common_defines.h"





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
  
  
  // Vertexhandle not present as child in this map must be original vertex
  //TODO refactor
  void Subdivider::execute(common::MyMesh& mesh, ChildrenParents& children_parents_map) const
  {
    ChildrenParents retval;
  
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
