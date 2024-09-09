#include "subdivider.h"
#include "common_defines.h"
#include "mesh.h"
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


  common::MyMesh::VertexHandle deb1;
  
  void Subdivider::processQuadFace(common::MyMesh::FaceHandle& face, common::MyMesh& mesh, ChildrenParents& children_parents_map, const std::map<common::MyMesh::EdgeHandle, common::MyMesh::VertexHandle>& halfpoint_map) const
  {
    std::vector<common::MyMesh::VertexHandle> border_vertices; /* layout: original-halfpoint-original...halfpoint */
    for(const common::MyMesh::HalfedgeHandle f_hedge : mesh.fh_range(face)) {
      border_vertices.push_back(mesh.from_vertex_handle(f_hedge));
      border_vertices.push_back(halfpoint_map.at(mesh.edge_handle(f_hedge)));
    }

    common::MyMesh::Point middle_point = (mesh.point(border_vertices[1]) + mesh.point(border_vertices[5])) / 2.0;
    common::MyMesh::VertexHandle middle_vertex = mesh.add_vertex(middle_point);

    deb1=middle_vertex;

    mesh.delete_face(face, false);

    mesh.add_face({border_vertices[7], border_vertices[0], border_vertices[1], middle_vertex});
    mesh.add_face({border_vertices[1], border_vertices[2], border_vertices[3], middle_vertex});
    mesh.add_face({border_vertices[3], border_vertices[4], border_vertices[5], middle_vertex});
    mesh.add_face({border_vertices[5], border_vertices[6], border_vertices[7], middle_vertex});


    mesh.children_parents_map.insert({middle_vertex, {border_vertices[1], border_vertices[5]}});
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
      children_parents_map.insert({vnew, {v1, v2}});
      mesh.edge_neighbor_map.insert({vnew, {v1, v2}});

      if(mesh.edge_neighbor_map.count(v1) && mesh.edge_neighbor_map.at(v1)[0] == v2) {
	mesh.edge_neighbor_map[v1][0] = vnew;
      }
      if(mesh.edge_neighbor_map.count(v1) && mesh.edge_neighbor_map.at(v1)[1] == v2) {
	mesh.edge_neighbor_map[v1][1] = vnew;
      }
      if(mesh.edge_neighbor_map.count(v2) && mesh.edge_neighbor_map.at(v2)[0] == v1) {
	mesh.edge_neighbor_map[v2][0] = vnew;
      }
      if(mesh.edge_neighbor_map.count(v2) && mesh.edge_neighbor_map.at(v2)[1] == v1) {
	mesh.edge_neighbor_map[v2][1] = vnew;
      }
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

    
    
    for (auto& de : original_edges) {

      mesh.delete_edge(de, false);
    }
    
    mesh.garbage_collection();
    
    return;
    
    
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


  void Subdivider::execute(common::MyMesh& mesh, const size_t iteration_count) const
  {
    for (size_t i = 0; i < iteration_count; i++) {
      execute(mesh, mesh.children_parents_map);
    }
  }


  
} //namespace core
