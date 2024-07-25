#include "discretefairer.h"
#include <OpenMesh/Tools/Subdivider/Uniform/MidpointT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "common_defines.h"
#include "curvaturecalculator.h"
#include "organize.hpp"
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <map>
#include <ostream>
#include <stdexcept>
#include <unordered_map>
#include <optional>
#include <utility>
#include <vector>
#include "curvaturebased.h"


namespace core {
  
  namespace {



    unsigned fact(unsigned v)
    {
      unsigned res = 1;
      while(v>0){
	res *=v;
	v--;
      }
      return res;
    }
    
    double bernstein(const unsigned k, const unsigned n, const double t)
    {
      return (fact(n) / (fact(k) * fact(n-k))) * std::pow(t, k) * std::pow(1-t, n-k);
    }

    double bernstein_interpol(const double t)
    {
      const unsigned k = 3;
      double res = 0;
      for(size_t i = 0; i <=k; i++) {
	res += bernstein(i, 2.0 * k + 1, t);
      }
      return 1.0 - res;
    }
   

    common::MyMesh::EdgeHandle findEdgeConnectingVertices(common::MyMesh& mesh, common::MyMesh::VertexHandle v1, common::MyMesh::VertexHandle v2) {
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
  } //namespace

  
    // Vertexhandle not present as child in this map must be original vertex
    //TODO refactor
  void DiscreteFairer::subdivide(common::MyMesh& mesh, ChildrenParents& child_parent_map) const
    {
  
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

	  child_parent_map.insert({new_vertex, {v1, v2}});
	

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
  
  void DiscreteFairer::getEffectorsHelper(const common::MyMesh::VertexHandle& to,
			    std::set<common::MyMesh::VertexHandle>& visited_vertices,
			    std::set<common::MyMesh::VertexHandle>& effectors,
			    const ChildrenParents& child_parents_map) const
    {
      if(child_parents_map.count(to) == 0) {
	// Original vertex
	effectors.insert(to);
	return;
      }

      visited_vertices.insert(to);
    
      const auto& parents = child_parents_map.at(to);

      for(size_t i = 0; i < 2; i++){
	if(std::find(visited_vertices.begin(), visited_vertices.end(), parents[i]) == visited_vertices.end()){
	  getEffectorsHelper(parents[i], visited_vertices, effectors, child_parents_map);
	}
      }    
    }

  std::set<common::MyMesh::VertexHandle> DiscreteFairer::getEffectors(const common::MyMesh::VertexHandle& to, const ChildrenParents& child_parents_map) const
    {
      std::set<common::MyMesh::VertexHandle> result;
      std::set<common::MyMesh::VertexHandle> tmp_visited_vertices;

      getEffectorsHelper(to, tmp_visited_vertices, result, child_parents_map);


      return result;
    }

  std::vector<std::pair<common::MyMesh::VertexHandle, double>> DiscreteFairer::getWeighedEffectors(const common::MyMesh::VertexHandle& to, const ChildrenParents& child_parents_map, const common::MyMesh& mesh) const
    {
      std::vector<std::pair<common::MyMesh::VertexHandle, double>> retval;
      const auto effectors = getEffectors(to, child_parents_map);

      if(effectors.size() == 2) {
	std::vector<common::MyMesh::VertexHandle> t_effectors;
	for(const auto& a : effectors) {
	  t_effectors.push_back(a); //TODO delete this
	}
	const auto p0 = mesh.point(to);
	const auto p1 = mesh.point(t_effectors[0]);
	const auto p2 = mesh.point(t_effectors[1]);
	auto d1 = (p0 - p1).norm();
	auto d2 = (p0 - p2).norm();
	const auto normalizer = d2 + d1;
	d1 = d1 / normalizer;
	d2 = d2 / normalizer;
	retval.push_back({t_effectors[0], d2});//TODO refactor
	retval.push_back({t_effectors[1], d1});
      }
      else if (effectors.size() == 3) {
	std::vector<common::MyMesh::VertexHandle> t_effectors;
	for(const auto& a : effectors) {
	  t_effectors.push_back(a);
	}
	const auto p0 = mesh.point(to);
	const auto p1 = mesh.point(t_effectors[0]);
	const auto p2 = mesh.point(t_effectors[1]);
	const auto p3 = mesh.point(t_effectors[2]);

	const auto v0 = p2 - p1;
	const auto v1 = p3 - p1;
	const auto v2 = p0 - p1;
	const float d00 = v0.dot(v0);
	const float d01 = v0.dot(v1);
	const float d11 = v1.dot(v1);
	const float d20 = v2.dot(v0);
	const float d21 = v2.dot(v1);
	const float denom = d00 * d11 - d01 * d01;
	const float v = (d11 * d20 - d01 * d21) / denom;
	const float w = (d00 * d21 - d01 * d20) / denom;
	const float u = 1.0f - v - w;
	retval.push_back({t_effectors[0], u});
	retval.push_back({t_effectors[1], v});
	retval.push_back({t_effectors[2], w});

      }

      return retval;
    }


  


  
  Eigen::Vector3d DiscreteFairer::iterateVertex(common::MyMesh& mesh, common::MyMesh::VertexHandle& iteratable, const DiscreteFairer::ExtendedVertexStaticInfo& extended_vertex_static_info) const
    {
      
      std::vector<common::MyMesh::VertexHandle> neighbors;
      for (common::MyMesh::ConstVertexVertexIter vv_it = mesh.cvv_begin(iteratable); vv_it != mesh.cvv_end(iteratable); ++vv_it)
	{
	  neighbors.push_back(*vv_it);
	}


    
      std::array<Eigen::Vector3d, 6> e_neighbors;
      for(size_t i = 0; i < 6; i++){
	const auto& p = mesh.point(neighbors[i]);
	e_neighbors[i] =  Eigen::Vector3d(p[0], p[1], p[2]);
      }

      CurvatureCalculator cc(mesh, true);
      cc.execute(iteratable);
  
      const auto& normal = cc.getNormal();
      const auto& fe = cc.getFundamentalElements();

      const auto H = calcTargetCurvature(extended_vertex_static_info.weighed_effectors);

      const auto Qm = mesh.point(iteratable);
      const Eigen::Vector3d Q(Qm[0], Qm[1], Qm[2]);
      
      //return DiscreteFairer::Q(e_neighbors, normal, H, fe, Q);
      return DiscreteFairer::Q2(e_neighbors, normal, H, fe, Q, cc.getLastM());
  
    }

  std::map<common::MyMesh::VertexHandle, DiscreteFairer::ExtendedVertexStaticInfo> DiscreteFairer::generateExtendedVertexSaticInfos(common::MyMesh& mesh, const ChildrenParents& child_parents_map) const
    {
      std::map<common::MyMesh::VertexHandle, DiscreteFairer::ExtendedVertexStaticInfo> retval;

      CurvatureCalculator cc(mesh);
      
      for (auto vh : mesh.vertices()) {
	
	DiscreteFairer::ExtendedVertexStaticInfo evsi;
	const auto weighed_effectors = getWeighedEffectors(vh, child_parents_map, mesh);
	evsi.is_original_vertex = weighed_effectors.size() < 2;

	if (!evsi.is_original_vertex) {
	  evsi.weighed_effectors = weighed_effectors;
	}
	retval.insert({vh, evsi});
      }
      
      return retval;
    }

        
  void DiscreteFairer::triangleExecuteDemo(common::MyMesh& mesh)
    {
      ChildrenParents cp;
      subdivide(mesh, cp);
      subdivide(mesh, cp);
      subdivide(mesh, cp);
      subdivide(mesh, cp);
      
      OpenMesh::VPropHandleT<double> demo_color;
      mesh.add_property(demo_color, "demo_color");

      const auto static_info = generateExtendedVertexSaticInfos(mesh, cp);

      for(const auto& v : mesh.vertices()) {
	switch(v.idx()) {
	case 0:
	  {
	    vertex_curvature_map[v] = 0.5;
	    mesh.property(demo_color, v) = 0.5;
	    break;
	  }
	case 1:
	  {
	    vertex_curvature_map[v] = 1.0;
	    mesh.property(demo_color, v) = 1.0;
	    break;
	  }
	case 2:
	  {
	    vertex_curvature_map[v] = 2.0;
	    mesh.property(demo_color, v) = 2.0;
	    break;
	  }
	default:
	  break;
	}
      }
      
	
      for (auto vh : mesh.vertices()) {
	if (!static_info.at(vh).is_original_vertex) {
	  //std::cout<<calcTargetCurvature(static_info.at(vh).weighed_effectors)<<std::endl;
	  mesh.property(demo_color, vh) = calcTargetCurvature(static_info.at(vh).weighed_effectors);
	}
      }

      
    }
    


  
  double DiscreteFairer::calcTargetCurvature(const std::vector<std::pair<common::MyMesh::VertexHandle, double>>& weighed_effectors) const
  {
#define CT 1
      
    if (CT ==1){
      double H = 0.0;
      for(const auto& weighed_effector : weighed_effectors) {
	H += vertex_curvature_map.at(weighed_effector.first) * weighed_effector.second;
      }
      return H;
    }
    //////////////
    if(CT ==2){
      double H = 0.0;
      double bern_sum = 0.0;
      for(const auto& weighed_effector : weighed_effectors) {
	bern_sum += bernstein_interpol(weighed_effector.second);
      }
      //bern_sum = 3.0/2.0; TODO
      //std::cout<<bern_sum<<std::endl;
      for(const auto& weighed_effector : weighed_effectors) {
	H += vertex_curvature_map.at(weighed_effector.first) * (bernstein_interpol(weighed_effector.second)/bern_sum);
      }
      return H;
    }

    if(CT ==3){
      double H = 0.0;
      const auto alpha = -0.5;
      for(const auto& weighed_effector : weighed_effectors) {
	H += std::pow(vertex_curvature_map.at(weighed_effector.first),-alpha) * weighed_effector.second;
      }
	
      //std::cout<<std::pow(H, -1.0 / alpha)<<std::endl;
	
      return std::pow(H, -1.0 / alpha);
      return H;
    }
      
  }

  void DiscreteFairer::execute(common::MyMesh& mesh, size_t face_split_count, size_t iteration_count)
  {
    if (mesh.n_vertices() == 3) {
      return triangleExecuteDemo(mesh);
    }

    
    ChildrenParents child_parents_map;
  
    // Subdivide the base mesh
    for(size_t i = 0; i < face_split_count; i++) {
      subdivide(mesh, child_parents_map);
    }

    extended_vertex_static_infos = generateExtendedVertexSaticInfos(mesh, child_parents_map);
    //std::cout<<"DiscreteFairer: ExtendedVertexStaticInfos are generated."<<std::endl;

    metrics::CurvatureBased cb(mesh, *this);
    cb.startSession();

    for(size_t i = 0; i < iteration_count ; i++){
      // Calculate the curvature for each vertex at the beginning of each iteration
      CurvatureCalculator mcc(mesh, true);
      
      //TODO range operator (smarthandle....)
      for(common::MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
	auto vh = *v_it;
	if(extended_vertex_static_infos.at(vh).is_original_vertex){
	  mcc.execute(vh);
	  vertex_curvature_map[vh] =  mcc.getMeanCurvature();
	  //std::cout<<"Curvature: "<< curvature<<std::endl;
	}
      }

      cb.postIteration();

      
      // Calculate the new position of each new vertex
      std::vector<std::pair<common::MyMesh::VertexHandle, Eigen::Vector3d>> new_vertex_positions;
      for(common::MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
	auto vh = *v_it;
	if (!extended_vertex_static_infos.at(vh).is_original_vertex) {
	  if(/*vh.idx()>=29 && vh.idx()<=3100*/ vh.idx()==0 || true){
	    new_vertex_positions.emplace_back(vh, iterateVertex(mesh, vh, extended_vertex_static_infos.at(vh)));
	  }
	}
      }
      //continue;
      size_t tmp_id = 0;
      // Replace each vertex to its new position
      for (auto& vertex_with_new_pos : new_vertex_positions) {
	// TODO: conversion in common (new file for all of these)
	const auto& new_pos_e = vertex_with_new_pos.second;
	common::MyMesh::Point new_pos(new_pos_e[0], new_pos_e[1], new_pos_e[2]);
	
	  mesh.point(vertex_with_new_pos.first) = new_pos;
      }

      std::cout<<"-----Iter--------"<<std::endl;
      for(auto vh : mesh.vertices()){
	if(extended_vertex_static_infos.at(vh).is_original_vertex){
	  CurvatureCalculator cc(mesh);
	  cc.execute(vh);
	  //std::cout<<cc.getMeanCurvature()<<std::endl;
	}
      }

    }
    cb.endSession();

    OpenMesh::IO::write_mesh(mesh, "result.obj");


    



  }

  Eigen::Vector3d DiscreteFairer::Q(const std::array<Eigen::Vector3d, 6>& p,const Eigen::Vector3d& normal, double H,  const CurvatureCalculator::FundamentalElements& fe, const Eigen::Vector3d& Q)
  {

    const auto p_k = common::average({p[0], p[1], p[2], p[3], p[4], p[5]});

    const auto t = (Q-p_k).dot(normal) + (fe.E * fe.N - 2.0 * fe.M * fe.F + fe.G * fe.L - 2.0 * H * (fe.E * fe.G - fe.F * fe.F)) / (2.0 * (fe.E + fe.G));
    const auto t2 = ((-2.0 * H * (fe.E * fe.G - fe.F * fe.F)) +
      (1/3.0)* fe.E * (-p[0] + 2 * p[1] + 2 * p[2] - p[3] + 2 * p[4] + 2 * p[5]).dot(normal) -
      2 / std::sqrt(3) * (p[1] - p[2] + p[4] - p[5]).dot(normal) +
		    fe.G * (p[0]+p[3]).dot(normal)) / (2 * (fe.E + fe.G))
      -p_k.dot(normal);
    
    return p_k + normal * t;
  }


  Eigen::Vector3d DiscreteFairer::Q2(const std::array<Eigen::Vector3d, 6>& p,const Eigen::Vector3d& normal, double H,  const CurvatureCalculator::FundamentalElements& fe, const Eigen::Vector3d& Q, const Eigen::Matrix<double, 5 , 6>& M)
  {

     const auto p_k = common::average({p[0], p[1], p[2], p[3], p[4], p[5]});


    Eigen::RowVectorXd row3 = M.row(2);
    double dotsum_m3 = 0.0;
    for(size_t i = 0; i< 6; i++) {
      dotsum_m3 += (row3(i) * p[i]).dot(normal);
    }

    Eigen::RowVectorXd row4 = M.row(3);
    double dotsum_m4 = 0.0;
    for(size_t i = 0; i< 6; i++) {
      dotsum_m4 += (row4(i) * p[i]).dot(normal);
    }

    Eigen::RowVectorXd row5 = M.row(4);
    double dotsum_m5 = 0.0;
    for(size_t i = 0; i< 6; i++) {
      dotsum_m5 += (row5(i) * p[i]).dot(normal);
    }

    const double t = (-2.0 * H*(fe.E * fe.G - fe.F * fe.F) + fe.E * dotsum_m5 - 2.0 * fe.F * dotsum_m4 + fe.G * dotsum_m3) / ((fe.E * M.row(4).sum())-2.0 * fe.F * M.row(3).sum() + fe.G * M.row(2).sum()) - p_k.dot(normal);
    return p_k + normal * t;
    

  }




}
