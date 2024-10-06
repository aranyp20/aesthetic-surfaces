#include "discretefairer.h"
#include <OpenMesh/Tools/Subdivider/Uniform/MidpointT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "common_defines.h"
#include "curvaturecalculator.h"
#include "mainwindow.h"
#include "mesh.h"
#include "organize.hpp"
#include <QtCore/qmetatype.h>
#include <QtCore/qnamespace.h>
#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <stdexcept>
#include <unordered_map>
#include <optional>
#include <utility>
#include <vector>
#include "curvaturebased.h"
#include "subdivider.h"
#include "settings.h"
#include "type_converter.h"
#include "mymath.h"


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


    std::array<double, 4> barycentricCoordinates(const Eigen::Vector3d& P, 
						 const Eigen::Vector3d& P1, 
						 const Eigen::Vector3d& P2, 
						 const Eigen::Vector3d& P3, 
						 const Eigen::Vector3d& P4) {

      Eigen::Vector4d v1(1, P1[0], P1[1], P1[2]);
      Eigen::Vector4d v2(1, P2[0], P2[1], P2[2]);
      Eigen::Vector4d v3(1, P3[0], P3[1], P3[2]);
      Eigen::Vector4d v4(1, P4[0], P4[1], P4[2]);
      
      Eigen::Matrix4d matrix;
      matrix.col(0) = v1;
      matrix.col(1) = v2;
      matrix.col(2) = v3;
      matrix.col(3) = v4;


      Eigen::Vector4d v(1, P[0], P[1], P[2]);
      
      Eigen::Vector4d lambdas;
      try {
        lambdas = matrix.colPivHouseholderQr().solve(v);
      } catch (const std::exception& e) {
        throw std::runtime_error("Barycentric matrix inversion failed...");
      }


      return {lambdas[0], lambdas[1], lambdas[2], lambdas[3]};
    }

    std::array<double, 4> barycentricCoordinatesImproved(const Eigen::Vector3d& P, 
							 const Eigen::Vector3d& P1, 
							 const Eigen::Vector3d& P2, 
							 const Eigen::Vector3d& P3, 
							 const Eigen::Vector3d& P4) {
      
      const auto bar_of_vertex = [P](const Eigen::Vector3d& pi, const Eigen::Vector3d& pi0, const Eigen::Vector3d& pi2) {

	return common::math::heron(pi, pi0, pi2) / (common::math::heron(pi, pi0, P) * common::math::heron(pi, pi2, P));
	
      };

      auto l1 = bar_of_vertex(P1, P2, P4);
      auto l2 = bar_of_vertex(P2, P1, P3);
      auto l3 = bar_of_vertex(P3, P2, P4);
      auto l4 = bar_of_vertex(P4, P3, P1);
      const auto normalizer = l1 + l2 + l3 + l4;
      l1 /= normalizer;
      l2 /= normalizer;
      l3 /= normalizer;
      l4 /= normalizer;
      
      return {l1, l2, l3, l4};
    }

    std::array<double, 3> barycentricCoordinatesImproved(const Eigen::Vector3d& P, 
							 const Eigen::Vector3d& P1, 
							 const Eigen::Vector3d& P2, 
							 const Eigen::Vector3d& P3) {
      
      const auto bar_of_vertex = [P](const Eigen::Vector3d& pi, const Eigen::Vector3d& pi0, const Eigen::Vector3d& pi2) {

	return common::math::heron(pi, pi0, pi2) / (common::math::heron(pi, pi0, P) * common::math::heron(pi, pi2, P));
	
      };

      auto l1 = bar_of_vertex(P1, P2, P3);
      auto l2 = bar_of_vertex(P2, P1, P3);
      auto l3 = bar_of_vertex(P3, P1, P2);
      const auto normalizer = l1 + l2 + l3;
      l1 /= normalizer;
      l2 /= normalizer;
      l3 /= normalizer;
      
      return {l1, l2, l3};
    }

 
  } //namespace




  
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

    std::vector<common::MyMesh::VertexHandle> t_effectors;
    for(const auto& a : effectors) {
      t_effectors.push_back(a); //TODO delete this
    }

    if(effectors.size() == 1) {
      return retval;
    }

    sortEffectors(mesh.original_state, t_effectors); // todo this shouldnt be performed for every vertex separately

	
    if(effectors.size() == 2) {
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
      //retval.push_back({t_effectors[0], u});
      //retval.push_back({t_effectors[1], v});
      //retval.push_back({t_effectors[2], w});

      const auto barys = barycentricCoordinatesImproved(common::converter::meshPointToEigen(mesh.point(to)),
							common::converter::meshPointToEigen(mesh.point(t_effectors[0])),
							common::converter::meshPointToEigen(mesh.point(t_effectors[1])),
							common::converter::meshPointToEigen(mesh.point(t_effectors[2])));

      retval.push_back({t_effectors[0], barys[0]});
      retval.push_back({t_effectors[1], barys[1]});
      retval.push_back({t_effectors[2], barys[2]});

    }
    else if (effectors.size() == 4) { // todo sima else es a generalized mehet n-re

      const auto barys = barycentricCoordinatesImproved(common::converter::meshPointToEigen(mesh.point(to)),
							common::converter::meshPointToEigen(mesh.point(t_effectors[0])),
							common::converter::meshPointToEigen(mesh.point(t_effectors[1])),
							common::converter::meshPointToEigen(mesh.point(t_effectors[2])),
							common::converter::meshPointToEigen(mesh.point(t_effectors[3])));

      retval.push_back({t_effectors[0], barys[0]});
      retval.push_back({t_effectors[1], barys[1]});
      retval.push_back({t_effectors[2], barys[2]});
      retval.push_back({t_effectors[3], barys[3]});
	
    }
    else {
      throw std::runtime_error("cant get weighed effectors...: " + std::to_string(effectors.size()));
    }

    return retval;
  }


  


  
  Eigen::Vector3d DiscreteFairer::iterateVertex(common::MyMesh& mesh, common::MyMesh::VertexHandle& iteratable, const DiscreteFairer::ExtendedVertexStaticInfo& extended_vertex_static_info) const
  {
      
    std::vector<common::MyMesh::VertexHandle> neighbors;
    for (common::MyMesh::ConstVertexVertexIter vv_it = mesh.cvv_begin(iteratable); vv_it != mesh.cvv_end(iteratable); ++vv_it) {
      neighbors.push_back(*vv_it);
    }

    std::vector<Eigen::Vector3d> e_neighbors;
      for(const auto& neighbor : neighbors) {
	const auto& p = mesh.point(neighbor);
	e_neighbors.emplace_back(p[0], p[1], p[2]);
      }

      CurvatureCalculator cc(mesh, true);
      cc.execute(iteratable);
  
      const auto& normal = cc.getNormal();
      const auto& fe = cc.getFundamentalElements();

      std::vector<EffectorExtra> effectors_extra; //todo ezt a konverziot elobb megcsinalni
      for(const auto& effector : extended_vertex_static_info.weighed_effectors) {
	effectors_extra.push_back({effector.second, vertex_curvature_map.at(effector.first), common::converter::meshPointToEigen(mesh.point(effector.first)), common::converter::meshPointToEigen(mesh.point(iteratable)), effector.first.idx()});
      }

      const auto H = calcTargetCurvature(effectors_extra);
      
      const auto Qm = mesh.point(iteratable);
      const Eigen::Vector3d Q(Qm[0], Qm[1], Qm[2]);

      

      Eigen::Vector3d Q0;
      
      if(extended_vertex_static_info.weighed_effectors.size() == 2){
	const auto edge_neighbors = mesh.edge_neighbor_map.at(iteratable);
	Q0 = common::converter::meshPointToEigen((mesh.point(edge_neighbors[0]) + mesh.point(edge_neighbors[1]))/ 2);
      } else {
	Q0 = common::average(e_neighbors);
      }

      
      Eigen::Vector3d new_pos = Q;
      if(common::settings::selected_curvature == common::settings::GAUSSIAN) {
	new_pos = DiscreteFairer::Q_Gaussian(e_neighbors, normal, H, fe, Q, Q0, cc.getLastM());
      }
      else if(common::settings::selected_curvature == common::settings::MEAN) {
	if(iteratable.idx()==50 || true){
	  new_pos = DiscreteFairer::Q2(e_neighbors, normal, H.main, fe, Q,Q0, cc.getLastM());
	}
      }
      return new_pos;
    }

  std::map<common::MyMesh::VertexHandle, DiscreteFairer::ExtendedVertexStaticInfo> DiscreteFairer::generateExtendedVertexSaticInfos(common::MyMesh& mesh, const ChildrenParents& child_parents_map) const
  {
    std::map<common::MyMesh::VertexHandle, DiscreteFairer::ExtendedVertexStaticInfo> retval;

    CurvatureCalculator cc(mesh, true);
    
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

  void DiscreteFairer::sortEffectors(const common::MyMesh& mesh, std::vector<common::MyMesh::VertexHandle>& effectors) const
  {
    const auto are_neighbors = [mesh](const common::MyMesh::VertexHandle& vh1, const common::MyMesh::VertexHandle& vh2)
    {
      return mesh.find_halfedge(vh1, vh2).is_valid() || mesh.find_halfedge(vh2, vh1).is_valid();
    };

    std::vector<common::MyMesh::VertexHandle> sorted_effectors(effectors.size());
    std::vector<size_t> used_indicies(effectors.size());
    sorted_effectors[0] = effectors[0];
    used_indicies[0] = 0;

    for(size_t i = 1; i < effectors.size(); i++) {
      bool found = false;
      for(size_t j = 1; j < effectors.size(); j++) {

	if(std::find(used_indicies.begin(), used_indicies.end(), j) != used_indicies.end()) {
	  continue;
	}

	if (are_neighbors(sorted_effectors[i-1], effectors[j])) {
	  sorted_effectors[i] = effectors[j];
	  used_indicies[i] = j;
	  found = true;
	  break;
	}
      }
      if (!found) {
	throw std::runtime_error("sortEffectors: cant sort effectors.");
      }
    }

    effectors = sorted_effectors;
  }

        
  void DiscreteFairer::triangleExecuteDemo(common::MyMesh& mesh)
  {
      
      Subdivider subdivider;
      subdivider.execute(mesh, 5);
      
      OpenMesh::VPropHandleT<double> demo_color;
      mesh.add_property(demo_color, "demo_color");

      const auto static_info = generateExtendedVertexSaticInfos(mesh, mesh.children_parents_map);

      for(const auto& v : mesh.vertices()) {
	switch(v.idx()) {
	case 0:
	  {
	    vertex_curvature_map[v] = 4;
	    mesh.property(demo_color, v) = .5;
	    break;
	  }
	case 4:
	  {
	    vertex_curvature_map[v] = -5;
	    mesh.property(demo_color, v) = -.2;
	    break;
	  }
	case 5:
	  {
	    vertex_curvature_map[v] = 7;
	    mesh.property(demo_color, v) = .3;
	    break;
	  }
	default:
	  break;
	}
      }
      

      for(size_t i = 0; i < 10 ; i++){
	for (auto vh : mesh.vertices()) {
	  bool in_triangle = true;
	  for(const auto& effector : static_info.at(vh).weighed_effectors) {
	    if(effector.first.idx() != 0 && effector.first.idx() != 4 && effector.first.idx() != 5) {
	      in_triangle = false;
	    }
	  }
	  
	  if (in_triangle && !static_info.at(vh).is_original_vertex) {
	    const auto new_pos = iterateVertex(mesh, vh, static_info.at(vh));
	    common::MyMesh::Point new_posa(new_pos[0], new_pos[1], new_pos[2]);
	    mesh.point(vh) = new_posa;
	  }
	}
      }

      
    }

    void DiscreteFairer::pentaExecuteDemo(common::MyMesh& mesh)
  {
      
      Subdivider subdivider;
      subdivider.execute(mesh, 4);
      
      OpenMesh::VPropHandleT<double> demo_color;
      mesh.add_property(demo_color, "demo_color");

      const auto static_info = generateExtendedVertexSaticInfos(mesh, mesh.children_parents_map);

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
	case 3:
	  {
	    vertex_curvature_map[v] = 4.0;
	    mesh.property(demo_color, v) = 4.0;
	    break;
	  }
	case 4:
	  {
	    vertex_curvature_map[v] = 5.0;
	    mesh.property(demo_color, v) = 5.0;
	    break;
	  }
	default:
	  break;
	}
      }
      
	
      for (auto vh : mesh.vertices()) {

	
	if (!static_info.at(vh).is_original_vertex) {
	  throw std::runtime_error("not implemented");
	  //mesh.property(demo_color, vh) = calcTargetCurvature(static_info.at(vh).weighed_effectors);
	}
      }

      
    }


  double calcSignedGaussianCurvature(const CurvatureCalculator& cc)
  {

    
    const auto principle_curvatures = cc.getPrincipleCurvatures();
    //const auto chosen_curvature (std::fabs(principle_curvatures.min_val) < std::fabs(principle_curvatures.max_val) ? principle_curvatures.min_val : principle_curvatures.max_val);
    //return chosen_curvature < 0 ? -std::fabs(cc.getGaussianCurvature()) : std::fabs(cc.getGaussianCurvature());
    if(cc.getGaussianCurvature() > 0) {
      return principle_curvatures.min_val < 0 ? -cc.getGaussianCurvature() : cc.getGaussianCurvature();
    }

    
    return 0;
  }

  void DiscreteFairer::triangleGaussExecuteDemo(common::MyMesh& mesh)
  {

    Subdivider subdivider;
    subdivider.execute(mesh, 3);

    extended_vertex_static_infos = generateExtendedVertexSaticInfos(mesh, mesh.children_parents_map); //todo ezt mindig ujra kell generalni?


    //CurvatureCalculator cc(mesh, true);
    std::cout << "-----------" << std::endl;
    //for(common::MyMesh::VertexHandle vert : mesh.vertices()) {
    //  if(extended_vertex_static_infos.at(vert).is_original_vertex) {
    //cc.execute(vert);
	//vertex_curvature_map[vert] =  getCurvature(cc);
    //}
    //}
    
    auto vh1 = mesh.vertex_handle(0);
    vertex_curvature_map.insert({vh1, 3});
    auto vh2 = mesh.vertex_handle(2);
    vertex_curvature_map.insert({vh2, 5});
    

    std::vector<int> av{72, 18, 70, 7, 87, 21, 89};
    //std::vector<int> av{8};

    for(size_t i = 0; i < 20; i++) {
      
    
      std::vector<std::pair<common::MyMesh::VertexHandle, Eigen::Vector3d>> new_vertex_positions;
      for(common::MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
	auto vh = *v_it;
	if(std::find(av.begin(), av.end(),vh.idx()) == av.end()) {
	  continue;
	}
	if (!extended_vertex_static_infos.at(vh).is_original_vertex) {
	  new_vertex_positions.emplace_back(vh, iterateVertex(mesh, vh, extended_vertex_static_infos.at(vh)));
	}
      }
      
      // Replace each vertex to its new position
      for (auto& vertex_with_new_pos : new_vertex_positions) {
	// TODO: conversion in common (new file for all of these)
	const auto& new_pos_e = vertex_with_new_pos.second;
	common::MyMesh::Point new_pos(new_pos_e[0], new_pos_e[1], new_pos_e[2]);
	
	mesh.point(vertex_with_new_pos.first) = new_pos;
      }

    }
    /*
    OpenMesh::VPropHandleT<double> demo_color;
    mesh.add_property(demo_color, "demo_color");
    */
    for(auto vh : mesh.vertices()){
      if(std::find(av.begin(), av.end(),vh.idx()) == av.end()) {
	continue;
      }
      CurvatureCalculator cc2(mesh, true);
      cc2.execute(vh);
      //std::cout <<"res: "<< getCurvature(cc2) << std::endl;
      //mesh.property(demo_color, vh) = cc2.getGaussianCurvature();
    }
    
  }


  DiscreteFairer::TargetCurvature DiscreteFairer::logAestheticTargetCurvatureCore(const std::vector<EffectorExtra>& effectors, const double alpha) const
  {
    TargetCurvature retval;
    const double sign = effectors[0].H < 0 ? -1 : 1;
    for(const auto& effector : effectors) {
      retval.main += std::powf(effector.H * sign, -alpha) * effector.weight;
    }
      
    retval.main = std::powf(retval.main, -1.0 / alpha) * sign;

    return retval;
  }

  struct SlicePoints
  {
    Eigen::Vector3d start;

    Eigen::Vector3d neg_side_sp;
    double neg_side_H;
    
    Eigen::Vector3d pos_side_sp;
    double pos_side_H;
    
    Eigen::Vector3d end;
  };
  
  DiscreteFairer::TargetCurvature DiscreteFairer::logAestheticTargetCurvatureMixed(const std::vector<EffectorExtra>& effectors) const
  {
    const auto zero_point = [](const EffectorExtra& e_neg, const EffectorExtra& e_pos)
    {
      return (e_pos.H * e_neg.pos + std::fabs(e_neg.H) * e_pos.pos) / (e_pos.H + std::fabs(e_neg.H));
    };

    const auto trislice_points = [](const EffectorExtra& e_neg, const EffectorExtra& e_pos)
    {
      constexpr double middle_ratio = 0.2;
      
      SlicePoints retval;
      retval.start = e_neg.pos;
      retval.end = e_pos.pos;

      const auto total_value_interval = std::fabs(e_neg.H) + e_pos.H;
      const auto middle_value_interval = total_value_interval * middle_ratio;
      
      const auto zero_t = std::fabs(e_neg.H) / total_value_interval;

      const Eigen::Vector3d dir_vec = e_pos.pos - e_neg.pos;

      const Eigen::Vector3d zero_p = e_neg.pos + zero_t * dir_vec;

      const auto neg_side_t = zero_t - middle_ratio / 2.0;
      const auto pos_side_t = zero_t + middle_ratio / 2.0;

      const auto H_interval = e_neg.H + e_pos.H;
      retval.neg_side_H = e_neg.H + neg_side_t * H_interval;
      retval.pos_side_H = e_neg.H + pos_side_t * H_interval;
      
      if (neg_side_t > 0) { 
	const Eigen::Vector3d neg_side_point = e_neg.pos + neg_side_t * dir_vec;
	retval.neg_side_sp = neg_side_point;
      }
      else {
	retval.neg_side_sp = common::average({zero_p, e_neg.pos});
      }

      if (pos_side_t < 1) {
	const Eigen::Vector3d pos_side_point = e_neg.pos + pos_side_t * dir_vec;
	retval.pos_side_sp = pos_side_point;
      }
      else {
	retval.neg_side_sp = common::average({zero_p, e_pos.pos});
      }

      
      
      return retval;
    };


    std::vector<EffectorExtra> adjusted_effectors;
    double alpha = common::settings::log_aesthetic_alpha;
    
    if(effectors.size() == 2) {
      if (common::settings::log_aesthetic_alpha < 0) {
	const auto zp = zero_point(effectors[0], effectors[1]);

	const auto dist_to_e0 = (effectors[0].pos - effectors[0].subject_pos).norm();
	const auto dist_to_e1 = (effectors[1].pos - effectors[0].subject_pos).norm();
	const auto dist_to_zp = (zp - effectors[0].subject_pos).norm();
	if((effectors[0].pos - zp).norm() > dist_to_e0) {
	  /* Subject is closer to effectors[0] than the zero point -> part of the effectors[0]--zp curve. */
	  adjusted_effectors.push_back({dist_to_zp / (dist_to_e0 + dist_to_zp),effectors[0].H,{}, {} });
	}
	else {
	  /* Subject is closer to effectors[1] than the zero point -> part of the effectors[1]--zp curve. */
	  adjusted_effectors.push_back({dist_to_zp / (dist_to_e1 + dist_to_zp),effectors[1].H});
	}
      }
      else {
	const auto sps = trislice_points(effectors[0], effectors[1]);
	const auto& rp = effectors[0].subject_pos;
	const double dist_to_start = (rp - effectors[0].pos).norm();
	const double neg_side_to_start = (sps.neg_side_sp - sps.start).norm();
	const double pos_side_to_start = (sps.pos_side_sp - sps.start).norm();
	
	if(dist_to_start < neg_side_to_start) {
	  adjusted_effectors.push_back({1.0 - (dist_to_start / neg_side_to_start),effectors[0].H});
	  adjusted_effectors.push_back({dist_to_start / neg_side_to_start, sps.neg_side_H});
	}
	else if (dist_to_start < pos_side_to_start) {
	  const double middle_end = pos_side_to_start - neg_side_to_start;
	  const double rp_corrigated = dist_to_start - neg_side_to_start;
	  adjusted_effectors.push_back({rp_corrigated / middle_end, sps.pos_side_H});
	  adjusted_effectors.push_back({1.0 - rp_corrigated / middle_end, sps.neg_side_H});
	  alpha = -1;
	}
	else {
	  const double interval_end = 1.0 - pos_side_to_start;
	  const double rp_corrigated = dist_to_start - pos_side_to_start;
	  adjusted_effectors.push_back({1.0 - rp_corrigated / interval_end, sps.pos_side_H});
	  adjusted_effectors.push_back({rp_corrigated / interval_end, effectors[1].H});
	}
      }
      
    }
    else {
      const Eigen::Vector3d normal = (effectors[0].pos - effectors[1].pos).cross(effectors[2].pos - effectors[1].pos);
      if (common::settings::log_aesthetic_alpha < 0) { //todo
      
	if (effectors[1].H > 0) {
	  /* N P P */
	  const auto zp1 = zero_point(effectors[0], effectors[1]);
	  const auto zp2 = zero_point(effectors[0], effectors[2]);

	  const Eigen::Vector3d towards_one_side = normal.cross(zp1 - zp2);

	  const bool neg_side = towards_one_side.dot(effectors[0].pos - zp2) > 0;
	  const bool examined_side = towards_one_side.dot(effectors[0].subject_pos - zp2) > 0;

	  const bool on_neg_side = neg_side == examined_side;

	
	  if(on_neg_side) {
	    const auto bary_coords = barycentricCoordinatesImproved(effectors[0].subject_pos, effectors[0].pos, zp1, zp2);
	    adjusted_effectors.push_back({bary_coords[0], effectors[0].H});
	  } else {
	    const auto bary_coords = barycentricCoordinatesImproved(effectors[0].subject_pos, effectors[1].pos, zp1, zp2, effectors[2].pos);
	    adjusted_effectors.push_back({bary_coords[0], effectors[1].H});
	    adjusted_effectors.push_back({bary_coords[3], effectors[2].H});
	  }
	}
	else {
	  /* N N P */

	  const auto zp1 = zero_point(effectors[0], effectors[2]);
	  const auto zp2 = zero_point(effectors[1], effectors[2]);

	  const Eigen::Vector3d towards_one_side = normal.cross(zp1 - zp2);

	  const bool pos_side = towards_one_side.dot(effectors[2].pos - zp2) > 0;
	  const bool examined_side = towards_one_side.dot(effectors[0].subject_pos - zp2) > 0;

	  const bool on_pos_side = pos_side == examined_side;
	
	  if(on_pos_side) {
	    const auto bary_coords = barycentricCoordinatesImproved(effectors[0].subject_pos, effectors[2].pos, zp1, zp2);
	    adjusted_effectors.push_back({bary_coords[0], effectors[2].H});
	  } else {
	    const auto bary_coords = barycentricCoordinatesImproved(effectors[0].subject_pos, effectors[0].pos, zp1, zp2, effectors[1].pos);
	    adjusted_effectors.push_back({bary_coords[0], effectors[0].H});
	    adjusted_effectors.push_back({bary_coords[3], effectors[1].H});
	  }
	
	}
      }
      else {
	constexpr double eps = 1e-10;


	if (effectors[1].H > 0) {
	  /* N P P */
	  const auto slice_points1 = trislice_points(effectors[0], effectors[1]);
	  const auto slice_points2 = trislice_points(effectors[0], effectors[2]);

	  const Eigen::Vector3d sep_line_1 = slice_points1.neg_side_sp - slice_points2.neg_side_sp;
	  if(sep_line_1.norm() > eps) {
	    const Eigen::Vector3d towards_one_side = normal.cross(sep_line_1);
	    const bool neg_side = towards_one_side.dot(effectors[0].pos - slice_points2.neg_side_sp) > 0;
	    const bool examined_side = towards_one_side.dot(effectors[0].subject_pos - slice_points2.neg_side_sp) > 0;
	    const bool on_neg_side = neg_side == examined_side;

	    if (on_neg_side) {
	      const auto bary_coords = barycentricCoordinatesImproved(effectors[0].subject_pos,
								      effectors[0].pos,
								      slice_points1.neg_side_sp,
								      slice_points2.neg_side_sp);
	      
	      adjusted_effectors.push_back({bary_coords[0], effectors[0].H});
	      adjusted_effectors.push_back({bary_coords[1], slice_points1.neg_side_H});
	      adjusted_effectors.push_back({bary_coords[2], slice_points2.neg_side_H});

	      return logAestheticTargetCurvatureCore(adjusted_effectors, alpha);
	    }
	  } 
	  /* Either there is no neg side or subject is not on it */
	  const Eigen::Vector3d sep_line_2 = slice_points1.pos_side_sp - slice_points2.pos_side_sp;
	  const Eigen::Vector3d towards_one_side = normal.cross(sep_line_2);
	  const bool middle_side = towards_one_side.dot(effectors[0].pos - slice_points2.pos_side_sp) > 0;
	  const bool examined_side = towards_one_side.dot(effectors[0].subject_pos - slice_points2.pos_side_sp) > 0;
	  const bool on_middle_side = middle_side == examined_side;

	  if (on_middle_side) {
	    const auto bary_coords = barycentricCoordinatesImproved(effectors[0].subject_pos,
								    slice_points1.pos_side_sp,
								    slice_points2.pos_side_sp,
								    slice_points2.neg_side_sp,
								    slice_points1.neg_side_sp);
	    alpha = -1.0;

	    adjusted_effectors.push_back({bary_coords[0], slice_points1.pos_side_H});
	    adjusted_effectors.push_back({bary_coords[1], slice_points2.pos_side_H});
	    adjusted_effectors.push_back({bary_coords[2], slice_points2.neg_side_H});
	    adjusted_effectors.push_back({bary_coords[3], slice_points1.neg_side_H});
	    
	    return logAestheticTargetCurvatureCore(adjusted_effectors, alpha);
	  }

	  const auto bary_coords = barycentricCoordinatesImproved(effectors[0].subject_pos,
								  slice_points1.end,
								  slice_points1.pos_side_sp,
								  slice_points2.pos_side_sp,
								  slice_points2.end);

	  adjusted_effectors.push_back({bary_coords[0], effectors[1].H});
	  adjusted_effectors.push_back({bary_coords[1], slice_points1.pos_side_H});
	  adjusted_effectors.push_back({bary_coords[2], slice_points2.pos_side_H});
	  adjusted_effectors.push_back({bary_coords[3], effectors[2].H});
	  
	  return logAestheticTargetCurvatureCore(adjusted_effectors, alpha);
	  
	}
	else {
	  /* N N P */
	  const auto slice_points1 = trislice_points(effectors[0], effectors[2]);
	  const auto slice_points2 = trislice_points(effectors[1], effectors[2]);


	  const Eigen::Vector3d sep_line_1 = slice_points1.pos_side_sp - slice_points2.pos_side_sp;
	  if(sep_line_1.norm() > eps) {
	    const Eigen::Vector3d towards_one_side = normal.cross(sep_line_1);
	    const bool pos_side = towards_one_side.dot(effectors[2].pos - slice_points2.pos_side_sp) > 0;
	    const bool examined_side = towards_one_side.dot(effectors[0].subject_pos - slice_points2.pos_side_sp) > 0;
	    const bool on_pos_side = pos_side == examined_side;

	    if (on_pos_side) {
	      const auto bary_coords = barycentricCoordinatesImproved(effectors[0].subject_pos,
								      effectors[2].pos,
								      slice_points1.pos_side_sp,
								      slice_points2.pos_side_sp);
	      adjusted_effectors.push_back({bary_coords[0], effectors[2].H});
	      adjusted_effectors.push_back({bary_coords[1], slice_points1.pos_side_H});
	      adjusted_effectors.push_back({bary_coords[2], slice_points2.pos_side_H});

	      return logAestheticTargetCurvatureCore(adjusted_effectors, alpha);
	    }
	  } 
	  /* Either there is no pos side or subject is not on it */
	  const Eigen::Vector3d sep_line_2 = slice_points1.neg_side_sp - slice_points2.neg_side_sp;
	  const Eigen::Vector3d towards_one_side = normal.cross(sep_line_2);
	  const bool middle_side = towards_one_side.dot(effectors[2].pos - slice_points2.neg_side_sp) > 0;
	  const bool examined_side = towards_one_side.dot(effectors[0].subject_pos - slice_points2.neg_side_sp) > 0;
	  const bool on_middle_side = middle_side == examined_side;

	  if (on_middle_side) {
	    const auto bary_coords = barycentricCoordinatesImproved(effectors[0].subject_pos,
								    slice_points1.pos_side_sp,
								    slice_points2.pos_side_sp,
								    slice_points2.neg_side_sp,
								    slice_points1.neg_side_sp);
	    alpha = -1.0;

	    adjusted_effectors.push_back({bary_coords[0], slice_points1.pos_side_H});
	    adjusted_effectors.push_back({bary_coords[1], slice_points2.pos_side_H});
	    adjusted_effectors.push_back({bary_coords[2], slice_points2.neg_side_H});
	    adjusted_effectors.push_back({bary_coords[3], slice_points1.neg_side_H});
	    
	    return logAestheticTargetCurvatureCore(adjusted_effectors, alpha);
	  }

	  const auto bary_coords = barycentricCoordinatesImproved(effectors[0].subject_pos,
								  slice_points1.start,
								  slice_points1.neg_side_sp,
								  slice_points2.neg_side_sp,
								  slice_points2.start);

	  adjusted_effectors.push_back({bary_coords[0], effectors[0].H});
	  adjusted_effectors.push_back({bary_coords[1], slice_points1.neg_side_H});
	  adjusted_effectors.push_back({bary_coords[2], slice_points2.neg_side_H});
	  adjusted_effectors.push_back({bary_coords[3], effectors[1].H});
	  
	  return logAestheticTargetCurvatureCore(adjusted_effectors, alpha);
	}
      }

    }
    return logAestheticTargetCurvatureCore(adjusted_effectors, alpha);
  }
  
  DiscreteFairer::TargetCurvature DiscreteFairer::calcLogAestheticTargetCurvature(const std::vector<EffectorExtra>& weighed_effectors) const
  {
  
    auto w_effectors = weighed_effectors;

    
    if (w_effectors.size() > 3) {
      throw std::runtime_error("Log-aesthetic df only works for trimeshes.");
    }
    
    std::sort(w_effectors.begin(), w_effectors.end(), [this](const EffectorExtra& v1, const EffectorExtra& v2){
      return v1.H < v2.H;
    });

    if (w_effectors.front().H > 0 || w_effectors.back().H < 0) {
      /* For same signed curvatures simply use the core function. */
      return logAestheticTargetCurvatureCore(w_effectors);
    }
    return logAestheticTargetCurvatureMixed(w_effectors);
  }

  DiscreteFairer::TargetCurvature DiscreteFairer::calcLogAestheticTargetCurvature_offsetVersion(const std::vector<EffectorExtra>& weighed_effectors) const
  {
    TargetCurvature retval;
    
    auto min_curvature = std::numeric_limits<double>::min();
    for(const auto& weighed_effector : weighed_effectors) {
      if(weighed_effector.H < min_curvature) {
	min_curvature = weighed_effector.H;
      }
    }
    auto offseter = 0.0;
    if(min_curvature < 1.0) {
      offseter = 1.0;
      if(min_curvature < 0.0) {
	offseter += std::fabs(min_curvature);
      }
    }
      
    const auto alpha = common::settings::log_aesthetic_alpha;
    for(const auto& weighed_effector : weighed_effectors) {
      retval.main += std::pow(weighed_effector.H + offseter,-alpha) * weighed_effector.weight;
    }
      
    retval.main = std::pow(retval.main, -1.0 / alpha) - offseter;

      
    return retval;
  }
  
  
  DiscreteFairer::TargetCurvature DiscreteFairer::calcTargetCurvature(const std::vector<EffectorExtra>& weighed_effectors) const
  {

    TargetCurvature retval;
    switch(common::settings::selected_alg) {      
    case common::settings::Algorithm::BASIC: {
	for(const auto& weighed_effector : weighed_effectors) {
	  retval.main += weighed_effector.H * weighed_effector.weight;
	}
	return retval;
      }
      //////////////
    case common::settings::Algorithm::BEZIER: {
	double bern_sum = 0.0;
	for(const auto& weighed_effector : weighed_effectors) {
	  bern_sum += bernstein_interpol(weighed_effector.weight);
	}
	//bern_sum = 3.0/2.0; TODO
	for(const auto& weighed_effector : weighed_effectors) {
	  retval.main += weighed_effector.H * (bernstein_interpol(weighed_effector.weight) / bern_sum);
	}
	return retval;
      }

    case common::settings::Algorithm::LOG_AESTHETIC: {
      return calcLogAestheticTargetCurvature(weighed_effectors);
    }
    }
      
  }

  void DiscreteFairer::iterateVerticesAsync(common::MyMesh& mesh)
  {
    // Calculate the new position of each new vertex
    std::vector<std::pair<common::MyMesh::VertexHandle, Eigen::Vector3d>> new_vertex_positions;
    for(common::MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
      auto vh = *v_it;
      if (!extended_vertex_static_infos.at(vh).is_original_vertex) {
	new_vertex_positions.emplace_back(vh, iterateVertex(mesh, vh, extended_vertex_static_infos.at(vh)));
      }
    }
      
    // Replace each vertex to its new position
    for (auto& vertex_with_new_pos : new_vertex_positions) {
      // TODO: conversion in common (new file for all of these)
      const auto& new_pos_e = vertex_with_new_pos.second;
      common::MyMesh::Point new_pos(new_pos_e[0], new_pos_e[1], new_pos_e[2]);
      mesh.point(vertex_with_new_pos.first) = new_pos;
    }
  }

  void DiscreteFairer::iterateVerticesSync(common::MyMesh& mesh)
  {
    for(common::MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
      auto vh = *v_it;
      if (!extended_vertex_static_infos.at(vh).is_original_vertex) {
	const auto new_pos_e =  iterateVertex(mesh, vh, extended_vertex_static_infos.at(vh));
	common::MyMesh::Point new_pos(new_pos_e[0], new_pos_e[1], new_pos_e[2]);
	mesh.point(vh) = new_pos;
      }
    }
  }

  void DiscreteFairer::execute(common::MyMesh& mesh, size_t iteration_count, std::function<void(int)> cb)
  {
    
    //return triangleExecuteDemo(mesh);
    
    
    extended_vertex_static_infos = generateExtendedVertexSaticInfos(mesh, mesh.children_parents_map); //todo ezt mindig ujra kell generalni?

    
    //metrics::CurvatureBased cb(mesh, *this);
    //cb.startSession();

    OpenMesh::VPropHandleT<double> demo_color;
    mesh.add_property(demo_color, "demo_color");
    
    for(size_t i = 0; i < iteration_count ; i++){

     std::cout<<"-----Iter--------"<<std::endl;

      // Calculate the curvature for each vertex at the beginning of each iteration
      CurvatureCalculator mcc(mesh, true);
      
      //TODO range operator (smarthandle....)
      for(common::MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
	auto vh = *v_it;
	if(extended_vertex_static_infos.at(vh).is_original_vertex){
	  mcc.execute(vh);
	  vertex_curvature_map[vh] =  getCurvature(mcc);
	}
      }
      //cb.postIteration();
      
      if (common::settings::sync) {
	iterateVerticesSync(mesh);
      }
      else {
	iterateVerticesAsync(mesh);
      }
      

      for(auto vh : mesh.vertices()){
	CurvatureCalculator cc(mesh, true);
	cc.execute(vh);
	//std::cout << "ID: " << vh.idx() << std::endl;
	//std::cout << "exp: "<< calcTargetCurvature(extended_vertex_static_infos.at(vh).weighed_effectors).main << std::endl;
	//std::cout << "real: "<<cc.getMeanCurvature() << std::endl;
	mesh.property(demo_color, vh) = cc.getGaussianCurvature();

      }

      cb(static_cast<int>((i+1) * 100.0 / iteration_count));

    }
    //cb.endSession();



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


  Eigen::Vector3d DiscreteFairer::Q2(const std::vector<Eigen::Vector3d>& p,const Eigen::Vector3d& normal, double H,  const CurvatureCalculator::FundamentalElements& fe,
				     const Eigen::Vector3d& Q, const Eigen::Vector3d& Q0,
				     const Eigen::Matrix<double, 5 , Eigen::Dynamic>& M)
  {

    //const auto p_k = common::average(p);
    const auto neighbour_count = p.size();

    Eigen::RowVectorXd row3 = M.row(2);
    double dotsum_m3 = 0.0;
    for(size_t i = 0; i < neighbour_count; i++) {
      dotsum_m3 += (row3(i) * p[i]).dot(normal);
    }

    Eigen::RowVectorXd row4 = M.row(3);
    double dotsum_m4 = 0.0;
    for(size_t i = 0; i< neighbour_count; i++) {
      dotsum_m4 += (row4(i) * p[i]).dot(normal);
    }

    Eigen::RowVectorXd row5 = M.row(4);
    double dotsum_m5 = 0.0;
    for(size_t i = 0; i< neighbour_count; i++) {
      dotsum_m5 += (row5(i) * p[i]).dot(normal);
    }

    const double t = (-2.0 * H*(fe.E * fe.G - fe.F * fe.F) + fe.E * dotsum_m5 - 2.0 * fe.F * dotsum_m4 + fe.G * dotsum_m3) / ((fe.E * M.row(4).sum())-2.0 * fe.F * M.row(3).sum() + fe.G * M.row(2).sum()) - Q0.dot(normal);
    return Q0 + normal * t;
    

  }




  Eigen::Vector3d DiscreteFairer::Q_Gaussian(const std::vector<Eigen::Vector3d>& p,const Eigen::Vector3d& normal, TargetCurvature H,  const CurvatureCalculator::FundamentalElements& fe,
					     const Eigen::Vector3d& Q, const Eigen::Vector3d& Q0,
					     const Eigen::Matrix<double, 5 , Eigen::Dynamic>& M)
  {
    bool chose_max = true;
    if (H.main < 0) {
      chose_max = !chose_max;
      H.main *= -1;
    }


    //const auto p_k = common::average(p);
    const auto neighbour_count = p.size();


    const double row3sum = M.row(2).sum();
    const double row4sum = M.row(3).sum();
    const double row5sum = M.row(4).sum();

    Eigen::RowVectorXd row3 = M.row(2);
    double a3 = 0.0;
    for(size_t i = 0; i< neighbour_count; i++) {
      a3 += (row3(i) * p[i]).dot(normal);
    }

    Eigen::RowVectorXd row4 = M.row(3);
    double a4 = 0.0;
    for(size_t i = 0; i< neighbour_count; i++) {
      a4 += (row4(i) * p[i]).dot(normal);
    }

    Eigen::RowVectorXd row5 = M.row(4);
    double a5 = 0.0;
    for(size_t i = 0; i< neighbour_count; i++) {
      a5 += (row5(i) * p[i]).dot(normal);
    }
    const auto rat = H.main * (fe.E * fe.G - fe.F * fe.F);

    const auto pkndot = Q0.dot(normal);


    const auto a = row3sum * row5sum - row4sum * row4sum;
    const auto b = -a3 * row5sum -a5 * row3sum + 2 * row3sum * row5sum * pkndot + 2 * a4 * row4sum - 2 * pkndot * row4sum * row4sum;
    const auto c = a3 * a5 - a3 * pkndot * row5sum - a5 * pkndot * row3sum + row5sum * row3sum * pkndot * pkndot - a4 * a4 + 2 * a4 * pkndot * row4sum - pkndot * pkndot * row4sum * row4sum - rat;


    double t = 0;
    
    const auto determinant = b * b - 4 * a * c;
    if (determinant < -1e-5) {
      std::cout << H.main << std::endl;
      std::cout << "Gaussian Q det error. "<< determinant << std::endl;
      return Q0;
    }
    else if (determinant < 1e-5) {
      t = - b / (2 * a);
    }
    else {
      const auto t1 = (- b + sqrt(determinant)) / (2 * a);
      const auto t2 = (- b - sqrt(determinant)) / (2 * a);
      t = chose_max ? std::max(t1,t2) : std::min(t1, t2);
      //t = std::fabs(t1) < std::fabs(t2) ? t1 : t2;
      //t = std::min(t1,t2);
    }

    return Q0 + normal * t;
  }

  double DiscreteFairer::getCurvature(const CurvatureCalculator& cc) const
  {
    switch(common::settings::selected_curvature) {
    case common::settings::MEAN:
      return cc.getMeanCurvature();
    case common::settings::GAUSSIAN:
      return calcSignedGaussianCurvature(cc);
    }
  }

  

}
