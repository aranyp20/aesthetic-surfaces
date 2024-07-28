#pragma once

#include "curvaturecalculator.h"
#include "subdivider.h"


namespace core {


class DiscreteFairer
{
  std::map<common::MyMesh::VertexHandle, double> vertex_curvature_map;

  
public:



  
  struct ExtendedVertexStaticInfo
  {
    bool is_original_vertex = true;
    std::vector<std::pair<common::MyMesh::VertexHandle, double>> weighed_effectors;
  };

  double calcTargetCurvature(const std::vector<std::pair<common::MyMesh::VertexHandle, double>>& weighed_effectors) const;

  void triangleExecuteDemo(common::MyMesh& mesh);

  std::map<common::MyMesh::VertexHandle, DiscreteFairer::ExtendedVertexStaticInfo> generateExtendedVertexSaticInfos(common::MyMesh& mesh, const Subdivider::ChildrenParents& child_parents_map) const;

    void getEffectorsHelper(const common::MyMesh::VertexHandle& to,
			    std::set<common::MyMesh::VertexHandle>& visited_vertices,
			    std::set<common::MyMesh::VertexHandle>& effectors,
			    const Subdivider::ChildrenParents& child_parents_map) const;

  
  std::set<common::MyMesh::VertexHandle> getEffectors(const common::MyMesh::VertexHandle& to, const Subdivider::ChildrenParents& child_parents_map) const;

      std::vector<std::pair<common::MyMesh::VertexHandle, double>> getWeighedEffectors(const common::MyMesh::VertexHandle& to, const Subdivider::ChildrenParents& child_parents_map,const common::MyMesh& mesh) const;
  

  Eigen::Vector3d iterateVertex(common::MyMesh&, common::MyMesh::VertexHandle& iteratable, const DiscreteFairer::ExtendedVertexStaticInfo& extended_vertex_static_info) const;

//TODO replace
static Eigen::Vector3d Q(const std::array<Eigen::Vector3d, 6>& p,
			 const Eigen::Vector3d& normal, double H, const CurvatureCalculator::FundamentalElements& fe,
			 const Eigen::Vector3d& Q);


  static Eigen::Vector3d Q2(const std::array<Eigen::Vector3d, 6>& p,
			 const Eigen::Vector3d& normal, double H, const CurvatureCalculator::FundamentalElements& fe,
			    const Eigen::Vector3d& Q, const Eigen::Matrix<double, 5, 6>& M);

public:

void execute(common::MyMesh& mesh, size_t face_split_count, size_t iteration_count);

  const Subdivider::ChildrenParents& getChildParentsMap() const;
  

  std::map<common::MyMesh::VertexHandle, DiscreteFairer::ExtendedVertexStaticInfo> extended_vertex_static_infos;

};

}
//point-to-register, jump-to-register, c-x 0
