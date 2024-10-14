#pragma once

#include "curvaturecalculator.h"
#include "mesh.h"
#include "settings.h"
#include "subdivider.h"
#include <functional>


namespace core {


class DiscreteFairer
{
public:

  struct ExtendedCurvature
  {
    double val = 0;
    bool negative_positive = false;
    bool is_saddle = false;
  };

private:
  std::map<common::MyMesh::VertexHandle, ExtendedCurvature> vertex_curvature_map;

  
public:



  
  struct ExtendedVertexStaticInfo
  {
    bool is_original_vertex = true;
    std::vector<std::pair<common::MyMesh::VertexHandle, double>> weighed_effectors;
  };

  void sortEffectors(const common::MyMesh& mesh, std::vector<common::MyMesh::VertexHandle>& effectors) const;

  struct TargetCurvature
  {
    double main = 0;
    double k_max = 0;
    double k_min = 0;

    bool negative_positive = false;
    bool all_saddle = false;
  };

  struct EffectorExtra
  {
    double weight;
    double H;
    Eigen::Vector3d pos;

    Eigen::Vector3d subject_pos; //TODO out of here

    bool negative_positive = false;
    bool is_saddle = false;
    int d_id;
  };

  TargetCurvature calcTargetCurvature(const std::vector<EffectorExtra>& weighed_effectors) const;

  TargetCurvature calcGaussianTargetCurvature(const std::vector<EffectorExtra>& weighed_effectors) const;


  /* Precondition: weights are all same signed */
  TargetCurvature logAestheticTargetCurvatureCore(const std::vector<EffectorExtra>& effectors,
						  double alpha = common::settings::log_aesthetic_alpha,
						  bool enable_flipping = true) const;

  /* Precondition: negative and positive weight(s) are present as well.
     Precondition II: ordered by curvature. */
  TargetCurvature logAestheticTargetCurvatureMixed(const std::vector<EffectorExtra>& effectors) const;
  
  TargetCurvature calcLogAestheticTargetCurvature(const std::vector<EffectorExtra>& weighed_effectors) const;
  TargetCurvature calcLogAestheticTargetCurvature_offsetVersion(const std::vector<EffectorExtra>& weighed_effectors) const;


    
  
  
  void triangleExecuteDemo(common::MyMesh& mesh);
  void pentaExecuteDemo(common::MyMesh& mesh);
  void triangleGaussExecuteDemo(common::MyMesh& mesh);

  std::map<common::MyMesh::VertexHandle, DiscreteFairer::ExtendedVertexStaticInfo> generateExtendedVertexSaticInfos(common::MyMesh& mesh, const ChildrenParents& child_parents_map) const;

    void getEffectorsHelper(const common::MyMesh::VertexHandle& to,
			    std::set<common::MyMesh::VertexHandle>& visited_vertices,
			    std::set<common::MyMesh::VertexHandle>& effectors,
			    const ChildrenParents& child_parents_map) const;

  
  std::set<common::MyMesh::VertexHandle> getEffectors(const common::MyMesh::VertexHandle& to, const ChildrenParents& child_parents_map) const;

      std::vector<std::pair<common::MyMesh::VertexHandle, double>> getWeighedEffectors(const common::MyMesh::VertexHandle& to, const ChildrenParents& child_parents_map,const common::MyMesh& mesh) const;
  

  Eigen::Vector3d iterateVertex(common::MyMesh&, common::MyMesh::VertexHandle& iteratable, const DiscreteFairer::ExtendedVertexStaticInfo& extended_vertex_static_info) const;


  void iterateVerticesAsync(common::MyMesh& mesh);
  void iterateVerticesSync(common::MyMesh& mesh);


  ExtendedCurvature getCurvature(const CurvatureCalculator& cc) const;
  
//TODO replace
static Eigen::Vector3d Q(const std::array<Eigen::Vector3d, 6>& p,
			 const Eigen::Vector3d& normal, double H, const CurvatureCalculator::FundamentalElements& fe,
			 const Eigen::Vector3d& Q);


  static Eigen::Vector3d Q2(const std::vector<Eigen::Vector3d>& p,
			    const Eigen::Vector3d& normal, double H, const CurvatureCalculator::FundamentalElements& fe,
			    const Eigen::Vector3d& Q, const Eigen::Vector3d& Q0,
			    const Eigen::Matrix<double, 5, Eigen::Dynamic>& M);


   static Eigen::Vector3d Q_Gaussian(const std::vector<Eigen::Vector3d>& p,
				     const Eigen::Vector3d& normal, TargetCurvature H, const CurvatureCalculator::FundamentalElements& fe,
				     const Eigen::Vector3d& Q, const Eigen::Vector3d& Q0,
				     const Eigen::Matrix<double, 5, Eigen::Dynamic>& M,
				     bool negative_positive);
public:

  void execute(common::MyMesh& mesh, size_t iteration_count, std::function<void(int)> cb);

  const ChildrenParents& getChildParentsMap() const;
  

  std::map<common::MyMesh::VertexHandle, DiscreteFairer::ExtendedVertexStaticInfo> extended_vertex_static_infos;

};

}
//point-to-register, jump-to-register, c-x 0
