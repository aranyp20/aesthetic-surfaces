#pragma once

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include "common_defines.h"

namespace core
{
    class CurvatureCalculator
    {

        struct SurfaceParam
	{
            double u = 0;
            double v = 0;
        };

        struct InputPoints
        {
            Eigen::Vector3d center{0, 0, 0};
            std::vector<Eigen::Vector3d> P;

            InputPoints translatePoints() const;

            std::vector<SurfaceParam> calcUVs(bool use_adaptive_uvs) const;
        };

        struct DerResults
        {
            Eigen::Vector3d Su{0, 0, 0};
            Eigen::Vector3d Sv{0, 0, 0};
            Eigen::Vector3d Suu{0, 0, 0};
            Eigen::Vector3d Suv{0, 0, 0};
            Eigen::Vector3d Svv{0, 0, 0};
        } eq;

      DerResults calcDer(const InputPoints &ipp);
        Eigen::Vector3d S(const double u, const double v, const DerResults &Ss) const;
        void calcCurvature(const InputPoints &ipp); //TODO rename
        
        
        //void calcCurvatures(common::MyMesh &mesh) const;
      void tessellateSurface(const size_t resolution, DerResults Ss) const;


      
      Eigen::Matrix2d getShapeOperator() const;

    public:
        struct FundamentalElements
        {
            double E = 0;
            double F = 0;
            double G = 0;

            double L = 0;
            double M = 0;
            double N = 0;
        };

      struct PrincipleCurvatures
      {
	bool umbolic;
	Eigen::Vector3d max_dir;
	Eigen::Vector3d min_dir;
	double max_val;
	double min_val;
      };

      PrincipleCurvatures getPrincipleCurvatures() const;

      
      
        FundamentalElements getFundamentalElements() const;
        Eigen::Vector3d getNormal() const;
      //TODO: rename to mean curvature 
        double getMeanCurvature() const;
        double getGaussianCurvature() const;


        void execute(common::MyMesh::VertexHandle &vh);
        void execute(const Eigen::Vector3d vertex_pos, const std::vector<Eigen::Vector3d>& neighbors);


      CurvatureCalculator(common::MyMesh &mesh, bool _use_adaptive_uvs = true);


      Eigen::Matrix<double, 5, Eigen::Dynamic> getLastM() const;
    private:

        common::MyMesh& mesh;

        FundamentalElements fundamental_elements;
        Eigen::Vector3d normal;

      bool use_adaptive_uvs = false;

      Eigen::Matrix<double, 5, Eigen::Dynamic> last_M;
    };

}
