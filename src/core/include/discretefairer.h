#pragma once

#include "curvaturecalculator.h"


namespace core {


class DiscreteFairer
{
    public:
//TODO replace
static Eigen::Vector3d Q(const std::array<Eigen::Vector3d, 6>& p,
			 const Eigen::Vector3d& normal, double H, const CurvatureCalculator::FundamentalElements& fe,
			 const Eigen::Vector3d& Q);


  static Eigen::Vector3d Q2(const std::array<Eigen::Vector3d, 6>& p,
			 const Eigen::Vector3d& normal, double H, const CurvatureCalculator::FundamentalElements& fe,
			    const Eigen::Vector3d& Q, const Eigen::Matrix<double, 5, 6>& M);

public:

void execute(common::MyMesh& mesh, size_t face_split_count, size_t iteration_count);



};

}
//point-to-register, jump-to-register, c-x 0
