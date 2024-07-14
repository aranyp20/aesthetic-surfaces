#pragma once

#include <eigen3/Eigen/Dense>


namespace LSQPlane {

  //using namespace Geometry;
  typedef Eigen::Vector3d Point3D;
  typedef Eigen::Vector3d Vector3D;
  typedef std::vector<Eigen::Vector3d> PointVector;
  typedef std::vector<Eigen::Vector2d> Point2DVector;

struct Plane {
  Point3D p;
  Vector3D u, v, n;
};

Plane fitPlane(const PointVector &pv);

Point2DVector projectToBestFitPlane(const PointVector &pv);

}
