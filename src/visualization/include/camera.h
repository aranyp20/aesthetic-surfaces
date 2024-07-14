#pragma once


#include <eigen3/Eigen/Dense>


class Camera
{

    Eigen::Vector3d w_eye, w_up, w_lookat;
    double fov, aspect, fp, bp;

public:

    Camera();

    Eigen::Matrix4d V() const;
    Eigen::Matrix4d P() const;

};