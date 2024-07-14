#include "camera.h"



Camera::Camera()
{
    w_lookat = {0, 0, 0};
    w_up = {0, 0, 1};
    w_eye = {-3, 0, 0};

    fp = 1;
    bp = 20;

    aspect = 1.0;

    fov = 70.0 * M_PI / 180.0;
}

Eigen::Matrix4d Camera::V() const
{
    const auto w = (w_eye - w_lookat).normalized();
    const auto u = w_up.cross(w).normalized();
    const auto v = w.cross(u); // normalized by default

    Eigen::MatrixXd translate_to_origin = Eigen::Matrix4d::Identity();
    translate_to_origin.block<3, 1>(0, 3) = Eigen::Vector3d(-w_eye);

    Eigen::Matrix4d rotate_to_axis;
    rotate_to_axis << u[0], u[1], u[2], 0,
                      v[0], v[1], v[2], 0,
                      w[0], w[1], w[2], 0,
                      0,    0,    0,    1;

    return rotate_to_axis * translate_to_origin;
}

Eigen::Matrix4d Camera::P() const
{
    Eigen::Matrix4d retval;
    retval << 1 / (tan(fov / 2) * aspect), 0, 0, 0,
              0, 1 / tan(fov / 2), 0, 0,
              0, 0, -(fp + bp) / (bp - fp), -2 * fp * bp / (bp - fp),
              0, 0, -1, 0;
            
    return retval;
}