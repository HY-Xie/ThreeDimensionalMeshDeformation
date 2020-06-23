#pragma once
#include "Differential_Coordinates.hpp"

void isEqual(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B);
void validate_differential_coordiate(const Eigen::MatrixXd &D, TriMesh &mesh);
void computeOriginalMesh(const Eigen::MatrixXd &D,  const TriMesh &mesh);
