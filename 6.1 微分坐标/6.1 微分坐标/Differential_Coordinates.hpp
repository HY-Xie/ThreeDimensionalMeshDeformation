#pragma once

#include <Helper/Helper.hpp>

Eigen::MatrixXd buildDifferentialCoord_fromDefinition(const TriMesh &mesh);

Eigen::MatrixXd buildDifferentialCoord_laplace(const TriMesh &mesh);
