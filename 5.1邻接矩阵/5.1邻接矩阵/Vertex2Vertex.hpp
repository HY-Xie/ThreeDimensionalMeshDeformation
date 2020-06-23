#pragma once
#include <Helper/Helper.hpp>

Eigen::MatrixXd buildAdjacentMatrixVV_Dense(const TriMesh &mesh);

SpMat buildAdjacentMatrix_VV(const TriMesh &mesh);

SpMat buildMatrixDegree(const TriMesh &mesh);

SpMat buildAdjacentMatrix_FV(const TriMesh &mesh);

SpMat buildAdjacentMatrix_FF(const TriMesh &mesh);

SpMat buildAdjacentMatrix_VE(const TriMesh &mesh);

SpMat buildAdjacentMatrix_VF(const TriMesh &mesh); // transpose of FV

SpMat buildAdjacentMatrix_EE(const TriMesh &mesh);	// 总结

SpMat buildAdjacentMatrix_EV(const TriMesh &mesh); // trranspose of VE

SpMat buildAdjacentMatrix_EF(const TriMesh &mesh);	// 总结

SpMat buildAdjacentMatrix_FE(const TriMesh &mesh);	// transpose of EF