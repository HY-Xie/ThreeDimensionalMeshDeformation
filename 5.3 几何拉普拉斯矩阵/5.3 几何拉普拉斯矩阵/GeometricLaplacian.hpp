#pragma once

#include <Helper/Helper.hpp>

void laplacianValidate(const SpMat & L);

// 5.3.1 余切拉普拉斯
SpMat buildLaplacianCot(const TriMesh &mesh);


// 5.3.2 面积余切拉普拉斯
SpMat buildLaplacianCotArea(const TriMesh &mesh);

// 5.3.3 中值拉普拉斯
SpMat buildLaplacianMeanValue(const TriMesh &mesh);
