#pragma once
//#include <Helper/Helper.hpp>
#include <Vertex2Vertex.hpp>
// 验证
void laplacianValidate(const SpMat & L);

// 5.1.1 度数拉普拉斯矩阵
SpMat buildLaplacianDegree(const TriMesh &mesh);

// 5.1.2 Tutte拉普拉斯
SpMat buildLaplacianTutte(const TriMesh &mesh);


// 5.1.3 归一化图形拉普拉斯
SpMat builLaplacianNormalizedGraph(const TriMesh &mesh);





