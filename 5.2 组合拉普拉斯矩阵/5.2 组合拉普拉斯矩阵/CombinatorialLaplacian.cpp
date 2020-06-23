#include "CombinatorialLaplacian.hpp"
#include <cmath>
const double EPS = 1e-9;

SpMat buildLaplacianDegree(const TriMesh &mesh)
{
	SpMat laplacian;
	SpMat adjVertices = buildAdjacentMatrix_VV(mesh);	// W
	SpMat adjDegree = buildMatrixDegree(mesh);	// D
	laplacian = adjDegree - adjVertices; // L = D-W

	return laplacian;
}

SpMat buildLaplacianTutte(const TriMesh &mesh)
{
	SpMat tutte; tutte.resize(mesh.n_vertices(), mesh.n_vertices());
	std::vector<T> coeffs; coeffs.clear();
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		int i = v_it->idx();
		int valance = mesh.valence(*v_it);
		for (auto vv_it = mesh.cvv_iter(*v_it); vv_it.is_valid(); ++vv_it)
		{
			int j = vv_it->idx();
			coeffs.push_back(T(i, j, -1.0 / valance));
		}
	}
	tutte.setFromTriplets(coeffs.begin(), coeffs.end());

	for (int i = 0; i < mesh.n_vertices(); ++i)
		tutte.coeffRef(i, i) = 1.0;

	return tutte;
}

// 5.1.3 归一化图形拉普拉斯
SpMat builLaplacianNormalizedGraph(const TriMesh &mesh)
{
	int n = mesh.n_vertices();
	SpMat normalizeGraph; normalizeGraph.resize(n, n);
	std::vector<T> coeffs;

	for (int i = 0; i < n; ++i)
		normalizeGraph.coeffRef(i, i) = 1.0;

	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		int i = v_it->idx();
		int i_valance = mesh.valence(*v_it);
		for (auto vv_it = mesh.cvv_iter(*v_it); vv_it.is_valid(); ++vv_it)
		{
			int j = vv_it->idx();
			int j_valance = mesh.valence(*vv_it);
			coeffs.push_back(T(i, j, -1.0 / sqrtf(i_valance * j_valance)));
		}
	}

	normalizeGraph.setFromTriplets(coeffs.begin(), coeffs.end());

	return normalizeGraph;
}


bool isZero(const SpMat &L)
{
	int numRows = L.rows();
	int numCols = L.cols();
	// 每行和为0
	//std::cout << "The sum of each row is zero." << std::endl;
	for (int i = 0; i < numRows; ++i) {
		double sum = 0;
		for (int j = 0; j < numCols; ++j) {
			sum += L.coeff(i, j);
		}
		if (fabs(sum) > EPS)
			return false;
	}
	return true;
}


bool isSymmetric(const SpMat &L)
{
	int numRows = L.rows();
	int numCols = L.cols();
	// 对称： L矩阵本身不是对称矩阵
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++)
		{
			if (fabs(L.coeff(i, j) - L.coeff(j, i)) > EPS)
				return false;
		}
	}
	return true;
}


void laplacianValidate(const SpMat & L)
{
	
	if (isZero(L))
		std::cout << "每行和为零." << std::endl;
	else
		std::cout << "不是每行和为零." << std::endl;

	if (isSymmetric(L))
		std::cout << "对称矩阵" << std::endl;
	else
		std::cout << "非对称" << std::endl;
}



