#include "Differential_Coordinates.hpp"
#include <cstdlib>
#include <CombinatorialLaplacian.hpp>

Eigen::MatrixXd buildDifferentialCoord_fromDefinition(const TriMesh &mesh)
{
	size_t numVer = mesh.n_vertices();
	Eigen::MatrixXd diffCoord(numVer, 3);
	diffCoord.setZero();
	int i = 0;
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		
		TriMesh::Point neighbor_sum(0.f, 0.f, 0.f);
		int numNeighbors = 0;
		for (auto vv_it = mesh.cvv_iter(*v_it); vv_it.is_valid(); ++vv_it)
		{
			neighbor_sum += mesh.point(*vv_it);
			numNeighbors++;
		}
		neighbor_sum /= numNeighbors;
		TriMesh::Point currentDiffCoord_openmesh= mesh.point(*v_it) - neighbor_sum;
		Eigen::Vector3d currentDiffCoord(currentDiffCoord_openmesh.data()[0], currentDiffCoord_openmesh.data()[1], currentDiffCoord_openmesh.data()[2]);
		diffCoord.row(i) = currentDiffCoord;
		++i;
	}

	return diffCoord;
}

Eigen::MatrixXd buildDifferentialCoord_laplace(const TriMesh &mesh)
{
	// 1. 构建欧式坐标矩阵 P
	size_t n = static_cast<int>(mesh.n_vertices());
	Eigen::MatrixXd P(n, 3); P.setZero();

	int i = 0;
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		Eigen::Vector3d vertex(mesh.point(*v_it).data()[0], mesh.point(*v_it).data()[1], mesh.point(*v_it).data()[2]);
		P.row(i) = vertex;
		i++;
	}

	// 2. 构建矩阵M
	SpMat M = buildLaplacianTutte(mesh);

	// D = M * P
	Eigen::MatrixXd D = M * P;
	
	return D;
}