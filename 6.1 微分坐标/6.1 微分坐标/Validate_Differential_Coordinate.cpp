#include "Validate_Differential_Coordinate.hpp"
#include "CombinatorialLaplacian.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <fstream>
const double EPS = 1e-6;
bool matrixIsEqual(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B)
{
	assert(A.size() == B.size());
	int numRows = A.rows(); 
	int numCols = A.cols();
	for(int i = 0; i < numRows; ++i)
		for (int j = 0; j < numCols; ++j)
		{
			if (fabs(A.coeff(i, j) - B.coeff(i, j)) > EPS)
				return false;
		}
	return true;
}

void isEqual(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B)
{
	if (matrixIsEqual(A, B))
		std::cout << "矩阵相等" << std::endl;
	else
		std::cout << "矩阵不等" << std::endl;
}
// 测试微分坐标的方向是否与顶点法向方向相同
bool validate_normal_direction(const Eigen::MatrixXd &D, TriMesh &mesh)
{
	if (!mesh.has_face_normals())
		mesh.request_face_normals();
	if (!mesh.has_vertex_normals())
		mesh.request_vertex_normals();
	mesh.update_normals();
	

	int n = mesh.n_vertices();
	Eigen::MatrixXd normals(n, 3);
	int i = 0;
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		TriMesh::Point normal_openmesh = mesh.normal(*v_it);
		Eigen::Vector3d normal(normal_openmesh.data()[0], normal_openmesh.data()[1], normal_openmesh.data()[2]);
		normals.row(i) = normal;
		i++;
	}


	i = 0;
	for (; i < n; ++i)
	{
		Eigen::Vector3d vd = D.row(i);
		Eigen::Vector3d vn = normals.row(i);
		double cos_theta = (vd.dot(vn) / (vd.lpNorm<2>() * vn.lpNorm<2>()));
		if (cos_theta < 0.8)
			return false;
	}
	return true;
}


void validate_differential_coordiate(const Eigen::MatrixXd &D, TriMesh &mesh)
{
	// 1. 测试方向是否与法向相同
	if (validate_normal_direction(D, mesh))
		std::cout << "与法向方向一致---YES" << std::endl;
	else
		std::cout << "与法向不一致---NO" << std::endl;



}

void computeOriginalMesh(const Eigen::MatrixXd &D,  const TriMesh &mesh)
{
	SpMat M = buildLaplacianTutte(mesh);
	Eigen::MatrixXd M_dense = M;
	Eigen::MatrixXd V = M_dense.inverse() * D;
	// 写入obj
	std::ofstream of("outMesh.obj");
	if (of.is_open())
	{
		std::cout << "开始写入" << std::endl;
		of << "# Vertices: " << M_dense.rows() << " Faces: " << mesh.n_faces() << std::endl;
		for (int i = 0; i < V.rows(); ++i)
		{
			of << "v " << V.row(i) << std::endl;
		}
		for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
		{
			of << "f ";
			for (auto fv_it = mesh.cfv_ccwiter(*f_it); fv_it.is_valid(); ++fv_it)
			{
				of << fv_it->idx() + 1 << " ";
			}
			of << std::endl;
		}
		of.close();
	}
	std::cout << "写入完成" << std::endl;
}

