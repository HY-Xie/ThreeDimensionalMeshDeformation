#include "GeometricLaplacian.hpp"
const double EPS = 1e-6;
#include <limits>
#include <iomanip>

void computeAreaMix(const TriMesh &mesh, std::vector<double> &areas)  {
	areas.clear();
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		double area = 0.0;
		
		for (auto vf_it = mesh.cvf_iter(*v_it); vf_it.is_valid(); ++vf_it)
		{
			Eigen::Vector3d v[3];
			int i = 0;
			for (auto vfv_it = mesh.cfv_iter(*vf_it); vfv_it.is_valid(); ++vfv_it)
			{
				v[i] = Eigen::Vector3d(mesh.point(*vfv_it).data()[0], mesh.point(*vfv_it).data()[1], mesh.point(*vfv_it).data()[1]);
				++i;
			}
			area += 0.5 * (v[2] - v[0]).cross(v[1] - v[0]).norm() * 0.33333;
		}
		areas.push_back(area);
	}
}

SpMat buildLaplacianCot(const TriMesh &mesh)
{
	int n = mesh.n_vertices();
	SpMat cotLaplacian; cotLaplacian.resize(n, n);

	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		Eigen::Vector3d v[3] = { Eigen::Vector3d(0.0, 0.0, 0.0) };
		int c[3] = { 0 };
		double cotx[3] = { 0 };
		int i = 0;
		for (auto fv_it = mesh.cfv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			c[i] = fv_it->idx();
			v[i] = Eigen::Vector3d(mesh.point(*fv_it).data()[0], mesh.point(*fv_it).data()[1], mesh.point(*fv_it).data()[1]);
			++i;
		}
		cotx[0] = (v[1] - v[0]).dot(v[2] - v[0]) / (v[1] - v[0]).cross(v[2] - v[0]).lpNorm<2>();
		cotx[1] = (v[2] - v[1]).dot(v[0] - v[1]) / (v[2] - v[1]).cross(v[0] - v[1]).lpNorm<2>();
		cotx[2] = (v[1] - v[2]).dot(v[0] - v[2]) / (v[1] - v[2]).cross(v[0] - v[2]).lpNorm<2>();

		cotLaplacian.coeffRef(c[0], c[1]) += -cotx[2] / 2; cotLaplacian.coeffRef(c[1], c[0]) += -cotx[2] / 2;
		cotLaplacian.coeffRef(c[1], c[2]) += -cotx[0] / 2; cotLaplacian.coeffRef(c[2], c[1]) += -cotx[0] / 2;
		cotLaplacian.coeffRef(c[0], c[2]) += -cotx[1] / 2; cotLaplacian.coeffRef(c[2], c[0]) += -cotx[1] / 2;
	}
	for (int i = 0; i < n; ++i)
	{
		double sum = 0;
		sum = -cotLaplacian.row(i).sum();
		cotLaplacian.coeffRef(i, i) += sum;
	}
	
	return cotLaplacian;
}


SpMat buildLaplacianCotArea(const TriMesh &mesh)
{
	int n = mesh.n_vertices();
	std::vector<double> areas;
	computeAreaMix(mesh, areas);
	SpMat cotAreaLaplacian; cotAreaLaplacian.resize(n, n);
	cotAreaLaplacian = buildLaplacianCot(mesh);
	for (int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		cotAreaLaplacian.coeffRef(i,j) = cotAreaLaplacian.coeff(i,j) / areas[i];
	}

	return cotAreaLaplacian;

}


SpMat buildLaplacianMeanValue(const TriMesh &mesh)
{
	std::setprecision(10);
	int n = mesh.n_vertices();
	SpMat L; L.resize(n, n);
	std::vector<T> coeffs; coeffs.clear();

	for (auto h_it = mesh.halfedges_begin(); h_it != mesh.halfedges_end(); ++h_it)
	{
		TriMesh::VertexHandle fromJ = mesh.from_vertex_handle(*h_it);
		TriMesh::VertexHandle toI = mesh.to_vertex_handle(*h_it);
		Eigen::Vector3d J = Eigen::Vector3d(mesh.point(fromJ).data()[0], mesh.point(fromJ).data()[1], mesh.point(fromJ).data()[2]);
		Eigen::Vector3d I = Eigen::Vector3d(mesh.point(toI).data()[0], mesh.point(toI).data()[1], mesh.point(toI).data()[2]);
		Eigen::Vector3d jToi = (I-J).normalized();

		TriMesh::VertexHandle toT = mesh.to_vertex_handle(mesh.next_halfedge_handle(*h_it));
		Eigen::Vector3d M = Eigen::Vector3d(mesh.point(toT).data()[0], mesh.point(toT).data()[1], mesh.point(toT).data()[2]);
		Eigen::Vector3d alphaToI = (M - I).normalized();

		TriMesh::VertexHandle toK = mesh.from_vertex_handle(mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(*h_it)));
		Eigen::Vector3d K = Eigen::Vector3d(mesh.point(toK).data()[0], mesh.point(toK).data()[1], mesh.point(toK).data()[2]);
		Eigen::Vector3d betaTOI = (K - I).normalized();


		double cosGammaIJ = jToi.dot(alphaToI) / (jToi.norm() * alphaToI.norm());
		double cosThetaIJ = jToi.dot(betaTOI) / (jToi.norm() * betaTOI.norm());
		double angleGamma = acosf(cosGammaIJ);
		double angleTheta = acosf(cosThetaIJ);
		double wij = (tan(angleGamma / 2) + tan(angleTheta / 2)) / (I - J).norm();
		int idxI = toI.idx();
		int idxJ = fromJ.idx();

		coeffs.push_back(T(idxI, idxJ, wij));
	}

	L.setFromTriplets(coeffs.begin(), coeffs.end());


	for (int i = 0; i < n; ++i)
	{
		double sum = L.row(i).sum();
		L.coeffRef(i, i) = -sum;
	}

	return L;
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
		{
			std::cout <<"第 " << i << " 行: " <<  sum << std::endl;
			return false;
		}
			
	}
	return true;
}

bool allPositive(const SpMat &L)
{
	int numRows = L.rows();
	int numCols = L.cols();
	for (int i = 0; i < numRows; ++i) {
		for (int j = 0; j < numCols; ++j) {
			if (L.coeff(i, j) < 0)
			{
				return false;
			}
		}

	}
	return true;

}

bool nonDiagPositive(const SpMat &L)
{
	int numRows = L.rows();
	int numCols = L.cols();
	for (int i = 0; i < numRows; ++i) {
		for (int j = 0; j < numCols; ++j) {
			if (i != j && L.coeff(i, j) < 0)
				return false;
		}
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

	if (allPositive(L))
		std::cout << "全部正正值" << std::endl;
	else
		std::cout << "有负值" << std::endl;

	if (nonDiagPositive(L))
		std::cout << "非对角线元素全为正" << std::endl;
	else
		std::cout << "非对角线元素有负值" << std::endl;
}