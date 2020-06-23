#include "Vertex2Vertex.hpp"

Eigen::MatrixXd buildAdjacentMatrixVV_Dense(const TriMesh &mesh)
{
	const int numOfVertices = static_cast<int>(mesh.n_vertices());
	
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> v2vMatrix;
		v2vMatrix.resize(numOfVertices, numOfVertices);
		v2vMatrix.setZero();

		// const default
		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
		{
			int firstIdx = v_it->idx();
			for (auto vv_it = mesh.cvv_iter(*v_it); vv_it.is_valid(); ++vv_it)
			{
				int secIdx = vv_it->idx();
				v2vMatrix(firstIdx, secIdx) = 1.0;
			}
		}

		return v2vMatrix;
}

SpMat buildAdjacentMatrix_VV(const TriMesh &mesh)
{
	const int numOfVertices = static_cast<int>(mesh.n_vertices());
	SpMat v2vMatrix;
	v2vMatrix.resize(numOfVertices, numOfVertices);
	std::vector<T> coeffs;
	// const default
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		int firstIdx = v_it->idx();
		for (auto vv_it = mesh.cvv_iter(*v_it); vv_it.is_valid(); ++vv_it)
		{
			int secIdx = vv_it->idx();
			coeffs.push_back(T(firstIdx, secIdx, 1));
		}
	}

	v2vMatrix.setFromTriplets(coeffs.begin(), coeffs.end());

	return v2vMatrix;
}

SpMat buildMatrixDegree(const TriMesh &mesh)
{
	const int numOfVertices = static_cast<int>(mesh.n_vertices());
	SpMat degreeMatrix;
	degreeMatrix.resize(numOfVertices, numOfVertices);
	std::vector<T> coeffs; coeffs.clear();

	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		int idx = v_it->idx();
		int valence = mesh.valence(*v_it);
		coeffs.push_back(T(idx, idx, valence));
	}
	degreeMatrix.setFromTriplets(coeffs.begin(), coeffs.end());
	return degreeMatrix;
}

SpMat buildAdjacentMatrix_FV(const TriMesh &mesh)
{
	const int numOfVertices = static_cast<int>(mesh.n_vertices());
	const int numOfFaces = static_cast<int>(mesh.n_faces());
	SpMat face_vertex_matrix; face_vertex_matrix.resize(numOfFaces, numOfVertices);
	std::vector<T> coeffs; coeffs.clear();

	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		int firstIdx = f_it->idx();
		for (auto fv_it = mesh.cfv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			int secIdx = fv_it->idx();
			coeffs.push_back(T(firstIdx, secIdx, 1.0));
		}
	}
	face_vertex_matrix.setFromTriplets(coeffs.begin(), coeffs.end());

	return face_vertex_matrix;
}

SpMat buildAdjacentMatrix_FF(const TriMesh &mesh)
{
	int numOfFaces = static_cast<int>(mesh.n_faces());
	SpMat adjMat; adjMat.resize(numOfFaces, numOfFaces);
	std::vector<T> coeffs; coeffs.clear();

	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		int i = f_it->idx();
		for (auto ff_it = mesh.cff_iter(*f_it); ff_it.is_valid(); ++ff_it)
		{
			int j = ff_it->idx();
			coeffs.push_back(T(i, j, 1));
		}
	}

	adjMat.setFromTriplets(coeffs.begin(), coeffs.end());
	return adjMat;
}

SpMat buildAdjacentMatrix_VE(const TriMesh &mesh)
{
	const int numOfVertices = static_cast<int>(mesh.n_vertices());
	const int numOfEdges = static_cast<int>(mesh.n_edges());
	SpMat adjMat; adjMat.resize(numOfVertices, numOfEdges);
	std::vector<T> coeffs; coeffs.clear();

	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		int i = v_it->idx();
		for (auto ve_it = mesh.cve_iter(*v_it); ve_it.is_valid(); ++ve_it)
		{
			int j = ve_it->idx();
			coeffs.push_back(T(i, j, 1));
		}
	}

	adjMat.setFromTriplets(coeffs.begin(), coeffs.end());
	return adjMat;
}


SpMat buildAdjacentMatrix_VF(const TriMesh &mesh)
{
	// The commented code is used for testing.
	/*const int numOfVertices = static_cast<int>(mesh.n_vertices());
	int numOfFaces = static_cast<int>(mesh.n_faces());
	SpMat adjMat; adjMat.resize(numOfVertices, numOfFaces);
	std::vector<T> coeffs; coeffs.clear();

	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		int i = v_it->idx();
		for (auto vf_it = mesh.cvf_iter(*v_it); vf_it.is_valid(); ++vf_it)
		{
			int j = vf_it->idx();
			coeffs.push_back(T(i, j, 1));
		}
	}
	adjMat.setFromTriplets(coeffs.begin(), coeffs.end());
	return adjMat;*/

	return buildAdjacentMatrix_FV(mesh).transpose();
}


SpMat buildAdjacentMatrix_EE(const TriMesh &mesh)
{
	const int numOfEdges = static_cast<int>(mesh.n_edges());
	SpMat adjMat; adjMat.resize(numOfEdges, numOfEdges);
	std::vector<T> coeffs; coeffs.clear();

	for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		int i = e_it->idx();
		std::vector<OpenMesh::VertexHandle> vertices; vertices.clear();
		OpenMesh::VertexHandle v1 = mesh.to_vertex_handle(mesh.halfedge_handle(*e_it, 0));
		OpenMesh::VertexHandle v2 = mesh.from_vertex_handle(mesh.halfedge_handle(*e_it, 0));
		vertices.push_back(v1); vertices.push_back(v2);

		for (auto it = vertices.begin(); it != vertices.end(); ++it)
		{
			for (auto vh_it = mesh.cvoh_iter(*it); vh_it.is_valid(); ++vh_it)
			{
				OpenMesh::EdgeHandle edge = mesh.edge_handle(*vh_it);
				int j = edge.idx();
				if (i != j)
				{
					coeffs.push_back(T(i, j, 1.0));
				}
				
			}
		}
	}

	adjMat.setFromTriplets(coeffs.begin(), coeffs.end());
	return adjMat;
}


SpMat buildAdjacentMatrix_EV(const TriMesh &mesh)
{
	/*const int numOfEdges = static_cast<int>(mesh.n_edges());
	const int numOfVertices = static_cast<int>(mesh.n_vertices());
	SpMat adjMat; adjMat.resize(numOfEdges, numOfVertices);
	std::vector<T> coeffs; coeffs.clear();

	for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		int i = e_it->idx();
		OpenMesh::VertexHandle v1 = mesh.to_vertex_handle(mesh.halfedge_handle(*e_it, 0));
		OpenMesh::VertexHandle v2 = mesh.from_vertex_handle(mesh.halfedge_handle(*e_it, 0));
		int v1Idx = v1.idx();
		int v2Idx = v2.idx();

		coeffs.push_back(T(i, v1Idx, 1.0));
		coeffs.push_back(T(i, v2Idx, 1.0));
	}
	adjMat.setFromTriplets(coeffs.begin(), coeffs.end());
	return adjMat;*/

	return buildAdjacentMatrix_VE(mesh).transpose();
}

SpMat buildAdjacentMatrix_EF(const TriMesh &mesh)
{
	const int numOfEdges = static_cast<int>(mesh.n_edges());
	const int numOfFaces = static_cast<int>(mesh.n_faces());
	SpMat adjMat; adjMat.resize(numOfEdges, numOfFaces);
	std::vector<T> coeffs; coeffs.clear();

	for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		int i = e_it->idx();
		
		OpenMesh::HalfedgeHandle hf1 = mesh.halfedge_handle(*e_it, 0);
		OpenMesh::FaceHandle f1 = mesh.face_handle(hf1);
		int f1Idx = f1.idx();
		coeffs.push_back(T(i, f1Idx, 1.0));

		if (!mesh.is_boundary(*e_it))	// boundry has only one halfedge
		{
			OpenMesh::HalfedgeHandle hf2 = mesh.halfedge_handle(*e_it, 1);
			OpenMesh::FaceHandle f2 = mesh.face_handle(hf2);
			int f2Idx = f2.idx();
			coeffs.push_back(T(i, f2Idx, 1.0));
		}
		
	}
	adjMat.setFromTriplets(coeffs.begin(), coeffs.end());
	return adjMat;
}

SpMat buildAdjacentMatrix_FE(const TriMesh &mesh)
{
	/*const int numOfEdges = static_cast<int>(mesh.n_edges());
	const int numOfFaces = static_cast<int>(mesh.n_faces());
	SpMat adjMat; adjMat.resize(numOfFaces, numOfEdges);
	std::vector<T> coeffs; coeffs.clear();

	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		int i = f_it->idx();
		for (auto fe_it = mesh.cfe_iter(*f_it); fe_it.is_valid(); ++fe_it)
		{
			int j = fe_it->idx();
			coeffs.push_back(T(i, j, 1.0));
		}
	}

	adjMat.setFromTriplets(coeffs.begin(), coeffs.end());
	return adjMat;*/

	return buildAdjacentMatrix_EF(mesh).transpose();

}