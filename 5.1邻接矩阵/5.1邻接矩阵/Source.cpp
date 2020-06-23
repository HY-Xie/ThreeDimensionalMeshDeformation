#include <iostream>
#include "Vertex2Vertex.hpp"

int main(void) 
{
	TriMesh myMesh;
	if (!OpenMesh::IO::read_mesh(myMesh, "../../data/graph.off"))
	//if (!OpenMesh::IO::read_mesh(myMesh, "../../data/simple_square.obj"))
		std::cerr << "Mesh Load Failed." << std::endl;
	
	std::cout << "----------5.1.1 顶点和顶点的邻接矩阵------------" << std::endl;
	Eigen::MatrixXd v2vMatrix_Dense;
	v2vMatrix_Dense = buildAdjacentMatrixVV_Dense(myMesh);
	std::cout << "Vertex to Vertex Adj Matrix(Dense): " << std::endl << v2vMatrix_Dense << std::endl;

	std::cout << "----------5.1.1 顶点和顶点的邻接矩阵（稀疏）------------" << std::endl;
	SpMat v2vMatrix_Sparse;
	v2vMatrix_Sparse = buildAdjacentMatrix_VV(myMesh);
	std::cout << "Vertex to Vertex Adj Matrix(Sparse): " << std::endl << v2vMatrix_Sparse << std::endl;

	std::cout << "----------5.1.2 度数矩阵（稀疏）------------" << std::endl;
	SpMat degreeMatrix;
	degreeMatrix = buildMatrixDegree(myMesh);
	std::cout << "The Degree Matrix: " << std::endl << degreeMatrix << std::endl;

	std::cout << "----------5.1.3 面和顶点的邻接矩阵（稀疏）------------" << std::endl;
	SpMat f2vMatrix;
	f2vMatrix = buildAdjacentMatrix_FV(myMesh);
	std::cout << "Face to Vertex Adj Matrix(Sparse): " << std::endl << f2vMatrix << std::endl;

	std::cout << "----------5.1.4 面和面的邻接矩阵（稀疏）------------" << std::endl;
	SpMat f2fMatrix;
	f2fMatrix = buildAdjacentMatrix_FF(myMesh);
	std::cout << "Face to Face Adj Matrix(Sparse): " << std::endl << f2fMatrix << std::endl;

	std::cout << "----------5.1.5 顶点和边的邻接矩阵（稀疏）------------" << std::endl;
	SpMat v2eMatrix;
	v2eMatrix = buildAdjacentMatrix_VE(myMesh);
	std::cout << "Vertex to Edge Adj Matrix(Sparse): " << std::endl << v2eMatrix << std::endl;

	std::cout << "----------5.1.6 顶点和面的邻接矩阵（稀疏）------------" << std::endl;
	SpMat v2fMatrix;
	v2fMatrix = buildAdjacentMatrix_VF(myMesh);
	std::cout << "Vertex to Face Adj Matrix(Sparse): " << std::endl << v2fMatrix << std::endl;


	std::cout << "----------5.1.7 边和边的邻接矩阵（稀疏）------------" << std::endl;
	SpMat e2eMatrix;
	e2eMatrix = buildAdjacentMatrix_EE(myMesh);
	std::cout << "Edge to Edge Adj Matrix(Sparse): " << std::endl << e2eMatrix << std::endl;

	std::cout << "----------5.1.8 边和顶点的邻接矩阵（稀疏）------------" << std::endl;
	SpMat e2vMatrix;
	e2vMatrix = buildAdjacentMatrix_EV(myMesh);
	std::cout << "Edge to Vertex Adj Matrix(Sparse): " << std::endl << e2vMatrix << std::endl;

	std::cout << "----------5.1.9 边和面的邻接矩阵（稀疏）------------" << std::endl;
	SpMat e2fMatrix;
	e2fMatrix = buildAdjacentMatrix_EF(myMesh);
	std::cout << "Edge to Face Adj Matrix(Sparse): " << std::endl << e2fMatrix << std::endl;

	std::cout << "----------5.1.10 面和边的邻接矩阵（稀疏）------------" << std::endl;
	SpMat f2efMatrix;
	f2efMatrix = buildAdjacentMatrix_FE(myMesh);
	std::cout << "Face to Edge Adj Matrix(Sparse): " << std::endl << f2efMatrix << std::endl;

	/*SpMat m = f2efMatrix.transpose();
	for (int i = 0; i < e2fMatrix.rows(); ++i)
	{
		for (int j = 0; j < e2fMatrix.cols(); ++j)
		{
			if (e2fMatrix.coeffRef(i, j) != m.coeffRef(i,j))
			{
				std::cout << "False" << std::endl;
			}
		}
	}*/

	system("pause");
	return 0;
}

// 输出结果
/*
----------5.1.1 顶点和顶点的邻接矩阵------------
Vertex to Vertex Adj Matrix(Dense):
0 1 1 0 0 0
1 0 1 1 0 0
1 1 0 1 0 1
0 1 1 0 1 1
0 0 0 1 0 1
0 0 1 1 1 0
----------5.1.1 顶点和顶点的邻接矩阵（稀疏）------------
Vertex to Vertex Adj Matrix(Sparse):
Nonzero entries:
(1,1) (1,2) (1,0) (1,2) (1,3) (1,0) (1,1) (1,3) (1,5) (1,1) (1,2) (1,4) (1,5) (1,3) (1,5) (1,2) (1,3) (1,4)

Outer pointers:
0 2 5 9 13 15  $

0 1 1 0 0 0
1 0 1 1 0 0
1 1 0 1 0 1
0 1 1 0 1 1
0 0 0 1 0 1
0 0 1 1 1 0

----------5.1.2 度数矩阵（稀疏）------------
The Degree Matrix:
Nonzero entries:
(2,0) (3,1) (4,2) (4,3) (2,4) (3,5)

Outer pointers:
0 1 2 3 4 5  $

2 0 0 0 0 0
0 3 0 0 0 0
0 0 4 0 0 0
0 0 0 4 0 0
0 0 0 0 2 0
0 0 0 0 0 3

----------5.1.3 面和顶点的邻接矩阵（稀疏）------------
Face to Vertex Adj Matrix(Sparse):
Nonzero entries:
(1,0) (1,0) (1,1) (1,0) (1,1) (1,2) (1,1) (1,2) (1,3) (1,3) (1,2) (1,3)

Outer pointers:
0 1 3 6 9 10  $

1 1 1 0 0 0
0 1 1 1 0 0
0 0 1 1 0 1
0 0 0 1 1 1

----------5.1.4 面和面的邻接矩阵（稀疏）------------
Face to Face Adj Matrix(Sparse):
Nonzero entries:
(1,1) (1,0) (1,2) (1,1) (1,3) (1,2)

Outer pointers:
0 1 3 5  $

0 1 0 0
1 0 1 0
0 1 0 1
0 0 1 0

----------5.1.5 顶点和边的邻接矩阵（稀疏）------------
Vertex to Edge Adj Matrix(Sparse):
Nonzero entries:
(1,0) (1,1) (1,1) (1,2) (1,0) (1,2) (1,1) (1,3) (1,2) (1,3) (1,3) (1,5) (1,2) (1,5) (1,3) (1,4) (1,4) (1,5)

Outer pointers:
0 2 4 6 8 10 12 14 16  $

1 0 1 0 0 0 0 0 0
1 1 0 1 0 0 0 0 0
0 1 1 0 1 0 1 0 0
0 0 0 1 1 1 0 1 0
0 0 0 0 0 0 0 1 1
0 0 0 0 0 1 1 0 1

----------5.1.6 顶点和面的邻接矩阵（稀疏）------------
Vertex to Face Adj Matrix(Sparse):
Nonzero entries:
(1,0) (1,1) (1,2) (1,1) (1,2) (1,3) (1,2) (1,3) (1,5) (1,3) (1,4) (1,5)

Outer pointers:
0 3 6 9  $

1 0 0 0
1 1 0 0
1 1 1 0
0 1 1 1
0 0 0 1
0 0 1 1

----------5.1.7 边和边的邻接矩阵（稀疏）------------
Edge to Edge Adj Matrix(Sparse):
Nonzero entries:
(1,1) (1,2) (1,3) (1,0) (1,2) (1,3) (1,4) (1,6) (1,0) (1,1) (1,4) (1,6) (1,0) (1,1) (1,4)
(1,5) (1,7) (1,1) (1,2) (1,3) (1,5) (1,6) (1,7) (1,3) (1,4) (1,6) (1,7) (1,8) (1,1) (1,2)
(1,4) (1,5) (1,8) (1,3) (1,4) (1,5) (1,8) (1,5) (1,6) (1,7)

Outer pointers:
0 3 8 12 17 23 28 33 37  $

0 1 1 1 0 0 0 0 0
1 0 1 1 1 0 1 0 0
1 1 0 0 1 0 1 0 0
1 1 0 0 1 1 0 1 0
0 1 1 1 0 1 1 1 0
0 0 0 1 1 0 1 1 1
0 1 1 0 1 1 0 0 1
0 0 0 1 1 1 0 0 1
0 0 0 0 0 1 1 1 0

----------5.1.8 边和顶点的邻接矩阵（稀疏）------------
Edge to Vertex Adj Matrix(Sparse):
Nonzero entries:
(1,0) (1,2) (1,0) (1,1) (1,3) (1,1) (1,2) (1,4) (1,6) (1,3) (1,4) (1,5) (1,7) (1,7) (1,8) (1,5) (1,6) (1,8)

Outer pointers:
0 2 5 9 13 15  $

1 1 0 0 0 0
0 1 1 0 0 0
1 0 1 0 0 0
0 1 0 1 0 0
0 0 1 1 0 0
0 0 0 1 0 1
0 0 1 0 0 1
0 0 0 1 1 0
0 0 0 0 1 1

----------5.1.9 边和面的邻接矩阵（稀疏）------------
Edge to Face Adj Matrix(Sparse):
Nonzero entries:
(1,0) (1,1) (1,2) (1,1) (1,3) (1,4) (1,4) (1,5) (1,6) (1,5) (1,7) (1,8)

Outer pointers:
0 3 6 9  $

1 0 0 0
1 1 0 0
1 0 0 0
0 1 0 0
0 1 1 0
0 0 1 1
0 0 1 0
0 0 0 1
0 0 0 1

----------5.1.10 面和边的邻接矩阵（稀疏）------------
Face to Edge Adj Matrix(Sparse):
Nonzero entries:
(1,0) (1,0) (1,1) (1,0) (1,1) (1,1) (1,2) (1,2) (1,3) (1,2) (1,3) (1,3)

Outer pointers:
0 1 3 4 5 7 9 10 11  $

1 1 1 0 0 0 0 0 0
0 1 0 1 1 0 0 0 0
0 0 0 0 1 1 1 0 0
0 0 0 0 0 1 0 1 1

Press any key to continue . . .


*/