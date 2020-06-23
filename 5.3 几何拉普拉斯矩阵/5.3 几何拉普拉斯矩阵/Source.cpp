#include <iostream>
#include "GeometricLaplacian.hpp"

int main(void)
{
	TriMesh mesh; 
	if (OpenMesh::IO::read_mesh(mesh, "../../data/graph.off"))
	//if (OpenMesh::IO::read_mesh(mesh, "../../data/sphere_small.obj"))
	//if (OpenMesh::IO::read_mesh(mesh, "../../data/SPRING0001.obj"))
		std::cout << "Mesh Load Success" << std::endl;
	else {
		std::cout << "Mesh Load Fail." << std::endl;
		system("pause");
		return -1;
	}
	
	std::cout << "--------------------5.3.1 余切拉普拉斯--------------------" << std::endl;
	SpMat cotLaplacian = buildLaplacianCot(mesh);
	//std::cout << cotLaplacian << std::endl;
	laplacianValidate(cotLaplacian);

	std::cout << "--------------------5.3.2 面积余切拉普拉斯--------------------" << std::endl;
	SpMat cotAreaLaplacian = buildLaplacianCotArea(mesh);
	//std::cout << cotAreaLaplacian << std::endl;
	laplacianValidate(cotAreaLaplacian);

	std::cout << "--------------------5.3.3 中值拉普拉斯--------------------" << std::endl;
	SpMat meanValueAreaLaplacian = buildLaplacianMeanValue(mesh);
	std::cout << meanValueAreaLaplacian << std::endl;
	laplacianValidate(meanValueAreaLaplacian);
	





	system("pause");
	return 0;
}