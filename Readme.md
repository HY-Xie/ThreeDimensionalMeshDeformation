## 说明
1. 全部使用相对路径，内含所有依赖(OpenMesh, Eigen), 直接复制到其它主机即可运行；
2. 本文件夹内包含《三维网格变形》一书中的部分学习代码，子项目标题对应该书章节；

---
### 5.1 邻接矩阵
1. 涉及到边操作的时候往往比较复杂，可以通过边对应的半边作为跳板再进行后续操作，比如`SpMat buildAdjacentMatrix_EE(const TriMesh &mesh)`和`SpMat buildAdjacentMatrix_EV(const TriMesh &mesh)`函数中：   
```
OpenMesh::VertexHandle v1 = mesh.to_vertex_handle(mesh.halfedge_handle(*e_it, 0));

OpenMesh::VertexHandle v2 = mesh.from_vertex_handle(mesh.halfedge_handle(*e_it, 0));
```
2. 处理边时，需考虑是否为边界，进而特殊处理：
```
if (!mesh.is_boundary(*e_it))	// boundry has only one halfedge
{
    OpenMesh::HalfedgeHandle hf2 =      mesh.halfedge_handle(*e_it, 1);
    OpenMesh::FaceHandle f2 = mesh.face_handle(hf2);
    int f2Idx = f2.idx();
    coeffs.push_back(T(i, f2Idx, 1.0));
}
```
3. X和Y邻接矩阵与Y和X邻接矩阵互为转置，其余相同。
---

5.2 组合拉普拉斯矩阵
1. 一般情况下，拉普拉斯矩阵不是对称矩阵(也有对称的情况)
2. 一般情况下，拉普拉斯矩阵每行的和为0(也有不为0的情况，如5.1.3归一化图形矩阵，每行和不为0)

