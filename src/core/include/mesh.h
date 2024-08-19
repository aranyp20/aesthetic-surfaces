#pragma once

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>



namespace common {


struct MyTraits : public OpenMesh::DefaultTraits
{
  FaceAttributes(OpenMesh::Attributes::Status);
  EdgeAttributes(OpenMesh::Attributes::Status);
};
  typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> BaseMesh;



  class MyMesh : public BaseMesh
  {
  public:
    std::map<common::MyMesh::VertexHandle, std::array<common::MyMesh::VertexHandle,2>> children_parents_map;
    inline void reset()
    {
      
    }
    
  };
}

        typedef std::map<common::MyMesh::VertexHandle, std::array<common::MyMesh::VertexHandle,2>> ChildrenParents;

