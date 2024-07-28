#pragma once

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

namespace common {

struct MyTraits : public OpenMesh::DefaultTraits
{
    FaceAttributes(OpenMesh::Attributes::Status);
};
  //typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

  typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> MyMesh;
}
