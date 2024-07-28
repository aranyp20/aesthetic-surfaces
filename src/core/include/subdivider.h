#pragma once

#include "common_defines.h"


namespace core {

  class Subdivider
  {
    common::MyMesh::EdgeHandle findEdgeConnectingVertices(common::MyMesh& mesh, common::MyMesh::VertexHandle v1, common::MyMesh::VertexHandle v2) const;
    
  public:
    typedef std::map<common::MyMesh::VertexHandle, std::array<common::MyMesh::VertexHandle,2>> ChildrenParents;

    ChildrenParents execute(common::MyMesh& mesh, size_t iteration_count = 1) const;
    
  private:
    void execute(common::MyMesh& mesh, ChildrenParents& children_parents_map) const;
    

  };


  
}

