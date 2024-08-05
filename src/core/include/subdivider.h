#pragma once

#include "common_defines.h"


namespace core {

  class Subdivider
  {
    common::MyMesh::EdgeHandle findEdgeConnectingVertices(common::MyMesh& mesh, common::MyMesh::VertexHandle v1, common::MyMesh::VertexHandle v2) const;


    void execute(common::MyMesh& mesh, ChildrenParents& children_parents_map) const;
    
  public:

    void execute(common::MyMesh& mesh, size_t iteration_count = 1) const;
    
  private:

    void processQuadFace(common::MyMesh::FaceHandle& face, common::MyMesh& mesh, ChildrenParents& children_parents_map, const std::map<common::MyMesh::EdgeHandle, common::MyMesh::VertexHandle>& halfpoint_map) const;

    
    void processTriFace(common::MyMesh::FaceHandle& face, common::MyMesh& mesh, const std::map<common::MyMesh::EdgeHandle, common::MyMesh::VertexHandle>& halfpoint_map) const;

  };


  
}

