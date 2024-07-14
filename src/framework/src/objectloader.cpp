#include "objectloader.h"
#include <iostream>
#include <filesystem>

#include "common_defines.h"
#include "discretefairer.h"

namespace framework {

void ObjectLoader::modelToShow(common::MyMesh& mesh) const
{
  if (mesh.vertices_empty()) {
    throw std::runtime_error("Cannot preprocess empty mesh.");
  }

  double max_x = std::numeric_limits<double>::min();
  double max_y = std::numeric_limits<double>::min();
  double max_z = std::numeric_limits<double>::min();
  double min_x = std::numeric_limits<double>::max();
  double min_y = std::numeric_limits<double>::max();
  double min_z = std::numeric_limits<double>::max();

  for (common::MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
  {
    const auto& vh = *v_it;
    const auto& vertexPos = mesh.point(vh);

    max_x = std::max(max_x, static_cast<double>(vertexPos[0]));
    max_y = std::max(max_y, static_cast<double>(vertexPos[1]));
    max_z = std::max(max_z, static_cast<double>(vertexPos[2]));

    min_x = std::min(min_x, static_cast<double>(vertexPos[0]));
    min_y = std::min(min_y, static_cast<double>(vertexPos[1]));
    min_z = std::min(min_z, static_cast<double>(vertexPos[2]));
  }

  Eigen::Matrix4d translate_to_origin = Eigen::Matrix4d::Identity();
  
  const Eigen::Vector3d original_center((max_x + min_x)/2, (max_y + min_y)/2, (max_z + min_z)/2);
  translate_to_origin.block<3,1>(0,3) = -original_center;

  double normalize_scaler = 2 / std::max(max_x - min_x, std::max(max_y - min_y, max_z - min_z));

  Eigen::Matrix4d scale_to_normalized = Eigen::Matrix4d::Identity();

  scale_to_normalized(0,0) = normalize_scaler;
  scale_to_normalized(1,1) = normalize_scaler;  
  scale_to_normalized(2,2) = normalize_scaler;  
 

  const Eigen::Matrix4d T_show_model = scale_to_normalized * translate_to_origin;


  for (common::MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
  {
    auto& vh = *v_it;
    auto& vertexPos = mesh.point(vh);

    vertexPos = common::toOpenMeshV3(common::transformAsPoint(T_show_model, common::toEigenV3(vertexPos)));
  }

}

  #include "curvaturecalculator.h"
void ObjectLoader::preprocessMesh(common::MyMesh& mesh) const
{
  //modelToShow(mesh); //TODOcheck

}

void aaaaa(common::MyMesh& mesh)
{

  common::MyMesh::VertexHandle v0 = mesh.add_vertex(common::MyMesh::Point(0, 0, 0)); // A
  common::MyMesh::VertexHandle v1 = mesh.add_vertex(common::MyMesh::Point(1, 0, 0)); // B
  common::MyMesh::VertexHandle v2 = mesh.add_vertex(common::MyMesh::Point(0.5, 1, 0)); // C

  common::MyMesh::FaceHandle f0 = mesh.add_face(v0, v1, v2);

}


std::shared_ptr<common::MyMesh> ObjectLoader::loadFromFile(const std::string& path) const
{  
  common::MyMesh input_mesh;

  if (!OpenMesh::IO::read_mesh(input_mesh, path))
  {
      std::cout << "Error: Cannot read mesh from file." << std::endl;
      return nullptr;
  }
  else {

    common::MyMesh mesh;

    aaaaa(mesh);

    //    preprocessMesh(input_mesh);

    return std::make_shared<common::MyMesh>(std::move(input_mesh));
  }

  return nullptr;
}



  std::vector<std::string> ObjectLoader::loadFileOptions()
  {
    std::vector<std::string> retval;
    
    for (const auto& entry : std::filesystem::directory_iterator("./")) {
      if (entry.path().extension() == ".obj" || entry.path().extension() == ".stl") {
        std::cout << entry.path() << std::endl;
	retval.push_back(entry.path());
      }
    }

    file_options = retval;
    return retval;
  }



std::shared_ptr<common::MyMesh> ObjectLoader::loadFromFile(const size_t file_option_index) const
{
  return loadFromFile(file_options[file_option_index]);
}


}
