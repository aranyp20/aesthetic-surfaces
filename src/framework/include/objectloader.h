#pragma once

#include <OpenMesh/Core/IO/MeshIO.hh>
#include "common_defines.h"
#include <memory>
#include "organize.hpp"


namespace framework {

class ObjectLoader  
{
  std::vector<std::string> file_options;


    void modelToShow(common::MyMesh& mesh) const;
    void preprocessMesh(common::MyMesh& mesh) const;

  std::shared_ptr<common::MyMesh> loadFromFile(const std::string& path) const;
  
public:
  std::vector<std::string> loadFileOptions();


  std::shared_ptr<common::MyMesh> loadFromFile(const size_t file_option_index) const;


};
}
