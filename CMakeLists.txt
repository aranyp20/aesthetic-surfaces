cmake_minimum_required(VERSION 3.14)

project(a_fit1 VERSION 1.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)





find_package(OpenGL REQUIRED)



find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)




include_directories(${CMAKE_CURRENT_BINARY_DIR})



qt5_wrap_ui(UI_HEADERS src/visualization/qt/mainwindow.ui)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

find_package(glfw3 REQUIRED)
find_package(OpenMesh REQUIRED)




add_executable(${PROJECT_NAME}
src/visualization/qt/main.cpp
src/visualization/qt/mainwindow.cpp
src/visualization/qt/mainwindow.ui
src/visualization/src/canvas.cpp
src/visualization/src/camera.cpp
src/visualization/include/canvas.h
${UI_HEADERS}    
src/core/src/curvaturecalculator.cpp
src/framework/src/objectloader.cpp
src/common/include/organize.hpp
src/core/src/discretefairer.cpp
src/metrics/curvaturebased.cpp
src/core/src/subdivider.cpp
src/core/src/parametricsurface.cpp
src/core/src/parametriclogaesthetic1.cpp
src/core/src/parametric_spline_3d.cpp
src/core/src/logaesthetic_spline.cpp
src/core/src/planar_la.cpp
dependencies/lsq-plane.cc
)

target_include_directories(${PROJECT_NAME} PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/src/visualization/qt/
${CMAKE_CURRENT_SOURCE_DIR}/src/
${CMAKE_CURRENT_SOURCE_DIR}/src/core/include/
${CMAKE_CURRENT_SOURCE_DIR}/src/common/include/
${CMAKE_CURRENT_SOURCE_DIR}/src/visualization/include/
${CMAKE_CURRENT_SOURCE_DIR}/src/framework/include
${CMAKE_CURRENT_SOURCE_DIR}/src/metrics/
)



target_include_directories(${PROJECT_NAME} PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/dependencies/
)

link_directories(dependencies)


target_link_libraries(${PROJECT_NAME} PRIVATE OpenGL::GL Qt5::Widgets ${OPENMESH_LIBRARIES})
