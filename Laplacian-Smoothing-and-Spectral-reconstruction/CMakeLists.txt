cmake_minimum_required(VERSION 3.1)
project(Laplacian_processing_viewer)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use IMGUI"          ON)
#option(LIBIGL_WITH_PNG "Use PNG"          ON)

find_package(LIBIGL REQUIRED QUIET)

# ANN
include_directories(${PROJECT_SOURCE_DIR}/../spectra-0.7.0/include)

# Add your project files

set(SRC_FILES ${PROJECT_SOURCE_DIR}/main.cpp
              )
              

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(${PROJECT_NAME}_bin ${SRC_FILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw igl::opengl_glfw_imgui)
