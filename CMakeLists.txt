cmake_minimum_required(VERSION 3.24)

project(GMD VERSION 0.1.0 LANGUAGES CXX)

option(GMD_ENABLE_CUDA "Enable CUDA language support" OFF)
option(GMD_BUILD_PYTHON "Build Python bindings" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CompilerOptions)

if(GMD_ENABLE_CUDA)
    enable_language(CUDA)
    include(CUDAOptions)
endif()

add_library(gmd_headers INTERFACE)
target_include_directories(gmd_headers
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

gmd_apply_default_warnings(gmd_headers)

add_executable(gmd app/gmd_main.cpp)
target_link_libraries(gmd PRIVATE gmd_headers)
