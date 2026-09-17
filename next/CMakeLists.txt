cmake_minimum_required(VERSION 3.24)

# This project intentionally has no connection to the parent build graph.
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(FATAL_ERROR "GMD Next must be configured independently: cmake -S next -B next/build")
endif()
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "Use an out-of-source build, for example -B next/build")
endif()

project(GMDNext VERSION 0.1.0 LANGUAGES CXX)

option(GMD_NEXT_BUILD_TESTS "Build CPU mathematical contract checks" ON)
option(GMD_NEXT_ENABLE_CUDA "Build the CUDA Runtime/CUB environment probe" OFF)
option(GMD_NEXT_TEST_CUDA "Register the probe as a test; requires an NVIDIA GPU" OFF)

if(GMD_NEXT_TEST_CUDA AND NOT GMD_NEXT_ENABLE_CUDA)
    message(FATAL_ERROR "GMD_NEXT_TEST_CUDA requires GMD_NEXT_ENABLE_CUDA=ON")
endif()

function(gmd_next_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/W4>)
    else()
        target_compile_options(${target} PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>)
    endif()
endfunction()

add_library(gmd_next_reference STATIC src/reference.cpp)
add_library(GMDNext::reference ALIAS gmd_next_reference)
target_compile_features(gmd_next_reference PUBLIC cxx_std_20)
set_target_properties(gmd_next_reference PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(gmd_next_reference PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
gmd_next_warnings(gmd_next_reference)

add_executable(gmd_next_reference_demo app/reference_demo.cpp)
target_link_libraries(gmd_next_reference_demo PRIVATE GMDNext::reference)
gmd_next_warnings(gmd_next_reference_demo)

if(GMD_NEXT_BUILD_TESTS OR GMD_NEXT_TEST_CUDA)
    enable_testing()
endif()

if(GMD_NEXT_BUILD_TESTS)
    foreach(test_name IN ITEMS operators lj)
        add_executable(gmd_next_${test_name}_tests tests/${test_name}_tests.cpp)
        target_link_libraries(gmd_next_${test_name}_tests PRIVATE GMDNext::reference)
        gmd_next_warnings(gmd_next_${test_name}_tests)
        add_test(NAME gmd_next.${test_name} COMMAND gmd_next_${test_name}_tests)
    endforeach()
endif()

if(GMD_NEXT_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(NOT CMAKE_CUDA_COMPILER)
        message(FATAL_ERROR
            "CUDA was requested, but no CUDA compiler was found. Configure on an NVIDIA "
            "CUDA development host, or use GMD_NEXT_ENABLE_CUDA=OFF for CPU reference checks.")
    endif()
    enable_language(CUDA)
    find_package(CUDAToolkit 12.0 REQUIRED)
    add_subdirectory(backends/cuda)
endif()
