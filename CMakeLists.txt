cmake_minimum_required(VERSION 3.24)

project(GMD VERSION 0.1.0 LANGUAGES CXX)

include(CTest)

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

add_library(gmd_core STATIC
    src/boundary/minimum_image.cpp
    src/boundary/periodic_boundary.cpp
    src/core/simulation.cpp
    src/force/bonded_force_provider.cpp
    src/force/classical_force_provider.cpp
    src/force/ewald_force_provider.cpp
    src/force/pme_force_provider.cpp
    src/integrator/thermostat.cpp
    src/integrator/velocity_rescaling_thermostat.cpp
    src/integrator/nose_hoover_thermostat.cpp
    src/integrator/berendsen_barostat.cpp
    src/integrator/mc_barostat.cpp
    src/integrator/velocity_verlet_integrator.cpp
    src/io/config_loader.cpp
    src/io/trajectory_writer.cpp
    src/ml/ml_force_provider.cpp
    src/ml/model_runtime_adapter.cpp
    src/neighbor/verlet_neighbor_builder.cpp
    src/system/initializer.cpp
)
target_link_libraries(gmd_core PUBLIC gmd_headers)

gmd_apply_default_warnings(gmd_core)

add_executable(gmd app/gmd_main.cpp)
target_link_libraries(gmd PRIVATE gmd_core)

if(BUILD_TESTING)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests)
    configure_file(tests/smoke_lj.xyz ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_xyz.in COPYONLY)
    configure_file(tests/smoke_lj.run ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_run.in COPYONLY)

    add_test(
        NAME gmd_smoke_inline_lj
        COMMAND gmd
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_xyz.in
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_run.in
    )
    set_tests_properties(gmd_smoke_inline_lj PROPERTIES LABELS "smoke;integration")
endif()
