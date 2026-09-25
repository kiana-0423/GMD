add_executable(gmd_next_cuda_probe probe.cu)
target_link_libraries(gmd_next_cuda_probe PRIVATE CUDA::cudart)
set_target_properties(gmd_next_cuda_probe PROPERTIES
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_EXTENSIONS OFF
)

# P1.2 execution resources. Host C++20 against the CUDA Runtime API; kernels
# that use it are CUDA C++20 so they can share the contracts headers.
add_library(gmd_next_cuda_runtime STATIC
    runtime/buffer.cpp
    runtime/context.cpp)
add_library(GMDNext::cuda_runtime ALIAS gmd_next_cuda_runtime)
target_link_libraries(gmd_next_cuda_runtime PUBLIC GMDNext::contracts CUDA::cudart)
target_include_directories(gmd_next_cuda_runtime
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(gmd_next_cuda_runtime PROPERTIES CXX_EXTENSIONS OFF)
gmd_next_warnings(gmd_next_cuda_runtime)

if(GMD_NEXT_TEST_CUDA)
    add_test(NAME gmd_next.cuda_environment COMMAND gmd_next_cuda_probe)
    set_tests_properties(gmd_next.cuda_environment PROPERTIES LABELS "cuda;hardware")

    add_executable(gmd_next_cuda_runtime_tests tests/runtime_tests.cu)
    target_link_libraries(gmd_next_cuda_runtime_tests PRIVATE GMDNext::cuda_runtime)
    target_include_directories(gmd_next_cuda_runtime_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/tests)
    set_target_properties(gmd_next_cuda_runtime_tests PROPERTIES
        CUDA_STANDARD 20
        CUDA_STANDARD_REQUIRED ON
        CUDA_EXTENSIONS OFF)
    add_test(NAME gmd_next.cuda_runtime COMMAND gmd_next_cuda_runtime_tests)
    set_tests_properties(gmd_next.cuda_runtime PROPERTIES LABELS "cuda;hardware")
endif()
