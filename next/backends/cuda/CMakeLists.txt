add_executable(gmd_next_cuda_probe probe.cu)
target_link_libraries(gmd_next_cuda_probe PRIVATE CUDA::cudart)
set_target_properties(gmd_next_cuda_probe PROPERTIES
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_EXTENSIONS OFF
)

if(GMD_NEXT_TEST_CUDA)
    add_test(NAME gmd_next.cuda_environment COMMAND gmd_next_cuda_probe)
    set_tests_properties(gmd_next.cuda_environment PROPERTIES LABELS "cuda;hardware")
endif()
