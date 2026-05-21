include_guard(GLOBAL)

function(gmd_configure_mpi_target target_name)
    target_link_libraries(${target_name} PUBLIC MPI::MPI_CXX)
    target_compile_definitions(${target_name} PUBLIC GMD_ENABLE_MPI)
endfunction()
