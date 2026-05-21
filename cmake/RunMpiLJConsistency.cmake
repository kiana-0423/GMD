foreach(required_var
        GMD_EXECUTABLE
        GMD_COMPARE_EXECUTABLE
        MPIEXEC_EXECUTABLE
        MPIEXEC_NUMPROC_FLAG
        MPI_LJ_XYZ
        MPI_LJ_RUN
        MPI_LJ_WORK_ROOT)
    if(NOT DEFINED ${required_var})
        message(FATAL_ERROR "RunMpiLJConsistency.cmake missing ${required_var}")
    endif()
endforeach()

set(serial_dir "${MPI_LJ_WORK_ROOT}/serial")
set(mpi_dir "${MPI_LJ_WORK_ROOT}/mpi")
file(REMOVE_RECURSE "${serial_dir}" "${mpi_dir}")
file(MAKE_DIRECTORY "${serial_dir}" "${mpi_dir}")

file(COPY_FILE "${MPI_LJ_XYZ}" "${serial_dir}/smoke_mpi_lj.xyz")
file(COPY_FILE "${MPI_LJ_RUN}" "${serial_dir}/smoke_mpi_lj.run")
file(COPY_FILE "${MPI_LJ_XYZ}" "${mpi_dir}/smoke_mpi_lj.xyz")
file(COPY_FILE "${MPI_LJ_RUN}" "${mpi_dir}/smoke_mpi_lj.run")

execute_process(
    COMMAND "${GMD_EXECUTABLE}"
            "${serial_dir}/smoke_mpi_lj.xyz"
            "${serial_dir}/smoke_mpi_lj.run"
            --np 1
    RESULT_VARIABLE serial_result
    OUTPUT_VARIABLE serial_stdout
    ERROR_VARIABLE serial_stderr
)
if(NOT serial_result EQUAL 0)
    message(FATAL_ERROR
        "Serial LJ reference run failed (${serial_result})\n"
        "stdout:\n${serial_stdout}\n"
        "stderr:\n${serial_stderr}")
endif()

execute_process(
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
            "${GMD_EXECUTABLE}"
            "${mpi_dir}/smoke_mpi_lj.xyz"
            "${mpi_dir}/smoke_mpi_lj.run"
            --np 2
    RESULT_VARIABLE mpi_result
    OUTPUT_VARIABLE mpi_stdout
    ERROR_VARIABLE mpi_stderr
)
if(NOT mpi_result EQUAL 0)
    message(FATAL_ERROR
        "Two-rank MPI LJ run failed (${mpi_result})\n"
        "stdout:\n${mpi_stdout}\n"
        "stderr:\n${mpi_stderr}")
endif()

execute_process(
    COMMAND "${GMD_COMPARE_EXECUTABLE}"
            "${serial_dir}/output.log"
            "${mpi_dir}/output.log"
            1e-6
            1e-3
            1e-6
    RESULT_VARIABLE compare_result
    OUTPUT_VARIABLE compare_stdout
    ERROR_VARIABLE compare_stderr
)
if(NOT compare_result EQUAL 0)
    message(FATAL_ERROR
        "MPI LJ numeric comparison failed (${compare_result})\n"
        "stdout:\n${compare_stdout}\n"
        "stderr:\n${compare_stderr}\n"
        "serial run:\n${serial_stdout}\n"
        "mpi run:\n${mpi_stdout}")
endif()
