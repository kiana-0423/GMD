foreach(required_var
        GMD_EXECUTABLE
        GMD_COMPARE_EXECUTABLE
        MPIEXEC_EXECUTABLE
        MPIEXEC_NUMPROC_FLAG
        MPI_CASE_NAME
        MPI_PROCESS_COUNT
        MPI_XYZ
        MPI_RUN
        MPI_WORK_ROOT)
    if(NOT DEFINED ${required_var})
        message(FATAL_ERROR "RunMpiConsistency.cmake missing ${required_var}")
    endif()
endforeach()

if(NOT DEFINED MPI_PE_TOLERANCE)
    set(MPI_PE_TOLERANCE 1e-6)
endif()
if(NOT DEFINED MPI_TEMPERATURE_TOLERANCE)
    set(MPI_TEMPERATURE_TOLERANCE 1e-3)
endif()
if(NOT DEFINED MPI_DRIFT_TOLERANCE)
    set(MPI_DRIFT_TOLERANCE 1e-6)
endif()

set(serial_dir "${MPI_WORK_ROOT}/${MPI_CASE_NAME}/serial")
set(mpi_dir "${MPI_WORK_ROOT}/${MPI_CASE_NAME}/mpi")
file(REMOVE_RECURSE "${serial_dir}" "${mpi_dir}")
file(MAKE_DIRECTORY "${serial_dir}" "${mpi_dir}")

get_filename_component(xyz_name "${MPI_XYZ}" NAME)
get_filename_component(run_name "${MPI_RUN}" NAME)
file(COPY_FILE "${MPI_XYZ}" "${serial_dir}/${xyz_name}")
file(COPY_FILE "${MPI_RUN}" "${serial_dir}/${run_name}")
file(COPY_FILE "${MPI_XYZ}" "${mpi_dir}/${xyz_name}")
file(COPY_FILE "${MPI_RUN}" "${mpi_dir}/${run_name}")

execute_process(
    COMMAND "${GMD_EXECUTABLE}"
            "${serial_dir}/${xyz_name}"
            "${serial_dir}/${run_name}"
            --np 1
    RESULT_VARIABLE serial_result
    OUTPUT_VARIABLE serial_stdout
    ERROR_VARIABLE serial_stderr
)
if(NOT serial_result EQUAL 0)
    message(FATAL_ERROR
        "Serial ${MPI_CASE_NAME} reference run failed (${serial_result})\n"
        "stdout:\n${serial_stdout}\n"
        "stderr:\n${serial_stderr}")
endif()

set(mpi_grid_args)
if(DEFINED MPI_PROC_GRID_X AND DEFINED MPI_PROC_GRID_Y AND DEFINED MPI_PROC_GRID_Z)
    list(APPEND mpi_grid_args
         --proc-grid "${MPI_PROC_GRID_X}" "${MPI_PROC_GRID_Y}" "${MPI_PROC_GRID_Z}")
endif()

execute_process(
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${MPI_PROCESS_COUNT}"
            "${GMD_EXECUTABLE}"
            "${mpi_dir}/${xyz_name}"
            "${mpi_dir}/${run_name}"
            --np "${MPI_PROCESS_COUNT}"
            ${mpi_grid_args}
    RESULT_VARIABLE mpi_result
    OUTPUT_VARIABLE mpi_stdout
    ERROR_VARIABLE mpi_stderr
)
if(NOT mpi_result EQUAL 0)
    message(FATAL_ERROR
        "${MPI_PROCESS_COUNT}-rank MPI ${MPI_CASE_NAME} run failed (${mpi_result})\n"
        "stdout:\n${mpi_stdout}\n"
        "stderr:\n${mpi_stderr}")
endif()

execute_process(
    COMMAND "${GMD_COMPARE_EXECUTABLE}"
            "${serial_dir}/output.log"
            "${mpi_dir}/output.log"
            "${MPI_PE_TOLERANCE}"
            "${MPI_TEMPERATURE_TOLERANCE}"
            "${MPI_DRIFT_TOLERANCE}"
    RESULT_VARIABLE compare_result
    OUTPUT_VARIABLE compare_stdout
    ERROR_VARIABLE compare_stderr
)
if(NOT compare_result EQUAL 0)
    message(FATAL_ERROR
        "MPI ${MPI_CASE_NAME} numeric comparison failed (${compare_result})\n"
        "stdout:\n${compare_stdout}\n"
        "stderr:\n${compare_stderr}\n"
        "serial run:\n${serial_stdout}\n"
        "mpi run:\n${mpi_stdout}")
endif()
