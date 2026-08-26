cmake_minimum_required(VERSION 3.24)

project(GMD VERSION 2.4.0 LANGUAGES CXX)

include(CTest)

option(GMD_ENABLE_CUDA "Enable CUDA language support" OFF)
option(GMD_BUILD_PYTHON "Build Python bindings" OFF)
option(GMD_ENABLE_TORCH "Enable TorchScript ML force provider (requires LibTorch)" OFF)
option(GMD_ENABLE_MPI "Enable MPI parallel support" OFF)
option(GMD_ENABLE_DISTRIBUTED_PME
    "Enable optional distributed PME/FFT support (requires MPI and FFTW-MPI)"
    OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CompilerOptions)

if(GMD_ENABLE_CUDA)
    enable_language(CUDA)
    include(CUDAOptions)
endif()

if(GMD_ENABLE_TORCH)
    find_package(Torch REQUIRED)
endif()

if(GMD_ENABLE_MPI)
    find_package(MPI REQUIRED)
    include(MPIOptions)
endif()

if(GMD_ENABLE_DISTRIBUTED_PME)
    if(NOT GMD_ENABLE_MPI)
        message(FATAL_ERROR "GMD_ENABLE_DISTRIBUTED_PME requires GMD_ENABLE_MPI=ON")
    endif()
    find_path(FFTW3_MPI_INCLUDE_DIR fftw3-mpi.h)
    find_library(FFTW3_MPI_LIBRARY NAMES fftw3_mpi)
    find_library(FFTW3_LIBRARY NAMES fftw3)
    if(NOT FFTW3_MPI_INCLUDE_DIR OR NOT FFTW3_MPI_LIBRARY OR NOT FFTW3_LIBRARY)
        message(FATAL_ERROR
            "GMD_ENABLE_DISTRIBUTED_PME=ON requires FFTW-MPI. "
            "Could not find fftw3-mpi.h, libfftw3_mpi, and libfftw3. "
            "Set FFTW3_MPI_INCLUDE_DIR, FFTW3_MPI_LIBRARY, and FFTW3_LIBRARY "
            "or disable GMD_ENABLE_DISTRIBUTED_PME.")
    endif()
endif()

add_library(gmd_headers INTERFACE)
target_include_directories(gmd_headers
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

gmd_apply_default_warnings(gmd_headers)

add_library(gmd_core STATIC
    src/core/runtime_context.cpp
    src/core/simulation.cpp
    src/force/bonded_force_provider.cpp
    src/force/classical_force_provider.cpp
    src/force/ewald_force_provider.cpp
    src/force/ml_force_provider.cpp
    src/force/model_runtime_adapter.cpp
    src/force/pme_force_provider.cpp
    src/integrator/constraint_solver.cpp
    src/integrator/thermostat.cpp
    src/integrator/velocity_rescaling_thermostat.cpp
    src/integrator/nose_hoover_thermostat.cpp
    src/integrator/berendsen_barostat.cpp
    src/integrator/mc_barostat.cpp
    src/integrator/velocity_verlet_integrator.cpp
    src/io/checkpoint.cpp
    src/io/config_loader.cpp
    src/io/trajectory_writer.cpp
    src/parallel/domain_decomposition.cpp
    src/parallel/mpi_communicator.cpp
    src/parallel/mpi_environment.cpp
    src/parallel/pme_parallel.cpp
    src/system/initializer.cpp
    src/system/minimum_image.cpp
    src/system/periodic_boundary.cpp
    src/system/verlet_neighbor_builder.cpp
)
target_link_libraries(gmd_core PUBLIC gmd_headers)

if(GMD_ENABLE_TORCH)
    target_sources(gmd_core PRIVATE src/force/torchscript_adapter.cpp)
    target_link_libraries(gmd_core PUBLIC ${TORCH_LIBRARIES})
    target_compile_definitions(gmd_core PUBLIC GMD_ENABLE_TORCH)
endif()

if(GMD_ENABLE_MPI)
    gmd_configure_mpi_target(gmd_core)
endif()

if(GMD_ENABLE_DISTRIBUTED_PME)
    target_compile_definitions(gmd_core PUBLIC GMD_ENABLE_DISTRIBUTED_PME)
    target_include_directories(gmd_core PUBLIC ${FFTW3_MPI_INCLUDE_DIR})
    target_link_libraries(gmd_core PUBLIC ${FFTW3_MPI_LIBRARY} ${FFTW3_LIBRARY})
endif()

gmd_apply_default_warnings(gmd_core)

add_executable(gmd app/gmd_main.cpp)
target_link_libraries(gmd PRIVATE gmd_core)

add_executable(gmd_validate app/gmd_validate_main.cpp)
target_link_libraries(gmd_validate PRIVATE gmd_core)

if(BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)

    add_executable(gmd_special_pair_tests tests/special_pair_tests.cpp)
    target_link_libraries(gmd_special_pair_tests PRIVATE gmd_core)
    add_executable(gmd_constraint_solver_tests tests/constraint_solver_tests.cpp)
    target_link_libraries(gmd_constraint_solver_tests PRIVATE gmd_core)
    add_executable(gmd_checkpoint_restart_tests tests/checkpoint_restart_tests.cpp)
    target_link_libraries(gmd_checkpoint_restart_tests PRIVATE gmd_core)
    add_executable(gmd_pme_mode_tests tests/pme_mode_tests.cpp)
    target_link_libraries(gmd_pme_mode_tests PRIVATE gmd_core)
    add_executable(gmd_degrees_of_freedom_tests tests/degrees_of_freedom_tests.cpp)
    target_link_libraries(gmd_degrees_of_freedom_tests PRIVATE gmd_core)
    add_executable(gmd_box_image_flag_tests tests/box_image_flag_tests.cpp)
    target_link_libraries(gmd_box_image_flag_tests PRIVATE gmd_core)
    add_executable(gmd_barostat_force_refresh_tests tests/barostat_force_refresh_tests.cpp)
    target_link_libraries(gmd_barostat_force_refresh_tests PRIVATE gmd_core)
    add_executable(gmd_virial_finite_difference_tests tests/virial_finite_difference_tests.cpp)
    target_link_libraries(gmd_virial_finite_difference_tests PRIVATE gmd_core)
    add_executable(gmd_nose_hoover_restart_tests tests/nose_hoover_restart_tests.cpp)
    target_link_libraries(gmd_nose_hoover_restart_tests PRIVATE gmd_core)
    add_executable(gmd_constraint_validation_tests tests/constraint_validation_tests.cpp)
    target_link_libraries(gmd_constraint_validation_tests PRIVATE gmd_core)
    add_executable(gmd_bonded_force_gradient_tests tests/bonded_force_gradient_tests.cpp)
    target_link_libraries(gmd_bonded_force_gradient_tests PRIVATE gmd_core)
    add_executable(gmd_constraint_virial_tests tests/constraint_virial_tests.cpp)
    target_link_libraries(gmd_constraint_virial_tests PRIVATE gmd_core)

    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests)
    configure_file(tests/smoke_lj.xyz         ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_xyz.in       COPYONLY)
    configure_file(tests/smoke_lj.run         ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_run.in       COPYONLY)
    configure_file(tests/smoke_mpi_lj.xyz     ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_lj.xyz   COPYONLY)
    configure_file(tests/smoke_mpi_ewald.xyz  ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_ewald.xyz COPYONLY)
    configure_file(tests/smoke_mpi_molecular.xyz ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_molecular.xyz COPYONLY)
    configure_file(tests/smoke_mpi_lj.run     ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_lj.run   COPYONLY)
    configure_file(tests/smoke_mpi_ewald_consistency.run ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_ewald_consistency.run COPYONLY)
    configure_file(tests/smoke_mpi_pme_consistency.run ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_pme_consistency.run COPYONLY)
 configure_file(tests/smoke_mpi_periodic_lj.xyz ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_periodic_lj.xyz COPYONLY)
 configure_file(tests/smoke_mpi_periodic_lj.run ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_periodic_lj.run COPYONLY)
 configure_file(tests/smoke_mpi_boundary_migration.xyz ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_boundary_migration.xyz COPYONLY)
 configure_file(tests/smoke_mpi_boundary_migration.run ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_boundary_migration.run COPYONLY)
    configure_file(tests/smoke_ewald.xyz      ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_ewald.xyz    COPYONLY)
    configure_file(tests/smoke_ewald.run      ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_ewald.run    COPYONLY)
    configure_file(tests/smoke_pme.run        ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_pme.run      COPYONLY)
    configure_file(tests/smoke_mc_barostat.xyz ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mc_barostat.xyz COPYONLY)
    configure_file(tests/smoke_mc_barostat.run ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mc_barostat.run COPYONLY)
    configure_file(tests/smoke_molecular.xyz  ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_molecular.xyz COPYONLY)
    configure_file(tests/smoke_molecular.ff   ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_molecular.ff  COPYONLY)
    configure_file(tests/smoke_molecular.top  ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_molecular.top COPYONLY)
    configure_file(tests/smoke_molecular.run  ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_molecular.run COPYONLY)

    add_test(
        NAME gmd_smoke_inline_lj
        COMMAND gmd
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_xyz.in
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_run.in
    )
    set_tests_properties(gmd_smoke_inline_lj PROPERTIES LABELS "smoke;integration")

    add_test(
        NAME gmd_smoke_ewald
        COMMAND gmd
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_ewald.xyz
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_ewald.run
    )
    set_tests_properties(gmd_smoke_ewald PROPERTIES LABELS "smoke;integration")

    add_test(
        NAME gmd_smoke_pme
        COMMAND gmd
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_ewald.xyz
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_pme.run
    )
    set_tests_properties(gmd_smoke_pme PROPERTIES LABELS "smoke;integration")

    add_test(
        NAME gmd_smoke_mc_barostat
        COMMAND gmd
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mc_barostat.xyz
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mc_barostat.run
    )
    set_tests_properties(gmd_smoke_mc_barostat PROPERTIES LABELS "smoke;integration")

    # Molecular FF smoke: run from the binary test dir so relative paths in .run resolve.
    add_test(
        NAME gmd_smoke_molecular
        COMMAND gmd
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_molecular.xyz
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_molecular.run
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_molecular.ff
                ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_molecular.top
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(gmd_smoke_molecular PROPERTIES LABELS "smoke;integration")

    add_test(NAME gmd_special_pair_rules COMMAND gmd_special_pair_tests)
    set_tests_properties(gmd_special_pair_rules PROPERTIES LABELS "unit;force")
    add_test(NAME gmd_constraint_solver COMMAND gmd_constraint_solver_tests)
    set_tests_properties(gmd_constraint_solver PROPERTIES LABELS "unit;integrator")
    add_test(NAME gmd_checkpoint_restart COMMAND gmd_checkpoint_restart_tests)
    set_tests_properties(gmd_checkpoint_restart PROPERTIES LABELS "unit;io;restart")
    add_test(
        NAME gmd_restart_continuity_cli
        COMMAND Python3::Interpreter
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/restart_continuity.py
                --gmd $<TARGET_FILE:gmd>
                --work-root ${CMAKE_CURRENT_BINARY_DIR}/tests/restart_continuity_serial
                --np 1
                --position-tolerance 1e-10
                --velocity-tolerance 1e-10
                --energy-tolerance 1e-6
    )
    set_tests_properties(gmd_restart_continuity_cli
        PROPERTIES LABELS "integration;restart")
    add_test(
        NAME gmd_nose_hoover_restart_continuity_cli
        COMMAND Python3::Interpreter
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/nose_hoover_restart_continuity.py
                --gmd $<TARGET_FILE:gmd>
                --gmd-validate $<TARGET_FILE:gmd_validate>
                --work-root ${CMAKE_CURRENT_BINARY_DIR}/tests/nh_restart_continuity_serial
                --np 1
    )
    set_tests_properties(gmd_nose_hoover_restart_continuity_cli
        PROPERTIES LABELS "integration;restart;integrator")
    add_test(NAME gmd_pme_modes COMMAND gmd_pme_mode_tests)
    add_test(NAME gmd_degrees_of_freedom COMMAND gmd_degrees_of_freedom_tests)
    set_tests_properties(gmd_degrees_of_freedom PROPERTIES LABELS "unit;integrator")
    add_test(NAME gmd_box_image_flags COMMAND gmd_box_image_flag_tests)
    set_tests_properties(gmd_box_image_flags PROPERTIES LABELS "unit;force")
    add_test(NAME gmd_barostat_force_refresh COMMAND gmd_barostat_force_refresh_tests)
    set_tests_properties(gmd_barostat_force_refresh PROPERTIES LABELS "unit;integrator")
    add_test(NAME gmd_virial_finite_difference COMMAND gmd_virial_finite_difference_tests)
    set_tests_properties(gmd_virial_finite_difference PROPERTIES LABELS "unit;force;validation")
    add_test(NAME gmd_nose_hoover_restart COMMAND gmd_nose_hoover_restart_tests)
    set_tests_properties(gmd_nose_hoover_restart PROPERTIES LABELS "unit;integrator;restart")
    add_test(NAME gmd_constraint_validation COMMAND gmd_constraint_validation_tests)
    set_tests_properties(gmd_constraint_validation PROPERTIES LABELS "unit;integrator")
    add_test(NAME gmd_bonded_force_gradient COMMAND gmd_bonded_force_gradient_tests)
    set_tests_properties(gmd_bonded_force_gradient PROPERTIES LABELS "unit;force;validation")
    add_test(NAME gmd_constraint_virial COMMAND gmd_constraint_virial_tests)
    set_tests_properties(gmd_constraint_virial PROPERTIES LABELS "unit;integrator;force")
    set_tests_properties(gmd_pme_modes PROPERTIES LABELS "unit;force;pme")

    add_test(
        NAME gmd_validation_static_lj_cluster
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/validation/run_validation.py
                --case static_lj_cluster
                --gmd $<TARGET_FILE:gmd>
                --gmd-validate $<TARGET_FILE:gmd_validate>
                --work-root ${CMAKE_CURRENT_BINARY_DIR}/validation
    )
    set_tests_properties(gmd_validation_static_lj_cluster
        PROPERTIES LABELS "validation;unit")

    add_test(
        NAME gmd_validation_static_special_pairs
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/validation/run_validation.py
                --case static_special_pairs
                --gmd $<TARGET_FILE:gmd>
                --gmd-validate $<TARGET_FILE:gmd_validate>
                --work-root ${CMAKE_CURRENT_BINARY_DIR}/validation
    )
    set_tests_properties(gmd_validation_static_special_pairs
        PROPERTIES LABELS "validation;unit")

    add_test(
        NAME gmd_validation_static_coulomb
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/validation/run_validation.py
                --case static_coulomb
                --gmd $<TARGET_FILE:gmd>
                --gmd-validate $<TARGET_FILE:gmd_validate>
                --work-root ${CMAKE_CURRENT_BINARY_DIR}/validation
    )
    set_tests_properties(gmd_validation_static_coulomb
        PROPERTIES LABELS "validation;unit")

    add_custom_target(gmd_validation_short
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/validation/run_validation.py
                --case short
                --gmd $<TARGET_FILE:gmd>
                --gmd-validate $<TARGET_FILE:gmd_validate>
                --work-root ${CMAKE_CURRENT_BINARY_DIR}/validation
        DEPENDS gmd gmd_validate
        USES_TERMINAL
    )

    add_custom_target(gmd_validation_long
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/validation/run_validation.py
                --case long
                --gmd $<TARGET_FILE:gmd>
                --gmd-validate $<TARGET_FILE:gmd_validate>
                --work-root ${CMAKE_CURRENT_BINARY_DIR}/validation
        DEPENDS gmd gmd_validate
        USES_TERMINAL
    )

    if(GMD_ENABLE_MPI)
        add_executable(gmd_mpi_periodic_1d_tests tests/mpi_periodic_1d.cpp)
        target_link_libraries(gmd_mpi_periodic_1d_tests PRIVATE gmd_core)
        add_executable(gmd_mpi_domain_decomposition_3d_tests
            tests/mpi_domain_decomposition_3d.cpp)
        target_link_libraries(gmd_mpi_domain_decomposition_3d_tests PRIVATE gmd_core)
        add_executable(gmd_mpi_ghost_exchange_3d_tests
            tests/mpi_ghost_exchange_3d.cpp)
        target_link_libraries(gmd_mpi_ghost_exchange_3d_tests PRIVATE gmd_core)
        add_executable(gmd_mpi_special_pair_tests tests/mpi_special_pair.cpp)
        target_link_libraries(gmd_mpi_special_pair_tests PRIVATE gmd_core)
        add_executable(gmd_mpi_constraint_solver_tests tests/mpi_constraint_solver.cpp)
        target_link_libraries(gmd_mpi_constraint_solver_tests PRIVATE gmd_core)
        add_executable(gmd_mpi_checkpoint_restart_tests tests/mpi_checkpoint_restart.cpp)
        target_link_libraries(gmd_mpi_checkpoint_restart_tests PRIVATE gmd_core)
        add_executable(gmd_mpi_neighbor_list_tests tests/mpi_neighbor_list.cpp)
        target_link_libraries(gmd_mpi_neighbor_list_tests PRIVATE gmd_core)
        add_executable(gmd_mpi_dof_virial_tests tests/mpi_dof_virial.cpp)
        target_link_libraries(gmd_mpi_dof_virial_tests PRIVATE gmd_core)
        add_executable(gmd_mpi_constraint_virial_tests tests/mpi_constraint_virial.cpp)
        target_link_libraries(gmd_mpi_constraint_virial_tests PRIVATE gmd_core)

        # Same fixture at 1, 2 and 4 ranks: the constraint virial must come out
        # identical, which is what catches a dropped or rank-multiplied value.
        foreach(constraint_virial_np 1 2 4)
            add_test(
                NAME gmd_mpi_constraint_virial_${constraint_virial_np}proc
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${constraint_virial_np}
                        $<TARGET_FILE:gmd_mpi_constraint_virial_tests>
            )
            set_tests_properties(gmd_mpi_constraint_virial_${constraint_virial_np}proc
                PROPERTIES
                    LABELS "mpi;parallel;integrator;force"
                    PROCESSORS ${constraint_virial_np})
        endforeach()

        add_test(
            NAME gmd_mpi_degrees_of_freedom_4proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:gmd_mpi_dof_virial_tests> dof
        )
        set_tests_properties(gmd_mpi_degrees_of_freedom_4proc
            PROPERTIES
                LABELS "mpi;parallel;integrator"
                PROCESSORS 4)

        add_test(
            NAME gmd_mpi_virial_reduction_4proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:gmd_mpi_dof_virial_tests> virial
        )
        set_tests_properties(gmd_mpi_virial_reduction_4proc
            PROPERTIES
                LABELS "mpi;parallel;force"
                PROCESSORS 4)

        add_test(
            NAME gmd_mpi_neighbor_rebuild_collective_4proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:gmd_mpi_neighbor_list_tests> rebuild
        )
        set_tests_properties(gmd_mpi_neighbor_rebuild_collective_4proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 4)

        add_test(
            NAME gmd_mpi_ghost_image_flags_4proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:gmd_mpi_neighbor_list_tests> image
        )
        set_tests_properties(gmd_mpi_ghost_image_flags_4proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 4)

        add_test(
            NAME gmd_mpi_periodic_1d_ghost_exchange_2proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                    $<TARGET_FILE:gmd_mpi_periodic_1d_tests> ghost
        )
        set_tests_properties(gmd_mpi_periodic_1d_ghost_exchange_2proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 2)

        add_test(
            NAME gmd_mpi_periodic_1d_force_consistency_2proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                    $<TARGET_FILE:gmd_mpi_periodic_1d_tests> force
        )
        set_tests_properties(gmd_mpi_periodic_1d_force_consistency_2proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 2)

        add_test(
            NAME gmd_mpi_periodic_1d_migration_2proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                    $<TARGET_FILE:gmd_mpi_periodic_1d_tests> migration
        )
        set_tests_properties(gmd_mpi_periodic_1d_migration_2proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 2)

        add_test(
            NAME gmd_mpi_special_pair_consistency_2proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                    $<TARGET_FILE:gmd_mpi_special_pair_tests>
        )
        set_tests_properties(gmd_mpi_special_pair_consistency_2proc
            PROPERTIES
                LABELS "mpi;parallel;force"
                PROCESSORS 2)

        add_test(
            NAME gmd_mpi_constraint_solver_2proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                    $<TARGET_FILE:gmd_mpi_constraint_solver_tests>
        )
        set_tests_properties(gmd_mpi_constraint_solver_2proc
            PROPERTIES
                LABELS "mpi;parallel;integrator"
                PROCESSORS 2)

        add_test(
            NAME gmd_mpi_checkpoint_restart_2proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                    $<TARGET_FILE:gmd_mpi_checkpoint_restart_tests>
        )
        set_tests_properties(gmd_mpi_checkpoint_restart_2proc
            PROPERTIES
                LABELS "mpi;parallel;restart"
                PROCESSORS 2)

        add_test(
            NAME gmd_mpi_nose_hoover_restart_continuity_cli_2proc
            COMMAND Python3::Interpreter
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/nose_hoover_restart_continuity.py
                    --gmd $<TARGET_FILE:gmd>
                    --gmd-validate $<TARGET_FILE:gmd_validate>
                    --work-root ${CMAKE_CURRENT_BINARY_DIR}/tests/nh_restart_continuity_mpi_2
                    --np 2
                    --mpiexec ${MPIEXEC_EXECUTABLE}
                    --mpiexec-np-flag=${MPIEXEC_NUMPROC_FLAG}
        )
        set_tests_properties(gmd_mpi_nose_hoover_restart_continuity_cli_2proc
            PROPERTIES
                LABELS "mpi;parallel;integration;restart;integrator"
                PROCESSORS 2)

        add_test(
            NAME gmd_mpi_restart_continuity_cli_2proc
            COMMAND Python3::Interpreter
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/restart_continuity.py
                    --gmd $<TARGET_FILE:gmd>
                    --work-root ${CMAKE_CURRENT_BINARY_DIR}/tests/restart_continuity_mpi_2
                    --np 2
                    --mpiexec ${MPIEXEC_EXECUTABLE}
                    --mpiexec-np-flag=${MPIEXEC_NUMPROC_FLAG}
                    --position-tolerance 1e-9
                    --velocity-tolerance 1e-9
                    --energy-tolerance 1e-6
        )
        set_tests_properties(gmd_mpi_restart_continuity_cli_2proc
            PROPERTIES
                LABELS "mpi;parallel;restart;integration"
                PROCESSORS 2)

        add_test(
            NAME gmd_mpi_restart_continuity_cli_4proc
            COMMAND Python3::Interpreter
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/restart_continuity.py
                    --gmd $<TARGET_FILE:gmd>
                    --work-root ${CMAKE_CURRENT_BINARY_DIR}/tests/restart_continuity_mpi_4
                    --np 4
                    --mpiexec ${MPIEXEC_EXECUTABLE}
                    --mpiexec-np-flag=${MPIEXEC_NUMPROC_FLAG}
                    --position-tolerance 1e-9
                    --velocity-tolerance 1e-9
                    --energy-tolerance 1e-6
        )
        set_tests_properties(gmd_mpi_restart_continuity_cli_4proc
            PROPERTIES
                LABELS "mpi;parallel;restart;integration"
                PROCESSORS 4)

        add_test(
            NAME gmd_mpi_domain_decomposition_3d_grid_4proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:gmd_mpi_domain_decomposition_3d_tests> grid4
        )
        set_tests_properties(gmd_mpi_domain_decomposition_3d_grid_4proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 4)

        add_test(
            NAME gmd_mpi_domain_decomposition_3d_migration_8proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:gmd_mpi_domain_decomposition_3d_tests> migrate8
        )
        set_tests_properties(gmd_mpi_domain_decomposition_3d_migration_8proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 8)

        add_test(
            NAME gmd_mpi_ghost_exchange_3d_face_8proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:gmd_mpi_ghost_exchange_3d_tests> face
        )
        set_tests_properties(gmd_mpi_ghost_exchange_3d_face_8proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 8)

        add_test(
            NAME gmd_mpi_ghost_exchange_3d_edge_8proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:gmd_mpi_ghost_exchange_3d_tests> edge
        )
        set_tests_properties(gmd_mpi_ghost_exchange_3d_edge_8proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 8)

        add_test(
            NAME gmd_mpi_ghost_exchange_3d_corner_8proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:gmd_mpi_ghost_exchange_3d_tests> corner
        )
        set_tests_properties(gmd_mpi_ghost_exchange_3d_corner_8proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 8)

        add_test(
            NAME gmd_mpi_ghost_exchange_3d_periodic_8proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:gmd_mpi_ghost_exchange_3d_tests> periodic
        )
        set_tests_properties(gmd_mpi_ghost_exchange_3d_periodic_8proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 8)

        add_test(
            NAME gmd_mpi_reverse_force_3d_face_8proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:gmd_mpi_ghost_exchange_3d_tests> reverse-face
        )
        set_tests_properties(gmd_mpi_reverse_force_3d_face_8proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 8)

        add_test(
            NAME gmd_mpi_reverse_force_3d_edge_8proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:gmd_mpi_ghost_exchange_3d_tests> reverse-edge
        )
        set_tests_properties(gmd_mpi_reverse_force_3d_edge_8proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 8)

        add_test(
            NAME gmd_mpi_reverse_force_3d_corner_8proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:gmd_mpi_ghost_exchange_3d_tests> reverse-corner
        )
        set_tests_properties(gmd_mpi_reverse_force_3d_corner_8proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 8)

        add_test(
            NAME gmd_mpi_lj_3d_periodic_consistency_8proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:gmd_mpi_ghost_exchange_3d_tests> lj
        )
        set_tests_properties(gmd_mpi_lj_3d_periodic_consistency_8proc
            PROPERTIES
                LABELS "mpi;parallel"
                PROCESSORS 8)

        add_test(
            NAME gmd_smoke_mpi_lj_2proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                    $<TARGET_FILE:gmd>
                    ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_lj.xyz
                    ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_lj.run
        )
        set_tests_properties(gmd_smoke_mpi_lj_2proc
            PROPERTIES
                LABELS "smoke;integration;mpi"
                PROCESSORS 2)

        add_test(
            NAME gmd_smoke_mpi_ewald_2proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                    $<TARGET_FILE:gmd>
                    ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_ewald.xyz
                    ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_ewald.run
        )
        set_tests_properties(gmd_smoke_mpi_ewald_2proc
            PROPERTIES
                LABELS "smoke;integration;mpi"
                PROCESSORS 2)

        add_test(
            NAME gmd_smoke_mpi_pme_2proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                    $<TARGET_FILE:gmd>
                    ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_ewald.xyz
                    ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_pme.run
        )
        set_tests_properties(gmd_smoke_mpi_pme_2proc
            PROPERTIES
                LABELS "smoke;integration;mpi"
                PROCESSORS 2)

        add_test(
            NAME gmd_smoke_mpi_molecular_2proc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                    $<TARGET_FILE:gmd>
                    ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_molecular.xyz
                    ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_molecular.run
                    ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_molecular.ff
                    ${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_molecular.top
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
        set_tests_properties(gmd_smoke_mpi_molecular_2proc
            PROPERTIES
                LABELS "smoke;integration;mpi"
                PROCESSORS 2)

        add_executable(gmd_compare_energy_logs tests/compare_energy_logs.cpp)
        add_test(
            NAME gmd_smoke_mpi_lj_consistency
            COMMAND ${CMAKE_COMMAND}
                    -DGMD_EXECUTABLE=$<TARGET_FILE:gmd>
                    -DGMD_COMPARE_EXECUTABLE=$<TARGET_FILE:gmd_compare_energy_logs>
                    -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
                    -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
                    -DMPI_LJ_XYZ=${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_lj.xyz
                    -DMPI_LJ_RUN=${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_lj.run
                    -DMPI_LJ_WORK_ROOT=${CMAKE_CURRENT_BINARY_DIR}/tests/mpi_lj_consistency
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunMpiLJConsistency.cmake
        )
        set_tests_properties(gmd_smoke_mpi_lj_consistency
            PROPERTIES
                LABELS "smoke;integration;mpi"
                PROCESSORS 2)

        add_test(
            NAME gmd_smoke_mpi_ewald_consistency_8proc
            COMMAND ${CMAKE_COMMAND}
                    -DGMD_EXECUTABLE=$<TARGET_FILE:gmd>
                    -DGMD_COMPARE_EXECUTABLE=$<TARGET_FILE:gmd_compare_energy_logs>
                    -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
                    -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
                    -DMPI_CASE_NAME=ewald
                    -DMPI_PROCESS_COUNT=8
                    -DMPI_PROC_GRID_X=2
                    -DMPI_PROC_GRID_Y=2
                    -DMPI_PROC_GRID_Z=2
                    -DMPI_XYZ=${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_lj.xyz
                    -DMPI_RUN=${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_ewald_consistency.run
                    -DMPI_WORK_ROOT=${CMAKE_CURRENT_BINARY_DIR}/tests/mpi_consistency
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunMpiConsistency.cmake
        )
        set_tests_properties(gmd_smoke_mpi_ewald_consistency_8proc
            PROPERTIES
                LABELS "smoke;integration;mpi"
                PROCESSORS 8)

        add_test(
            NAME gmd_smoke_mpi_pme_consistency_8proc
            COMMAND ${CMAKE_COMMAND}
                    -DGMD_EXECUTABLE=$<TARGET_FILE:gmd>
                    -DGMD_COMPARE_EXECUTABLE=$<TARGET_FILE:gmd_compare_energy_logs>
                    -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
                    -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
                    -DMPI_CASE_NAME=pme
                    -DMPI_PROCESS_COUNT=8
                    -DMPI_PROC_GRID_X=2
                    -DMPI_PROC_GRID_Y=2
                    -DMPI_PROC_GRID_Z=2
                    -DMPI_XYZ=${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_lj.xyz
                    -DMPI_RUN=${CMAKE_CURRENT_BINARY_DIR}/tests/smoke_mpi_pme_consistency.run
                    -DMPI_WORK_ROOT=${CMAKE_CURRENT_BINARY_DIR}/tests/mpi_consistency
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunMpiConsistency.cmake
        )
        set_tests_properties(gmd_smoke_mpi_pme_consistency_8proc
            PROPERTIES
                LABELS "smoke;integration;mpi"
                PROCESSORS 8)
    endif()
endif()
