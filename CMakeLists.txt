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
    add_executable(gmd_constraint_independence_tests tests/constraint_independence_tests.cpp)
    target_link_libraries(gmd_constraint_independence_tests PRIVATE gmd_core)
    add_executable(gmd_pme_bspline_tests tests/pme_bspline_tests.cpp)
    target_link_libraries(gmd_pme_bspline_tests PRIVATE gmd_core)
    add_executable(gmd_electrostatic_constant_tests tests/electrostatic_constant_tests.cpp)
    target_link_libraries(gmd_electrostatic_constant_tests PRIVATE gmd_core)
    add_executable(gmd_boltzmann_constant_tests tests/boltzmann_constant_tests.cpp)
    target_link_libraries(gmd_boltzmann_constant_tests PRIVATE gmd_core)
    add_executable(gmd_pressure_unit_tests tests/pressure_unit_tests.cpp)
    target_link_libraries(gmd_pressure_unit_tests PRIVATE gmd_core)
    add_executable(gmd_pme_energy_force_tests tests/pme_energy_force_tests.cpp)
    target_link_libraries(gmd_pme_energy_force_tests PRIVATE gmd_core)
    add_executable(gmd_ml_virial_contract_tests tests/ml_virial_contract_tests.cpp)
    target_link_libraries(gmd_ml_virial_contract_tests PRIVATE gmd_core)

    # The TorchScript edge_shift test observes the tensor from inside a real
    # scripted model, so it can only exist when LibTorch is available. It builds
    # its own .pt fixtures at run time into the test working directory; nothing
    # is checked in and no Python is required.
    if(GMD_ENABLE_TORCH)
        add_executable(gmd_torchscript_edge_shift_tests
            tests/torchscript_edge_shift_tests.cpp)
        target_link_libraries(gmd_torchscript_edge_shift_tests PRIVATE gmd_core)
    endif()
    add_executable(gmd_virial_source_inventory_tests tests/virial_source_inventory_tests.cpp)
    target_link_libraries(gmd_virial_source_inventory_tests PRIVATE gmd_core)
    add_executable(gmd_pme_reciprocal_virial_tests tests/pme_reciprocal_virial_tests.cpp)
    target_link_libraries(gmd_pme_reciprocal_virial_tests PRIVATE gmd_core)

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
    add_test(NAME gmd_constraint_independence COMMAND gmd_constraint_independence_tests)
    set_tests_properties(gmd_constraint_independence PROPERTIES LABELS "unit;integrator")

    # The production cardinal B-spline kernel, evaluated directly. Registered
    # as an ordinary single-process unit test in both the serial and the
    # MPI-enabled build: it exercises no MPI and its result must not depend on
    # whether MPI support was compiled in.
    add_test(NAME gmd_boltzmann_constant COMMAND gmd_boltzmann_constant_tests)
    set_tests_properties(gmd_boltzmann_constant
        PROPERTIES LABELS "unit;integrator;validation")
    add_test(NAME gmd_electrostatic_constant COMMAND gmd_electrostatic_constant_tests)
    set_tests_properties(gmd_electrostatic_constant
        PROPERTIES LABELS "unit;force;validation")
    # The bar <-> eV/A^3 conversion, measured back out of every production path
    # that uses it rather than read from a production symbol.
    add_test(NAME gmd_pressure_units COMMAND gmd_pressure_unit_tests)
    set_tests_properties(gmd_pressure_units
        PROPERTIES LABELS "unit;integrator;io;validation")
    add_test(NAME gmd_pme_bspline COMMAND gmd_pme_bspline_tests)
    set_tests_properties(gmd_pme_bspline PROPERTIES LABELS "unit;force;pme")

    # The production PME energy and force path: the transform normalisation
    # against an explicit direct DFT, and F = -dU/dr by central differences.
    add_test(NAME gmd_pme_energy_force COMMAND gmd_pme_energy_force_tests)
    set_tests_properties(gmd_pme_energy_force
        PROPERTIES LABELS "unit;force;pme;validation")

    # MLForceProvider reports no virial, and nothing downstream may consume one.
    add_test(NAME gmd_ml_virial_contract COMMAND gmd_ml_virial_contract_tests)
    set_tests_properties(gmd_ml_virial_contract PROPERTIES LABELS "unit;force")

    if(GMD_ENABLE_TORCH)
        # WORKING_DIRECTORY keeps the generated .pt fixtures in the build tree.
        add_test(NAME gmd_torchscript_edge_shift
                 COMMAND gmd_torchscript_edge_shift_tests
                 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests)
        set_tests_properties(gmd_torchscript_edge_shift
            PROPERTIES LABELS "unit;force;ml;torch")
    else()
        message(STATUS
            "GMD_ENABLE_TORCH is OFF: gmd_torchscript_edge_shift is NOT registered. "
            "The TorchScript edge_shift contract is unverified in this build.")
    endif()

    # Every virial source held against an independent reference, all nine
    # components. The orthorhombic `Box` cannot express a shear strain, so the
    # off-diagonal claims rest on analytical pair references, independent
    # bonded force moments, and a test-only triclinic reciprocal-space
    # reference that CAN be sheared -- see tests/virial_reference.hpp.
    add_test(NAME gmd_virial_source_inventory COMMAND gmd_virial_source_inventory_tests)
    set_tests_properties(gmd_virial_source_inventory
        PROPERTIES LABELS "unit;force;validation")
    add_test(NAME gmd_pme_reciprocal_virial COMMAND gmd_pme_reciprocal_virial_tests)
    set_tests_properties(gmd_pme_reciprocal_virial
        PROPERTIES LABELS "unit;force;validation;pme")

    # Constrained NVE energy conservation and constrained Berendsen NPT pressure
    # reporting, through the real CLI: many consecutive rescaling steps, and the
    # reported pressure column must never silently become zeros.
    add_test(
        NAME gmd_constrained_pressure_reporting
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/constrained_pressure_reporting.py
                --gmd $<TARGET_FILE:gmd>
                --work-root ${CMAKE_CURRENT_BINARY_DIR}/tests/constrained_pressure
    )
    set_tests_properties(gmd_constrained_pressure_reporting
        PROPERTIES LABELS "integration;integrator;force")
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

    # Compares GMD PME against a checked-in OpenMM PME result, cross-checked
    # against LAMMPS PPPM and LAMMPS Ewald. Regenerating the reference needs
    # OpenMM and LAMMPS; this comparison needs neither, so it runs in normal CI.
    add_test(
        NAME gmd_validation_pme_external
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/validation/run_validation.py
                --case pme_external
                --gmd $<TARGET_FILE:gmd>
                --gmd-validate $<TARGET_FILE:gmd_validate>
                --work-root ${CMAKE_CURRENT_BINARY_DIR}/validation
    )
    set_tests_properties(gmd_validation_pme_external
        PROPERTIES LABELS "validation;unit;pme;external")

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
        add_executable(gmd_mpi_constraint_independence_tests tests/mpi_constraint_independence.cpp)
        target_link_libraries(gmd_mpi_constraint_independence_tests PRIVATE gmd_core)
        add_executable(gmd_mpi_virial_sources_tests tests/mpi_virial_sources.cpp)
        target_link_libraries(gmd_mpi_virial_sources_tests PRIVATE gmd_core)
        add_executable(gmd_mpi_boltzmann_constant_tests
            tests/mpi_boltzmann_constant.cpp)
        target_link_libraries(gmd_mpi_boltzmann_constant_tests PRIVATE gmd_core)

        # Both halves of a temperature are reduced: compute_twice_ke() allreduces
        # the kinetic energy and the degrees of freedom come from a reduced global
        # atom count. A per-rank energy sum and a local DOF count cancel exactly at
        # np=1, so the constant is measured at 1, 2 and 4 ranks. At np=4 one rank
        # deliberately owns no atoms and must still enter every collective.
        foreach(kb_np 1 2 4)
            add_test(
                NAME gmd_mpi_boltzmann_constant_${kb_np}proc
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${kb_np}
                        $<TARGET_FILE:gmd_mpi_boltzmann_constant_tests>
            )
            set_tests_properties(gmd_mpi_boltzmann_constant_${kb_np}proc
                PROPERTIES
                    LABELS "mpi;parallel;integrator;validation"
                    PROCESSORS ${kb_np})
        endforeach()

        add_executable(gmd_mpi_electrostatic_constant_tests
            tests/mpi_electrostatic_constant.cpp)
        target_link_libraries(gmd_mpi_electrostatic_constant_tests PRIVATE gmd_core)

        add_executable(gmd_mpi_pressure_unit_tests tests/mpi_pressure_units.cpp)
        target_link_libraries(gmd_mpi_pressure_unit_tests PRIVATE gmd_core)

        # The bar conversion itself is a compile-time constant, but the pressure
        # it converts is not: 2K is allreduced and the virial is separately
        # reduced. A term summed once per rank, or a conversion applied per rank
        # on the way into a sum, multiplies the reported pressure by the rank
        # count and is an identity at np=1. These runs require the same bar value
        # at 1, 2 and 4 ranks and on every rank, with one empty domain at np=4.
        foreach(pressure_units_np 1 2 4)
            add_test(
                NAME gmd_mpi_pressure_units_${pressure_units_np}proc
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${pressure_units_np}
                        $<TARGET_FILE:gmd_mpi_pressure_unit_tests>
            )
            set_tests_properties(gmd_mpi_pressure_units_${pressure_units_np}proc
                PROPERTIES
                    LABELS "mpi;parallel;integrator;io;validation"
                    PROCESSORS ${pressure_units_np})
        endforeach()

        # The serial audit proves each electrostatic path carries exactly one
        # factor of k_e. Under MPI the providers allreduce their own energies,
        # so a term applied once per rank instead of once per pair would stay
        # self-consistent within a rank while scaling with the rank count.
        # These runs measure the constant at 1, 2 and 4 ranks and require the
        # same number every time and on every rank.
        foreach(ke_np 1 2 4)
            add_test(
                NAME gmd_mpi_electrostatic_constant_${ke_np}proc
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ke_np}
                        $<TARGET_FILE:gmd_mpi_electrostatic_constant_tests>
            )
            set_tests_properties(gmd_mpi_electrostatic_constant_${ke_np}proc
                PROPERTIES
                    LABELS "mpi;parallel;force;validation"
                    PROCESSORS ${ke_np})
        endforeach()

        add_executable(gmd_mpi_ml_edge_shift_tests tests/mpi_ml_edge_shift.cpp)
        target_link_libraries(gmd_mpi_ml_edge_shift_tests PRIVATE gmd_core)

        # The ML path refuses to run under MPI domain decomposition, so
        # edge_shift never reaches a model there. These runs prove that refusal
        # -- at np=1 the provider must work, at np=2 and np=4 both entry points
        # must throw on every rank and the adapter must never be reached. They
        # need no LibTorch: the refusal happens before any adapter is touched.
        foreach(ml_edge_shift_np 1 2 4)
            add_test(
                NAME gmd_mpi_ml_edge_shift_${ml_edge_shift_np}proc
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ml_edge_shift_np}
                        $<TARGET_FILE:gmd_mpi_ml_edge_shift_tests>
            )
            set_tests_properties(gmd_mpi_ml_edge_shift_${ml_edge_shift_np}proc
                PROPERTIES
                    LABELS "mpi;parallel;force;ml"
                    PROCESSORS ${ml_edge_shift_np})
        endforeach()

        # Every virial source, all nine components, at 1, 2 and 4 ranks, held
        # against the same rank-count-independent references the serial test
        # uses. Catches a tensor reduced twice, never reduced, or a pair /
        # bonded term / special-pair correction counted on more than one rank.
        foreach(virial_sources_np 1 2 4)
            add_test(
                NAME gmd_mpi_virial_sources_${virial_sources_np}proc
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${virial_sources_np}
                        $<TARGET_FILE:gmd_mpi_virial_sources_tests>
            )
            set_tests_properties(gmd_mpi_virial_sources_${virial_sources_np}proc
                PROPERTIES
                    LABELS "mpi;parallel;force;validation"
                    PROCESSORS ${virial_sources_np})
        endforeach()

        # The independence verdict must not depend on how atoms are distributed.
        foreach(independence_np 1 2 4)
            add_test(
                NAME gmd_mpi_constraint_independence_${independence_np}proc
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${independence_np}
                        $<TARGET_FILE:gmd_mpi_constraint_independence_tests>
            )
            set_tests_properties(gmd_mpi_constraint_independence_${independence_np}proc
                PROPERTIES
                    LABELS "mpi;parallel;integrator"
                    PROCESSORS ${independence_np})
        endforeach()

        # Same fixture at 1, 2 and 4 ranks: the constraint virial must come out
        # identical, which is what catches a dropped or rank-multiplied value.
        # The MPI CLI output path: write_global_frame() copies frame state onto a
        # separate output System, so only a real run through the executable
        # covers it. Each of these also runs the same fixture serially and
        # compares the two logs column by column.
        foreach(pressure_np 2 4)
            add_test(
                NAME gmd_mpi_constrained_pressure_reporting_${pressure_np}proc
                COMMAND ${Python3_EXECUTABLE}
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/constrained_pressure_reporting.py
                        --gmd $<TARGET_FILE:gmd>
                        --work-root ${CMAKE_CURRENT_BINARY_DIR}/tests/constrained_pressure_mpi_${pressure_np}
                        --np ${pressure_np}
                        --mpiexec ${MPIEXEC_EXECUTABLE}
                        --mpiexec-np-flag=${MPIEXEC_NUMPROC_FLAG}
            )
            set_tests_properties(gmd_mpi_constrained_pressure_reporting_${pressure_np}proc
                PROPERTIES
                    LABELS "mpi;parallel;integration;integrator;force"
                    PROCESSORS ${pressure_np})
        endforeach()

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
            NAME gmd_mpi_nose_hoover_restart_continuity_cli_4proc
            COMMAND Python3::Interpreter
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/nose_hoover_restart_continuity.py
                    --gmd $<TARGET_FILE:gmd>
                    --gmd-validate $<TARGET_FILE:gmd_validate>
                    --work-root ${CMAKE_CURRENT_BINARY_DIR}/tests/nh_restart_continuity_mpi_4
                    --np 4
                    --mpiexec ${MPIEXEC_EXECUTABLE}
                    --mpiexec-np-flag=${MPIEXEC_NUMPROC_FLAG}
        )
        set_tests_properties(gmd_mpi_nose_hoover_restart_continuity_cli_4proc
            PROPERTIES
                LABELS "mpi;parallel;integration;restart;integrator"
                PROCESSORS 4)

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
