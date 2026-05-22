#include "gmd/integrator/thermostat.hpp"

#include "gmd/system/system.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {

double compute_twice_ke(const System& system) noexcept {
    const auto masses    = system.masses();
    const auto velocities = system.velocities();
    double twice_ke = 0.0;
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        const double m = masses[i];
        const auto&  v = velocities[i];
        twice_ke += m * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }

#ifdef GMD_ENABLE_MPI
    int is_initialized = 0;
    int is_finalized = 0;
    MPI_Initialized(&is_initialized);
    MPI_Finalized(&is_finalized);
    if (is_initialized != 0 && is_finalized == 0) {
        double global_twice_ke = 0.0;
        MPI_Allreduce(&twice_ke, &global_twice_ke, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        twice_ke = global_twice_ke;
    }
#endif

    return twice_ke;
}

std::size_t global_atom_count(const System& system) noexcept {
    std::size_t atom_count = system.num_local_atoms();

#ifdef GMD_ENABLE_MPI
    int is_initialized = 0;
    int is_finalized = 0;
    MPI_Initialized(&is_initialized);
    MPI_Finalized(&is_finalized);
    if (is_initialized != 0 && is_finalized == 0) {
        unsigned long long local_count =
            static_cast<unsigned long long>(system.num_local_atoms());
        unsigned long long global_count = 0;
        MPI_Allreduce(&local_count,
                      &global_count,
                      1,
                      MPI_UNSIGNED_LONG_LONG,
                      MPI_SUM,
                      MPI_COMM_WORLD);
        atom_count = static_cast<std::size_t>(global_count);
    }
#endif

    return atom_count;
}

double temperature_from_twice_ke(double twice_ke, std::size_t dof) noexcept {
    if (dof == 0) return 0.0;
    return twice_ke / (static_cast<double>(dof) * kBoltzmann);
}

}  // namespace gmd
