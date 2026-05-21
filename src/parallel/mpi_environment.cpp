#include "gmd/parallel/mpi_environment.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {

MpiEnvironment::MpiEnvironment(int& argc, char**& argv) {
#ifdef GMD_ENABLE_MPI
    if (!initialized()) {
        MPI_Init(&argc, &argv);
        owns_mpi_ = true;
    }
#else
    (void)argc;
    (void)argv;
#endif
}

MpiEnvironment::~MpiEnvironment() {
#ifdef GMD_ENABLE_MPI
    if (owns_mpi_) {
        finalize();
    }
#endif
}

bool MpiEnvironment::initialized() {
#ifdef GMD_ENABLE_MPI
    int is_initialized = 0;
    MPI_Initialized(&is_initialized);
    return is_initialized != 0;
#else
    return false;
#endif
}

void MpiEnvironment::finalize() {
#ifdef GMD_ENABLE_MPI
    int is_finalized = 0;
    MPI_Finalized(&is_finalized);
    if (!initialized() || is_finalized != 0) {
        return;
    }

    MPI_Finalize();
#endif
}

}  // namespace gmd
