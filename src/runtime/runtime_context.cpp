#include "gmd/runtime/runtime_context.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {
namespace {

#ifdef GMD_ENABLE_MPI
bool mpi_is_available() {
    int is_initialized = 0;
    int is_finalized = 0;
    MPI_Initialized(&is_initialized);
    MPI_Finalized(&is_finalized);
    return is_initialized != 0 && is_finalized == 0;
}

MPI_Comm runtime_comm(void* mpicom) {
    if (mpicom == nullptr) {
        return MPI_COMM_WORLD;
    }

    return *static_cast<MPI_Comm*>(mpicom);
}
#endif

}  // namespace

int RuntimeContext::rank() const {
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        int mpi_rank = 0;
        MPI_Comm_rank(runtime_comm(mpicom), &mpi_rank);
        return mpi_rank;
    }
#endif

    return 0;
}

int RuntimeContext::size() const {
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        int mpi_size = 1;
        MPI_Comm_size(runtime_comm(mpicom), &mpi_size);
        return mpi_size;
    }
#endif

    return 1;
}

}  // namespace gmd
