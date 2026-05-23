#pragma once

namespace gmd {

class RuntimeContext {
public:
    RuntimeContext() = default;
    ~RuntimeContext() = default;

    int rank() const;
    int size() const;

    // Optional pointer to an MPI_Comm. Null selects MPI_COMM_WORLD.
    void* mpicom{nullptr};
};

}  // namespace gmd
