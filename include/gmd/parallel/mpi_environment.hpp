#pragma once

namespace gmd {

class MpiEnvironment {
public:
    MpiEnvironment(int& argc, char**& argv);
    ~MpiEnvironment();

    MpiEnvironment(const MpiEnvironment&) = delete;
    MpiEnvironment& operator=(const MpiEnvironment&) = delete;
    MpiEnvironment(MpiEnvironment&&) = delete;
    MpiEnvironment& operator=(MpiEnvironment&&) = delete;

    static bool initialized();
    static void finalize();

private:
#ifdef GMD_ENABLE_MPI
    bool owns_mpi_{false};
#endif
};

}  // namespace gmd
