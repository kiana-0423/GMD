#include "gmd/force/bonded_force_provider.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <string>
#include <stdexcept>

#include "gmd/boundary/minimum_image.hpp"
#include "gmd/runtime/runtime_context.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
namespace {

using Vec3 = std::array<double, 3>;

inline double dot(const Vec3& a, const Vec3& b) noexcept {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline Vec3 cross(const Vec3& a, const Vec3& b) noexcept {
    return {a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]};
}

inline double norm2(const Vec3& v) noexcept {
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

inline double norm(const Vec3& v) noexcept {
    return std::sqrt(norm2(v));
}

inline Vec3 scale(double s, const Vec3& v) noexcept {
    return {s*v[0], s*v[1], s*v[2]};
}

inline Vec3 add(const Vec3& a, const Vec3& b) noexcept {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

inline Vec3 sub(const Vec3& a, const Vec3& b) noexcept {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

// Apply minimum-image on a displacement vector in-place.
inline Vec3 min_image(Vec3 dr, const Box& box) noexcept {
    apply_minimum_image(dr, box);
    return dr;
}

// Accumulate a force vector onto atom `idx` in the result forces array.
inline void accum_force(ForceResult& result, int idx, const Vec3& f) noexcept {
    result.forces[static_cast<std::size_t>(idx)][0] += f[0];
    result.forces[static_cast<std::size_t>(idx)][1] += f[1];
    result.forces[static_cast<std::size_t>(idx)][2] += f[2];
}

inline std::array<int, 2> atom_indices(const BondTerm& term) noexcept {
    return {term.i, term.j};
}

inline std::array<int, 3> atom_indices(const AngleTerm& term) noexcept {
    return {term.i, term.j, term.k};
}

inline std::array<int, 4> atom_indices(const DihedralTerm& term) noexcept {
    return {term.i, term.j, term.k, term.l};
}

inline std::array<int, 4> atom_indices(const ImproperTerm& term) noexcept {
    return {term.i, term.j, term.k, term.l};
}

template <typename Term>
void validate_nonnegative_atom_indices(const Term& term, const char* term_name) {
    for (int atom_index : atom_indices(term)) {
        if (atom_index < 0) {
            throw std::runtime_error(std::string("BondedForceProvider: ") +
                                     term_name + " atom index out of range");
        }
    }
}

template <typename Term>
void validate_atom_indices(const Term& term,
                           std::size_t atom_count,
                           const char* term_name) {
    validate_nonnegative_atom_indices(term, term_name);
    for (int atom_index : atom_indices(term)) {
        if (static_cast<std::size_t>(atom_index) >= atom_count) {
            throw std::runtime_error(std::string("BondedForceProvider: ") +
                                     term_name + " atom index out of range");
        }
    }
}

template <typename Term>
void validate_atom_indices(const std::vector<Term>& terms,
                           std::size_t atom_count,
                           const char* term_name) {
    for (const auto& term : terms) {
        validate_atom_indices(term, atom_count, term_name);
    }
}

template <typename Term>
void validate_nonnegative_atom_indices(const std::vector<Term>& terms,
                                       const char* term_name) {
    for (const auto& term : terms) {
        validate_nonnegative_atom_indices(term, term_name);
    }
}

template <typename Term>
bool has_local_atom(const System& system, const Term& term) noexcept {
    for (int atom_index : atom_indices(term)) {
        if (system.is_local_atom(static_cast<std::size_t>(atom_index))) {
            return true;
        }
    }
    return false;
}

template <typename Term>
bool all_atoms_local(const System& system, const Term& term) noexcept {
    for (int atom_index : atom_indices(term)) {
        if (!system.is_local_atom(static_cast<std::size_t>(atom_index))) {
            return false;
        }
    }
    return true;
}

template <typename Term>
int min_tag_owner(const System& system, const Term& term, int local_rank) noexcept {
    int owner = -1;
    int min_tag = std::numeric_limits<int>::max();
    for (int atom_index : atom_indices(term)) {
        const auto system_index = static_cast<std::size_t>(atom_index);
        const int tag = system.atom_tag(system_index);
        if (tag >= min_tag) {
            continue;
        }

        min_tag = tag;
        owner = system.is_local_atom(system_index)
            ? local_rank
            : system.atom_owner(system_index);
    }
    return owner;
}

#ifdef GMD_ENABLE_MPI
std::size_t infer_global_atom_count(const Topology& topology) noexcept {
    int max_index = -1;
    for (const auto& term : topology.bonds) {
        max_index = std::max(max_index, std::max(term.i, term.j));
    }
    for (const auto& term : topology.angles) {
        max_index = std::max(max_index, std::max({term.i, term.j, term.k}));
    }
    for (const auto& term : topology.dihedrals) {
        max_index = std::max(max_index, std::max({term.i, term.j, term.k, term.l}));
    }
    for (const auto& term : topology.impropers) {
        max_index = std::max(max_index, std::max({term.i, term.j, term.k, term.l}));
    }
    return max_index >= 0 ? static_cast<std::size_t>(max_index) + 1 : 0;
}
#endif

// -----------------------------------------------------------------------
// Dihedral geometry: given four positions, compute the torsion angle φ and
// distribute the generalized force dV/dφ (scalar) onto the four atoms.
//
// Convention (Blondel & Karplus, J. Comput. Chem. 17, 1996):
//   b1 = rj - ri   (i→j)
//   b2 = rk - rj   (j→k)
//   b3 = rl - rk   (k→l)
//   m  = b1 × b2   (normal to plane i-j-k)
//   n  = b2 × b3   (normal to plane j-k-l)
//   φ  = atan2(sin_φ, cos_φ)
//
// The gradient distribution:
//   fi = (-dV/dφ) *  (|b2|/|m|²) * m
//   fl = (-dV/dφ) * -(|b2|/|n|²) * n
//   p  = (b1·b2) / |b2|²
//   q  =  (b3·b2) / |b2|²
//   fj = (p-1)*fi - q*fl
//   fk = -(fi + fj + fl)
// -----------------------------------------------------------------------
inline double dihedral_angle(const Vec3& ri, const Vec3& rj,
                              const Vec3& rk, const Vec3& rl,
                              const Box& box) noexcept {
    Vec3 b1 = min_image(sub(rj, ri), box);
    Vec3 b2 = min_image(sub(rk, rj), box);
    Vec3 b3 = min_image(sub(rl, rk), box);

    Vec3 m  = cross(b1, b2);
    Vec3 n  = cross(b2, b3);

    const double len_m = norm(m);
    const double len_n = norm(n);
    if (len_m < 1.0e-12 || len_n < 1.0e-12) return 0.0;

    // cos φ and sin φ (from triple product)
    const double cos_phi = dot(m, n) / (len_m * len_n);
    const double b2_len  = norm(b2);
    const Vec3   cb      = cross(m, n);          // m × n direction for sin sign
    const double sin_phi = dot(cb, b2) / (len_m * len_n * b2_len);

    return std::atan2(sin_phi, cos_phi);
}

inline void apply_dihedral_forces(ForceResult& result,
                                  int i, int j, int k, int l,
                                  const Vec3& ri, const Vec3& rj,
                                  const Vec3& rk, const Vec3& rl,
                                  double dV_dphi,   // dV/dφ (scalar)
                                  const Box& box) noexcept {
    Vec3 b1 = min_image(sub(rj, ri), box);
    Vec3 b2 = min_image(sub(rk, rj), box);
    Vec3 b3 = min_image(sub(rl, rk), box);

    Vec3 m   = cross(b1, b2);
    Vec3 n   = cross(b2, b3);

    const double m2    = norm2(m);
    const double n2    = norm2(n);
    const double b2len = norm(b2);

    if (m2 < 1.0e-20 || n2 < 1.0e-20 || b2len < 1.0e-12) return;

    // F = -dV/dφ
    const double F   = -dV_dphi;
    const double s   = b2len;

    // Force on i and l
    Vec3 fi = scale( F * s / m2, m);
    Vec3 fl = scale(-F * s / n2, n);

    // Projection coefficients
    const double b2len2 = b2len * b2len;
    const double p = dot(b1, b2) / b2len2;
    const double q = dot(b3, b2) / b2len2;

    // Force on j and k (via constraint ∑ F = 0)
    Vec3 fj = add(scale(p - 1.0, fi), scale(-q, fl));
    Vec3 fk = {-(fi[0]+fj[0]+fl[0]), -(fi[1]+fj[1]+fl[1]), -(fi[2]+fj[2]+fl[2])};

    if (i >= 0) accum_force(result, i, fi);
    if (j >= 0) accum_force(result, j, fj);
    if (k >= 0) accum_force(result, k, fk);
    if (l >= 0) accum_force(result, l, fl);
}

}  // anonymous namespace

// ---------------------------------------------------------------------------
// BondedForceProvider – public API
// ---------------------------------------------------------------------------

int BondedForceProvider::add_bond_type(BondParams p) {
    bond_types_.push_back(p);
    return static_cast<int>(bond_types_.size()) - 1;
}

int BondedForceProvider::add_angle_type(AngleParams p) {
    angle_types_.push_back(p);
    return static_cast<int>(angle_types_.size()) - 1;
}

int BondedForceProvider::add_dihedral_type(DihedralParams p) {
    dihedral_types_.push_back(p);
    return static_cast<int>(dihedral_types_.size()) - 1;
}

int BondedForceProvider::add_improper_type(ImproperParams p) {
    improper_types_.push_back(p);
    return static_cast<int>(improper_types_.size()) - 1;
}

void BondedForceProvider::set_topology(std::shared_ptr<Topology> topology) {
    topology_ = std::move(topology);
}

// ---------------------------------------------------------------------------
// ForceProvider interface
// ---------------------------------------------------------------------------

void BondedForceProvider::initialize(RuntimeContext& /*runtime*/) {
    if (!topology_) {
        throw std::runtime_error("BondedForceProvider: no topology set before initialize()");
    }
    // Validate type indices at initialisation time.
    for (const auto& b : topology_->bonds) {
        if (b.type_idx < 0 || b.type_idx >= static_cast<int>(bond_types_.size()))
            throw std::runtime_error("BondedForceProvider: bond type_idx out of range");
    }
    for (const auto& a : topology_->angles) {
        if (a.type_idx < 0 || a.type_idx >= static_cast<int>(angle_types_.size()))
            throw std::runtime_error("BondedForceProvider: angle type_idx out of range");
    }
    for (const auto& d : topology_->dihedrals) {
        if (d.type_idx < 0 || d.type_idx >= static_cast<int>(dihedral_types_.size()))
            throw std::runtime_error("BondedForceProvider: dihedral type_idx out of range");
    }
    for (const auto& ip : topology_->impropers) {
        if (ip.type_idx < 0 || ip.type_idx >= static_cast<int>(improper_types_.size()))
            throw std::runtime_error("BondedForceProvider: improper type_idx out of range");
    }

    // The System storage is supplied to compute(), so initialize() can catch
    // malformed negative topology indices here and compute() verifies their
    // upper bound against the current local + ghost atom storage.
    validate_nonnegative_atom_indices(topology_->bonds, "bond");
    validate_nonnegative_atom_indices(topology_->angles, "angle");
    validate_nonnegative_atom_indices(topology_->dihedrals, "dihedral");
    validate_nonnegative_atom_indices(topology_->impropers, "improper");
}

void BondedForceProvider::compute(const ForceRequest& request,
                                   ForceResult& result,
                                   RuntimeContext& runtime) {
    const std::size_t n = request.coordinates.size();
    result.success = true;
    result.potential_energy = 0.0;
    result.forces.assign(n, Force3D{0.0, 0.0, 0.0});
    result.virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    result.virial_valid = true;

    if (topology_) {
        compute_system_ = request.system;
        compute_rank_ = runtime.rank();
        compute_use_global_coordinates_ = false;
        compute_global_atom_count_ = 0;
        compute_global_coordinates_ = nullptr;
        compute_owner_by_tag_ = nullptr;
        compute_local_index_by_tag_ = nullptr;

        std::vector<Coordinate3D> global_coordinates;
        std::vector<int> owner_by_tag;
        std::vector<int> local_index_by_tag;

#ifdef GMD_ENABLE_MPI
        if (request.system != nullptr && runtime.size() > 1) {
            constexpr int coordinate_record_width = 5;
            compute_use_global_coordinates_ = true;
            compute_global_atom_count_ = infer_global_atom_count(*topology_);
            global_coordinates.assign(compute_global_atom_count_, Coordinate3D{0.0, 0.0, 0.0});
            owner_by_tag.assign(compute_global_atom_count_, -1);
            local_index_by_tag.assign(compute_global_atom_count_, -1);

            std::vector<double> local_coordinates;
            local_coordinates.reserve(
                request.system->num_local_atoms() * coordinate_record_width);
            for (std::size_t atom_index = 0;
                 atom_index < request.system->num_local_atoms();
                 ++atom_index) {
                const int tag = request.system->atom_tag(atom_index);
                if (tag < 0 ||
                    static_cast<std::size_t>(tag) >= compute_global_atom_count_) {
                    throw std::runtime_error(
                        "BondedForceProvider: local atom tag is outside the bonded topology");
                }

                local_index_by_tag[static_cast<std::size_t>(tag)] =
                    static_cast<int>(atom_index);
                local_coordinates.push_back(static_cast<double>(tag));
                local_coordinates.push_back(request.coordinates[atom_index][0]);
                local_coordinates.push_back(request.coordinates[atom_index][1]);
                local_coordinates.push_back(request.coordinates[atom_index][2]);
                local_coordinates.push_back(static_cast<double>(compute_rank_));
            }

            const int send_count = static_cast<int>(local_coordinates.size());
            std::vector<int> recv_counts(static_cast<std::size_t>(runtime.size()), 0);
            MPI_Allgather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

            std::vector<int> displs(static_cast<std::size_t>(runtime.size()), 0);
            int total_count = 0;
            for (int index = 0; index < runtime.size(); ++index) {
                displs[static_cast<std::size_t>(index)] = total_count;
                total_count += recv_counts[static_cast<std::size_t>(index)];
            }
            if (total_count < 0 || total_count % coordinate_record_width != 0) {
                throw std::runtime_error(
                    "BondedForceProvider: malformed MPI coordinate gather");
            }

            std::vector<double> gathered_coordinates(static_cast<std::size_t>(total_count));
            MPI_Allgatherv(local_coordinates.data(),
                           send_count,
                           MPI_DOUBLE,
                           gathered_coordinates.data(),
                           recv_counts.data(),
                           displs.data(),
                           MPI_DOUBLE,
                           MPI_COMM_WORLD);

            for (std::size_t offset = 0;
                 offset < gathered_coordinates.size();
                 offset += coordinate_record_width) {
                const int tag = static_cast<int>(gathered_coordinates[offset]);
                if (tag < 0 ||
                    static_cast<std::size_t>(tag) >= compute_global_atom_count_) {
                    throw std::runtime_error(
                        "BondedForceProvider: gathered atom tag is outside the bonded topology");
                }

                global_coordinates[static_cast<std::size_t>(tag)] = {
                    gathered_coordinates[offset + 1],
                    gathered_coordinates[offset + 2],
                    gathered_coordinates[offset + 3],
                };
                owner_by_tag[static_cast<std::size_t>(tag)] =
                    static_cast<int>(gathered_coordinates[offset + 4]);
            }

            compute_global_coordinates_ = &global_coordinates;
            compute_owner_by_tag_ = &owner_by_tag;
            compute_local_index_by_tag_ = &local_index_by_tag;
        }
#endif

        const std::size_t atom_count_for_validation =
            compute_use_global_coordinates_ ? compute_global_atom_count_ : n;
        validate_atom_indices(topology_->bonds, atom_count_for_validation, "bond");
        validate_atom_indices(topology_->angles, atom_count_for_validation, "angle");
        validate_atom_indices(topology_->dihedrals, atom_count_for_validation, "dihedral");
        validate_atom_indices(topology_->impropers, atom_count_for_validation, "improper");

        compute_bonds    (request, result);
        compute_angles   (request, result);
        compute_dihedrals(request, result);
        compute_impropers(request, result);

        const bool use_global_coordinates = compute_use_global_coordinates_;

        // Bonded forces are already accumulated per atom. For internal forces with
        // zero net translation, the configurational virial can be formed from the
        // outer product r_i ⊗ F_i and summed over all atoms.
        const std::size_t virial_atom_count = use_global_coordinates
            ? request.system->num_local_atoms()
            : n;
        for (std::size_t i = 0; i < virial_atom_count; ++i) {
            const auto& r = request.coordinates[i];
            const auto& f = result.forces[i];
            result.virial[0] += r[0] * f[0];
            result.virial[1] += r[0] * f[1];
            result.virial[2] += r[0] * f[2];
            result.virial[3] += r[1] * f[0];
            result.virial[4] += r[1] * f[1];
            result.virial[5] += r[1] * f[2];
            result.virial[6] += r[2] * f[0];
            result.virial[7] += r[2] * f[1];
            result.virial[8] += r[2] * f[2];
        }

        compute_system_ = nullptr;
        compute_rank_ = 0;
        compute_use_global_coordinates_ = false;
        compute_global_atom_count_ = 0;
        compute_global_coordinates_ = nullptr;
        compute_owner_by_tag_ = nullptr;
        compute_local_index_by_tag_ = nullptr;
    }

#ifdef GMD_ENABLE_MPI
    int is_initialized = 0;
    int is_finalized = 0;
    MPI_Initialized(&is_initialized);
    MPI_Finalized(&is_finalized);
    if (is_initialized != 0 && is_finalized == 0) {
        auto local_virial = result.virial;
        MPI_Allreduce(local_virial.data(),
                      result.virial.data(),
                      static_cast<int>(result.virial.size()),
                      MPI_DOUBLE,
                      MPI_SUM,
                      MPI_COMM_WORLD);
    }
#endif
}

void BondedForceProvider::finalize(RuntimeContext& /*runtime*/) {}

int BondedForceProvider::compute_owner(const BondTerm& term) const {
    if (compute_owner_by_tag_ != nullptr) {
        const auto& owner_by_tag = *compute_owner_by_tag_;
        return owner_by_tag[static_cast<std::size_t>(std::min(term.i, term.j))];
    }
    return compute_system_ != nullptr
        ? min_tag_owner(*compute_system_, term, compute_rank_)
        : compute_rank_;
}

int BondedForceProvider::compute_owner(const AngleTerm& term) const {
    if (compute_owner_by_tag_ != nullptr) {
        const auto min_tag = std::min({term.i, term.j, term.k});
        return (*compute_owner_by_tag_)[static_cast<std::size_t>(min_tag)];
    }
    return compute_system_ != nullptr
        ? min_tag_owner(*compute_system_, term, compute_rank_)
        : compute_rank_;
}

int BondedForceProvider::compute_owner(const DihedralTerm& term) const {
    if (compute_owner_by_tag_ != nullptr) {
        const auto min_tag = std::min({term.i, term.j, term.k, term.l});
        return (*compute_owner_by_tag_)[static_cast<std::size_t>(min_tag)];
    }
    return compute_system_ != nullptr
        ? min_tag_owner(*compute_system_, term, compute_rank_)
        : compute_rank_;
}

int BondedForceProvider::compute_owner(const ImproperTerm& term) const {
    if (compute_owner_by_tag_ != nullptr) {
        const auto min_tag = std::min({term.i, term.j, term.k, term.l});
        return (*compute_owner_by_tag_)[static_cast<std::size_t>(min_tag)];
    }
    return compute_system_ != nullptr
        ? min_tag_owner(*compute_system_, term, compute_rank_)
        : compute_rank_;
}

bool BondedForceProvider::should_compute(const BondTerm& term) const {
    if (compute_local_index_by_tag_ != nullptr) {
        const auto& local_index_by_tag = *compute_local_index_by_tag_;
        return local_index_by_tag[static_cast<std::size_t>(term.i)] >= 0 ||
               local_index_by_tag[static_cast<std::size_t>(term.j)] >= 0;
    }
    if (compute_system_ == nullptr || all_atoms_local(*compute_system_, term)) {
        return true;
    }
    return has_local_atom(*compute_system_, term) && compute_owner(term) == compute_rank_;
}

bool BondedForceProvider::should_compute(const AngleTerm& term) const {
    if (compute_local_index_by_tag_ != nullptr) {
        const auto& local_index_by_tag = *compute_local_index_by_tag_;
        return local_index_by_tag[static_cast<std::size_t>(term.i)] >= 0 ||
               local_index_by_tag[static_cast<std::size_t>(term.j)] >= 0 ||
               local_index_by_tag[static_cast<std::size_t>(term.k)] >= 0;
    }
    if (compute_system_ == nullptr || all_atoms_local(*compute_system_, term)) {
        return true;
    }
    return has_local_atom(*compute_system_, term) && compute_owner(term) == compute_rank_;
}

bool BondedForceProvider::should_compute(const DihedralTerm& term) const {
    if (compute_local_index_by_tag_ != nullptr) {
        const auto& local_index_by_tag = *compute_local_index_by_tag_;
        return local_index_by_tag[static_cast<std::size_t>(term.i)] >= 0 ||
               local_index_by_tag[static_cast<std::size_t>(term.j)] >= 0 ||
               local_index_by_tag[static_cast<std::size_t>(term.k)] >= 0 ||
               local_index_by_tag[static_cast<std::size_t>(term.l)] >= 0;
    }
    if (compute_system_ == nullptr || all_atoms_local(*compute_system_, term)) {
        return true;
    }
    return has_local_atom(*compute_system_, term) && compute_owner(term) == compute_rank_;
}

bool BondedForceProvider::should_compute(const ImproperTerm& term) const {
    if (compute_local_index_by_tag_ != nullptr) {
        const auto& local_index_by_tag = *compute_local_index_by_tag_;
        return local_index_by_tag[static_cast<std::size_t>(term.i)] >= 0 ||
               local_index_by_tag[static_cast<std::size_t>(term.j)] >= 0 ||
               local_index_by_tag[static_cast<std::size_t>(term.k)] >= 0 ||
               local_index_by_tag[static_cast<std::size_t>(term.l)] >= 0;
    }
    if (compute_system_ == nullptr || all_atoms_local(*compute_system_, term)) {
        return true;
    }
    return has_local_atom(*compute_system_, term) && compute_owner(term) == compute_rank_;
}

bool BondedForceProvider::counts_global_energy(const BondTerm& term) const {
    return !compute_use_global_coordinates_ || compute_owner(term) == compute_rank_;
}

bool BondedForceProvider::counts_global_energy(const AngleTerm& term) const {
    return !compute_use_global_coordinates_ || compute_owner(term) == compute_rank_;
}

bool BondedForceProvider::counts_global_energy(const DihedralTerm& term) const {
    return !compute_use_global_coordinates_ || compute_owner(term) == compute_rank_;
}

bool BondedForceProvider::counts_global_energy(const ImproperTerm& term) const {
    return !compute_use_global_coordinates_ || compute_owner(term) == compute_rank_;
}

const Coordinate3D& BondedForceProvider::coordinate_for_tag(
        std::span<const Coordinate3D> local_coordinates,
        int tag) const {
    if (compute_global_coordinates_ != nullptr) {
        return (*compute_global_coordinates_)[static_cast<std::size_t>(tag)];
    }
    return local_coordinates[static_cast<std::size_t>(tag)];
}

int BondedForceProvider::force_index_for_tag(int tag) const noexcept {
    if (compute_local_index_by_tag_ != nullptr) {
        return (*compute_local_index_by_tag_)[static_cast<std::size_t>(tag)];
    }
    return tag;
}

// ---------------------------------------------------------------------------
// Kernel: harmonic bonds
//   V   = k (r - r0)^2
//   dV/dr = 2k (r - r0)
//   F_i = -dV/dr * r_hat (in the direction from j to i)
//   F_j = -F_i
// ---------------------------------------------------------------------------
void BondedForceProvider::compute_bonds(const ForceRequest& req,
                                         ForceResult& result) const {
    const auto& coords = req.coordinates;
    const Box&  box    = *req.box;

    for (const auto& b : topology_->bonds) {
        if (!should_compute(b)) {
            continue;
        }

        const auto& pi = coordinate_for_tag(coords, b.i);
        const auto& pj = coordinate_for_tag(coords, b.j);

        Vec3 dr = min_image(sub(pj, pi), box);  // vector i→j
        const double r = norm(dr);
        if (r < 1.0e-12) continue;

        const auto& p   = bond_types_[static_cast<std::size_t>(b.type_idx)];
        const double dev = r - p.r0;

        // Energy
        if (counts_global_energy(b)) {
            result.potential_energy += p.k * dev * dev;
        }

        // Force magnitude along i→j
        const double f_mag = -2.0 * p.k * dev / r;  // dF/dr component; applied toward i

        Vec3 fi = scale(-f_mag, dr);  // force on i (toward j when dev > 0, i.e. pull)
        Vec3 fj = scale( f_mag, dr);  // Newton's 3rd law

        const int force_i = force_index_for_tag(b.i);
        const int force_j = force_index_for_tag(b.j);
        if (force_i >= 0) accum_force(result, force_i, fi);
        if (force_j >= 0) accum_force(result, force_j, fj);
    }
}

// ---------------------------------------------------------------------------
// Kernel: harmonic angles
//   Vertex atom is j.
//   V    = k (θ - θ0)^2
//   dV/dθ = 2k (θ - θ0)
//
//   Gradient of θ w.r.t. atom i:
//     ∂θ/∂ri =  1/(|b_ji| sinθ) (cosθ e_ji - e_jk)
//   w.r.t. atom k:
//     ∂θ/∂rk =  1/(|b_jk| sinθ) (cosθ e_jk - e_ji)
//   w.r.t. atom j:
//     ∂θ/∂rj = -(∂θ/∂ri + ∂θ/∂rk)
//   Forces: F = -(dV/dθ) * ∂θ/∂r
// ---------------------------------------------------------------------------
void BondedForceProvider::compute_angles(const ForceRequest& req,
                                          ForceResult& result) const {
    const auto& coords = req.coordinates;
    const Box&  box    = *req.box;

    for (const auto& a : topology_->angles) {
        if (!should_compute(a)) {
            continue;
        }

        const auto& pi = coordinate_for_tag(coords, a.i);
        const auto& pj = coordinate_for_tag(coords, a.j);
        const auto& pk = coordinate_for_tag(coords, a.k);

        Vec3 b_ji = min_image(sub(pi, pj), box);  // j→i
        Vec3 b_jk = min_image(sub(pk, pj), box);  // j→k

        const double len_ji = norm(b_ji);
        const double len_jk = norm(b_jk);
        if (len_ji < 1.0e-12 || len_jk < 1.0e-12) continue;

        // Unit vectors
        Vec3 e_ji = scale(1.0 / len_ji, b_ji);
        Vec3 e_jk = scale(1.0 / len_jk, b_jk);

        double cos_theta = dot(e_ji, e_jk);
        // Clamp to avoid acos domain error
        if (cos_theta >  1.0) cos_theta =  1.0;
        if (cos_theta < -1.0) cos_theta = -1.0;

        const double theta = std::acos(cos_theta);
        const double sin_theta = std::sin(theta);

        const auto& p   = angle_types_[static_cast<std::size_t>(a.type_idx)];
        const double dev = theta - p.theta0;

        // Energy
        if (counts_global_energy(a)) {
            result.potential_energy += p.k * dev * dev;
        }

        if (sin_theta < 1.0e-12) continue;  // near-linear, skip force

        const double dV_dtheta = 2.0 * p.k * dev;
        const double inv_sin   = 1.0 / sin_theta;

        // ∂θ/∂ri direction: (cosθ e_ji - e_jk) / (|b_ji| sinθ)
        Vec3 grad_i = scale(inv_sin / len_ji,
                             sub(scale(cos_theta, e_ji), e_jk));
        // ∂θ/∂rk direction: (cosθ e_jk - e_ji) / (|b_jk| sinθ)
        Vec3 grad_k = scale(inv_sin / len_jk,
                             sub(scale(cos_theta, e_jk), e_ji));

        Vec3 fi = scale(-dV_dtheta, grad_i);
        Vec3 fk = scale(-dV_dtheta, grad_k);
        Vec3 fj = {-(fi[0]+fk[0]), -(fi[1]+fk[1]), -(fi[2]+fk[2])};

        const int force_i = force_index_for_tag(a.i);
        const int force_j = force_index_for_tag(a.j);
        const int force_k = force_index_for_tag(a.k);
        if (force_i >= 0) accum_force(result, force_i, fi);
        if (force_j >= 0) accum_force(result, force_j, fj);
        if (force_k >= 0) accum_force(result, force_k, fk);
    }
}

// ---------------------------------------------------------------------------
// Kernel: periodic proper dihedrals
//   V    = k [1 + cos(n φ - δ)]
//   dV/dφ = -k n sin(n φ - δ)
// ---------------------------------------------------------------------------
void BondedForceProvider::compute_dihedrals(const ForceRequest& req,
                                             ForceResult& result) const {
    const auto& coords = req.coordinates;
    const Box&  box    = *req.box;

    for (const auto& d : topology_->dihedrals) {
        if (!should_compute(d)) {
            continue;
        }

        const auto& pi = coordinate_for_tag(coords, d.i);
        const auto& pj = coordinate_for_tag(coords, d.j);
        const auto& pk = coordinate_for_tag(coords, d.k);
        const auto& pl = coordinate_for_tag(coords, d.l);

        const double phi = dihedral_angle(pi, pj, pk, pl, box);

        const auto& p = dihedral_types_[static_cast<std::size_t>(d.type_idx)];
        const double n_phi_delta = static_cast<double>(p.n) * phi - p.delta;

        // Energy
        if (counts_global_energy(d)) {
            result.potential_energy += p.k * (1.0 + std::cos(n_phi_delta));
        }

        // dV/dφ = -k n sin(nφ - δ)
        const double dV_dphi = -p.k * static_cast<double>(p.n) * std::sin(n_phi_delta);

        apply_dihedral_forces(result,
                              force_index_for_tag(d.i),
                              force_index_for_tag(d.j),
                              force_index_for_tag(d.k),
                              force_index_for_tag(d.l),
                              pi, pj, pk, pl,
                              dV_dphi, box);
    }
}

// ---------------------------------------------------------------------------
// Kernel: harmonic improper dihedrals
//   V    = k (φ - φ0)^2
//   dV/dφ = 2k (φ - φ0)
// Same geometry as proper dihedral; different energy function.
// ---------------------------------------------------------------------------
void BondedForceProvider::compute_impropers(const ForceRequest& req,
                                             ForceResult& result) const {
    const auto& coords = req.coordinates;
    const Box&  box    = *req.box;

    for (const auto& ip : topology_->impropers) {
        if (!should_compute(ip)) {
            continue;
        }

        const auto& pi = coordinate_for_tag(coords, ip.i);
        const auto& pj = coordinate_for_tag(coords, ip.j);
        const auto& pk = coordinate_for_tag(coords, ip.k);
        const auto& pl = coordinate_for_tag(coords, ip.l);

        const double phi = dihedral_angle(pi, pj, pk, pl, box);

        const auto& p = improper_types_[static_cast<std::size_t>(ip.type_idx)];
        const double dev = phi - p.phi0;

        // Energy
        if (counts_global_energy(ip)) {
            result.potential_energy += p.k * dev * dev;
        }

        // dV/dφ = 2k (φ - φ0)
        const double dV_dphi = 2.0 * p.k * dev;

        apply_dihedral_forces(result,
                              force_index_for_tag(ip.i),
                              force_index_for_tag(ip.j),
                              force_index_for_tag(ip.k),
                              force_index_for_tag(ip.l),
                              pi, pj, pk, pl,
                              dV_dphi, box);
    }
}

}  // namespace gmd
