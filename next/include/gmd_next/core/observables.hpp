#pragma once

#include "gmd_next/core/numeric.hpp"
#include "gmd_next/core/status.hpp"
#include "gmd_next/core/units.hpp"
#include "gmd_next/core/version.hpp"

#include <array>
#include <bit>
#include <cstdint>
#include <initializer_list>
#include <string_view>

namespace gmd_next::core {

// The quantities a potential evaluation can produce. Forces stay on the device;
// a host record only carries their availability.
enum class Observable : std::uint8_t { forces, potential_energy, virial };

inline constexpr std::array<Observable, 3> kAllObservables{
    Observable::forces, Observable::potential_energy, Observable::virial};

std::string_view name(Observable observable);

class ObservableSet {
public:
    constexpr ObservableSet() = default;
    constexpr ObservableSet(std::initializer_list<Observable> observables) {
        for (const auto observable : observables) insert(observable);
    }

    constexpr ObservableSet& insert(Observable observable) {
        bits_ = static_cast<std::uint8_t>(bits_ | bit(observable));
        return *this;
    }
    constexpr bool contains(Observable observable) const { return (bits_ & bit(observable)) != 0; }
    constexpr bool empty() const { return bits_ == 0; }
    constexpr std::size_t size() const { return static_cast<std::size_t>(std::popcount(bits_)); }
    constexpr ObservableSet operator|(ObservableSet other) const {
        ObservableSet result;
        result.bits_ = static_cast<std::uint8_t>(bits_ | other.bits_);
        return result;
    }
    constexpr ObservableSet operator&(ObservableSet other) const {
        ObservableSet result;
        result.bits_ = static_cast<std::uint8_t>(bits_ & other.bits_);
        return result;
    }
    friend constexpr bool operator==(ObservableSet, ObservableSet) = default;

private:
    static constexpr std::uint8_t bit(Observable observable) {
        return static_cast<std::uint8_t>(1u << static_cast<unsigned>(observable));
    }
    std::uint8_t bits_ = 0;
};

// Where inside a step a sample is taken. Two samples with the same step index
// but different stages describe different states and must not be merged.
//   initial_state            before the first step, at the input geometry
//   after_force_evaluation   post-drift geometry, velocities at the half step
//   completed_step           after the second kick, the state a log line reports
enum class SamplingStage : std::uint8_t { initial_state, after_force_evaluation, completed_step };

std::string_view name(SamplingStage stage);

// How an evaluation leaves the force array it writes into.
//   overwrite   the evaluation owns the array and initializes it
//   accumulate  contributions are added to what is already there
enum class Accumulation : std::uint8_t { overwrite, accumulate };

std::string_view name(Accumulation accumulation);

// One evaluation's requirements, with dynamics kept apart from output so that
// turning output off cannot drop work the integrator needs, and so that an
// output-only quantity is computed only when a subscription is due.
class ObservableRequest {
public:
    ObservableRequest() = default;
    ObservableRequest(ObservableSet dynamics, ObservableSet output, SamplingStage stage,
                      Accumulation accumulation = Accumulation::overwrite);

    // The dynamics of every supported integrator needs the force and nothing else.
    static ObservableRequest for_dynamics_step(ObservableSet output, SamplingStage stage);

    ObservableSet dynamics() const { return dynamics_; }
    ObservableSet output() const { return output_; }
    // What must actually be computed: the union of both needs.
    ObservableSet required() const { return dynamics_ | output_; }
    SamplingStage stage() const { return stage_; }
    Accumulation accumulation() const { return accumulation_; }

    bool needs(Observable observable) const { return required().contains(observable); }
    bool is_output_only(Observable observable) const {
        return output_.contains(observable) && !dynamics_.contains(observable);
    }

    // Subscriptions that share a stage and an accumulation mode share one
    // evaluation; anything else is a different capture point and is rejected.
    ObservableRequest merged_with(const ObservableRequest& other) const;

private:
    ObservableSet dynamics_;
    ObservableSet output_;
    SamplingStage stage_ = SamplingStage::after_force_evaluation;
    Accumulation accumulation_ = Accumulation::overwrite;
};

// Host-side result of one evaluation. It is created from the request, so a
// quantity nobody asked for stays unavailable and can never be read as a zero.
// The state versions and the sampling stage travel with it, which is what later
// output records need in order to prove which state a line describes.
class ObservableRecord {
public:
    ObservableRecord(const ObservableRequest& request, const VersionStamp& versions,
                     const StepStamp& stamp);

    const ObservableRequest& request() const { return request_; }
    ObservableSet requested() const { return request_.required(); }
    SamplingStage stage() const { return request_.stage(); }
    const VersionStamp& versions() const { return versions_; }
    const StepStamp& stamp() const { return stamp_; }

    // Storing a quantity that was not requested is a planning error, not a
    // silently accepted extra. Non-finite values are stored as invalid.
    void store_potential_energy(double energy);
    void store_virial(const Tensor3& virial);
    void store_forces();
    void mark_invalid(Observable observable, ErrorCode code, std::string detail);

    const Quantity<double>& potential_energy() const { return potential_energy_; }
    const Quantity<Tensor3>& virial() const { return virial_; }
    const Presence& forces() const { return forces_; }
    Availability availability(Observable observable) const;
    ErrorCode error(Observable observable) const;

    // Every requested quantity carries a valid value.
    ValidationReport check_complete() const;
    // The recorded state is still the current one for the given domains.
    ValidationReport check_current(const VersionStamp& current,
                                   VersionDomainSet required = kCachedForceDomains) const;

private:
    void require_requested(Observable observable) const;

    ObservableRequest request_;
    VersionStamp versions_;
    StepStamp stamp_;
    Quantity<double> potential_energy_;
    Quantity<Tensor3> virial_;
    Presence forces_;
};

}  // namespace gmd_next::core
