#include "gmd_next/core/observables.hpp"

#include <string>

namespace gmd_next::core {

std::string_view name(Observable observable) {
    switch (observable) {
    case Observable::forces: return "forces";
    case Observable::potential_energy: return "potential_energy";
    case Observable::virial: return "virial";
    }
    return "unknown";
}

std::string_view name(SamplingStage stage) {
    switch (stage) {
    case SamplingStage::initial_state: return "initial_state";
    case SamplingStage::after_force_evaluation: return "after_force_evaluation";
    case SamplingStage::completed_step: return "completed_step";
    }
    return "unknown";
}

std::string_view name(Accumulation accumulation) {
    switch (accumulation) {
    case Accumulation::overwrite: return "overwrite";
    case Accumulation::accumulate: return "accumulate";
    }
    return "unknown";
}

ObservableRequest::ObservableRequest(ObservableSet dynamics, ObservableSet output,
                                     SamplingStage stage, Accumulation accumulation)
    : dynamics_(dynamics), output_(output), stage_(stage), accumulation_(accumulation) {}

ObservableRequest ObservableRequest::for_dynamics_step(ObservableSet output, SamplingStage stage) {
    return ObservableRequest{ObservableSet{Observable::forces}, output, stage};
}

ObservableRequest ObservableRequest::merged_with(const ObservableRequest& other) const {
    if (stage_ != other.stage_) {
        throw ContractError({ErrorCode::unsupported_configuration, "observable_request.stage",
                             "cannot merge samples taken at " + std::string(name(stage_)) +
                                 " and " + std::string(name(other.stage_))});
    }
    if (accumulation_ != other.accumulation_) {
        throw ContractError({ErrorCode::unsupported_configuration,
                             "observable_request.accumulation",
                             "cannot merge " + std::string(name(accumulation_)) + " and " +
                                 std::string(name(other.accumulation_)) + " evaluations"});
    }
    return ObservableRequest{dynamics_ | other.dynamics_, output_ | other.output_, stage_,
                             accumulation_};
}

ObservableRecord::ObservableRecord(const ObservableRequest& request, const VersionStamp& versions,
                                   const StepStamp& stamp)
    : request_(request), versions_(versions), stamp_(stamp) {
    validate_step_stamp(stamp, "record").require_ok();
    if (!versions.is_complete()) {
        throw ContractError({ErrorCode::version_mismatch, "record.versions",
                             "a record must carry a complete state version stamp"});
    }
    // Requested quantities start out unevaluated; unrequested ones stay
    // unavailable for the record's whole life.
    const auto requested = request_.required();
    if (requested.contains(Observable::potential_energy)) {
        potential_energy_ = Quantity<double>::missing("requested, not evaluated yet");
    }
    if (requested.contains(Observable::virial)) {
        virial_ = Quantity<Tensor3>::missing("requested, not evaluated yet");
    }
    if (requested.contains(Observable::forces)) {
        forces_ = Presence::missing("requested, not evaluated yet");
    }
}

void ObservableRecord::require_requested(Observable observable) const {
    if (!request_.needs(observable)) {
        throw ContractError({ErrorCode::not_requested, std::string(name(observable)),
                             "this evaluation did not request the quantity"});
    }
}

void ObservableRecord::store_potential_energy(double energy) {
    require_requested(Observable::potential_energy);
    potential_energy_ = Quantity<double>::from_value(energy);
}

void ObservableRecord::store_virial(const Tensor3& virial) {
    require_requested(Observable::virial);
    virial_ = Quantity<Tensor3>::from_value(virial);
}

void ObservableRecord::store_forces() {
    require_requested(Observable::forces);
    forces_ = present();
}

void ObservableRecord::mark_invalid(Observable observable, ErrorCode code, std::string detail) {
    require_requested(observable);
    switch (observable) {
    case Observable::potential_energy:
        potential_energy_ = Quantity<double>::failed(code, std::move(detail));
        return;
    case Observable::virial:
        virial_ = Quantity<Tensor3>::failed(code, std::move(detail));
        return;
    case Observable::forces:
        forces_ = Presence::failed(code, std::move(detail));
        return;
    }
}

Availability ObservableRecord::availability(Observable observable) const {
    switch (observable) {
    case Observable::potential_energy: return potential_energy_.availability();
    case Observable::virial: return virial_.availability();
    case Observable::forces: return forces_.availability();
    }
    return Availability::unavailable;
}

ErrorCode ObservableRecord::error(Observable observable) const {
    switch (observable) {
    case Observable::potential_energy: return potential_energy_.error();
    case Observable::virial: return virial_.error();
    case Observable::forces: return forces_.error();
    }
    return ErrorCode::none;
}

ValidationReport ObservableRecord::check_complete() const {
    ValidationReport report;
    for (const auto observable : kAllObservables) {
        const auto state = availability(observable);
        if (request_.needs(observable)) {
            if (state == Availability::valid) continue;
            report.add(error(observable), std::string(name(observable)),
                       "requested quantity is " + std::string(name(state)));
        } else if (state != Availability::unavailable) {
            report.add(ErrorCode::not_requested, std::string(name(observable)),
                       "unrequested quantity must stay unavailable");
        }
    }
    return report;
}

ValidationReport ObservableRecord::check_current(const VersionStamp& current,
                                                 VersionDomainSet required) const {
    return check_versions(versions_, current, required);
}

}  // namespace gmd_next::core
