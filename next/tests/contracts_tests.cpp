#include "gmd_next/core/identity.hpp"
#include "gmd_next/core/observables.hpp"
#include "gmd_next/core/status.hpp"
#include "gmd_next/core/units.hpp"
#include "gmd_next/core/version.hpp"

#include "contract_support.hpp"
#include "test_support.hpp"

#include <cmath>
#include <exception>
#include <iostream>
#include <limits>

namespace {
using namespace gmd_next::core;

constexpr double kNan = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

void identity_and_capacity() {
    expect(!AtomId{}.is_valid(), "A default atom id must not be a usable identity");
    expect(AtomId{0}.is_valid(), "Zero is a legal stable id");
    expect(!AtomId{-1}.is_valid(), "Negative ids mark absence, not identity");
    expect(AtomId{4} < AtomId{9} && AtomId{4} == AtomId{4}, "Atom ids must order by value");

    expect(!LocalIndex{}.is_valid(), "A default local index must not address a slot");
    expect(LocalIndex::from_size(5).value() == 5, "Local index lost its slot");
    // The device index space is narrower than the host size_t that produced it.
    expect(fits_local_index_space(LocalIndex::kCapacityLimit - 1), "Capacity limit is off by one");
    expect(!fits_local_index_space(LocalIndex::kCapacityLimit), "Capacity limit is not enforced");
    throws_code<ContractError>(ErrorCode::capacity_overflow,
                               [] { LocalIndex::from_size(LocalIndex::kCapacityLimit); });
}

void version_rules() {
    StateVersions versions;
    const auto recorded = versions.stamp();
    expect(recorded.is_complete(), "Owner-side versions must start out set");
    expect(recorded.matches(versions.stamp()), "An unchanged state must match its own stamp");

    // An unset stamp certifies nothing, not even against another unset stamp.
    const VersionStamp unset;
    expect(!unset.is_complete(), "A default stamp must not be complete");
    expect(!unset.matches(unset), "Unset versions must never certify a cached result");
    expect(unset.mismatches(recorded).size() == 4, "Every unset domain must be reported");

    versions.advance(VersionDomain::coordinates);
    const auto stale = recorded.mismatches(versions.stamp());
    expect(stale.contains(VersionDomain::coordinates) && stale.size() == 1,
           "Advancing one domain must not disturb the others");

    // A neighbour list survives coordinate motion; that is what the skin is for.
    accepts(check_versions(recorded, versions.stamp(), kNeighborListDomains));
    const auto force_report = check_versions(recorded, versions.stamp(), kCachedForceDomains);
    rejects(force_report, ErrorCode::version_mismatch, "coordinates");
    expect(force_report.diagnostics().size() == 1, "Only the stale domain may be reported");

    // Device reordering invalidates both the list and any slot-indexed force.
    versions.advance(VersionDomain::permutation);
    rejects(check_versions(recorded, versions.stamp(), kNeighborListDomains),
            ErrorCode::version_mismatch, "permutation");
}

void validity_states() {
    const Quantity<double> unrequested;
    expect(unrequested.availability() == Availability::unavailable, "Default must be unavailable");
    expect(unrequested.error() == ErrorCode::not_requested, "Default must say why it is missing");
    expect(std::isnan(value_or_nan(unrequested)), "A missing scalar must read back as NaN");
    throws_code<ContractError>(ErrorCode::not_requested, [&] { unrequested.value(); });

    // Zero is a physical value and must stay distinguishable from "not sampled".
    const auto zero = Quantity<double>::from_value(0.0);
    expect(zero.is_valid() && zero.value() == 0.0, "Zero must be a valid quantity");
    expect(value_or_nan(zero) == 0.0, "A valid zero must not be reported as missing");

    const auto not_finite = Quantity<double>::from_value(kNan);
    expect(not_finite.availability() == Availability::invalid, "NaN must not be stored as valid");
    expect(not_finite.error() == ErrorCode::non_finite_value, "NaN must be classified");
    expect(not_finite.value_or(-1.0) == -1.0, "An invalid quantity must not hand out its payload");

    const auto stale = Quantity<Tensor3>::failed(ErrorCode::version_mismatch, "coordinates moved");
    throws_code<ContractError>(ErrorCode::version_mismatch, [&] { stale.value(); });
    expect(!Quantity<Tensor3>::from_value(Tensor3{0, 0, kInf, 0, 0, 0, 0, 0, 0}).is_valid(),
           "One non-finite component must invalidate a tensor");

    expect(present().is_valid(), "A produced device payload must report as present");
    expect(Presence{}.availability() == Availability::unavailable, "Default presence is missing");

    ValidationReport report;
    expect(report.ok(), "An empty report must be ok");
    report.add(ErrorCode::malformed_input, "field", "detail");
    report.add(ErrorCode::unit_mismatch, "other", "detail");
    expect(!report.ok() && report.contains(ErrorCode::unit_mismatch), "Report lost a diagnostic");
    throws_code<ContractError>(ErrorCode::malformed_input, [&] { report.require_ok(); });
}

void unit_conversions() {
    // Independent derivation, from SI definitions rather than from the header:
    // the internal time unit is A * sqrt(amu/eV), so in femtoseconds it is
    // 1e-10 * sqrt(m_u/e) / 1e-15 with m_u in kg and e in J per eV.
    constexpr double kAtomicMassUnitKg = 1.66053906892e-27;  // CODATA 2022
    constexpr double kJoulesPerEv = 1.602176634e-19;         // exact, SI 2019
    const double femtoseconds_per_internal = 1.0e5 * std::sqrt(kAtomicMassUnitKg / kJoulesPerEv);
    near(units::kFemtosecondsPerInternalTime, femtoseconds_per_internal, 0.0, 1e-12);
    near(to_internal_time(1.0, TimeUnit::femtoseconds, UnitSystem::metal),
         1.0 / femtoseconds_per_internal, 0.0, 1e-12);
    // One internal velocity unit in A/fs, sqrt(eV/amu) converted independently.
    near(1.0 / units::kFemtosecondsPerInternalTime,
         std::sqrt(kJoulesPerEv / kAtomicMassUnitKg) * 1.0e-5, 0.0, 1e-12);

    const double dt = 0.5;
    near(from_internal_time(to_internal_time(dt, TimeUnit::femtoseconds, UnitSystem::metal),
                            TimeUnit::femtoseconds, UnitSystem::metal),
         dt, 0.0, 1e-15);
    expect(to_internal_time(3.5, TimeUnit::internal, UnitSystem::reduced) == 3.5,
           "Internal time must pass through unchanged in any unit system");

    // reduced units carry no scale to seconds, pascal or kelvin.
    throws_code<ContractError>(ErrorCode::unit_mismatch, [] {
        to_internal_time(1.0, TimeUnit::femtoseconds, UnitSystem::reduced);
    });
    throws_code<ContractError>(ErrorCode::unit_mismatch, [] {
        from_internal_time(1.0, TimeUnit::femtoseconds, UnitSystem::reduced);
    });
    throws_code<ContractError>(ErrorCode::unit_mismatch,
                               [] { pressure_to_bar(1.0, UnitSystem::reduced); });
    throws_code<ContractError>(ErrorCode::non_finite_value, [] {
        to_internal_time(kInf, TimeUnit::internal, UnitSystem::metal);
    });

    // 1 eV/A^3 = e J per 1e-30 m^3, over 1e5 Pa per bar.
    const double bar_per_ev_per_cubic_angstrom = kJoulesPerEv / 1.0e-30 / 1.0e5;
    near(pressure_to_bar(1.0, UnitSystem::metal), bar_per_ev_per_cubic_angstrom, 0.0, 1e-15);
    near(pressure_from_bar(pressure_to_bar(2.75, UnitSystem::metal), UnitSystem::metal), 2.75, 0.0,
         1e-15);
    // k_B = 1.380649e-23 J/K over the same exact eV in joules.
    near(units::kBoltzmannEvPerKelvin, 1.380649e-23 / kJoulesPerEv, 0.0, 1e-15);

    accepts(validate_step_stamp(StepStamp{7, 1.5}, "record"));
    rejects(validate_step_stamp(StepStamp{-1, 0.0}, "record"), ErrorCode::value_out_of_range,
            "record.step");
    rejects(validate_step_stamp(StepStamp{0, kNan}, "record"), ErrorCode::non_finite_value,
            "record.time");
}

void observable_requests() {
    const ObservableSet none;
    expect(none.empty() && none.size() == 0, "An empty observable set must be empty");
    const ObservableSet scalars{Observable::potential_energy, Observable::virial};
    expect(scalars.size() == 2 && scalars.contains(Observable::virial), "Set lost an entry");
    expect((scalars | ObservableSet{Observable::forces}).size() == 3, "Union lost an entry");
    expect((scalars & ObservableSet{Observable::virial}) == ObservableSet{Observable::virial},
           "Intersection lost an entry");

    // Dynamics needs the force whether or not anything is being written out.
    const auto quiet = ObservableRequest::for_dynamics_step({}, SamplingStage::after_force_evaluation);
    expect(quiet.required() == ObservableSet{Observable::forces}, "NVE needs exactly the force");
    expect(!quiet.needs(Observable::potential_energy), "Nothing asked for the energy");

    const auto logging =
        ObservableRequest::for_dynamics_step(ObservableSet{Observable::potential_energy},
                                             SamplingStage::after_force_evaluation);
    expect(logging.required().size() == 2, "An output subscription must add to the work");
    expect(logging.is_output_only(Observable::potential_energy), "Energy is not needed by NVE");
    expect(!logging.is_output_only(Observable::forces), "The force is a dynamics requirement");

    const auto merged = quiet.merged_with(logging);
    expect(merged.required() == logging.required(), "Merging subscriptions must take the union");

    // Same step index, different state: these are different capture points.
    const ObservableRequest completed({}, ObservableSet{Observable::virial},
                                      SamplingStage::completed_step);
    throws_code<ContractError>(ErrorCode::unsupported_configuration,
                               [&] { (void)logging.merged_with(completed); });
    const ObservableRequest accumulating({}, ObservableSet{Observable::virial},
                                         SamplingStage::after_force_evaluation,
                                         Accumulation::accumulate);
    throws_code<ContractError>(ErrorCode::unsupported_configuration,
                               [&] { (void)logging.merged_with(accumulating); });
}

void observable_records() {
    StateVersions versions;
    const ObservableRequest request(ObservableSet{Observable::forces},
                                    ObservableSet{Observable::virial},
                                    SamplingStage::after_force_evaluation);
    ObservableRecord record(request, versions.stamp(), StepStamp{3, 0.25});

    // Nothing asked for the energy, so it can never read back as a number.
    expect(record.availability(Observable::potential_energy) == Availability::unavailable,
           "An unrequested quantity must stay unavailable");
    expect(std::isnan(value_or_nan(record.potential_energy())), "Unrequested must read as NaN");
    throws_code<ContractError>(ErrorCode::not_requested,
                               [&] { record.store_potential_energy(0.0); });

    const auto incomplete = record.check_complete();
    rejects(incomplete, ErrorCode::not_evaluated, "forces");
    rejects(incomplete, ErrorCode::not_evaluated, "virial");

    record.store_forces();
    record.store_virial(Tensor3{});
    accepts(record.check_complete());
    // A requested virial that happens to be zero is valid, not missing.
    expect(record.virial().is_valid() && record.virial().value()[0] == 0.0,
           "A zero virial must be a valid value");

    record.mark_invalid(Observable::forces, ErrorCode::capacity_overflow, "neighbour overflow");
    rejects(record.check_complete(), ErrorCode::capacity_overflow, "forces");

    ObservableRecord energy_record(
        ObservableRequest::for_dynamics_step(ObservableSet{Observable::potential_energy},
                                             SamplingStage::after_force_evaluation),
        versions.stamp(), StepStamp{3, 0.25});
    // Requested but not yet produced is also unavailable, never a default zero.
    expect(std::isnan(value_or_nan(energy_record.potential_energy())),
           "An unevaluated energy must read as NaN");
    rejects(energy_record.check_complete(), ErrorCode::not_evaluated, "potential_energy");
    energy_record.store_potential_energy(std::numeric_limits<double>::infinity());
    expect(energy_record.potential_energy().availability() == Availability::invalid,
           "A non-finite energy must not be stored as valid");
    rejects(energy_record.check_complete(), ErrorCode::non_finite_value, "potential_energy");

    // The record carries the state it describes, so a later consumer can prove it.
    accepts(record.check_current(versions.stamp()));
    versions.advance(VersionDomain::coordinates);
    rejects(record.check_current(versions.stamp()), ErrorCode::version_mismatch, "coordinates");
    accepts(record.check_current(versions.stamp(), kNeighborListDomains));
    expect(record.stage() == SamplingStage::after_force_evaluation, "Record lost its stage");
    expect(record.stamp() == (StepStamp{3, 0.25}), "Record lost its step and time");

    throws_code<ContractError>(ErrorCode::version_mismatch, [&] {
        ObservableRecord(request, VersionStamp{}, StepStamp{0, 0.0});
    });
    throws_code<ContractError>(ErrorCode::value_out_of_range, [&] {
        ObservableRecord(request, versions.stamp(), StepStamp{-1, 0.0});
    });
}

}  // namespace

int main() {
    try {
        identity_and_capacity();
        version_rules();
        validity_states();
        unit_conversions();
        observable_requests();
        observable_records();
        std::cout << "Identity, versions, validity, units and observable contracts passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
