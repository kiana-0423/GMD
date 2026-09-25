#include "gmd_next/core/version.hpp"

#include <string>

namespace gmd_next::core {
namespace {

std::uint64_t counter_of(const VersionStamp& stamp, VersionDomain domain) {
    switch (domain) {
    case VersionDomain::coordinates: return stamp.coordinates.counter();
    case VersionDomain::cell: return stamp.cell.counter();
    case VersionDomain::model: return stamp.model.counter();
    case VersionDomain::permutation: return stamp.permutation.counter();
    }
    return 0;
}

}  // namespace

std::string_view name(VersionDomain domain) {
    switch (domain) {
    case VersionDomain::coordinates: return "coordinates";
    case VersionDomain::cell: return "cell";
    case VersionDomain::model: return "model";
    case VersionDomain::permutation: return "permutation";
    }
    return "unknown";
}

bool VersionStamp::is_complete() const {
    for (const auto domain : kVersionDomains) {
        if (counter_of(*this, domain) == 0) return false;
    }
    return true;
}

VersionDomainSet VersionStamp::mismatches(const VersionStamp& other) const {
    VersionDomainSet stale;
    for (const auto domain : kVersionDomains) {
        const auto mine = counter_of(*this, domain);
        const auto theirs = counter_of(other, domain);
        // An unset counter never certifies anything, so it always mismatches.
        if (mine == 0 || theirs == 0 || mine != theirs) stale.insert(domain);
    }
    return stale;
}

bool VersionStamp::matches(const VersionStamp& other) const {
    return mismatches(other).empty();
}

void StateVersions::advance(VersionDomain domain) {
    counters_[static_cast<std::size_t>(domain)] += 1;
}

VersionStamp StateVersions::stamp() const {
    return VersionStamp{CoordinateVersion{counters_[0]}, CellVersion{counters_[1]},
                        ModelVersion{counters_[2]}, PermutationVersion{counters_[3]}};
}

ValidationReport check_versions(const VersionStamp& recorded, const VersionStamp& current,
                                VersionDomainSet required) {
    ValidationReport report;
    const auto stale = recorded.mismatches(current) & required;
    for (const auto domain : kVersionDomains) {
        if (!stale.contains(domain)) continue;
        report.add(ErrorCode::version_mismatch, std::string(name(domain)),
                   "recorded " + std::to_string(counter_of(recorded, domain)) + ", current " +
                       std::to_string(counter_of(current, domain)));
    }
    return report;
}

}  // namespace gmd_next::core
