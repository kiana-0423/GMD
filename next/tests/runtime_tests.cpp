#include "gmd_next/runtime/completion.hpp"
#include "gmd_next/runtime/device_view.hpp"

#include "contract_support.hpp"
#include "test_support.hpp"

#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>

// Host-side bookkeeping of the CUDA runtime. The CUDA calls themselves are
// covered by backends/cuda/tests and need a GPU; these checks do not.
namespace {
using namespace gmd_next;
using core::ContractError;
using core::ErrorCode;
using runtime::DeviceId;

constexpr std::size_t kMax = std::numeric_limits<std::size_t>::max();

void capacity_growth() {
    expect(runtime::grown_capacity(0, 10, 8) == 10, "First growth must cover the request");
    expect(runtime::grown_capacity(100, 101, 8) == 150, "Growth must be geometric");
    expect(runtime::grown_capacity(100, 400, 8) == 400, "A large request wins over x1.5");
    expect(runtime::grown_capacity(100, 60, 8) == 100, "Capacity never shrinks");
    // Near the top of size_t the geometric step clamps instead of wrapping.
    const std::size_t top = kMax / 8;
    expect(runtime::grown_capacity(top - 10, top - 5, 8) == top, "Growth must clamp, not wrap");
    throws_code<ContractError>(ErrorCode::capacity_overflow,
                               [&] { runtime::grown_capacity(0, top + 1, 8); });
    throws_code<ContractError>(ErrorCode::malformed_input, [] { runtime::grown_capacity(0, 1, 0); });
}

void views_and_generations() {
    const DeviceId gpu0{0};
    runtime::CapacityTracker tracker(sizeof(double));
    auto plan = tracker.plan(10);
    expect(plan.reallocate && plan.capacity == 10, "An empty buffer must allocate");
    tracker.commit(plan);
    const auto first = tracker.identity();
    accepts(tracker.check(first, 10, gpu0, gpu0));

    // Within capacity: no reallocation, but the size change still ends views.
    plan = tracker.plan(6);
    expect(!plan.reallocate && plan.capacity == 10, "Shrinking must not reallocate");
    tracker.commit(plan);
    rejects(tracker.check(first, 10, gpu0, gpu0), ErrorCode::stale_view, "view.generation");
    const auto second = tracker.identity();
    // A resize to the same size is a no-op for views: per-step calls stay cheap.
    tracker.commit(tracker.plan(6));
    accepts(tracker.check(second, 6, gpu0, gpu0));

    plan = tracker.plan(11);
    expect(plan.reallocate && plan.capacity == 15, "Growth past capacity must reallocate x1.5");
    tracker.commit(plan);
    rejects(tracker.check(second, 6, gpu0, gpu0), ErrorCode::stale_view, "view.generation");

    const auto current = tracker.identity();
    rejects(tracker.check(current, 11, DeviceId{1}, gpu0), ErrorCode::device_mismatch,
            "view.device");
    rejects(tracker.check(current, 12, gpu0, gpu0), ErrorCode::inconsistent_size, "view.size");
    runtime::CapacityTracker other(sizeof(double));
    rejects(other.check(current, 0, gpu0, gpu0), ErrorCode::stale_view, "view.buffer");

    tracker.release();
    expect(tracker.size() == 0 && tracker.capacity() == 0, "Release must drop storage");
    rejects(tracker.check(current, 0, gpu0, gpu0), ErrorCode::stale_view, "view.generation");
    throws_code<ContractError>(ErrorCode::malformed_input, [&] {
        tracker.commit(runtime::ResizePlan{false, 99, 1});
    });

    rejects(runtime::check_same_device(DeviceId{}, DeviceId{}, "view"), ErrorCode::device_mismatch,
            "view");
    const runtime::DeviceSpan<double> span(nullptr, 3, gpu0, current);
    expect(span.size_bytes() == 3 * sizeof(double), "View byte size is wrong");
}

void completions_and_faults() {
    accepts(runtime::decode_device_faults(0));
    const auto both = runtime::decode_device_faults(
        runtime::bit(runtime::DeviceFault::capacity_exceeded) |
        runtime::bit(runtime::DeviceFault::non_finite_result));
    rejects(both, ErrorCode::capacity_overflow, "device_status");
    rejects(both, ErrorCode::non_finite_value, "device_status");
    expect(both.diagnostics().size() == 2, "Each fault bit must be one diagnostic");
    rejects(runtime::decode_device_faults(runtime::bit(runtime::DeviceFault::invalid_argument)),
            ErrorCode::malformed_input, "device_status");
    // Bits nobody defined are reported, never dropped.
    rejects(runtime::decode_device_faults(1u << 12), ErrorCode::execution_failed, "device_status");

    runtime::CompletionLedger ledger;
    const auto first = ledger.begin();
    expect(ledger.is_outstanding(first), "A new ticket must be outstanding");
    throws_code<ContractError>(ErrorCode::resource_busy, [&] { ledger.begin(); });
    ledger.resolve(first);
    throws_code<ContractError>(ErrorCode::malformed_input, [&] { ledger.resolve(first); });
    const auto second = ledger.begin();
    expect(second != first, "Tickets must not be reused");
    throws_code<ContractError>(ErrorCode::malformed_input, [&] { ledger.resolve(first); });
    ledger.resolve(second);
}

}  // namespace

int main() {
    try {
        capacity_growth();
        views_and_generations();
        completions_and_faults();
        std::cout << "Runtime capacity, view, completion and fault bookkeeping passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
