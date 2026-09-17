#include "gmd_next/reference.hpp"
#include "test_support.hpp"

#include <array>
#include <exception>
#include <iostream>

int main() {
    using namespace gmd_next::reference;
    try {
        const std::array<Vec3, 3> x{{{0.1, 1.2, -2.0}, {3.0, -0.4, 0.5}, {-1.0, 2.0, 0.3}}};
        const std::array<std::size_t, 4> indices{2, 0, 2, 1};
        const std::array<Vec3, 4> z{{{1.0, 2.0, 3.0}, {-0.2, 0.3, 0.4},
                                    {0.5, -1.0, 0.7}, {0.6, 0.2, -0.9}}};
        const auto q = gather(x, indices);
        const auto assembled = assemble_sum(x.size(), indices, z);
        expect(q[0] == x[2] && q[2] == x[2], "Gather lost ordered repeated slots");
        for (std::size_t axis = 0; axis < 3; ++axis) {
            near(assembled[2][axis], z[0][axis] + z[2][axis]);
        }
        double lhs = 0.0;
        double rhs = 0.0;
        for (std::size_t slot = 0; slot < q.size(); ++slot) {
            for (std::size_t axis = 0; axis < 3; ++axis) lhs += q[slot][axis] * z[slot][axis];
        }
        for (std::size_t atom = 0; atom < x.size(); ++atom) {
            for (std::size_t axis = 0; axis < 3; ++axis) rhs += x[atom][axis] * assembled[atom][axis];
        }
        near(lhs, rhs);
        expect(gather(x, {}).empty(), "Empty gather must be empty");
        const auto empty = assemble_sum(3, {}, {});
        for (const auto& value : empty) expect(value == Vec3{}, "Empty assembly must be zero");
        const std::array<std::size_t, 1> invalid{3};
        throws<std::out_of_range>([&] { gather(x, invalid); });
        throws<std::out_of_range>([&] { assemble_sum(3, invalid, std::span(z).first(1)); });
        throws<std::invalid_argument>([&] { assemble_sum(3, indices, std::span(z).first(1)); });
        std::cout << "Gather/assembly contracts passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
