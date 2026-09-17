#include "gmd_next/reference.hpp"

#include <array>
#include <exception>
#include <iomanip>
#include <iostream>

int main() {
    using namespace gmd_next::reference;
    try {
        const std::array<Vec3, 2> positions{{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}};
        const std::array<Pair, 1> pairs{{{0, 1}}};
        const auto result = evaluate_lj(positions, pairs, LennardJones{1.0, 1.0});
        std::cout << std::setprecision(17)
                  << "GMD Next / CPU mathematical reference (not a GPU simulation)\n"
                  << "LJ epsilon=1, sigma=1, r=1, no cutoff\n"
                  << "U = " << result.energy << '\n'
                  << "F[0].x = " << result.forces[0][0] << '\n'
                  << "F[1].x = " << result.forces[1][0] << '\n'
                  << "W_xx = " << result.virial[0] << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
