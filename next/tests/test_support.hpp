#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

inline void expect(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

inline void near(double actual, double expected, double atol = 1e-12, double rtol = 1e-12) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > atol + rtol * std::abs(expected)) {
        throw std::runtime_error("Expected " + std::to_string(expected) +
                                 ", received " + std::to_string(actual));
    }
}

// Exception type plus the code it must carry, for APIs that classify failures.
template<class Exception, class Code, class Function>
void throws_code(Code code, Function&& function) {
    try {
        std::forward<Function>(function)();
    } catch (const Exception& error) {
        if (error.code() != code) {
            throw std::runtime_error("Expected error code " + std::string(name(code)) +
                                     ", received " + std::string(name(error.code())));
        }
        return;
    }
    throw std::runtime_error("Expected exception was not raised");
}

template<class Exception, class Function>
void throws(Function&& function) {
    try {
        std::forward<Function>(function)();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error("Expected exception was not raised");
}
