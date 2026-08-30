#!/usr/bin/env python3
"""Independent reimplementation of gmd/core/keyed_random.hpp.

The reference vectors in tests/keyed_random_tests.cpp came from here, not from
capturing the C++ output. This is a different language with arbitrary-precision
integers and shares no code with the header, so agreement between the two is
evidence that the algorithm is implemented as documented rather than that it is
self-consistent.

Run it to regenerate or re-check the vectors:

    python3 tests/keyed_random_reference.py
"""

import math

MASK64 = (1 << 64) - 1
GAMMA = 0x9E3779B97F4A7C15


def splitmix64_mix(z: int) -> int:
    z &= MASK64
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
    return z ^ (z >> 31)


def key_absorb(key: int, value: int) -> int:
    return splitmix64_mix((key + GAMMA + value) & MASK64)


def random_key(seed: int, stream: int, identity: int, component: int) -> int:
    key = 0x243F6A8885A308D3          # fractional bits of pi
    for value in (seed, stream, identity, component):
        key = key_absorb(key, value)
    return key


def random_bits(key: int, index: int) -> int:
    return splitmix64_mix((key + (index + 1) * GAMMA) & MASK64)


def uniform_open_unit(bits: int) -> float:
    # Twelve bits discarded, not eleven: with 53 the top of the range rounds up
    # to exactly 1.0 and log(u) becomes exactly 0.
    return ((bits >> 12) + 0.5) * (1.0 / 2.0 ** 52)


def standard_normal(seed: int, stream: int, identity: int, component: int) -> float:
    key = random_key(seed, stream, identity, component)
    u1 = uniform_open_unit(random_bits(key, 0))
    u2 = uniform_open_unit(random_bits(key, 1))
    return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


VELOCITY_STREAM = 1

if __name__ == "__main__":
    print("keys:")
    for args in [(0, 0, 0, 0), (1, 1, 0, 0), (20260830, 1, 0, 0), (20260830, 1, 0, 1),
                 (20260830, 1, 0, 2), (20260830, 1, 1, 0),
                 (20260830, 1, 4294967295, 2), (4294967295, 1, 123456789, 1)]:
        print(f"  random_key{args} = 0x{random_key(*args):016X}")
    key = random_key(20260830, 1, 0, 0)
    print("stream:")
    print(f"  random_bits(key, 0) = 0x{random_bits(key, 0):016X}")
    print(f"  random_bits(key, 1) = 0x{random_bits(key, 1):016X}")
    print("uniform endpoints:")
    print(f"  uniform_open_unit(0)          = {uniform_open_unit(0)!r}")
    print(f"  uniform_open_unit(2**64 - 1)  = {uniform_open_unit(MASK64)!r}")
    print("normals (seed 20260830, velocity stream):")
    for identity, component in [(0, 0), (0, 1), (0, 2), (1, 0), (7, 2), (23, 1)]:
        value = standard_normal(20260830, VELOCITY_STREAM, identity, component)
        print(f"  identity {identity} component {component} = {value!r}")
