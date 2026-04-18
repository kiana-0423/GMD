#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${ROOT_DIR}/build_demo"
BIN_PATH="${BUILD_DIR}/ethane_demo"

mkdir -p "${BUILD_DIR}"

c++ -std=c++20 -I"${ROOT_DIR}/include" \
  "${ROOT_DIR}/examples/ethane_demo/ethane_demo.cpp" \
  "${ROOT_DIR}/src/boundary/minimum_image.cpp" \
  "${ROOT_DIR}/src/boundary/periodic_boundary.cpp" \
  "${ROOT_DIR}/src/core/simulation.cpp" \
  "${ROOT_DIR}/src/force/bonded_force_provider.cpp" \
  "${ROOT_DIR}/src/force/classical_force_provider.cpp" \
  "${ROOT_DIR}/src/force/ewald_force_provider.cpp" \
  "${ROOT_DIR}/src/force/pme_force_provider.cpp" \
  "${ROOT_DIR}/src/integrator/thermostat.cpp" \
  "${ROOT_DIR}/src/integrator/velocity_rescaling_thermostat.cpp" \
  "${ROOT_DIR}/src/integrator/nose_hoover_thermostat.cpp" \
  "${ROOT_DIR}/src/integrator/berendsen_barostat.cpp" \
  "${ROOT_DIR}/src/integrator/mc_barostat.cpp" \
  "${ROOT_DIR}/src/integrator/velocity_verlet_integrator.cpp" \
  "${ROOT_DIR}/src/io/config_loader.cpp" \
  "${ROOT_DIR}/src/io/trajectory_writer.cpp" \
  "${ROOT_DIR}/src/ml/ml_force_provider.cpp" \
  "${ROOT_DIR}/src/ml/model_runtime_adapter.cpp" \
  "${ROOT_DIR}/src/neighbor/verlet_neighbor_builder.cpp" \
  "${ROOT_DIR}/src/system/initializer.cpp" \
  -o "${BIN_PATH}"

cp "${ROOT_DIR}/examples/ethane_demo/ethane.xyz" "${BUILD_DIR}/ethane.xyz"
cp "${ROOT_DIR}/examples/ethane_demo/ethane.run" "${BUILD_DIR}/ethane.run"
cp "${ROOT_DIR}/examples/ethane_demo/ethane.ff" "${BUILD_DIR}/ethane.ff"
cp "${ROOT_DIR}/examples/ethane_demo/ethane.top" "${BUILD_DIR}/ethane.top"

cd "${BUILD_DIR}"
"${BIN_PATH}" ethane.xyz ethane.run ethane.ff ethane.top ethane_out
