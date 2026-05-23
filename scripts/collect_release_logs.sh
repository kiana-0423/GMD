#!/usr/bin/env bash
set -u
set -o pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date +%Y%m%d-%H%M%S)"
OUT_DIR="${1:-/tmp/gmd-release-logs-${STAMP}}"
SERIAL_BUILD="${OUT_DIR}/build-serial"
MPI_BUILD="${OUT_DIR}/build-mpi"

mkdir -p "${OUT_DIR}/validation_summaries"

run_and_log() {
  local log_file="$1"
  shift
  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*" | tee "${log_file}"
  "$@" 2>&1 | tee -a "${log_file}"
  local status=${PIPESTATUS[0]}
  echo "exit_code=${status}" | tee -a "${log_file}"
  return "${status}"
}

{
  echo "GMD release log collection"
  echo "timestamp_utc $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "root ${ROOT_DIR}"
  echo "output ${OUT_DIR}"
  echo
  echo "[git]"
  git -C "${ROOT_DIR}" rev-parse --short HEAD 2>/dev/null || true
  git -C "${ROOT_DIR}" status --short 2>/dev/null || true
  echo
  echo "[system]"
  uname -a || true
  sw_vers 2>/dev/null || true
  lscpu 2>/dev/null || true
  sysctl -n machdep.cpu.brand_string 2>/dev/null || true
  sysctl -n hw.ncpu 2>/dev/null || true
  echo
  echo "[tools]"
  cmake --version || true
  c++ --version 2>/dev/null || c++ -v 2>&1 || true
  mpiexec --version 2>/dev/null || true
} > "${OUT_DIR}/environment.txt"

status=0

run_and_log "${OUT_DIR}/serial_configure.log" \
  cmake -S "${ROOT_DIR}" -B "${SERIAL_BUILD}" -DCMAKE_BUILD_TYPE=Release || status=1

run_and_log "${OUT_DIR}/serial_build.log" \
  cmake --build "${SERIAL_BUILD}" -j || status=1

run_and_log "${OUT_DIR}/serial_ctest.log" \
  ctest --test-dir "${SERIAL_BUILD}" --output-on-failure || status=1

run_and_log "${OUT_DIR}/mpi_configure.log" \
  cmake -S "${ROOT_DIR}" -B "${MPI_BUILD}" -DCMAKE_BUILD_TYPE=Release -DGMD_ENABLE_MPI=ON || status=1

run_and_log "${OUT_DIR}/mpi_build.log" \
  cmake --build "${MPI_BUILD}" -j || status=1

run_and_log "${OUT_DIR}/mpi_ctest.log" \
  ctest --test-dir "${MPI_BUILD}" --output-on-failure || status=1

find "${SERIAL_BUILD}/validation" "${MPI_BUILD}/validation" \
  \( -name summary.json -o -name result.json \) 2>/dev/null | while read -r summary; do
    rel="${summary#${SERIAL_BUILD}/validation/}"
    rel="${rel#${MPI_BUILD}/validation/}"
    dest="${OUT_DIR}/validation_summaries/${rel//\//__}"
    cp "${summary}" "${dest}"
  done

echo "release_logs_output=${OUT_DIR}"
echo "overall_exit_code=${status}"
exit "${status}"
