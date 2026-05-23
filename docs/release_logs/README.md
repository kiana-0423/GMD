# Release Logs

This directory is reserved for release evidence, not for generated build trees.

For `v2.4`, release evidence should include:

- `environment.txt`: OS, CPU, compiler, CMake, MPI, and git revision.
- `serial_configure.log`
- `serial_build.log`
- `serial_ctest.log`
- `mpi_configure.log`
- `mpi_build.log`
- `mpi_ctest.log`
- `validation_summaries/`: copied JSON summaries from short validation cases.

Use the helper script from the repository root:

```bash
scripts/collect_release_logs.sh /tmp/gmd-v2.4-release-logs
```

The script intentionally builds under the chosen output directory and does not
delete or reuse an existing `build/` or `build-mpi/` directory in the working
tree. MPI tests may require running outside restricted sandboxes because
OpenMPI/PRTE needs local socket access.

Do not treat a copied `provisional_gmd_baseline` validation result as an
independent scientific reference. It is release regression evidence only.
