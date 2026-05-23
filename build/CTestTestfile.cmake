# CMake generated Testfile for 
# Source directory: /Users/guo/Documents/workspace/GMD
# Build directory: /Users/guo/Documents/workspace/GMD/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[gmd_smoke_inline_lj]=] "/Users/guo/Documents/workspace/GMD/build/gmd" "/Users/guo/Documents/workspace/GMD/build/tests/smoke_xyz.in" "/Users/guo/Documents/workspace/GMD/build/tests/smoke_run.in")
set_tests_properties([=[gmd_smoke_inline_lj]=] PROPERTIES  LABELS "smoke;integration" _BACKTRACE_TRIPLES "/Users/guo/Documents/workspace/GMD/CMakeLists.txt;109;add_test;/Users/guo/Documents/workspace/GMD/CMakeLists.txt;0;")
add_test([=[gmd_smoke_ewald]=] "/Users/guo/Documents/workspace/GMD/build/gmd" "/Users/guo/Documents/workspace/GMD/build/tests/smoke_ewald.xyz" "/Users/guo/Documents/workspace/GMD/build/tests/smoke_ewald.run")
set_tests_properties([=[gmd_smoke_ewald]=] PROPERTIES  LABELS "smoke;integration" _BACKTRACE_TRIPLES "/Users/guo/Documents/workspace/GMD/CMakeLists.txt;117;add_test;/Users/guo/Documents/workspace/GMD/CMakeLists.txt;0;")
add_test([=[gmd_smoke_pme]=] "/Users/guo/Documents/workspace/GMD/build/gmd" "/Users/guo/Documents/workspace/GMD/build/tests/smoke_ewald.xyz" "/Users/guo/Documents/workspace/GMD/build/tests/smoke_pme.run")
set_tests_properties([=[gmd_smoke_pme]=] PROPERTIES  LABELS "smoke;integration" _BACKTRACE_TRIPLES "/Users/guo/Documents/workspace/GMD/CMakeLists.txt;125;add_test;/Users/guo/Documents/workspace/GMD/CMakeLists.txt;0;")
add_test([=[gmd_smoke_mc_barostat]=] "/Users/guo/Documents/workspace/GMD/build/gmd" "/Users/guo/Documents/workspace/GMD/build/tests/smoke_mc_barostat.xyz" "/Users/guo/Documents/workspace/GMD/build/tests/smoke_mc_barostat.run")
set_tests_properties([=[gmd_smoke_mc_barostat]=] PROPERTIES  LABELS "smoke;integration" _BACKTRACE_TRIPLES "/Users/guo/Documents/workspace/GMD/CMakeLists.txt;133;add_test;/Users/guo/Documents/workspace/GMD/CMakeLists.txt;0;")
add_test([=[gmd_smoke_molecular]=] "/Users/guo/Documents/workspace/GMD/build/gmd" "/Users/guo/Documents/workspace/GMD/build/tests/smoke_molecular.xyz" "/Users/guo/Documents/workspace/GMD/build/tests/smoke_molecular.run" "/Users/guo/Documents/workspace/GMD/build/tests/smoke_molecular.ff" "/Users/guo/Documents/workspace/GMD/build/tests/smoke_molecular.top")
set_tests_properties([=[gmd_smoke_molecular]=] PROPERTIES  LABELS "smoke;integration" WORKING_DIRECTORY "/Users/guo/Documents/workspace/GMD/build" _BACKTRACE_TRIPLES "/Users/guo/Documents/workspace/GMD/CMakeLists.txt;142;add_test;/Users/guo/Documents/workspace/GMD/CMakeLists.txt;0;")
