# CMake generated Testfile for 
# Source directory: /home/guozy/workspace/GMD
# Build directory: /home/guozy/workspace/GMD/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[gmd_smoke_inline_lj]=] "/home/guozy/workspace/GMD/build/gmd" "/home/guozy/workspace/GMD/build/tests/smoke_xyz.in" "/home/guozy/workspace/GMD/build/tests/smoke_run.in")
set_tests_properties([=[gmd_smoke_inline_lj]=] PROPERTIES  LABELS "smoke;integration" _BACKTRACE_TRIPLES "/home/guozy/workspace/GMD/CMakeLists.txt;103;add_test;/home/guozy/workspace/GMD/CMakeLists.txt;0;")
add_test([=[gmd_smoke_ewald]=] "/home/guozy/workspace/GMD/build/gmd" "/home/guozy/workspace/GMD/build/tests/smoke_ewald.xyz" "/home/guozy/workspace/GMD/build/tests/smoke_ewald.run")
set_tests_properties([=[gmd_smoke_ewald]=] PROPERTIES  LABELS "smoke;integration" _BACKTRACE_TRIPLES "/home/guozy/workspace/GMD/CMakeLists.txt;111;add_test;/home/guozy/workspace/GMD/CMakeLists.txt;0;")
add_test([=[gmd_smoke_pme]=] "/home/guozy/workspace/GMD/build/gmd" "/home/guozy/workspace/GMD/build/tests/smoke_ewald.xyz" "/home/guozy/workspace/GMD/build/tests/smoke_pme.run")
set_tests_properties([=[gmd_smoke_pme]=] PROPERTIES  LABELS "smoke;integration" _BACKTRACE_TRIPLES "/home/guozy/workspace/GMD/CMakeLists.txt;119;add_test;/home/guozy/workspace/GMD/CMakeLists.txt;0;")
add_test([=[gmd_smoke_mc_barostat]=] "/home/guozy/workspace/GMD/build/gmd" "/home/guozy/workspace/GMD/build/tests/smoke_mc_barostat.xyz" "/home/guozy/workspace/GMD/build/tests/smoke_mc_barostat.run")
set_tests_properties([=[gmd_smoke_mc_barostat]=] PROPERTIES  LABELS "smoke;integration" _BACKTRACE_TRIPLES "/home/guozy/workspace/GMD/CMakeLists.txt;127;add_test;/home/guozy/workspace/GMD/CMakeLists.txt;0;")
add_test([=[gmd_smoke_molecular]=] "/home/guozy/workspace/GMD/build/gmd" "/home/guozy/workspace/GMD/build/tests/smoke_molecular.xyz" "/home/guozy/workspace/GMD/build/tests/smoke_molecular.run" "/home/guozy/workspace/GMD/build/tests/smoke_molecular.ff" "/home/guozy/workspace/GMD/build/tests/smoke_molecular.top")
set_tests_properties([=[gmd_smoke_molecular]=] PROPERTIES  LABELS "smoke;integration" WORKING_DIRECTORY "/home/guozy/workspace/GMD/build" _BACKTRACE_TRIPLES "/home/guozy/workspace/GMD/CMakeLists.txt;136;add_test;/home/guozy/workspace/GMD/CMakeLists.txt;0;")
