cmake_minimum_required(VERSION 3.15)
# cmake -S src/isosurface -B build/isosurface works without MFA, MPI or Torch.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(adaptive_isosurface LANGUAGES CXX)
    include(CTest)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE HINTS $ENV{EIGEN_PATH})
    set(isosurface_install_dir bin)
else()
    set(isosurface_install_dir src/isosurface)
endif()
add_executable(root_finding_explicit root_finding_explicit.cpp)
add_executable(isosurface_explicit isosurface_explicit.cpp)
add_executable(critical_curves_explicit critical_curves_explicit.cpp)
foreach(target root_finding_explicit isosurface_explicit critical_curves_explicit)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_link_libraries(${target} PRIVATE Eigen3::Eigen)
endforeach()
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    add_executable(root_finding_mfa root_finding_mfa.cpp)
    add_executable(isosurface_mfa isosurface_mfa.cpp)
    add_executable(critical_values_mfa critical_values_mfa.cpp)
    add_executable(critical_curves_mfa critical_curves_mfa.cpp)
    foreach(target root_finding_mfa isosurface_mfa critical_values_mfa critical_curves_mfa)
        target_compile_features(${target} PRIVATE cxx_std_17)
        target_link_libraries(${target} PRIVATE Eigen3::Eigen ${libraries})
    endforeach()
    install(TARGETS root_finding_mfa isosurface_mfa critical_values_mfa critical_curves_mfa
        RUNTIME DESTINATION ${isosurface_install_dir})
endif()
if(BUILD_TESTING)
    add_executable(critical_curves_tests critical_curves_tests.cpp)
    target_compile_features(critical_curves_tests PRIVATE cxx_std_17)
    target_link_libraries(critical_curves_tests PRIVATE Eigen3::Eigen)
    add_test(NAME critical_curves_explicit COMMAND critical_curves_tests)
    set(critical_curve_cli_args)
    if(TARGET isosurface_mfa)
        add_executable(critical_curves_mfa_tests critical_curves_mfa_tests.cpp)
        target_compile_features(critical_curves_mfa_tests PRIVATE cxx_std_17)
        target_link_libraries(critical_curves_mfa_tests PRIVATE Eigen3::Eigen ${libraries})
        add_test(NAME critical_curves_mfa_adapter COMMAND critical_curves_mfa_tests)
        list(APPEND critical_curve_cli_args
            -DMFA_VALUES=$<TARGET_FILE:critical_values_mfa>
            -DMFA_CURVES=$<TARGET_FILE:critical_curves_mfa>
            -DMFA_TEST=$<TARGET_FILE:critical_curves_mfa_tests>)
        add_executable(isosurface_mfa_tests mfa_tests.cpp)
        target_link_libraries(isosurface_mfa_tests PRIVATE Eigen3::Eigen ${libraries})
        add_test(NAME isosurface_mfa_adapter COMMAND isosurface_mfa_tests)
    endif()
    add_test(NAME critical_curves_cli COMMAND ${CMAKE_COMMAND}
        -DEXPLICIT=$<TARGET_FILE:critical_curves_explicit>
        -DWORK=${CMAKE_CURRENT_BINARY_DIR}/critical_curve_cli_test
        ${critical_curve_cli_args}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/critical_curve_cli.cmake)
    add_executable(isosurface_tests tests.cpp)
    target_compile_features(isosurface_tests PRIVATE cxx_std_17)
    target_link_libraries(isosurface_tests PRIVATE Eigen3::Eigen)
    foreach(case geometry sphere ellipsoid torus components scale quality open)
        add_test(NAME isosurface_${case} COMMAND isosurface_tests ${case})
    endforeach()
    add_test(NAME isosurface_cli COMMAND ${CMAKE_COMMAND}
        -DROOT_FINDER=$<TARGET_FILE:root_finding_explicit>
        -DMESHER=$<TARGET_FILE:isosurface_explicit>
        -DWORK=${CMAKE_CURRENT_BINARY_DIR}/cli_test
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_test.cmake)
endif()
install(TARGETS root_finding_explicit isosurface_explicit critical_curves_explicit
    RUNTIME DESTINATION ${isosurface_install_dir})
