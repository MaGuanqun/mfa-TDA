cmake_minimum_required(VERSION 3.15)
# cmake -S src/isosurface -B build/isosurface works without MFA, MPI or Torch.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(adaptive_isosurface LANGUAGES CXX)
    include(CTest)
endif()
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
add_executable(root_finding_explicit root_finding_explicit.cpp)
add_executable(isosurface_explicit isosurface_explicit.cpp)
foreach(target root_finding_explicit isosurface_explicit)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_link_libraries(${target} PRIVATE Eigen3::Eigen)
endforeach()
if(BUILD_TESTING)
    add_executable(isosurface_tests tests.cpp)
    target_compile_features(isosurface_tests PRIVATE cxx_std_17)
    target_link_libraries(isosurface_tests PRIVATE Eigen3::Eigen)
    foreach(case geometry sphere ellipsoid torus components scale quality)
        add_test(NAME isosurface_${case} COMMAND isosurface_tests ${case})
    endforeach()
    add_test(NAME isosurface_cli COMMAND ${CMAKE_COMMAND}
        -DROOT_FINDER=$<TARGET_FILE:root_finding_explicit>
        -DMESHER=$<TARGET_FILE:isosurface_explicit>
        -DWORK=${CMAKE_CURRENT_BINARY_DIR}/cli_test
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_test.cmake)
endif()
install(TARGETS root_finding_explicit isosurface_explicit RUNTIME DESTINATION bin)
