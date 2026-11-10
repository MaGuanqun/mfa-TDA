add_executable(sample_local_field sample_local_field.cpp)
target_include_directories(sample_local_field PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(sample_local_field PRIVATE ${libraries} Eigen3::Eigen Torch::Torch)
set_property(TARGET sample_local_field PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS sample_local_field RUNTIME DESTINATION src/critical_point_tracking/local_reference)

if(BUILD_TESTING)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_executable(local_sampling_grid_tests grid_tests.cpp)
    add_test(NAME local_sampling_grid COMMAND local_sampling_grid_tests)
    add_executable(local_sampling_fixtures sampling_fixtures.cpp)
    target_include_directories(local_sampling_fixtures PRIVATE ${CMAKE_SOURCE_DIR}/src)
    target_link_libraries(local_sampling_fixtures PRIVATE ${libraries} Eigen3::Eigen Torch::Torch)
    add_test(NAME local_sampling_models COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_sampling.py
        --sampler $<TARGET_FILE:sample_local_field>
        --fixtures $<TARGET_FILE:local_sampling_fixtures>
        --work ${CMAKE_CURRENT_BINARY_DIR}/sampling_test)
endif()
