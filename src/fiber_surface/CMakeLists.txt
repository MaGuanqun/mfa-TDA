add_library(fiber_surface_core STATIC find_seeds.cpp trace_endpoint_fibers.cpp
    march_surface.cpp stitch_surface.cpp evaluate_surface.cpp)
target_link_libraries(fiber_surface_core PUBLIC Eigen3::Eigen)
target_compile_features(fiber_surface_core PUBLIC cxx_std_17)
add_library(fiber_surface_pipeline STATIC pipeline_io.cpp model_io.cpp pipeline_cli.cpp)
target_link_libraries(fiber_surface_pipeline PUBLIC fiber_surface_core ${libraries})
set(fiber_stage_targets fiber_find_seeds fiber_trace_endpoint_fibers fiber_march_surface
    fiber_stitch_surface fiber_evaluate_surface)
set(stage_number 0)
foreach(target IN LISTS fiber_stage_targets)
    math(EXPR stage_number "${stage_number} + 1")
    add_executable(${target} stage_main.cpp)
    target_compile_definitions(${target} PRIVATE FIBER_STAGE=${stage_number})
    target_link_libraries(${target} PRIVATE fiber_surface_pipeline)
endforeach()
install(TARGETS ${fiber_stage_targets} RUNTIME DESTINATION src/fiber_surface)
add_executable(fiber_surface_explicit fiber_surface_explicit.cpp)
add_executable(fiber_surface_mfa fiber_surface_mfa.cpp)
target_link_libraries(fiber_surface_explicit PRIVATE fiber_surface_core)
target_link_libraries(fiber_surface_mfa PRIVATE fiber_surface_core ${libraries})
install(TARGETS fiber_surface_explicit fiber_surface_mfa RUNTIME DESTINATION src/fiber_surface)
if(BUILD_TESTING)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_executable(fiber_surface_mfa_tests mfa_tests.cpp)
    target_link_libraries(fiber_surface_mfa_tests PRIVATE fiber_surface_core ${libraries})
    add_test(NAME fiber_surface_mfa_adapter COMMAND fiber_surface_mfa_tests)
    set_tests_properties(fiber_surface_mfa_adapter PROPERTIES FIXTURES_SETUP fiber_surface_mfa_model)
    add_executable(fiber_surface_tests tests.cpp)
    target_link_libraries(fiber_surface_tests PRIVATE fiber_surface_core)
    foreach(case derivatives sphere closed-sphere strip disconnected)
        add_test(NAME fiber_surface_${case} COMMAND fiber_surface_tests ${case})
    endforeach()
    add_test(NAME fiber_surface_cli COMMAND ${CMAKE_COMMAND}
        -DEXPLICIT=$<TARGET_FILE:fiber_surface_explicit>
        -DWORK=${CMAKE_CURRENT_BINARY_DIR}/cli_test
        -P ${CMAKE_CURRENT_SOURCE_DIR}/cli_test.cmake)
    add_test(NAME fiber_surface_pipeline COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/pipeline_test.py
        --bin $<TARGET_FILE_DIR:fiber_find_seeds>
        --work ${CMAKE_CURRENT_BINARY_DIR}/pipeline_test)
    add_test(NAME fiber_surface_pipeline_mfa COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/pipeline_test.py
        --bin $<TARGET_FILE_DIR:fiber_find_seeds>
        --work ${CMAKE_CURRENT_BINARY_DIR}/pipeline_mfa_test
        --fixture ${CMAKE_CURRENT_BINARY_DIR}/fiber_surface_mfa_fixture/sphere.mfa)
    set_tests_properties(fiber_surface_pipeline_mfa PROPERTIES FIXTURES_REQUIRED fiber_surface_mfa_model)
endif()
