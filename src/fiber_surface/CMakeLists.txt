add_library(fiber_surface_core STATIC find_seeds.cpp trace_endpoint_fibers.cpp recover_fiber.cpp
    march_surface.cpp stitch_surface.cpp curve_surface.cpp evaluate_surface.cpp)
target_link_libraries(fiber_surface_core PUBLIC Eigen3::Eigen TBB::tbb)
target_compile_features(fiber_surface_core PUBLIC cxx_std_17)
add_library(fiber_coordnet_field STATIC coordnet_field.cpp)
target_link_libraries(fiber_coordnet_field PUBLIC fiber_surface_core PRIVATE Torch::Torch)
add_library(fiber_surface_pipeline STATIC pipeline_io.cpp model_io.cpp pipeline_cli.cpp)
target_link_libraries(fiber_surface_pipeline PUBLIC fiber_surface_core ${libraries} PRIVATE fiber_coordnet_field)
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
option(FIBER_ENABLE_VOLUME "Build filled-range-triangle extraction (TetGen and earcut)" ON)
if(FIBER_ENABLE_VOLUME)
    include(volume_dependencies/CMakeLists.txt)
    add_library(fiber_volume_core STATIC volume_boundary.cpp volume_validation.cpp volume_tetgen.cpp volume_io.cpp
        volume_ply.cpp)
    target_link_libraries(fiber_volume_core PUBLIC fiber_surface_core PRIVATE fiber_tetgen)
    target_include_directories(fiber_volume_core SYSTEM PRIVATE ${FIBER_EARCUT_INCLUDE_DIR})
    target_compile_definitions(fiber_volume_core PRIVATE TETLIBRARY)
    add_executable(fiber_volume_explicit fiber_volume_explicit.cpp)
    target_link_libraries(fiber_volume_explicit PRIVATE fiber_volume_core)
    add_executable(fiber_volume_mfa fiber_volume_mfa.cpp)
    target_link_libraries(fiber_volume_mfa PRIVATE fiber_volume_core ${libraries})
    add_executable(fiber_volume_coordnet fiber_volume_coordnet.cpp)
    target_link_libraries(fiber_volume_coordnet PRIVATE fiber_volume_core fiber_coordnet_field)
    add_executable(fiber_volume_from_patches fiber_volume_from_patches.cpp)
    target_link_libraries(fiber_volume_from_patches PRIVATE fiber_volume_core fiber_surface_pipeline)
    add_executable(fiber_refine_saved_patch experiment/refine_saved_patch.cpp)
    target_link_libraries(fiber_refine_saved_patch PRIVATE fiber_volume_core fiber_surface_pipeline)
    add_executable(fiber_repair_boundary experiment/repair_boundary.cpp)
    target_link_libraries(fiber_repair_boundary PRIVATE fiber_volume_core fiber_surface_pipeline)
    add_executable(fiber_grid_recovery experiment/grid_recovery.cpp)
    target_link_libraries(fiber_grid_recovery PRIVATE fiber_volume_core fiber_surface_pipeline)
    install(TARGETS fiber_volume_explicit fiber_volume_mfa fiber_volume_coordnet fiber_volume_from_patches RUNTIME DESTINATION src/fiber_surface)
endif()
if(BUILD_TESTING)
    if(FIBER_ENABLE_VOLUME)
        add_executable(fiber_volume_tests volume_tests.cpp)
        target_link_libraries(fiber_volume_tests PRIVATE fiber_volume_core)
        add_executable(fiber_volume_review_tests volume_review_tests.cpp)
        target_link_libraries(fiber_volume_review_tests PRIVATE fiber_volume_core)
        foreach(case invalid-patches domain-guard uncached-boundary invalid-cells cap-constraints scale-small scale-tiny translated-box)
            add_test(NAME fiber_volume_${case} COMMAND fiber_volume_review_tests ${case})
        endforeach()
        foreach(case prism curved-prism box hollow-cylinder hollow-sphere sphere disconnected validation)
            add_test(NAME fiber_volume_${case} COMMAND fiber_volume_tests ${case})
        endforeach()
        add_executable(fiber_volume_mfa_tests volume_mfa_tests.cpp)
        target_link_libraries(fiber_volume_mfa_tests PRIVATE fiber_volume_core ${libraries})
        add_test(NAME fiber_volume_mfa_adapter COMMAND fiber_volume_mfa_tests)
    endif()
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    if(FIBER_ENABLE_VOLUME)
        add_test(NAME fiber_surface_boundary_repair COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/experiment/test_boundary_repair.py
            $<TARGET_FILE:fiber_find_seeds> $<TARGET_FILE:fiber_grid_recovery> $<TARGET_FILE:fiber_repair_boundary>
            $<TARGET_FILE:fiber_refine_saved_patch>)
        add_test(NAME fiber_surface_grid_recovery COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/experiment/test_grid_recovery.py
            $<TARGET_FILE:fiber_find_seeds> $<TARGET_FILE:fiber_grid_recovery>)
        add_test(NAME fiber_volume_cli COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/volume_cli_test.py
            --exe $<TARGET_FILE:fiber_volume_explicit>
            --work ${CMAKE_CURRENT_BINARY_DIR}/volume_cli_test)
        add_test(NAME fiber_volume_patch_restart COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/volume_patch_restart_test.py
            $<TARGET_FILE_DIR:fiber_volume_from_patches> ${CMAKE_CURRENT_BINARY_DIR}/volume_patch_restart)
    endif()
    add_executable(fiber_surface_mfa_tests mfa_tests.cpp)
    target_link_libraries(fiber_surface_mfa_tests PRIVATE fiber_surface_core ${libraries})
    add_test(NAME fiber_surface_mfa_adapter COMMAND fiber_surface_mfa_tests)
    set_tests_properties(fiber_surface_mfa_adapter PROPERTIES FIXTURES_SETUP fiber_surface_mfa_model)
    add_executable(fiber_surface_coordnet_tests coordnet_tests.cpp)
    target_link_libraries(fiber_surface_coordnet_tests PRIVATE fiber_surface_pipeline Torch::Torch)
    add_test(NAME fiber_surface_coordnet_adapter COMMAND fiber_surface_coordnet_tests
        ${CMAKE_CURRENT_BINARY_DIR}/coordnet_fixture)
    if(FIBER_ENABLE_VOLUME)
        set_tests_properties(fiber_surface_coordnet_adapter PROPERTIES FIXTURES_SETUP fiber_coordnet_models)
        add_test(NAME fiber_volume_coordnet_cli COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/volume_coordnet_test.py
            $<TARGET_FILE:fiber_volume_coordnet> ${CMAKE_CURRENT_BINARY_DIR}/coordnet_fixture
            ${CMAKE_CURRENT_BINARY_DIR}/volume_coordnet_test)
        set_tests_properties(fiber_volume_coordnet_cli PROPERTIES FIXTURES_REQUIRED fiber_coordnet_models)
    endif()
    add_executable(fiber_surface_tests tests.cpp)
    target_link_libraries(fiber_surface_tests PRIVATE fiber_surface_core)
    foreach(case derivatives segment-sampling fiber-tolerances seed-coverage seeds contacts box-sliver polylines sphere closed-sphere strip disconnected)
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
