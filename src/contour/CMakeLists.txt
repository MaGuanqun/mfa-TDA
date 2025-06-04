include_directories         (${CMAKE_SOURCE_DIR}/src)          # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/src/contour)          # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/src/utility)  

add_executable              (compute_critical_point_jacobi_set        compute_critical_point_jacobi_set.cpp)
add_executable              (contour        contour.cpp)
add_executable              (jacobi_set        jacobi_set.cpp)
add_executable              (ridge_valley_graph        ridge_valley_graph.cpp)
add_executable              (compute_function_error compute_function_error.cpp)

target_link_libraries         (compute_critical_point_jacobi_set ${libraries})
target_link_libraries       (contour        ${libraries})
target_link_libraries       (jacobi_set        ${libraries})
target_link_libraries       (ridge_valley_graph        ${libraries})
target_link_libraries       (compute_function_error ${libraries})

set_property                (TARGET compute_critical_point_jacobi_set jacobi_set ridge_valley_graph compute_function_error contour PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install                     (TARGETS compute_critical_point_jacobi_set jacobi_set ridge_valley_graph compute_function_error contour
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/src/contour/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)