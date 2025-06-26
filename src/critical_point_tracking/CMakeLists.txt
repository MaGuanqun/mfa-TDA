include_directories         (${CMAKE_SOURCE_DIR}/src)          # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/src/utility)              # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/src/critical_point_tracking)          # needed for block.hpp and opts.h


add_executable              (critical_point_tracking        critical_point_tracking.cpp)
add_executable              (degenerate_case        degenerate_case.cpp)

target_link_libraries       (critical_point_tracking        ${libraries})
target_link_libraries       (degenerate_case        ${libraries})

set_property                (TARGET critical_point_tracking degenerate_case PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install                     (TARGETS critical_point_tracking degenerate_case
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/src/critical_point_tracking/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)
