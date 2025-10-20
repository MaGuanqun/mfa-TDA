include_directories         (${CMAKE_SOURCE_DIR}/src)          # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/src/utility)              # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/src/critical_point_tracking)          # needed for block.hpp and opts.h


add_executable              (critical_point_tracking        critical_point_tracking.cpp)
add_executable              (critical_point_tracking_explicit        critical_point_tracking_explicit.cpp)
add_executable              (degenerate_case        degenerate_case.cpp)
add_executable              (degenerate_case_explicit        degenerate_case_explicit.cpp)
add_executable              (degenerate_case_INR        degenerate_case_INR.cpp)

add_executable              (control_point_smoothing        control_point_smoothing.cpp)

target_link_libraries       (critical_point_tracking        ${libraries} Torch::Torch)
target_link_libraries       (critical_point_tracking_explicit        ${libraries} Torch::Torch)
target_link_libraries       (degenerate_case        ${libraries} Torch::Torch)
target_link_libraries       (degenerate_case_explicit        ${libraries} Torch::Torch)
target_link_libraries       (degenerate_case_INR        ${libraries} Torch::Torch)
target_link_libraries       (control_point_smoothing        ${libraries} Torch::Torch)

set(_rpaths "$ENV{TORCH_PATH}/lib;$ENV{TBB_PATH}/lib;/usr/local/cuda-12.6/lib64")

set_property(TARGET
  critical_point_tracking
  critical_point_tracking_explicit
  degenerate_case
  degenerate_case_explicit
  degenerate_case_INR
  control_point_smoothing
  PROPERTY INSTALL_RPATH "${_rpaths}"
)

set_property                (TARGET critical_point_tracking critical_point_tracking_explicit degenerate_case degenerate_case_explicit degenerate_case_INR control_point_smoothing PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install                     (TARGETS critical_point_tracking critical_point_tracking_explicit degenerate_case degenerate_case_explicit degenerate_case_INR control_point_smoothing
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/src/critical_point_tracking/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)
