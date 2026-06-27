include_directories         (${CMAKE_SOURCE_DIR}/src)          # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/src/utility)
include_directories         (${CMAKE_SOURCE_DIR}/src/critical_point_tracking)  # degenerate_case_tracing.h, CP_Trace.h
include_directories         (${CMAKE_SOURCE_DIR}/src/feature_flow_fields)      # feature_flow_fields.h

add_executable              (feature_flow_fields        feature_flow_fields.cpp)
add_executable              (stable_feature_flow_fields stable_feature_flow_fields.cpp)
add_executable              (gradient_norm              gradient_norm.cpp)
add_executable              (hessian_determinant        hessian_determinant.cpp)
add_executable              (spatial_acceleration       spatial_acceleration.cpp)

target_link_libraries       (feature_flow_fields        ${libraries} Torch::Torch)
target_link_libraries       (stable_feature_flow_fields ${libraries} Torch::Torch)
target_link_libraries       (gradient_norm              ${libraries} Torch::Torch)
target_link_libraries       (hessian_determinant        ${libraries} Torch::Torch)
target_link_libraries       (spatial_acceleration       ${libraries} Torch::Torch)

set(_rpaths "$ENV{TORCH_PATH}/lib;$ENV{TBB_PATH}/lib;/usr/local/cuda-12.6/lib64")

set_property                (TARGET feature_flow_fields PROPERTY INSTALL_RPATH "${_rpaths}")
set_property                (TARGET feature_flow_fields PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
set_property                (TARGET stable_feature_flow_fields PROPERTY INSTALL_RPATH "${_rpaths}")
set_property                (TARGET stable_feature_flow_fields PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
set_property                (TARGET gradient_norm PROPERTY INSTALL_RPATH "${_rpaths}")
set_property                (TARGET gradient_norm PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
set_property                (TARGET hessian_determinant PROPERTY INSTALL_RPATH "${_rpaths}")
set_property                (TARGET hessian_determinant PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
set_property                (TARGET spatial_acceleration PROPERTY INSTALL_RPATH "${_rpaths}")
set_property                (TARGET spatial_acceleration PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

install                     (TARGETS feature_flow_fields stable_feature_flow_fields gradient_norm hessian_determinant spatial_acceleration
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/src/feature_flow_fields/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)
