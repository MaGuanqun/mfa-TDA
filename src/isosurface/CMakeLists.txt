include_directories         (${CMAKE_SOURCE_DIR}/src)
include_directories         (${CMAKE_SOURCE_DIR}/src/utility)
include_directories         (${CMAKE_SOURCE_DIR}/src/isosurface)
include_directories         (${CMAKE_SOURCE_DIR}/src/critical_point_tracking)

add_executable              (root_finding_explicit        root_finding_explicit.cpp)
add_executable              (isosurface_explicit          isosurface_explicit.cpp)

target_link_libraries       (root_finding_explicit        ${libraries} Torch::Torch)
target_link_libraries       (isosurface_explicit          ${libraries} Torch::Torch)

set(_rpaths "$ENV{TORCH_PATH}/lib;$ENV{TBB_PATH}/lib;/usr/local/cuda-12.6/lib64")

set_property(TARGET
  root_finding_explicit isosurface_explicit
  PROPERTY INSTALL_RPATH "${_rpaths}"
)

set_property                (TARGET root_finding_explicit isosurface_explicit
PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

install                     (TARGETS root_finding_explicit isosurface_explicit
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/src/isosurface/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)
