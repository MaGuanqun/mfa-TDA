include_directories         (${CMAKE_SOURCE_DIR}/src)          # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/src/utility)              # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/src/isosurface)         

add_executable              (isosurface        isosurface.cpp)


target_link_libraries       (isosurface        ${libraries} Torch::Torch)


set(_rpaths "$ENV{TORCH_PATH}/lib;$ENV{TBB_PATH}/lib;/usr/local/cuda-12.6/lib64")

set_property(TARGET
  isosurface
  PROPERTY INSTALL_RPATH "${_rpaths}"
)

set_property                (TARGET isosurface
PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install                     (TARGETS isosurface
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/src/isosurface/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)
