include_directories         (${CMAKE_SOURCE_DIR}/src)          # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/src/utility)

add_executable              (write_vtk              write_vtk.cpp writer.cpp)
add_executable              (write_gradient              write_gradient.cpp)
add_executable              (write_two_value              write_two_value.cpp writer.cpp)
add_executable              (write_h              write_h.cpp writer.cpp)
add_executable              (write_bin_select_slices              write_bin_select_slices.cpp writer.cpp)

target_link_libraries       (write_vtk              ${libraries})
target_link_libraries       (write_gradient              ${libraries} Torch::Torch)
target_link_libraries       (write_two_value              ${libraries})
target_link_libraries       (write_h              ${libraries})
target_link_libraries       (write_bin_select_slices              ${libraries})

# write_gradient now uses INRModel.h (LibTorch); give it the runtime rpaths
set(_wg_rpaths "$ENV{TORCH_PATH}/lib;$ENV{TBB_PATH}/lib;/usr/local/cuda-12.6/lib64")
set_property                (TARGET write_gradient PROPERTY INSTALL_RPATH "${_wg_rpaths}")

set_property                (TARGET write_vtk write_gradient write_two_value write_h write_bin_select_slices           PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

install                     (TARGETS write_vtk write_gradient write_two_value write_h write_bin_select_slices
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/src/convert/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)
