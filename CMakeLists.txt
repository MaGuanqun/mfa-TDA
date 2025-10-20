project                     (mfa-TDA)
cmake_minimum_required      (VERSION 3.15)

set                         (CMAKE_CXX_STANDARD 17)
set                         (CMAKE_CXX_STANDARD_REQUIRED ON)
set                         (CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# if                          (NOT DEFINED CMAKE_CXX_STANDARD)
#     set                     (CMAKE_CXX_STANDARD 17)
# endif                       ()   

message(STATUS "${CMAKE_CXX_STANDARD}")


if                          (NOT CMAKE_BUILD_TYPE)
    set                     (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property            (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif                       (NOT CMAKE_BUILD_TYPE)

option                      (mfa_thread     "tbb, kokkos, sycl, or serial (default tbb)"                                 "tbb"  )
set(mfa_thread "tbb")




# MPI
find_package                (MPI REQUIRED)
set                         (libraries MPI::MPI_C)

list(APPEND libraries ${CMAKE_DL_LIBS} ${CMAKE_THREAD_LIBS_INIT})

# MFA threading model
if                          (mfa_thread MATCHES "tbb")
    message                 ("Using TBB threading")

    add_definitions         (-DMFA_TBB)
    
    find_package(TBB REQUIRED)

    
    if (TBB_FOUND)
        message("Found TBB using CMake: includes in ${TBB_INCLUDE_DIRS}, libraries = TBB::tbb")
        list(APPEND libraries TBB::tbb)
        # include_directories(${TBB_INCLUDE_DIRS})   # Optional; usually not necessary
        # set(libraries ${libraries} TBB::tbb)
    else()
        message(FATAL_ERROR "Could not find TBB")
    endif()
    # find_package            (TBB QUIET)
    # find_path               (TBB_INCLUDE_DIR        tbb/tbb.h)
    # find_library            (TBB_LIBRARY NAMES      tbb)
    # if                      (TBB_INCLUDE_DIR AND TBB_LIBRARY)
    #     include_directories ("${TBB_INCLUDE_DIR}")
    #     set                 (libraries ${libraries} ${TBB_LIBRARY})
    #     message             ("Found TBB in ${TBB_INCLUDE_DIR} and in ${TBB_LIBRARY}")
    # else                    (TBB_INCLUDE_DIR AND TBB_LIBRARY)
    #     message             (FATAL_ERROR "Could not find TBB")
    # endif                   ()
else                        ()
    message                 ("Using no threading (serial)")
    add_definitions         (-DMFA_SERIAL)
endif                       ()

# -----------------------------------------------------------------------------
# Find LibTorch (installed via Spack)
# Make sure Torch_DIR points to the correct cmake package directory, e.g.:
# export Torch_DIR=/home/guanqunma/spack/opt/spack/linux-x86_64_v4/libtorch-2.8.0_cu126-xxxx/share/cmake/Torch
# -----------------------------------------------------------------------------
find_package(Torch REQUIRED)

set(Torch_DIR $ENV{Torch_DIR} CACHE PATH "Path to Torch CMake config" FORCE)
set(TBB_DIR $ENV{TBB_DIR} CACHE PATH "Path to TBB CMake config" FORCE)
set(CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH} CACHE PATH "Prefix path for CMake packages" FORCE)

message(STATUS "Found LibTorch at: ${Torch_DIR}")
message(STATUS "Torch include dirs: ${TORCH_INCLUDE_DIRS}")
message(STATUS "Torch libraries: ${TORCH_LIBRARIES}")

if(NOT TARGET Torch::Torch)
  add_library(Torch::Torch INTERFACE IMPORTED)
  set_property(TARGET Torch::Torch PROPERTY
    INTERFACE_LINK_LIBRARIES "${TORCH_LIBRARIES}")
  target_include_directories(Torch::Torch INTERFACE ${TORCH_INCLUDE_DIRS})
endif()

# Add Torch C++ flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
# list(APPEND libraries ${TORCH_LIBRARIES})   # do this BEFORE the target


# MFA
if                          (NOT DEFINED MFA_PATH AND DEFINED ENV{MFA_PATH})
    set                     (MFA_PATH $ENV{MFA_PATH})
endif                       ()
add_definitions             (-DMFA_NO_WEIGHTS)
# add_definitions             (-DMFA_SERIAL)

# Include dirs
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories         (${MFA_PATH}/include
                            ${MFA_PATH}/include/cppoptlib/include
                            ${MFA_PATH}/include/diy/include
                            ${MFA_PATH}/include/eigen3
                            ${HIGHFIVE_PATH}/include/
                            SYSTEM ${MPI_C_INCLUDE_PATH})

# ---- Prefer Spack libs at runtime (avoid system TBB) ----
# Requires env script to export TORCH_PATH and TBB_PATH (see below).
list(APPEND CMAKE_BUILD_RPATH
  $ENV{TORCH_PATH}/lib
  $ENV{TBB_PATH}/lib
  /usr/local/cuda-12.6/lib64
)

set(CMAKE_INSTALL_RPATH "${CMAKE_BUILD_RPATH}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

# Libraries
# set                         (libraries
                            #  ${libraries}
                            #  ${CMAKE_DL_LIBS} 
                            #  ${CMAKE_THREAD_LIBS_INIT} )

# Subdirectories
add_subdirectory            (src)

# # Install the headers
# file                        (GLOB DEPLOY_FILES_AND_DIRS "${PROJECT_SOURCE_DIR}/include/*")
# foreach                     (ITEM ${DEPLOY_FILES_AND_DIRS})
#    if                       (IS_DIRECTORY "${ITEM}")
#       list                  (APPEND DIRS_TO_DEPLOY "${ITEM}")
#    else                     ()
#       list                  (APPEND FILES_TO_DEPLOY "${ITEM}")
#    endif                    ()
# endforeach                  ()
# install                     (FILES ${FILES_TO_DEPLOY} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
# install                     (DIRECTORY ${DIRS_TO_DEPLOY} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

