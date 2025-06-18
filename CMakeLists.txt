project                     (mfa-TDA)
cmake_minimum_required      (VERSION 3.15)

# set                         (CMAKE_CXX_STANDARD 17)
# set                         (CMAKE_CXX_STANDARD_REQUIRED ON)
# set                         (CMAKE_CXX_EXTENSIONS OFF)

if                          (NOT DEFINED CMAKE_CXX_STANDARD)
    set                     (CMAKE_CXX_STANDARD 17)
endif                       ()   

message(STATUS "${CMAKE_CXX_STANDARD}")

# set                         (CMAKE_POSITION_INDEPENDENT_CODE ON)

# ignore variadic args warning (no other good solution)
# set                         (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-varargs")

# Default to Release
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




# Threads
# find_package                (Threads)

# # DIY
# option                      (build_examples "" OFF)
# option                      (build_tests "" OFF)
# add_subdirectory            (ext/diy)
# set                         (libraries ${libraries} diy)
# 
# # fmt
# option                      (FMT_INSTALL "" OFF)
# add_subdirectory            (ext/fmt)
# set                         (libraries ${libraries} fmt::fmt)

# # zlib
# find_package 		        (ZLIB)
# set 			            (libraries ${libraries} ${ZLIB_LIBRARIES})


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

