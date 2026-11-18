# Prefer the pinned Spack packages exported by load-env.sh. CMake does not
# provision the environment; source overrides remain available as a fallback.
include(FetchContent)
find_path(FIBER_TETGEN_INCLUDE_DIR tetgen.h
    HINTS ENV TETGEN_PATH PATH_SUFFIXES include)
find_library(FIBER_TETGEN_LIBRARY NAMES tet tetgen
    HINTS ENV TETGEN_PATH PATH_SUFFIXES lib lib64)
if(FIBER_TETGEN_INCLUDE_DIR AND FIBER_TETGEN_LIBRARY)
    message(STATUS "Using installed TetGen: ${FIBER_TETGEN_LIBRARY}")
    add_library(fiber_tetgen UNKNOWN IMPORTED)
    set_target_properties(fiber_tetgen PROPERTIES
        IMPORTED_LOCATION "${FIBER_TETGEN_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${FIBER_TETGEN_INCLUDE_DIR}")
else()
    FetchContent_Declare(fiber_tetgen
        URL https://codeload.github.com/TetGen/TetGen/tar.gz/refs/tags/v1.6.0
        URL_HASH SHA256=524bddc38afc205db075b74b67da424b87e676cfd12eccf7915d2e945afcac10)
    FetchContent_GetProperties(fiber_tetgen)
    if(NOT fiber_tetgen_POPULATED)
        FetchContent_Populate(fiber_tetgen)
    endif()
    add_library(fiber_tetgen STATIC ${fiber_tetgen_SOURCE_DIR}/tetgen.cxx ${fiber_tetgen_SOURCE_DIR}/predicates.cxx)
    target_include_directories(fiber_tetgen SYSTEM PUBLIC ${fiber_tetgen_SOURCE_DIR})
    target_compile_definitions(fiber_tetgen PRIVATE TETLIBRARY)
endif()
find_path(FIBER_EARCUT_INCLUDE_DIR mapbox/earcut.hpp
    HINTS ENV EARCUT_PATH PATH_SUFFIXES include)
if(NOT FIBER_EARCUT_INCLUDE_DIR)
    FetchContent_Declare(fiber_earcut
        URL https://codeload.github.com/mapbox/earcut.hpp/tar.gz/refs/tags/v2.2.4
        URL_HASH SHA256=fcfa6a47a52d4c94dc960bdb747f17e077609235517b0bb5ce8097d6b747695a)
    FetchContent_GetProperties(fiber_earcut)
    if(NOT fiber_earcut_POPULATED)
        FetchContent_Populate(fiber_earcut)
    endif()
    set(FIBER_EARCUT_INCLUDE_DIR ${fiber_earcut_SOURCE_DIR}/include)
endif()
message(STATUS "Using earcut.hpp headers: ${FIBER_EARCUT_INCLUDE_DIR}")
