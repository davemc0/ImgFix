# ImgFix Solution

cmake_minimum_required(VERSION 3.13)

project(ImgFix)

# CMake settings
set(PROJECT_LIBRARY_TYPE STATIC)
set(CMAKE_CXX_STANDARD 11 CACHE STRING "" FORCE)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(PROJECT_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# warning level and all warnings as errors
add_compile_options(/W3 /WX)
# build in parallel
add_compile_options(/MP)
# disable security warnings
add_compile_options(/sdl-)
# disable exception handler warnings
add_compile_options(/EHsc)

if (MSVC_VERSION GREATER_EQUAL "1920")
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("/std:c++latest" _cpp_latest_flag_supported)
    if (_cpp_latest_flag_supported)
        add_compile_options("/std:c++latest")
    endif()
endif()
set(__ARCH "x64")

set(SOURCES
    Driver.cpp
    FixOps.cpp
    FixOps.h
    Helpers.cpp
    Helpers.h
)

source_group("src"  FILES ${SOURCES})

add_subdirectory(${PROJECT_ROOT_DIR}/../DMcTools ${CMAKE_CURRENT_BINARY_DIR}/DMcTools)

add_executable(ImgFix ${SOURCES})

set_target_properties(ImgFix PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${PROJECT_ROOT_DIR})

target_include_directories(ImgFix
    PRIVATE "${PROJECT_ROOT_DIR}/../DMcTools"
    PRIVATE "${PROJECT_ROOT_DIR}/../DMcTools/deps"
)

target_link_libraries(ImgFix PRIVATE DMcTools)

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ImgFix)
