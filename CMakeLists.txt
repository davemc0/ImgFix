# ImgFix Solution

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

set(EXE_NAME ImgFix)

project(${EXE_NAME})

set(CMAKE_CXX_STANDARD 17)
set(PROJECT_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

add_compile_options(/MP) # build in parallel

set(SOURCES
    Driver.cpp
    FixOps.cpp
    FixOps.h
    Helpers.cpp
    Helpers.h
)

source_group("src"  FILES ${SOURCES})

add_subdirectory(${PROJECT_ROOT_DIR}/../DMcTools ${CMAKE_CURRENT_BINARY_DIR}/DMcTools)

add_executable(${EXE_NAME} ${SOURCES})

set_target_properties(${EXE_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${PROJECT_ROOT_DIR} )
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${EXE_NAME})

target_link_libraries(${EXE_NAME} PRIVATE DMcTools)
