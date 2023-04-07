set(CMAKE_LEGACY_CYGWIN_WIN32 0)
project(KSNP)
cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Detecting C++ default compiler
if (NOT UNIX)
    message(FATAL_ERROR "Non Unix system not supported yet")
endif()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "Default C++ compiler: GCC")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    message(STATUS "Default C++ compiler: LLVM Clang")
else ()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER_ID} currently not supported by CMake build")
endif()
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3")

# Source files
set(SOURCE_FILES
        realignment.h  realignment.cpp
        ksnp_reader.h  ksnp_reader.cpp
        snp_dbg.h      snp_dbg.cpp
        cut_bubble.h   cut_bubble.cpp
        phased_block.h phased_block.cpp
        time_stamp.h)

# Dependencies
if (HTSLIB)
    link_directories(${HTSLIB})
endif()
add_executable(ksnp ${SOURCE_FILES} main.cpp)
target_link_libraries(ksnp -lhts -lz)