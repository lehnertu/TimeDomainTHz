cmake_minimum_required (VERSION 2.6)
project (TimeDomainTHz)

set(CMAKE_VERBOSE_MAKEFILE OFF)

# Set default build type.
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "CMAKE_BUILD_TYPE not given; setting to 'Release'.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# HDF5 library
find_package(HDF5 COMPONENTS CXX REQUIRED)
if (HDF5_FOUND)
    message(STATUS "HDF5 includes ${HDF5_INCLUDE_DIRS}")
    message(STATUS "HDF5 libraries ${HDF5_LIBRARIES}")
endif (HDF5_FOUND)
find_package(ZLIB REQUIRED)

# where to find header files
include_directories( ${PROJECT_SOURCE_DIR}/lib_src 
                     ${HDF5_INCLUDE_DIRS} )

# we put all objects into a library that can be linked to the various executables
aux_source_directory(${PROJECT_SOURCE_DIR}/lib_src lib_sources )
message(STATUS "Found source files ${lib_sources}")
add_library(thzlib ${lib_sources})

# compiler flags
add_definitions(-std=c++11 -Wall -g0 -m64 -pthread)

add_executable(GaussianWavePacket ${PROJECT_SOURCE_DIR}/src/GaussianWavePacket.cpp)
target_link_libraries(GaussianWavePacket
        thzlib
        ${HDF5_LIBRARIES}
        ${ZLIB_LIBRARY_RELEASE} )

