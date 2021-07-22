cmake_minimum_required(VERSION 3.0)

project(bitpared)

add_executable(demo src/demo.cpp)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -mpopcnt -std=gnu++11")

# set the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

install(TARGETS demo DESTINATION bin)
add_subdirectory(unittest)
