cmake_minimum_required(VERSION 2.8.0)
project(tropicity)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

if(CMAKE_VERSION VERSION_LESS "3.12.0")
  # find_package(PythonInterp 3.0 REQUIRED)
  find_package(PythonLibs 3.0 REQUIRED)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
else()
  find_package(Python3 COMPONENTS Development REQUIRED)
  INCLUDE_DIRECTORIES(${Python3_INCLUDE_DIRS})
endif()

SET(CMAKE_CXX_FLAGS "-fPIC")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wunused -Wpedantic")

if(NOT DEFINED ${CMAKE_INSTALL_PREFIX})
  SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
endif()

add_subdirectory(src)

