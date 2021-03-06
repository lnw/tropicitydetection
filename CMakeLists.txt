cmake_minimum_required(VERSION 3.12)
project(tropicity)
enable_language(CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

set(CMAKE_VERBOSE_MAKEFILE FALSE)
# set(CMAKE_CXX_CLANG_TIDY clang-tidy)

option(ENABLE_PYTHON  "enable building the python interface"  OFF)
option(ENABLE_CUDA    "enable CUDA"                           OFF)
option(ENABLE_OMP     "enable OpenMP"                         OFF)

if(NOT DEFINED CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "release")
endif()

if(ENABLE_PYTHON)
  if(CMAKE_VERSION VERSION_LESS "3.12.0")
    # find_package(PythonInterp 3.0 REQUIRED)
    find_package(PythonLibs 3.0 REQUIRED)
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
  else()
    find_package(Python3 COMPONENTS Development REQUIRED)
    INCLUDE_DIRECTORIES(${Python3_INCLUDE_DIRS})
  endif()
endif()


if(ENABLE_CUDA)
  enable_language(CUDA)
  add_definitions(-DHAVE_CUDA)
  include_directories(${CUDA_INCLUDE_DIRS})
  message("Cuda include directories ${CUDA_INCLUDE_DIRS}")
  message("Cuda libraries ${CUDA_LIBRARIES}")
  set(LIBS ${LIBS} ${CUDA_LIBRARIES})
endif ()


set(CMAKE_CXX_FLAGS "-fPIC")
if(ENABLE_OMP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  add_definitions(-DHAVE_OMP)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wshadow -Wunused -fsanitize=undefined -fsanitize=leak")

if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
endif()

message(STATUS "enable python: " ${ENABLE_PYTHON} )
message(STATUS "enable CUDA: "   ${ENABLE_CUDA} )
message(STATUS "enable OpenMP: " ${ENABLE_OMP} )
message(STATUS "build type: "    ${CMAKE_BUILD_TYPE} )

add_subdirectory(src)
add_subdirectory(tests)

