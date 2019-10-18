cmake_minimum_required(VERSION 3.12.0)
project(tropicity)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
INCLUDE_DIRECTORIES(${Python3_INCLUDE_DIRS})

SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wpedantic")


add_subdirectory(src)
