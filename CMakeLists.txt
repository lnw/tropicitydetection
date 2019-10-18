cmake_minimum_required(VERSION 3.1.0)
project(tropicity)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)


SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wpedantic")


add_subdirectory(src)
