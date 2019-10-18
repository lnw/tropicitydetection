
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wpedantic")

add_executable(main main.cc geometry3.cc cube.cc trajectory.cc)

