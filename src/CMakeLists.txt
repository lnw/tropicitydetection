
add_library(trop SHARED geometry3.cc cube.cc trajectory.cc)

add_executable(main main.cc)
target_link_libraries(main trop)

add_library(tropicity SHARED interface.cc)
target_link_libraries(tropicity trop)

