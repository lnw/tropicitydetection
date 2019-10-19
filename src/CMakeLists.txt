
add_library(trop STATIC geometry3.cc cube.cc trajectory.cc)

add_executable(main main.cc)
target_link_libraries(main trop)

add_library(tropicity SHARED interface.cc)
target_link_libraries(tropicity trop)

install(TARGETS main DESTINATION bin)
install(TARGETS tropicity DESTINATION lib)
install(FILES splitgrid.py assign_plane.py DESTINATION bin)

