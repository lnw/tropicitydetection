
set(SOURCES
    geometry3.cc
    cube.cc
    trajectory.cc
   )

if(CUDA_FOUND)
    set(SOURCES ${SOURCES}
        gpu-info.cu
        cube.cu
       )
    cuda_add_library(trop STATIC ${SOURCES})
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS})
else()
    add_library(trop STATIC ${SOURCES})
endif()


add_executable(main main.cc)
target_link_libraries(main trop)


if(ENABLE_PYTHON)
  add_library(tropicity SHARED interface.cc)
  target_link_libraries(tropicity trop)
endif()

install(TARGETS main DESTINATION bin)
if(ENABLE_PYTHON)
  install(TARGETS tropicity DESTINATION lib)
  install(FILES splitgrid.py assign_plane.py DESTINATION bin)
endif()
