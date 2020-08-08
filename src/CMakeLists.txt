
set(SOURCES
  cube.cc
  cube.cu
  geometry3.cc
  gpu-info.cu
  trajectory.cc
)

add_library(trop SHARED ${SOURCES})

if(ENABLE_CUDA)
  set(CUDA_NVCC_ARCH "52;60") # build for compute capabilties 5.2 and 6.0
  foreach (ARCH ${CUDA_NVCC_ARCH})
    target_compile_options(trop PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: "SHELL:-gencode arch=compute_${ARCH},code=sm_${ARCH}" >)
  endforeach()

  target_compile_options(trop PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>) 

  if(CMAKE_BUILD_TYPE MATCHES DEBUG)
    target_compile_options(trop PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>) 
  endif()
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

