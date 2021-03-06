#include <cstdio>
#include <cstdlib>
#include <cuda.h>
#include <cuda_runtime_api.h>


void print_device_props_short() {
  const size_t kb = 1024;
  const size_t mb = kb * kb;

  int devCount;
  cudaGetDeviceCount(&devCount);
  printf("Found the following GPUs:\n");

  for (int i = 0; i < devCount; ++i) {
    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, i);
    printf("  | %d: %s\n", i, props.name);
    printf("  | arch version / compute capability: %d.%d\n", props.major, props.minor);
    printf("  | Global memory: %zu MB\n", props.totalGlobalMem / mb);
  }
}


void print_device_props_complete() {
  const int kb = 1024;
  const int mb = kb * kb;

  int devCount;
  cudaGetDeviceCount(&devCount);
  printf("Found the following GPUs:\n");

  for (int i = 0; i < devCount; ++i) {
    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, i);
    printf("  | %d: %s\n", i, props.name);
    printf("  | arch version / compute capability: %d.%d\n", props.major, props.minor);
    printf("  | global memory: %zu MB\n", props.totalGlobalMem / mb);
    printf("  | shared memory: %zu KB\n", props.sharedMemPerBlock / kb);
    printf("  | constant memory: %zu KB\n", props.totalConstMem / kb);
    printf("  | 32b-registers per block: %d\n", props.regsPerBlock);
    printf("  | warp size: %d\n", props.warpSize);
    printf("  | max pitch: %zu KB\n", props.memPitch / kb);
    printf("  | threads per block: %d\n", props.maxThreadsPerBlock);
    printf("  | max block dimensions: %d, %d, %d\n", props.maxThreadsDim[0], props.maxThreadsDim[1], props.maxThreadsDim[2]);
    printf("  | max grid dimensions: %d, %d, %d\n", props.maxGridSize[0], props.maxGridSize[1], props.maxGridSize[2]);
    printf("  | clock rate: %d Hz\n", props.clockRate);
  }
}


int number_cuda_devices() {
  int devCount;
  cudaGetDeviceCount(&devCount);
  return devCount;
}


int number_cuda_devices_minimum_mem_mb(size_t min_mem) {
  const size_t kb = 1024;
  const size_t mb = kb * kb;

  int devCount;
  int count = 0;
  cudaGetDeviceCount(&devCount);
  for (int i = 0; i < devCount; ++i) {
    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, i);
    if (props.totalGlobalMem / mb >= min_mem)
      count++;
  }

  return count;
}


int number_cuda_devices_minimum_cc(int min_cc) {
  int devCount;
  int count = 0;
  cudaGetDeviceCount(&devCount);
  for (int i = 0; i < devCount; ++i) {
    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, i);
    if (props.major * 10 + props.minor >= min_cc)
      count++;
  }

  return count;
}
