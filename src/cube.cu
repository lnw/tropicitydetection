
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "dir-enum.hh"
#include "geometry3.hh"
#include "plane.hh"
#include "trop-enum.hh"

#include <cuda.h>
#include <cuda_runtime_api.h>


void allocate() {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  // cudaMalloc((void**)&n_x_d, sizeof(int));
  // cudaMalloc((void**)&n_y_d, sizeof(int));
  // cudaMalloc((void**)&n_z_d, sizeof(int));
  // cudaMalloc((void**)&field_d, size() * sizeof(coord3d));
  // cudaMalloc((void**)&origin_d, sizeof(coord3d));
  // cudaMalloc((void**)&spacing_d, sizeof(coord3d));
}


void upload() {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  // cudaMemcpy(n_x_d, &n_x, sizeof(int), cudaMemcpyHostToDevice);
  // cudaMemcpy(n_y_d, &n_y, sizeof(int), cudaMemcpyHostToDevice);
  // cudaMemcpy(n_z_d, &n_z, sizeof(int), cudaMemcpyHostToDevice);
  // cudaMemcpy(field_d, field.data(), size() * sizeof(coord3d), cudaMemcpyHostToDevice);
  // cudaMemcpy(origin_d, &origin, sizeof(coord3d), cudaMemcpyHostToDevice);
  // cudaMemcpy(spacing_d, &spacing, sizeof(coord3d), cudaMemcpyHostToDevice);
}


void download() {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  // CudaMemcpy( n_x_d, &n_x, sizeof(int) , cudaMemcpyDeviceToHost);
  // CudaMemcpy( n_y_d, &n_y, sizeof(int) , cudaMemcpyDeviceToHost);
  // CudaMemcpy( n_z_d, &n_z, sizeof(int) , cudaMemcpyDeviceToHost);
  // CudaMemcpy( field_d, field.data(), size() * sizeof(coord3d) , cudaMemcpyDeviceToHost);
  // CudaMemcpy( origin_d, &origin,  sizeof(coord3d) , cudaMemcpyDeviceToHost);
  // CudaMemcpy( spacing_d, &spacing,  sizeof(coord3d) , cudaMemcpyDeviceToHost);
}

void deallocate() {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  //cudaFree(n_x_d);
  //cudaFree(n_y_d);
  //cudaFree(n_z_d);
  // cudaFree(field_d);
  //cudaFree(origin_d);
  //cudaFree(spacing_d);
}


__global__ void classify_points_kernel(Tropicity* res_d, coord3d* field_d, int n_x, int n_y, int n_z) {
  // std::cout << __PRETTY_FUNCTION__ << std::endl;
 printf("hello from the gpu\n");
}


std::vector<Tropicity> classify_points_cudax(const std::vector<coord3d>& field, int nx, int ny, int nz, coord3d origin, coord3d spacing,
                                             const std::vector<coord3d>& coords, Direction bfielddir) {
  int64_t n = coords.size();
  std::vector<Tropicity> res(n);
  coord3d* field_d;
  Tropicity* res_d;

  cudaMalloc((void**)&field_d, n * sizeof(coord3d));
  cudaMalloc((void**)&res_d, n * sizeof(Tropicity));

  cudaMemcpy(field_d, field.data(), n * sizeof(coord3d), cudaMemcpyHostToDevice);

  classify_points_kernel<<<3, 3>>>(res_d, field_d, nx, ny, nz);

  cudaMemcpy(res.data(), res_d, n * sizeof(Tropicity), cudaMemcpyDeviceToHost);

  cudaFree(field_d);
  cudaFree(res_d);
  return res;
}


