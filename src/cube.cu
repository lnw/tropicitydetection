#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "dir-enum.hh"
#include "geometry3_d.hh"
#include "plane.hh"
#include "trop-enum.hh"

#include <cuda.h>
#include <cuda_runtime_api.h>
// #include <vector_functions.h>


__device__ inline double3 operator*(const double d, const double3 d3) {
  const double xval = d * d3.x;
  const double yval = d * d3.y;
  const double zval = d * d3.z;
  return make_double3(xval, yval, zval);
}

__device__ inline double3 operator+(const double3 d3a, const double3 d3b) {
  const double xval = d3a.x * d3b.x;
  const double yval = d3a.y * d3b.y;
  const double zval = d3a.z * d3b.z;
  return make_double3(xval, yval, zval);
}


__device__ bool outofbounds(const coord3d_d& pos, int64_t nx, int64_t ny, int64_t nz) {
  if (pos.x[0] >= nx - 1 || pos.x[1] >= ny - 1 || pos.x[2] >= nz - 1 ||
      pos.x[0] < 0 || pos.x[1] < 0 || pos.x[2] < 0) {
    return true;
  }
  return false;
}


// trilinear interpolation
__device__ bool getvector_v1(const coord3d_d& pos, const cudaPitchedPtr field_d, const int64_t nx, const int64_t ny, const int64_t nz, coord3d_d& res_vec) {
  // int indx = threadIdx.x + blockIdx.x * blockDim.x;
  // if (indx < 2) {
  //  printf("x %f, y %f, z %f \n", pos[0], pos[1], pos[2]);
  // }

  char* ptr = (char*)field_d.ptr;
  // pitch being xdim * sizeof(thing) rounded up to a multiple of 32
  size_t pitch = field_d.pitch;

  if (outofbounds(pos, nx, ny, nz))
    return false;

  double x = pos[0];
  double y = pos[1];
  double z = pos[2];
  int x0 = int(floor(pos[0]));
  int x1 = x0 + 1;
  int y0 = int(floor(pos[1]));
  int y1 = y0 + 1;
  int z0 = int(floor(pos[2]));
  int z1 = z0 + 1;
  coord3d_d v000 = *(coord3d_d*)(ptr + (pitch * ny * z0 + pitch * y0 + x0 * sizeof(coord3d_d)));
  coord3d_d v001 = *(coord3d_d*)(ptr + (pitch * ny * z1 + pitch * y0 + x0 * sizeof(coord3d_d)));
  coord3d_d v010 = *(coord3d_d*)(ptr + (pitch * ny * z0 + pitch * y1 + x0 * sizeof(coord3d_d)));
  coord3d_d v011 = *(coord3d_d*)(ptr + (pitch * ny * z1 + pitch * y1 + x0 * sizeof(coord3d_d)));
  coord3d_d v100 = *(coord3d_d*)(ptr + (pitch * ny * z0 + pitch * y0 + x1 * sizeof(coord3d_d)));
  coord3d_d v101 = *(coord3d_d*)(ptr + (pitch * ny * z1 + pitch * y0 + x1 * sizeof(coord3d_d)));
  coord3d_d v110 = *(coord3d_d*)(ptr + (pitch * ny * z0 + pitch * y1 + x1 * sizeof(coord3d_d)));
  coord3d_d v111 = *(coord3d_d*)(ptr + (pitch * ny * z1 + pitch * y1 + x1 * sizeof(coord3d_d)));
  coord3d_d aux0 = (x1 - x) * v000 + (x - x0) * v100;
  coord3d_d aux1 = (x1 - x) * v010 + (x - x0) * v110;
  coord3d_d aux2 = (x1 - x) * v001 + (x - x0) * v101;
  coord3d_d aux3 = (x1 - x) * v011 + (x - x0) * v111;
  coord3d_d aux4 = (y1 - y) * aux0 + (y - y0) * aux1;
  coord3d_d aux5 = (y1 - y) * aux2 + (y - y0) * aux3;
  res_vec = (z1 - z) * aux4 + (z - z0) * aux5;
  return true;
}


__device__ inline coord3d_d get_val_v2(const int x, const int y, const int z,
                                       cudaTextureObject_t field_x, cudaTextureObject_t field_y, cudaTextureObject_t field_z) {
  int2 xvalint2 = tex3D<int2>(field_x, x, y, z);
  double xval = __hiloint2double(xvalint2.y, xvalint2.x);
  int2 yvalint2 = tex3D<int2>(field_y, x, y, z);
  double yval = __hiloint2double(yvalint2.y, yvalint2.x);
  int2 zvalint2 = tex3D<int2>(field_z, x, y, z);
  double zval = __hiloint2double(zvalint2.y, zvalint2.x);
  return coord3d_d(xval, yval, zval);
}


// trilinear interpolation
__device__ bool getvector_v2(const coord3d_d& pos,
                             cudaTextureObject_t field_x, cudaTextureObject_t field_y, cudaTextureObject_t field_z,
                             const int64_t nx, const int64_t ny, const int64_t nz, coord3d_d& res_vec) {
  // int indx = threadIdx.x + blockIdx.x * blockDim.x;
  // if (indx < 2) {
  //  printf("x %f, y %f, z %f \n", pos[0], pos[1], pos[2]);
  // }

  if (outofbounds(pos, nx, ny, nz))
    return false;

  double x = pos[0];
  double y = pos[1];
  double z = pos[2];
  int x0 = int(floor(pos[0]));
  int x1 = x0 + 1;
  int y0 = int(floor(pos[1]));
  int y1 = y0 + 1;
  int z0 = int(floor(pos[2]));
  int z1 = z0 + 1;
  coord3d_d v000 = get_val_v2(x0, y0, z0, field_x, field_y, field_z);
  coord3d_d v001 = get_val_v2(x0, y0, z1, field_x, field_y, field_z);
  coord3d_d v010 = get_val_v2(x0, y1, z0, field_x, field_y, field_z);
  coord3d_d v011 = get_val_v2(x0, y1, z1, field_x, field_y, field_z);
  coord3d_d v100 = get_val_v2(x1, y0, z0, field_x, field_y, field_z);
  coord3d_d v101 = get_val_v2(x1, y0, z1, field_x, field_y, field_z);
  coord3d_d v110 = get_val_v2(x1, y1, z0, field_x, field_y, field_z);
  coord3d_d v111 = get_val_v2(x1, y1, z1, field_x, field_y, field_z);
  coord3d_d aux0 = (x1 - x) * v000 + (x - x0) * v100;
  coord3d_d aux1 = (x1 - x) * v010 + (x - x0) * v110;
  coord3d_d aux2 = (x1 - x) * v001 + (x - x0) * v101;
  coord3d_d aux3 = (x1 - x) * v011 + (x - x0) * v111;
  coord3d_d aux4 = (y1 - y) * aux0 + (y - y0) * aux1;
  coord3d_d aux5 = (y1 - y) * aux2 + (y - y0) * aux3;
  res_vec = (z1 - z) * aux4 + (z - z0) * aux5;
  return true;
}


// trilinear interpolation
__device__ inline bool getvector_v3(const coord3d_d& pos, cudaTextureObject_t field_d,
                                    const int64_t nx, const int64_t ny, const int64_t nz, coord3d_d& res_vec) {
  // int indx = threadIdx.x + blockIdx.x * blockDim.x;
  // if (indx < 2) {
  //  printf("x %f, y %f, z %f \n", pos[0], pos[1], pos[2]);
  // }

  if (outofbounds(pos, nx, ny, nz))
    return false;

  double x = pos[0];
  double y = pos[1];
  double z = pos[2];

  float4 val = tex3D<float4>(field_d, x, y, z);
  res_vec = coord3d_d(val.x, val.y, val.z);
  return true;
}


// Runge-Kutta method, 4th order
// c --> positions, k --> vectors at c
__device__ bool extend_rungekutta_v1(const cudaPitchedPtr field_d, const int64_t nx, const int64_t ny, const int64_t nz,
                                     const coord3d_d& prevpos, float step_length, coord3d_d& newpos) {
  // int indx = threadIdx.x + blockIdx.x * blockDim.x;
  coord3d_d c0 = prevpos;
  coord3d_d k0;
  bool good = getvector_v1(c0, field_d, nx, ny, nz, k0);
  k0 = k0.normalised() * step_length;

  const coord3d_d c1 = c0 + k0 * 0.5;
  coord3d_d k1;
  good = getvector_v1(c1, field_d, nx, ny, nz, k1);
  if (!good)
    return false;
  k1 = k1.normalised() * step_length;

  const coord3d_d c2 = c0 + k1 * 0.5;
  coord3d_d k2;
  good = getvector_v1(c2, field_d, nx, ny, nz, k2);
  if (!good)
    return false;
  k2 = k2.normalised() * step_length;

  const coord3d_d c3 = c0 + k2;
  coord3d_d k3;
  good = getvector_v1(c3, field_d, nx, ny, nz, k3);
  if (!good)
    return false;
  k3 = k3.normalised() * step_length;

  const coord3d_d c4 = c0 + (k0 + k1 * 2.0 + k2 * 2.0 + k3) / 6.0;
  coord3d_d k4;
  good = getvector_v1(c4, field_d, nx, ny, nz, k4);
  if (!good)
    return false;
  newpos = c4;
  return true;
}


// Runge-Kutta method, 4th order
// c --> positions, k --> vectors at c
__device__ bool extend_rungekutta_v2(const cudaTextureObject_t field_x, const cudaTextureObject_t field_y, const cudaTextureObject_t field_z,
                                     const int64_t nx, const int64_t ny, const int64_t nz,
                                     const coord3d_d& prevpos, float step_length, coord3d_d& newpos) {
  // int indx = threadIdx.x + blockIdx.x * blockDim.x;
  coord3d_d c0 = prevpos;
  coord3d_d k0;
  bool good = getvector_v2(c0, field_x, field_y, field_z, nx, ny, nz, k0);
  k0 = k0.normalised() * step_length;

  const coord3d_d c1 = c0 + k0 * 0.5;
  coord3d_d k1;
  good = getvector_v2(c1, field_x, field_y, field_z, nx, ny, nz, k1);
  if (!good)
    return false;
  k1 = k1.normalised() * step_length;

  const coord3d_d c2 = c0 + k1 * 0.5;
  coord3d_d k2;
  good = getvector_v2(c2, field_x, field_y, field_z, nx, ny, nz, k2);
  if (!good)
    return false;
  k2 = k2.normalised() * step_length;

  const coord3d_d c3 = c0 + k2;
  coord3d_d k3;
  good = getvector_v2(c3, field_x, field_y, field_z, nx, ny, nz, k3);
  if (!good)
    return false;
  k3 = k3.normalised() * step_length;

  const coord3d_d c4 = c0 + (k0 + k1 * 2.0 + k2 * 2.0 + k3) / 6.0;
  coord3d_d k4;
  good = getvector_v2(c4, field_x, field_y, field_z, nx, ny, nz, k4);
  if (!good)
    return false;
  newpos = c4;
  return true;
}


// Runge-Kutta method, 4th order
// c --> positions, k --> vectors at c
__device__ bool extend_rungekutta_v3(const cudaTextureObject_t field_d,
                                     const int64_t nx, const int64_t ny, const int64_t nz,
                                     const coord3d_d& prevpos, float step_length, coord3d_d& newpos) {
  // int indx = threadIdx.x + blockIdx.x * blockDim.x;
  coord3d_d c0 = prevpos;
  coord3d_d k0;
  bool good = getvector_v3(c0, field_d, nx, ny, nz, k0);
  k0 = k0.normalised() * step_length;

  const coord3d_d c1 = c0 + k0 * 0.5;
  coord3d_d k1;
  good = getvector_v3(c1, field_d, nx, ny, nz, k1);
  if (!good)
    return false;
  k1 = k1.normalised() * step_length;

  const coord3d_d c2 = c0 + k1 * 0.5;
  coord3d_d k2;
  good = getvector_v3(c2, field_d, nx, ny, nz, k2);
  if (!good)
    return false;
  k2 = k2.normalised() * step_length;

  const coord3d_d c3 = c0 + k2;
  coord3d_d k3;
  good = getvector_v3(c3, field_d, nx, ny, nz, k3);
  if (!good)
    return false;
  k3 = k3.normalised() * step_length;

  const coord3d_d c4 = c0 + (k0 + k1 * 2.0 + k2 * 2.0 + k3) / 6.0;
  coord3d_d k4;
  good = getvector_v3(c4, field_d, nx, ny, nz, k4);
  if (!good)
    return false;
  newpos = c4;
  return true;
}


__device__ void complete_trajectory_v1(const cudaPitchedPtr field_d, const int nx, const int ny, const int nz,
                                       coord3d_d* __restrict__ positions, int& index, int max_points_traj,
                                       float return_ratio, float step_length, bool& out_of_bounds) {
  // const int indx = threadIdx.x + blockIdx.x * blockDim.x;
  out_of_bounds = false;
  double dist2farthest = -1; // if this is set at 0 at declaration, the following while loop will never run
  if (index > 0) {
    for (int i = 0; i <= index; i++)
      dist2farthest = std::max(dist2farthest, (positions[i] - positions[0]).norm());
  }

  // if we get to a point that is less than return_ratio of the longest distance in the trajectory
  while ((positions[index] - positions[0]).norm() > return_ratio * dist2farthest) {
    if (!extend_rungekutta_v1(field_d, nx, ny, nz,
                              positions[index], step_length, positions[index + 1])) {
      out_of_bounds = true;
      // printf("%d: %d oob\n", indx, index);
      return;
    }
    index++;

    dist2farthest = std::max(dist2farthest, (positions[index] - positions[0]).norm());

    if (index == max_points_traj - 2) {
      step_length *= 1.5;
      index = 0;
      dist2farthest = -1;
    }
  }
  // printf("%d: %d\n", indx, index);
}


__device__ void complete_trajectory_v2(const cudaTextureObject_t field_x, const cudaTextureObject_t field_y, const cudaTextureObject_t field_z,
                                       const int nx, const int ny, const int nz,
                                       coord3d_d* __restrict__ positions, int& index, int max_points_traj,
                                       float return_ratio, float step_length, bool& out_of_bounds) {
  // const int indx = threadIdx.x + blockIdx.x * blockDim.x;
  out_of_bounds = false;
  double dist2farthest = -1; // if this is set at 0 at declaration, the following while loop will never run
  if (index > 0) {
    for (int i = 0; i <= index; i++)
      dist2farthest = std::max(dist2farthest, (positions[i] - positions[0]).norm());
  }

  // if we get to a point that is less than return_ratio of the longest distance in the trajectory
  while ((positions[index] - positions[0]).norm() > return_ratio * dist2farthest) {
    if (!extend_rungekutta_v2(field_x, field_y, field_z, nx, ny, nz,
                              positions[index], step_length, positions[index + 1])) {
      out_of_bounds = true;
      // printf("%d: %d oob\n", indx, index);
      return;
    }
    index++;

    dist2farthest = std::max(dist2farthest, (positions[index] - positions[0]).norm());

    if (index == max_points_traj - 2) {
      step_length *= 1.5;
      index = 0;
      dist2farthest = -1;
    }
  }
  // printf("%d: %d\n", indx, index);
}


__device__ void complete_trajectory_v3(const cudaTextureObject_t field_d,
                                       const int nx, const int ny, const int nz,
                                       coord3d_d* __restrict__ positions, int& index, int max_points_traj,
                                       float return_ratio, float step_length, bool& out_of_bounds) {
  // const int indx = threadIdx.x + blockIdx.x * blockDim.x;
  out_of_bounds = false;
  double dist2farthest = -1; // if this is set at 0 at declaration, the following while loop will never run
  if (index > 0) {
    for (int i = 0; i <= index; i++)
      dist2farthest = std::max(dist2farthest, (positions[i] - positions[0]).norm());
  }

  // if we get to a point that is less than return_ratio of the longest distance in the trajectory
  while ((positions[index] - positions[0]).norm() > return_ratio * dist2farthest) {
    if (!extend_rungekutta_v3(field_d, nx, ny, nz,
                              positions[index], step_length, positions[index + 1])) {
      out_of_bounds = true;
      // printf("%d: %d oob\n", indx, index);
      return;
    }
    index++;

    dist2farthest = std::max(dist2farthest, (positions[index] - positions[0]).norm());

    if (index == max_points_traj - 2) {
      step_length *= 1.5;
      index = 0;
      dist2farthest = -1;
    }
  }
  // printf("%d: %d\n", indx, index);
}


__device__ Tropicity classify_trajectory(const coord3d_d* __restrict__ positions, int n_points_in_traj, Direction bfielddir, bool out_of_bounds) {
  // int indx = threadIdx.x + blockIdx.x * blockDim.x;
  // printf("p in traj %d: %d\n", indx, n_points_in_traj);
  coord3d_d bfield;
  switch (bfielddir) {
    case Direction::pos_x: {
      bfield = coord3d_d(1, 0, 0);
      break;
    }
    case Direction::neg_x: {
      bfield = coord3d_d(-1, 0, 0);
      break;
    }
    case Direction::pos_y: {
      bfield = coord3d_d(0, 1, 0);
      break;
    }
    case Direction::neg_y: {
      bfield = coord3d_d(0, -1, 0);
      break;
    }
    case Direction::pos_z: {
      bfield = coord3d_d(0, 0, 1);
      break;
    }
    case Direction::neg_z: {
      bfield = coord3d_d(0, 0, -1);
      break;
    }
    default: {
      return Tropicity::input_error;
    }
  }

  if (out_of_bounds)
    return Tropicity::outofbounds;

  coord3d_d crosssum(0, 0, 0);
  for (size_t i = 0; i < n_points_in_traj; i++) {
    crosssum += positions[(i - 1 + n_points_in_traj) % n_points_in_traj].cross(positions[i]);
  }
  // crossum += positions[positions.size()-1].cross(positions[0]);
  // if (indx < 40) printf("cross: %f/%f/%f\n", crosssum[0], crosssum[1], crosssum[2]);

  double dot_product = bfield.dot(crosssum);
  if (dot_product > 0)
    return Tropicity::paratropic;
  else if (dot_product < 0)
    return Tropicity::diatropic;
  else
    return Tropicity::unclassifyable;
}


__global__ void classify_points_kernel_v1(coord3d_d* __restrict__ points, int64_t n_points,
                                          const cudaPitchedPtr field_d, const int64_t nx, const int64_t ny, const int64_t nz,
                                          coord3d_d* __restrict__ trajectories_d, float step_length, int64_t max_points_traj,
                                          Direction bfielddir, Tropicity* __restrict__ tropicities_d) {
  const int32_t indx = threadIdx.x + blockIdx.x * blockDim.x;
  // if (indx < 2) printf("hello from the gpu: %d\n", indx);

  if (indx > n_points - 1)
    return;

  coord3d_d vec(0, 0, 0);
  // if (indx < 2) printf("pos %d %f/%f/%f\n", indx, points[indx][0], points[indx][1], points[indx][2]);
  bool good = getvector_v1(points[indx], field_d, nx, ny, nz, vec);
  // if (indx < 2) printf("found vec %d %d: %f/%f/%f\n", indx, good, vec[0], vec[1], vec[2]);
  if (!good) {
    tropicities_d[indx] = Tropicity::outofbounds;
    return;
  }

  bool out_of_bounds;
  int current_index_in_traj = 0;
  float return_ratio = 0.2;
  trajectories_d[max_points_traj * indx] = points[indx];
  complete_trajectory_v1(field_d, nx, ny, nz,
                         trajectories_d + max_points_traj * indx, current_index_in_traj, max_points_traj,
                         return_ratio, step_length, out_of_bounds);
  tropicities_d[indx] = classify_trajectory(trajectories_d + max_points_traj * indx, current_index_in_traj + 1, bfielddir, out_of_bounds);
}


__global__ void classify_points_kernel_v2(coord3d_d* __restrict__ points, int64_t n_points,
                                          cudaTextureObject_t field_x, cudaTextureObject_t field_y, cudaTextureObject_t field_z,
                                          const int64_t nx, const int64_t ny, const int64_t nz,
                                          coord3d_d* __restrict__ trajectories_d, float step_length, int64_t max_points_traj,
                                          Direction bfielddir, Tropicity* __restrict__ tropicities_d) {

  const int32_t indx = threadIdx.x + blockIdx.x * blockDim.x;
  // if (indx < 2) printf("hello from the gpu: %d\n", indx);

  if (indx > n_points - 1)
    return;

  coord3d_d vec(0, 0, 0);
  // if (indx < 2) printf("pos %d %f/%f/%f\n", indx, points[indx][0], points[indx][1], points[indx][2]);
  bool good = getvector_v2(points[indx], field_x, field_y, field_z, nx, ny, nz, vec);
  // if (indx < 2) printf("found vec %d %d: %f/%f/%f\n", indx, good, vec[0], vec[1], vec[2]);
  if (!good) {
    tropicities_d[indx] = Tropicity::outofbounds;
    return;
  }

  bool out_of_bounds;
  int current_index_in_traj = 0;
  float return_ratio = 0.2;
  trajectories_d[max_points_traj * indx] = points[indx];
  complete_trajectory_v2(field_x, field_y, field_z, nx, ny, nz,
                         trajectories_d + max_points_traj * indx, current_index_in_traj, max_points_traj,
                         return_ratio, step_length, out_of_bounds);
  tropicities_d[indx] = classify_trajectory(trajectories_d + max_points_traj * indx, current_index_in_traj + 1, bfielddir, out_of_bounds);
}


__global__ void classify_points_kernel_v3(coord3d_d* __restrict__ points, int64_t n_points,
                                          cudaTextureObject_t field_d,
                                          const int64_t nx, const int64_t ny, const int64_t nz,
                                          coord3d_d* __restrict__ trajectories_d, float step_length, int64_t max_points_traj,
                                          Direction bfielddir, Tropicity* __restrict__ tropicities_d) {
  const int32_t indx = threadIdx.x + blockIdx.x * blockDim.x;
  // if (indx < 2) printf("hello from the gpu: %d\n", indx);

  if (indx > n_points - 1) {
    return;
  }

  coord3d_d vec(0, 0, 0);
  // if (indx < 2) printf("pos %d %f/%f/%f\n", indx, points[indx][0], points[indx][1], points[indx][2]);
  bool good = getvector_v3(points[indx], field_d, nx, ny, nz, vec);
  // if (indx < 2) printf("found vec %d %d: %f/%f/%f\n", indx, good, vec[0], vec[1], vec[2]);
  if (!good) {
    tropicities_d[indx] = Tropicity::outofbounds;
    return;
  }

  bool out_of_bounds;
  int current_index_in_traj = 0;
  float return_ratio = 0.2;
  trajectories_d[max_points_traj * indx] = points[indx];
  complete_trajectory_v3(field_d, nx, ny, nz,
                         trajectories_d + max_points_traj * indx, current_index_in_traj, max_points_traj,
                         return_ratio, step_length, out_of_bounds);
  tropicities_d[indx] = classify_trajectory(trajectories_d + max_points_traj * indx, current_index_in_traj + 1, bfielddir, out_of_bounds);
}


std::vector<Tropicity> classify_points_cudax_v1(const double* field_a, const int64_t nx, const int64_t ny, const int64_t nz, double* origin_a, double* spacing_a,
                                                const double* start_points_a, int64_t n_points, Direction bfielddir) {
  // std::cout << __PRETTY_FUNCTION__ << std::endl;
#if 0
  float steplength = 0.01;
#else
  float step_length_ratio = 0.05;
  float step_length = step_length_ratio * spacing_a[0];
#endif
  int64_t max_points_traj = 10000;

  std::vector<Tropicity> res(n_points);
  coord3d_d* field = new coord3d_d[nx * ny * nz];
  coord3d_d* start_points = new coord3d_d[n_points];

  coord3d_d* start_points_d;
  Tropicity* res_d;
  coord3d_d* trajectories_d;

  for (int64_t i = 0; i < nx * ny * nz; i++)
    for (int64_t j = 0; j < 3; j++)
      field[i][j] = field_a[3 * i + j];
  for (int64_t i = 0; i < n_points; i++)
    for (int64_t j = 0; j < 3; j++)
      start_points[i][j] = start_points_a[3 * i + j];

  cudaPitchedPtr field_d;
  cudaExtent field_extent = make_cudaExtent(nx * sizeof(coord3d_d), ny, nz);
  cudaMalloc3D(&field_d, field_extent);

  cudaMemcpy3DParms memCopyParameters = {0};
  memCopyParameters.srcPtr = make_cudaPitchedPtr(field, nx * sizeof(coord3d_d), ny, nz);
  memCopyParameters.dstPtr = field_d;
  memCopyParameters.extent = field_extent;
  memCopyParameters.kind = cudaMemcpyHostToDevice;

  cudaMemcpy3DAsync(&memCopyParameters, 0);
  // printf("nx, ny, nz %lu, %lu, %lu\n", field_d.pitch, field_d.xsize, field_d.ysize);

  // alloc
  cudaMalloc((void**)&start_points_d, n_points * sizeof(coord3d_d));
  cudaMalloc((void**)&trajectories_d, n_points * max_points_traj * sizeof(coord3d_d));
  cudaMalloc((void**)&res_d, n_points * sizeof(Tropicity));

  // copy to device
  cudaMemcpy(start_points_d, start_points, n_points * sizeof(coord3d_d), cudaMemcpyHostToDevice);
  // cout << "e " << cudaGetLastError() << endl;

  int block_size = 256;
  int grid_size = n_points / block_size + (n_points % block_size != 0);
  // std::cout << "points / gridsize / blocksize: " << n_points << ", " << grid_size << ", " << block_size << std::endl;
  classify_points_kernel_v1<<<grid_size, block_size>>>(start_points_d, n_points,
                                                       field_d, nx, ny, nz,
                                                       trajectories_d, step_length, max_points_traj,
                                                       bfielddir, res_d);
  // cout << "e " << cudaGetLastError() << endl;

  // copy from device
  cudaMemcpy(res.data(), res_d, n_points * sizeof(Tropicity), cudaMemcpyDeviceToHost);

  // dealloc
  cudaFree(field_d.ptr);
  cudaFree(start_points_d);
  cudaFree(trajectories_d);
  cudaFree(res_d);
  delete[] field;
  delete[] start_points;
  return res;
}


std::vector<Tropicity> classify_points_cudax_v2(double* field_x_a, double* field_y_a, double* field_z_a,
                                                const int64_t nx, const int64_t ny, const int64_t nz, double* origin_a, double* spacing_a,
                                                const double* start_points_a, int64_t n_points, Direction bfielddir) {
  // std::cout << __PRETTY_FUNCTION__ << std::endl;
#if 0
  float steplength = 0.01;
#else
  float step_length_ratio = 0.05;
  float step_length = step_length_ratio * spacing_a[0];
#endif
  int64_t max_points_traj = 10000;

  std::vector<Tropicity> res(n_points);
  coord3d_d* start_points = new coord3d_d[n_points];

  coord3d_d* start_points_d;
  Tropicity* res_d;
  coord3d_d* trajectories_d;

  for (int64_t i = 0; i < n_points; i++)
    for (int64_t j = 0; j < 3; j++)
      start_points[i][j] = start_points_a[3 * i + j];

  cudaArray_t field_x_d, field_y_d, field_z_d;
  // cudaChannelFormatDesc desc = cudaCreateChannelDesc<double>();
  cudaChannelFormatDesc desc = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindSigned); // we pretend to store int2 instead of double
  cudaExtent field_extent = make_cudaExtent(nx, ny, nz);
  cudaMalloc3DArray(&field_x_d, &desc, field_extent);
  cudaMalloc3DArray(&field_y_d, &desc, field_extent);
  cudaMalloc3DArray(&field_z_d, &desc, field_extent);

  cudaMemcpy3DParms memCopyParametersX = {0};
  memCopyParametersX.srcPtr = make_cudaPitchedPtr(field_x_a, nx * sizeof(double), nx, ny);
  memCopyParametersX.dstArray = field_x_d;
  memCopyParametersX.extent = field_extent;
  memCopyParametersX.kind = cudaMemcpyHostToDevice;
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  cudaMemcpy3DParms memCopyParametersY = {0};
  memCopyParametersY.srcPtr = make_cudaPitchedPtr(field_y_a, nx * sizeof(double), nx, ny);
  memCopyParametersY.dstArray = field_y_d;
  memCopyParametersY.extent = field_extent;
  memCopyParametersY.kind = cudaMemcpyHostToDevice;
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  cudaMemcpy3DParms memCopyParametersZ = {0};
  memCopyParametersZ.srcPtr = make_cudaPitchedPtr(field_z_a, nx * sizeof(double), nx, ny);
  memCopyParametersZ.dstArray = field_z_d;
  memCopyParametersZ.extent = field_extent;
  memCopyParametersZ.kind = cudaMemcpyHostToDevice;
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  cudaMemcpy3DAsync(&memCopyParametersX, 0);
  cudaMemcpy3DAsync(&memCopyParametersY, 0);
  cudaMemcpy3DAsync(&memCopyParametersZ, 0);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;
  // printf("nx, ny, nz %lu, %lu, %lu\n", field_d.pitch, field_d.xsize, field_d.ysize);

  // prepare textures
  struct cudaResourceDesc fieldXResDesc;
  memset(&fieldXResDesc, 0, sizeof(fieldXResDesc));
  fieldXResDesc.resType = cudaResourceTypeArray;
  fieldXResDesc.res.array.array = field_x_d;
  struct cudaResourceDesc fieldYResDesc;
  memset(&fieldYResDesc, 0, sizeof(fieldYResDesc));
  fieldYResDesc.resType = cudaResourceTypeArray;
  fieldYResDesc.res.array.array = field_y_d;
  struct cudaResourceDesc fieldZResDesc;
  memset(&fieldZResDesc, 0, sizeof(fieldZResDesc));
  fieldZResDesc.resType = cudaResourceTypeArray;
  fieldZResDesc.res.array.array = field_z_d;

  struct cudaTextureDesc fieldTexDesc;
  memset(&fieldTexDesc, 0, sizeof(fieldTexDesc));
  fieldTexDesc.addressMode[0] = cudaAddressModeBorder; // alternatively: wrap, clamp, mirror
  fieldTexDesc.addressMode[1] = cudaAddressModeBorder; // alternatively: wrap, clamp, mirror
  fieldTexDesc.addressMode[2] = cudaAddressModeBorder; // alternatively: wrap, clamp, mirror
  fieldTexDesc.filterMode = cudaFilterModePoint;       // ie, interpolate linearly
  fieldTexDesc.readMode = cudaReadModeElementType;
  fieldTexDesc.normalizedCoords = 0;

  cudaTextureObject_t fieldXTexture = 0, fieldYTexture = 0, fieldZTexture = 0;
  cudaCreateTextureObject(&fieldXTexture, &fieldXResDesc, &fieldTexDesc, nullptr);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;
  cudaCreateTextureObject(&fieldYTexture, &fieldYResDesc, &fieldTexDesc, nullptr);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;
  cudaCreateTextureObject(&fieldZTexture, &fieldZResDesc, &fieldTexDesc, nullptr);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;


  // alloc
  cudaMalloc((void**)&start_points_d, n_points * sizeof(coord3d_d));
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;
  cudaMalloc((void**)&trajectories_d, n_points * max_points_traj * sizeof(coord3d_d));
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;
  cudaMalloc((void**)&res_d, n_points * sizeof(Tropicity));
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  // copy to device
  cudaMemcpy(start_points_d, start_points, n_points * sizeof(coord3d_d), cudaMemcpyHostToDevice);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  int block_size = 256;
  int grid_size = n_points / block_size + (n_points % block_size != 0);
  // std::cout << "points / gridsize / blocksize: " << n_points << ", " << grid_size << ", " << block_size << std::endl;
  classify_points_kernel_v2<<<grid_size, block_size>>>(start_points_d, n_points,
                                                       fieldXTexture, fieldYTexture, fieldZTexture, nx, ny, nz,
                                                       trajectories_d, step_length, max_points_traj,
                                                       bfielddir, res_d);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  // copy from device
  cudaMemcpy(res.data(), res_d, n_points * sizeof(Tropicity), cudaMemcpyDeviceToHost);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  cudaDestroyTextureObject(fieldXTexture);
  cudaDestroyTextureObject(fieldYTexture);
  cudaDestroyTextureObject(fieldZTexture);

  // dealloc
  cudaFree(start_points_d);
  cudaFree(trajectories_d);
  cudaFree(res_d);
  delete[] start_points;
  return res;
}


std::vector<Tropicity> classify_points_cudax_v3(float* field_x_a, float* field_y_a, float* field_z_a,
                                                const int64_t nx, const int64_t ny, const int64_t nz, double* origin_a, double* spacing_a,
                                                const double* start_points_a, int64_t n_points, Direction bfielddir) {
  // std::cout << __PRETTY_FUNCTION__ << std::endl;
#if 0
  float steplength = 0.01;
#else
  float step_length_ratio = 0.05;
  float step_length = step_length_ratio * spacing_a[0];
#endif
  int64_t max_points_traj = 10000;

  std::vector<Tropicity> res(n_points);
  coord3d_d* start_points = new coord3d_d[n_points];

  coord3d_d* start_points_d;
  Tropicity* res_d;
  coord3d_d* trajectories_d;

  for (int64_t i = 0; i < n_points; i++)
    for (int64_t j = 0; j < 3; j++)
      start_points[i][j] = start_points_a[3 * i + j];

  float4* field_float4 = new float4[nx * ny * nz];
  for (int64_t i = 0; i < nx * ny * nz; i++) {
    field_float4[i].x = field_x_a[i];
    field_float4[i].y = field_y_a[i];
    field_float4[i].z = field_z_a[i];
  }

  cudaArray_t field_d;
  cudaChannelFormatDesc desc = cudaCreateChannelDesc<float4>();
  // cudaChannelFormatDesc desc = cudaCreateChannelDesc(32, 32, 32, 0, cudaChannelFormatKindFloat);
  cudaExtent field_extent = make_cudaExtent(nx, ny, nz);
  cudaMalloc3DArray(&field_d, &desc, field_extent);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  cudaMemcpy3DParms memCopyParameters = {0};
  memCopyParameters.srcPtr = make_cudaPitchedPtr(field_float4, nx * sizeof(float4), nx, ny);
  memCopyParameters.dstArray = field_d;
  memCopyParameters.extent = field_extent;
  memCopyParameters.kind = cudaMemcpyHostToDevice;
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  cudaMemcpy3DAsync(&memCopyParameters, 0);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;
  // printf("nx, ny, nz %lu, %lu, %lu\n", field_d.pitch, field_d.xsize, field_d.ysize);

  // prepare textures
  struct cudaResourceDesc fieldResDesc;
  memset(&fieldResDesc, 0, sizeof(fieldResDesc));
  fieldResDesc.resType = cudaResourceTypeArray;
  fieldResDesc.res.array.array = field_d;

  struct cudaTextureDesc fieldTexDesc;
  memset(&fieldTexDesc, 0, sizeof(fieldTexDesc));
  fieldTexDesc.addressMode[0] = cudaAddressModeBorder; // alternatively: wrap, clamp, mirror
  fieldTexDesc.addressMode[1] = cudaAddressModeBorder; // alternatively: wrap, clamp, mirror
  fieldTexDesc.addressMode[2] = cudaAddressModeBorder; // alternatively: wrap, clamp, mirror
  fieldTexDesc.filterMode = cudaFilterModeLinear;      // ie, interpolate linearly
  fieldTexDesc.readMode = cudaReadModeElementType;
  fieldTexDesc.normalizedCoords = 0;

  cudaTextureObject_t fieldTexture = 0;
  cudaCreateTextureObject(&fieldTexture, &fieldResDesc, &fieldTexDesc, nullptr);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;


  // alloc
  cudaMalloc((void**)&start_points_d, n_points * sizeof(coord3d_d));
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;
  cudaMalloc((void**)&trajectories_d, n_points * max_points_traj * sizeof(coord3d_d));
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;
  cudaMalloc((void**)&res_d, n_points * sizeof(Tropicity));
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  // copy to device
  cudaMemcpy(start_points_d, start_points, n_points * sizeof(coord3d_d), cudaMemcpyHostToDevice);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  int block_size = 256;
  int grid_size = n_points / block_size + (n_points % block_size != 0);
  // std::cout << "points / gridsize / blocksize: " << n_points << ", " << grid_size << ", " << block_size << std::endl;
  classify_points_kernel_v3<<<grid_size, block_size>>>(start_points_d, n_points,
                                                       fieldTexture, nx, ny, nz,
                                                       trajectories_d, step_length, max_points_traj,
                                                       bfielddir, res_d);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  // copy from device
  cudaMemcpy(res.data(), res_d, n_points * sizeof(Tropicity), cudaMemcpyDeviceToHost);
  // cout << "e " << cudaGetErrorName(cudaGetLastError()) << endl;

  cudaDestroyTextureObject(fieldTexture);

  // dealloc
  cudaFree(start_points_d);
  cudaFree(trajectories_d);
  cudaFree(res_d);
  delete[] start_points;
  return res;
}
