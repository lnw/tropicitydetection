
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


__device__ bool outofbounds(const coord3d_d& pos, int n_x, int n_y, int n_z) {
  if (pos.x[0] >= n_x - 1 || pos.x[1] >= n_y - 1 || pos.x[2] >= n_z - 1 ||
      pos.x[0] < 0 || pos.x[1] < 0 || pos.x[2] < 0) {
    return true;
  }
  return false;
}


// trilinear interpolation
__device__ bool getvector(const coord3d_d& pos, coord3d_d* field, int n_x, int n_y, int n_z, coord3d_d& res_vec) {
  if (outofbounds(pos, n_x, n_y, n_z))
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
  coord3d_d v000 = field[n_x * n_y * z0 + n_x * y0 + x0];
  coord3d_d v001 = field[n_x * n_y * z1 + n_x * y0 + x0];
  coord3d_d v010 = field[n_x * n_y * z0 + n_x * y1 + x0];
  coord3d_d v011 = field[n_x * n_y * z1 + n_x * y1 + x0];
  coord3d_d v100 = field[n_x * n_y * z0 + n_x * y0 + x1];
  coord3d_d v101 = field[n_x * n_y * z1 + n_x * y0 + x1];
  coord3d_d v110 = field[n_x * n_y * z0 + n_x * y1 + x1];
  coord3d_d v111 = field[n_x * n_y * z1 + n_x * y1 + x1];
  coord3d_d aux0 = (x1 - x) * v000 + (x - x0) * v100;
  coord3d_d aux1 = (x1 - x) * v010 + (x - x0) * v110;
  coord3d_d aux2 = (x1 - x) * v001 + (x - x0) * v101;
  coord3d_d aux3 = (x1 - x) * v011 + (x - x0) * v111;
  coord3d_d aux4 = (y1 - y) * aux0 + (y - y0) * aux1;
  coord3d_d aux5 = (y1 - y) * aux2 + (y - y0) * aux3;
  coord3d_d res = (z1 - z) * aux4 + (z - z0) * aux5;
  return true;
}


// Runge-Kutta method, 4th order
// c --> positions, k --> vectors at c
__device__ bool extend_rungekutta(coord3d_d* field, int n_x, int n_y, int n_z,
                                  const coord3d_d& prevpos, double step_length, coord3d_d& newpos) {
  coord3d_d c0 = prevpos;

  coord3d_d k0;
  bool good = getvector(c0, field, n_x, n_y, n_z, k0);
  k0 = k0.normalised() * step_length;

  const coord3d_d c1 = c0 + k0 * 0.5;

  coord3d_d k1;
  good = getvector(c1, field, n_x, n_y, n_z, k1);
  if (!good)
    return false;
  k1 = k1.normalised() * step_length;

  const coord3d_d c2 = c0 + k1 * 0.5;

  coord3d_d k2;
  good = getvector(c2, field, n_x, n_y, n_z, k2);
  if (!good)
    return false;
  k2 = k2.normalised() * step_length;

  const coord3d_d c3 = c0 + k2;

  coord3d_d k3;
  good = getvector(c3, field, n_x, n_y, n_z, k3);
  if (!good)
    return false;
  k3 = k3.normalised() * step_length;

  const coord3d_d c4 = c0 + (k0 + k1 * 2.0 + k2 * 2.0 + k3) / 6.0;

  coord3d_d k4;
  good = getvector(c4, field, n_x, n_y, n_z, k4);
  if (!good)
    return false;
  newpos = c4;
  return true;
}


__device__ void complete_trajectory(coord3d_d* field_d, int n_x, int n_y, int n_z,
                                    coord3d_d* positions, int& index, int max_points_traj,
                                    float return_ratio, float step_length, bool& out_of_bounds) {

  double dist2farthest = -1; // if this is set at 0 at declaration, the following while loop will never run
  if (index > 0) {
    for (int i = 0; i <= index; i++)
      dist2farthest = std::max(dist2farthest, (positions[i] - positions[0]).norm());
  }

  // if we get to a point that is less return_ratio of the longest distance in the trajectory
  while ((positions[index] - positions[0]).norm() > return_ratio * dist2farthest) {
    if (!extend_rungekutta(field_d, n_x, n_y, n_z,
                           positions[index], step_length, positions[index + 1])) {
      out_of_bounds = false;
      return;
    }
    index++;

    dist2farthest = std::max(dist2farthest, (positions[index] - positions[0]).norm());

    if (index >= max_points_traj) {
      step_length *= 1.5;
      index = 0;
      dist2farthest = -1;
    }
  }
  out_of_bounds = true;
}


__device__ Tropicity classify_trajectory(coord3d_d* positions, int n_points, Direction bfielddir, bool out_of_bounds) {
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
      //       cerr << "bfielddir value wasn't 0-5.\n";
      return Tropicity::input_error;
    }
  }

  if (out_of_bounds)
    return Tropicity::outofbounds;

  coord3d_d crossum(0, 0, 0);
  for (size_t i = 0; i < n_points; i++) {
    crossum += positions[(i - 1 + n_points) % n_points].cross(positions[i]);
  }
  // crossum += positions[positions.size()-1].cross(positions[0]);

  double dot_product = bfield.dot(crossum);
  if (dot_product > 0)
    return Tropicity::paratropic;
  else if (dot_product < 0)
    return Tropicity::diatropic;
  else
    return Tropicity::unclassifyable;
}


__global__ void classify_points_kernel(coord3d_d* points, int32_t n_points,
                                       coord3d_d* field_d, int32_t n_x, int32_t n_y, int32_t n_z,
                                       coord3d_d* trajectories_d, float step_length, int max_points_traj,
                                       Direction bfielddir, Tropicity* tropicities_d) {
  printf("hello from the gpu\n");

  for (int64_t i = 0; i < n_points; i++) {
    coord3d_d vec;
    bool good = getvector(points[i], field_d, n_x, n_y, n_z, vec);
    if (!good) {
      tropicities_d[i] = Tropicity::outofbounds;
      continue;
    }

    bool out_of_bounds;
    int n_points_in_traj = 0;
    float return_ratio = 0.2;
    complete_trajectory(field_d, n_x, n_y, n_z,
                        trajectories_d + 3 * max_points_traj * i, n_points_in_traj, max_points_traj,
                        return_ratio, step_length, out_of_bounds);
    tropicities_d[i] = classify_trajectory(trajectories_d + 3 * max_points_traj * i, max_points_traj, bfielddir, out_of_bounds);
  }
}


std::vector<Tropicity> classify_points_cudax(const double* field_a, int64_t nx, int64_t ny, int64_t nz, double* origin_a, double* spacing_a,
                                             const double* coords_a, int n_points, Direction bfielddir) {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
#if 0
  float steplength = 0.01;
#else
  float step_length_ratio = 0.05;
  float step_length = step_length_ratio * spacing_a[0];
#endif
 int32_t max_points_traj = 10000;

  std::vector<Tropicity> res(n_points);
  coord3d_d* field = new coord3d_d[nx * ny * nz];

  coord3d_d* field_d;
  coord3d_d* points_d;
  Tropicity* res_d;
  coord3d_d* trajectories_d;

  for (int i = 0; i < n_points; i++)
    for (int j = 0; j < 3; j++)
      field[i][j] = field_a[3 * i + j];

  // alloc
  cudaMalloc((void**)&field_d, nx * ny * nz * sizeof(coord3d_d));
  cudaMalloc((void**)&points_d, n_points * sizeof(coord3d_d));
  cudaMalloc((void**)&trajectories_d, n_points * 10000 * 3 * sizeof(coord3d_d));
  cudaMalloc((void**)&res_d, n_points * sizeof(Tropicity));

  // copy to device
  cudaMemcpy(field_d, field, n_points * sizeof(coord3d_d), cudaMemcpyHostToDevice);

  classify_points_kernel<<<3, 3>>>(points_d, n_points,
                                   field_d, nx, ny, nz,
                                   trajectories_d, step_length, max_points_traj,
 bfielddir, res_d);

  // copy from device
  cudaMemcpy(res.data(), res_d, n_points * sizeof(Tropicity), cudaMemcpyDeviceToHost);

  // dealloc
  cudaFree(field_d);
  cudaFree(points_d);
  cudaFree(trajectories_d);
  cudaFree(res_d);
  delete[] field;
  return res;
}
