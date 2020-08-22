#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
// #include <optional>
#include <vector>

#include "dir-enum.hh"
#include "geometry3.hh"
#include "plane.hh"
#include "trop-enum.hh"


std::vector<Tropicity> classify_points_cudax_v1(const double* field, int64_t nx, int64_t ny, int64_t nz,
                                                double* origin, double* spacing, const double* coords, int64_t n_points, Direction bfielddir);
std::vector<Tropicity> classify_points_cudax_v2(double* field_x, double* field_y, double* field_z, int64_t nx, int64_t ny, int64_t nz,
                                                double* origin, double* spacing, const double* coords, int64_t n_points, Direction bfielddir);


class Cube {
private:
  int n_x;
  int n_y;
  int n_z;
  std::vector<coord3d> field;
  coord3d origin = coord3d(0, 0, 0);
  coord3d spacing = coord3d(1, 1, 1);

public:
  Cube(int x, int y, int z): n_x(x), n_y(y), n_z(z), field(n_x * n_y * n_z) {}
  Cube(std::string filename);

  const coord3d& operator()(const int x, const int y, const int z) const {
    return field[n_x * n_y * z + n_x * y + x];
  }
  coord3d& operator()(const int x, const int y, const int z) {
    return field[n_x * n_y * z + n_x * y + x];
  }

  constexpr int nx() const { return n_x; }
  constexpr int ny() const { return n_y; }
  constexpr int nz() const { return n_z; }
  size_t size() const { return field.size(); }

  constexpr coord3d get_origin() const { return origin; }
  constexpr coord3d get_spacing() const { return spacing; }

  // where position is in 'cube coordinates', ie, with origin 0, 0, 0, and unit-spacing
  // here we're excluding n_x-1 etc, even though we have values for those coordinates, because one cannot interpolate there, and making the destinction is expensive/not very useful
  constexpr bool outofbounds(const coord3d& pos) const {
    if (pos[0] >= n_x - 1 || pos[1] >= n_y - 1 || pos[2] >= n_z - 1 || pos[0] < 0 || pos[1] < 0 || pos[2] < 0) {
      return true;
    }
    return false;
  }

  std::vector<Tropicity> classify_points_cuda_v1(const std::vector<coord3d>& coords, Direction bfielddir) const {
    double* field_a = new double[3 * field.size()];
    for (size_t i = 0; i < field.size(); i++) {
      field_a[3 * i + 0] = field[i][0];
      field_a[3 * i + 1] = field[i][1];
      field_a[3 * i + 2] = field[i][2];
    }
    double* origin_a = new double[3];
    for (int i = 0; i < 3; i++)
      origin_a[i] = origin[i];
    double* spacing_a = new double[3];
    for (int i = 0; i < 3; i++)
      spacing_a[i] = spacing[i];
    double* coords_a = new double[3 * coords.size()];
    for (size_t i = 0; i < coords.size(); i++) {
      coords_a[3 * i + 0] = coords[i][0];
      coords_a[3 * i + 1] = coords[i][1];
      coords_a[3 * i + 2] = coords[i][2];
    }
    // std::cout << spacing[0] << ", " << spacing[1] << ", " << spacing[2] <<std::endl;
    // std::cout << spacing_a[0] << ", " << spacing_a[1] << ", " << spacing_a[2] <<std::endl;
    std::vector<Tropicity> res = classify_points_cudax_v1(field_a, n_x, n_y, n_z, origin_a, spacing_a, coords_a, coords.size(), bfielddir);
    delete[] field_a;
    delete[] origin_a;
    delete[] spacing_a;
    delete[] coords_a;

    // std::cout << as_underlying(res[0]) << ", " << as_underlying(res[1]) << std::endl;
    return res;
  }


  std::vector<Tropicity> classify_points_cuda_v2(const std::vector<coord3d>& coords, Direction bfielddir) const {
    // x, y, and z coordinates of the vector field
    double* field_x_a = new double[field.size()];
    double* field_y_a = new double[field.size()];
    double* field_z_a = new double[field.size()];
    for (size_t i = 0; i < field.size(); i++) {
      field_x_a[i] = field[i][0];
      field_y_a[i] = field[i][1];
      field_z_a[i] = field[i][2];
    }
    double* origin_a = new double[3];
    for (int i = 0; i < 3; i++)
      origin_a[i] = origin[i];
    double* spacing_a = new double[3];
    for (int i = 0; i < 3; i++)
      spacing_a[i] = spacing[i];
    double* coords_a = new double[3 * coords.size()];
    for (size_t i = 0; i < coords.size(); i++) {
      coords_a[3 * i + 0] = coords[i][0];
      coords_a[3 * i + 1] = coords[i][1];
      coords_a[3 * i + 2] = coords[i][2];
    }
    // std::cout << spacing[0] << ", " << spacing[1] << ", " << spacing[2] <<std::endl;
    // std::cout << spacing_a[0] << ", " << spacing_a[1] << ", " << spacing_a[2] <<std::endl;
    std::vector<Tropicity> res = classify_points_cudax_v2(field_x_a, field_y_a, field_z_a, n_x, n_y, n_z, origin_a, spacing_a, coords_a, coords.size(), bfielddir);
    delete[] field_x_a;
    delete[] field_y_a;
    delete[] field_z_a;
    delete[] origin_a;
    delete[] spacing_a;
    delete[] coords_a;

    // std::cout << as_underlying(res[0]) << ", " << as_underlying(res[1]) << std::endl;
    return res;
  }

  std::vector<Tropicity> classify_points_cpu(const std::vector<coord3d>& coords, Direction bfielddir) const;
  std::vector<Tropicity> classify_points(const std::vector<coord3d>& coords, Direction bfielddir) const {
#if HAVE_CUDA
#if 1
    cout << "manual interpolation on pitched_ptr of coord3d_d" << endl;
    return classify_points_cuda_v1(coords, bfielddir);
#else
    cout << "manual interpolation on 3*texture of int2/double" << endl;
    return classify_points_cuda_v2(coords, bfielddir);
#endif
#else
    return classify_points_cpu(coords, bfielddir);
#endif
  }

  Plane<Tropicity> gettropplaneZ(double zcoord, bool debug) const;
  Plane<Tropicity> gettropplane(Direction bfielddir, int fixeddir, double fixedcoord, bool debug) const;
  void splitgrid(std::string gridfile, std::string weightfile, Direction bfielddir) const;

  std::tuple<bool, coord3d> getvector(coord3d pos) const;
  std::tuple<bool, coord3d> getvector(double x, double y, double z) const {
    return getvector(coord3d(x, y, z));
  }

  coord3d getvector3(coord3d pos) const;

  void writetropplane(std::string filename, Plane<Tropicity> tropicities) const;
  void writecube(const std::string& filename) const;
};
