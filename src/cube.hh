#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <optional>
#include <vector>

#include "dir-enum.hh"
#include "geometry3.hh"
#include "plane.hh"
#include "trop-enum.hh"


class Cube {
  int n_x;
  int n_y;
  int n_z;
  std::vector<coord3d> field;
  coord3d origin = coord3d(0, 0, 0);
  coord3d spacing = coord3d(1, 1, 1);

public:
  Cube(int x, int y, int z): n_x(x), n_y(y), n_z(z), field(n_x * n_y * n_z){};
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

  constexpr coord3d get_origin() const { return origin; }
  constexpr coord3d get_spacing() const { return spacing; }

  // where position is in 'cube coordinates', ie, with origin 0, 0, 0, and unit-spacing
  // here we're excluding n_x-1 etc, even though we have values for those coordinates, because one cannot interpolate there, and making the destinction is expensive/not very useful
  constexpr bool outofbounds(const coord3d& position) const {
    if (position[0] >= n_x-1 || position[1] >= n_y-1 || position[2] >= n_z-1 || position[0] < 0 || position[1] < 0 || position[2] < 0) {
      return true;
    }
    return false;
  }

  std::vector<Tropicity> classify_points(const std::vector<coord3d>& coords, Direction bfielddir) const;

  Plane<Tropicity> gettropplaneZ(double zcoord, bool debug) const;
  Plane<Tropicity> gettropplane(Direction bfielddir, int fixeddir, double fixedcoord, bool debug) const;
  void splitgrid(std::string gridfile, std::string weightfile, Direction bfielddir) const;

  std::optional<coord3d> getvector(coord3d position) const;
  std::optional<coord3d> getvector(double x, double y, double z) const {
    return getvector(coord3d(x, y, z));
  }

  coord3d getvector3(coord3d position) const;

  void writetropplane(std::string filename, Plane<Tropicity> tropicities) const;
  void writecube(const std::string& filename) const;
};
