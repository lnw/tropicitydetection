#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <optional>
#include <vector>

#include "geometry3.hh"
#include "trop-enum.hh"

using namespace std;


class Cube {
  vector<coord3d> field;
  int n_x;
  int n_y;
  int n_z;
  coord3d origin;
  coord3d spacing;

public:
  Cube(string filename);

  const coord3d& operator()(const int z, const int y, const int x) const {
    return field[n_x * n_y * z + n_x * y + x];
  }
  coord3d& operator()(const int z, const int y, const int x) {
    return field[n_x * n_y * z + n_x * y + x];
  }

  constexpr coord3d get_origin() const { return origin; }
  constexpr coord3d get_spacing() const { return spacing; }
  constexpr bool outofbounds(const coord3d& position) const { // where position is in 'cube coordinates', ie, with origin 0,0,0, and unit-spacing
    if (position[0] > n_x || position[1] > n_y || position[2] > n_z || position[0] < 0 || position[1] < 0 || position[2] < 0) {
      return true;
    }
    return false;
  }

  vector<vector<Tropicity>> gettropplaneZ(double zcoord) const;
  vector<vector<Tropicity>> gettropplane(int bfielddir, int fixeddir, double fixedcoord) const;
  void splitgrid(string gridfile, string weightfile, int bfielddir) const;

  std::optional<coord3d> getvector(coord3d position) const;
  coord3d getvector3(coord3d position) const;

  void writetropplane(string filename, vector<vector<Tropicity>> tropicities) const;
  void writecube(const string& filename) const;
};

