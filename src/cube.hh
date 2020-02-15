#ifndef CUBE_HH
#define CUBE_HH

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
  vector<double> origin;
  vector<double> spacing;

public:
  Cube(string filename);

  const coord3d& operator()(const int z, const int y, const int x) const {
    return field[n_x * n_y * z + n_x * y + x];
  }
  coord3d& operator()(const int z, const int y, const int x) {
    return field[n_x * n_y * z + n_x * y + x];
  }

  vector<double> get_origin() const { return origin; }
  vector<double> get_spacing() const { return spacing; }
  bool outofbounds(coord3d position) const;

  vector<vector<Tropicity>> gettropplaneZ(double zcoord) const;
  vector<vector<Tropicity>> gettropplane(int bfielddir, int fixeddir, double fixedcoord) const;
  void splitgrid(string gridfile, string weightfile, int bfielddir) const;

  optional<coord3d> getvector(coord3d position) const;
  coord3d getvector3(coord3d position) const;

  void writetropplane(string filename, vector<vector<Tropicity>> tropicities) const;
  void writecube(const string& filename) const;
};

#endif
