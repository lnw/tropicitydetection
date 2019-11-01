#ifndef CUBE_HH
#define CUBE_HH

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "geometry3.hh"
#include "trop-enum.hh"

using namespace std;


class Cube {
  vector<coord3d> field;
  int xrange;
  int yrange;
  int zrange;
  vector<double> origin;
  vector<double> spacing;

public:
  Cube(string filename);

  vector<double> get_origin() const { return origin; }
  vector<double> get_spacing() const { return spacing; }
  bool outofbounds(coord3d position) const;

  vector<vector<Tropicity>> gettropplaneZ(double zcoord) const;
  vector<vector<Tropicity>> gettropplane(int bfielddir, int fixeddir, double fixedcoord) const;
  void splitgrid(string gridfile, string weightfile, int bfielddir) const;

  coord3d getvector(coord3d position) const;
  coord3d getvector3(coord3d position) const;

  void writetropplane(string filename, vector<vector<Tropicity>> tropicities) const;
  void writecube(const string& filename) const;
};

#endif
