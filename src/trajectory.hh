#ifndef TRAJECTORY_HH
#define TRAJECTORY_HH

#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

#include "auxiliary.hh"
#include "cube.hh"
#include "geometry3.hh"
#include "trop-enum.hh"

using namespace std;


class trajectory {
  vector<coord3d> positions;
  vector<coord3d> directions;
  double step_length;
  bool out_of_bounds = false;

public:
  trajectory(const coord3d& pos, const coord3d& dir, const double eps = 0.01): step_length(eps) { positions.push_back(pos), directions.push_back(dir); }
  trajectory(const trajectory& traj): positions(traj.get_positions()), directions(traj.get_directions()), step_length(traj.get_step_length()) {}

  vector<coord3d> get_positions() const { return positions; }
  vector<coord3d> get_directions() const { return directions; }
  double get_step_length() const { return step_length; }
  double get_length() const { return positions.size(); }

  void append(const coord3d& pos, const coord3d& dir) {
    positions.push_back(pos);
    directions.push_back(dir);
  }
  void append(const pair<coord3d, coord3d>& p) {
    positions.push_back(p.first);
    directions.push_back(p.second);
  }
  pair<coord3d, coord3d> get(const unsigned int i) const { return make_pair(positions[i], directions[i]); }

  // extend trajectory by one element using Euler
  bool extend_euler(const Cube& cube);
  // extend trajectory by one element using Runge-Kutta
  bool extend_rungekutta(const Cube& cube);
  // extend trajectory until some criterion is met
  void complete(const Cube& cube);
  // return -1 or +1 for B dot (\sum r_i cross (p_i+1 - p_i)) less/greater zero
  Tropicity classify(const Cube& cube, int bfielddir) const;

  void write2mathematicalist(string filename);

  bool to_mathematica(const trajectory& t, FILE* file);


  friend ostream& operator<<(ostream& s, const trajectory& T) {
    s << fixed << T.positions << "," << T.directions << "," << T.step_length;
    return s;
  }
};

#endif
