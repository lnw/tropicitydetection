#pragma once

#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

#include "auxiliary.hh"
#include "cube.hh"
#include "dir-enum.hh"
#include "geometry3.hh"
#include "trop-enum.hh"

using namespace std;


class Trajectory {
  vector<coord3d> positions;
  vector<coord3d> directions;
  double step_length;
  bool out_of_bounds = false;

public:
  Trajectory() = default;
  Trajectory(const coord3d& pos, const coord3d& dir, const double step = 0.01): positions(1, pos), directions(1, dir), step_length(step) {}
  Trajectory(const Trajectory& traj): positions(traj.get_positions()), directions(traj.get_directions()), step_length(traj.get_step_length()) {}

  const vector<coord3d>& get_positions() const { return positions; }
  const vector<coord3d>& get_directions() const { return directions; }
  double get_step_length() const { return step_length; }
  void set_step_length(double step) { step_length = step; }
  size_t size() const { return positions.size(); }

  void append(const coord3d& pos, const coord3d& dir) {
    positions.push_back(pos);
    directions.push_back(dir);
  }

  // extend trajectory by one element using Euler
  bool extend_euler(const Cube& cube);
  // extend trajectory by one element using Runge-Kutta
  bool extend_rungekutta(const Cube& cube);
  // extend trajectory until some criterion is met
  void complete(const Cube& cube, double return_ratio = 0.2);
  // return -1 or +1 for B dot (\sum r_i cross (p_i+1 - p_i)) less/greater zero
  Tropicity classify(Direction bfielddir) const;

  void write2mathematicalist(string filename);
  bool to_mathematica(const Trajectory& t, FILE* file);

  friend ostream& operator<<(ostream& s, const Trajectory& T) {
    s << fixed << T.positions << "," << T.directions << "," << T.step_length;
    return s;
  }
};

