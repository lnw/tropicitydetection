
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "auxiliary.hh"
#include "geometry3.hh"
#include "trajectory.hh"

using namespace std;


bool Trajectory::extend_euler(const Cube& cube) { // Euler
  const coord3d nextposition(positions.back() + directions.back().normalised() * step_length);
  auto optvect = cube.getvector(nextposition);
  if (!std::get<0>(optvect))
    return false;
  append(nextposition, std::get<1>(cube.getvector(nextposition)));
  return true;
}


// Runge-Kutta method, 4th order
bool Trajectory::extend_rungekutta(const Cube& cube) {
  const coord3d c0 = positions.back();
  const coord3d k0 = std::get<1>(cube.getvector(c0)).normalised() * step_length;
  const coord3d c1 = c0 + k0 * 0.5;
  auto v1 = cube.getvector(c1);
  if (!std::get<0>(v1))
    return false;
  const coord3d k1 = std::get<1>(v1).normalised() * step_length;
  const coord3d c2 = c0 + k1 * 0.5;
  auto v2 = cube.getvector(c2);
  if (!std::get<0>(v2))
    return false;
  const coord3d k2 = std::get<1>(v2).normalised() * step_length;
  const coord3d c3 = c0 + k2;
  auto v3 = cube.getvector(c3);
  if (!std::get<0>(v3))
    return false;
  const coord3d k3 = std::get<1>(v3).normalised() * step_length;
  const coord3d c4 = c0 + (k0 + k1 * 2.0 + k2 * 2.0 + k3) / 6.0;
  auto v4 = cube.getvector(c4);
  if (!std::get<0>(v4))
    return false;
  append(c4, std::get<1>(v4));
  return true;
}


void Trajectory::complete(const Cube& cube, double return_ratio) {
  // const double step_length_ratio = 0.05;
  // step_length = step_length_ratio * cube.get_spacing()[0];
  int max_steps = 10000;

  double dist2farthest = -1; // if this is set at 0 at declaration, the following while loop will never run
  if (positions.size() > 1) {
    for (auto pos: positions)
      dist2farthest = std::max(dist2farthest, (pos - positions[0]).norm());
  }

  int step = 0;
  // if we get to a point that is less return_ratio of the longest distance in the trajectory
  while ((positions.back() - positions[0]).norm() > return_ratio * dist2farthest) {
    if (!extend_rungekutta(cube)) {
      out_of_bounds = true;
      return;
    }
    step++;

    dist2farthest = std::max(dist2farthest, (positions.back() - positions[0]).norm());

    if (step > max_steps) {
      step = 0;
      step_length += 2;
      positions.clear();
      directions.clear();
      dist2farthest = -1;
    }
  }
}


Tropicity Trajectory::classify(Direction bfielddir) const {
  coord3d bfield;
  switch (bfielddir) {
    case Direction::pos_x: {
      bfield = coord3d(1, 0, 0);
      break;
    }
    case Direction::neg_x: {
      bfield = coord3d(-1, 0, 0);
      break;
    }
    case Direction::pos_y: {
      bfield = coord3d(0, 1, 0);
      break;
    }
    case Direction::neg_y: {
      bfield = coord3d(0, -1, 0);
      break;
    }
    case Direction::pos_z: {
      bfield = coord3d(0, 0, 1);
      break;
    }
    case Direction::neg_z: {
      bfield = coord3d(0, 0, -1);
      break;
    }
    default: {
      cerr << "bfielddir value wasn't 0-5.\n";
      return Tropicity::input_error;
    }
  }

  if (out_of_bounds)
    return Tropicity::outofbounds;

  coord3d crossum(0, 0, 0);
  const size_t polygon_size(positions.size());
  for (size_t i = 0; i < polygon_size; i++) {
    crossum += positions[(i - 1 + polygon_size) % polygon_size].cross(positions[i]);
  }
  // crossum += positions[positions.size()-1].cross(positions[0]);

  const double dot_product = bfield.dot(crossum);
  if (dot_product > 0)
    return Tropicity::paratropic;
  else if (dot_product < 0)
    return Tropicity::diatropic;
  else
    return Tropicity::unclassifyable;
}


void Trajectory::write2mathematicalist(string filename) {
  ofstream outputfile;
  outputfile.open(filename);
  outputfile << "traj = {{";
  for (size_t i = 0; i < positions.size(); i++) {
    outputfile << "{" << positions[i][0] << "," << positions[i][1] << "," << positions[i][2] << "}";
    if (i < positions.size() - 1) {
      outputfile << ",";
    }
  }
  outputfile << "}}";
}

bool Trajectory::to_mathematica(const Trajectory& t, FILE* file) {
  ostringstream s;
  s << "traj = {{" << fixed << positions << "}};" << endl;
  s << "Graphics3D[{Arrow@#} & /@ data]" << endl;
  fputs(s.str().c_str(), file);

  return ferror(file) == 0;
}
