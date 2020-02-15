
#include <cassert>
#include <iostream>
#include <vector>

#include "auxiliary.hh"
#include "geometry3.hh"
#include "trajectory.hh"

using namespace std;


bool trajectory::extend_euler(const Cube& cube) { //Euler
  const coord3d nextposition(positions.back() + directions.back().normalised() * step_length);
  auto optvect = cube.getvector(nextposition);
  if (!optvect)
    return false;
  append(nextposition, cube.getvector(nextposition).value());
  return true;
}


// the numbers in extend-rungekutta are not magic numbers. (see wikipedia article for runge-kutta method).
// any other numbers (like "10000" or "0.05" are probably magic numbers.
// beware
bool trajectory::extend_rungekutta(const Cube& cube) {
  const coord3d c0 = positions.back();
  const coord3d k0 = cube.getvector(c0).value().normalised() * step_length;
  const coord3d c1 = c0 + k0 * 0.5;
  auto v1 = cube.getvector(c1);
  if (!v1)
    return false;
  const coord3d k1 = v1.value().normalised() * step_length;
  const coord3d c2 = c0 + k1 * 0.5;
  auto v2 = cube.getvector(c2);
  if (!v2)
    return false;
  const coord3d k2 = v2.value().normalised() * step_length;
  const coord3d c3 = c0 + k2;
  auto v3 = cube.getvector(c3);
  if (!v3)
    return false;
  const coord3d k3 = v3.value().normalised() * step_length;
  const coord3d c4 = c0 + (k0 + k1 * 2.0 + k2 * 2.0 + k3) / 6.0;
  auto v4 = cube.getvector(c4);
  if (!v4)
    return false;
  append(c4, v4.value());
  return true;
}


void trajectory::complete(const Cube& cube) {
  //const double threshold = 1e-2;
  //if (directions[0].norm() < threshold) {out_of_bounds=true; return;} //if the intensity is vero low, don't bother completing. classify as "out_of_bounds"
  //the above is commented out because i didn't figure out what would be a good value for this threshold
  //if someone does, this would probably save some computational time
  const double step_length_ratio = 0.05;
  step_length = step_length_ratio * cube.get_spacing()[0];

  int step = 0;
  double dist2farthest = -1; //if this is set at 0 at declaration, the following while loop will never run

  const double return_ratio = 0.2;
  //if we get to a point that is less than SOME WELL-GUESSED FRACTION (1/5) of the longest distance in the trajectory
  while ((positions.back() - positions[0]).norm() > return_ratio * dist2farthest) {
    if (!extend_rungekutta(cube)) {
       out_of_bounds = true;
      return;
    }
    step++;

    if ((positions.back() - positions[0]).norm() > dist2farthest) {
      dist2farthest = (positions.back() - positions[0]).norm();
    }

    if (step > 10000) { //a single trajectory must not be more than this WELL-GUESSED number 10 000 of steps
      step = 0;
      step_length += 2;
      int size = positions.size();
      for (int a = 0; a < size - 1; a++) { //also this can't be very effective, there must be a way to wipe the positions and directions lists and create new blanks ones that is cheaper than this
        positions.pop_back();
        directions.pop_back();
      }
      dist2farthest = -1;
    }
  }
}


Tropicity trajectory::classify(const Cube& cube, int bfielddir) const {
  coord3d bfield;
  switch (bfielddir) {
    case 0: {
      bfield = coord3d(1, 0, 0);
      break;
    }
    case 1: {
      bfield = coord3d(-1, 0, 0);
      break;
    }
    case 2: {
      bfield = coord3d(0, 1, 0);
      break;
    }
    case 3: {
      bfield = coord3d(0, -1, 0);
      break;
    }
    case 4: {
      bfield = coord3d(0, 0, 1);
      break;
    }
    case 5: {
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
  const size_t polygon_size(directions.size());
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


void trajectory::write2mathematicalist(string filename) {
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

bool trajectory::to_mathematica(const trajectory& t, FILE* file) {
  ostringstream s;
  s << "traj = {{" << fixed << positions << "}};" << endl;
  s << "Graphics3D[{Arrow@#} & /@ data]" << endl;
  fputs(s.str().c_str(), file);

  return ferror(file) == 0;
}
