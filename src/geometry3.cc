
#include <cmath>
#include <iomanip>
#include <limits>

#include "geometry3.hh"


// calculation of the angle beta at b(0,0,0)
double coord3d::angle(const coord3d& a, const coord3d& c) {
  const double L2 = a.dot(a);
  const double R2 = c.dot(c);
  const double M2 = (c - a).dot(c - a);
  //if (abs(M2)<1.e-10) return 0;
  const double den = 2.0 * sqrt(L2 * R2);
  double arg = (L2 + R2 - M2) / den;
  if (arg > 1)
    return 0;
  if (arg < -1)
    return M_PI;
  return acos(arg);
}

double coord3d::dihedral(const coord3d& b, const coord3d& c, const coord3d& d) {
  const coord3d ab = b; // a=0
  const coord3d bc = c - b;
  const coord3d cd = d - c;

  const coord3d abc = ab.cross(bc);
  const coord3d bcd = bc.cross(cd);

  const coord3d bc1 = bc / bc.norm();
  const coord3d abc1 = abc / abc.norm();
  const coord3d bcd1 = bcd / bcd.norm();
  const coord3d aux = abc1.cross(bc1);

  const double x = abc1.dot(bcd1);
  const double y = aux.dot(bcd1);

  return atan2(y, x);
}

void coord3d::scale(const double f) {
  for (int i = 0; i < 3; i++) {
    x[i] *= f;
  }
}
