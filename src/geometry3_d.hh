#pragma once

#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "auxiliary.hh"

using namespace std;

struct coord3d;

struct coord3d_d {
  double x[3] = {0, 0, 0};

  __host__ __device__ constexpr coord3d_d(const double y[3]) {
    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
  }
  __host__ __device__ constexpr coord3d_d(const double x_, const double y, const double z) {
    x[0] = x_;
    x[1] = y;
    x[2] = z;
  }
  __host__ __device__ coord3d_d() {}
  __host__ __device__ inline coord3d_d(const double theta, const double phi) {
    x[0] = sin(theta) * cos(phi);
    x[1] = sin(theta) * sin(phi);
    x[2] = cos(theta);
  } // and r=1

  __host__ __device__ constexpr coord3d_d operator/(const double d) const { return coord3d_d(*this) /= d; }
  __host__ __device__ constexpr coord3d_d& operator/=(const double d) {
    x[0] /= d;
    x[1] /= d;
    x[2] /= d;
    return *this;
  }
  __device__ constexpr coord3d_d operator*(const double d) const { return coord3d_d(*this) *= d; }
  friend constexpr coord3d_d operator*(const double d, const coord3d_d& c3d) { return c3d * d; }
  friend constexpr coord3d_d operator*(const int i, const coord3d_d& c3d) { return c3d * i; }
  __device__ constexpr coord3d_d& operator*=(const double d) {
    x[0] *= d;
    x[1] *= d;
    x[2] *= d;
    return *this;
  }
  __device__ constexpr coord3d_d operator+(const coord3d_d& y) const { return coord3d_d(*this) += y; }
  __device__ constexpr coord3d_d& operator+=(const coord3d_d& y) {
    x[0] += y[0];
    x[1] += y[1];
    x[2] += y[2];
    return *this;
  }
  __device__ constexpr coord3d_d operator-(const coord3d_d& y) const { return coord3d_d(*this) -= y; }
  __device__ constexpr coord3d_d& operator-=(const coord3d_d& y) {
    x[0] -= y[0];
    x[1] -= y[1];
    x[2] -= y[2];
    return *this;
  }
  __device__ constexpr coord3d_d operator-() const {
    coord3d_d y(-x[0], -x[1], -x[2]);
    return y;
  }

  __device__ constexpr coord3d_d cross(const coord3d_d& y) const {
    return coord3d_d(x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]);
  }
  __device__ constexpr double dot(const coord3d_d& y) const { return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]; }
  __device__ inline double norm() const { return sqrt(dot(*this)); }

  __device__ inline coord3d_d normalised() const {
    const double normi = norm();
    if (normi == 0) {
      return coord3d_d(0, 0, 0);
    }
    else {
      return *this / normi;
    }
  }

  __device__ constexpr double& operator[](const unsigned int i) { return x[i]; }
  __device__ constexpr double operator[](const unsigned int i) const { return x[i]; }


  __device__ static double dist(const coord3d_d& x, const coord3d_d& y) { return (x - y).norm(); }
//   // d^2/(dx_i dx_j) ||x|| = -x_i x_j/||x||^3 + [i==j]/||x||
//   // calculation of the angle beta at b(0,0,0)
//   static double angle(const coord3d_d& a, const coord3d_d& c);
//   // calculation of the dihedral angle theta at a(0,0,0), b, c and d,  the result is an angle between -\pi and +\pi (in radians)
//   static double dihedral(const coord3d_d& b, const coord3d_d& c, const coord3d_d& d);


//   friend vector<coord3d_d>& operator-=(vector<coord3d_d>& xs, const coord3d_d& y) {
//     for (size_t i = 0; i < xs.size(); i++)
//       xs[i] -= y;
//     return xs;
//   }
// 
//   friend vector<coord3d_d>& operator*=(vector<coord3d_d>& xs, const double& y) {
//     for (size_t i = 0; i < xs.size(); i++)
//       xs[i] *= y;
//     return xs;
//   }

};
