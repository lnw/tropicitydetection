#pragma once

#include <cassert>
#include <vector>


template <typename T>
class Plane {
  int nx, ny;
  std::vector<T> data_;

public:
  Plane(int x, int y): nx(x), ny(y), data_(nx * ny) {}
  Plane(int x, int y, std::vector<T> d): nx(x), ny(y) {
    assert(nx * ny == d.size());
    data_ = d;
  }

  constexpr T operator[](int64_t n) const { return data_[n]; }
  constexpr T& operator()(int64_t n) { return data_[n]; }
  constexpr T operator()(int x, int y) const { return data_[y * nx + x]; }
  constexpr T& operator()(int x, int y) { return data_[y * nx + x]; }
  constexpr size_t size() const { return data_.size(); }

  constexpr const std::vector<T>& data() const { return data_; }
};
