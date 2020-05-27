#pragma once

#include <cassert>
#include <vector>


template <typename T>
class Plane {
  int nx_, ny_;
  std::vector<T> data_;

public:
  Plane(int x, int y): nx_(x), ny_(y), data_(nx_ * ny_) {}
  Plane(int x, int y, std::vector<T> d): nx_(x), ny_(y) {
    assert(static_cast<size_t>(nx_ * ny_) == d.size());
    data_ = d;
  }

  constexpr T operator[](int64_t n) const { return data_[n]; }
  constexpr T& operator()(int64_t n) { return data_[n]; }
  constexpr T operator()(int x, int y) const { return data_[y * nx_ + x]; }
  constexpr T& operator()(int x, int y) { return data_[y * nx_ + x]; }

  constexpr typename std::vector<T>::iterator begin() { return data_.begin(); }
  constexpr typename std::vector<T>::const_iterator begin() const { return data_.begin(); }
  constexpr typename std::vector<T>::iterator end() { return data_.end(); }
  constexpr typename std::vector<T>::const_iterator end() const { return data_.end(); }

  constexpr size_t size() const { return data_.size(); }
  constexpr int64_t signedSize() const { return data_.size(); }

  constexpr const std::vector<T>& data() const { return data_; }

  constexpr int nx() const { return nx_; }
  constexpr int ny() const { return ny_; }
};
