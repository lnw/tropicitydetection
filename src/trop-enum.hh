#pragma once

enum class Tropicity : short { diatropic = 0,
                               paratropic,
                               outofbounds,
                               unclassifyable,
                               input_error };

// print one Tropicity as an integer
template <typename Enum>
auto as_underlying(Enum const value) -> typename std::underlying_type<Enum>::type {
  return static_cast<typename std::underlying_type<Enum>::type>(value);
}

// print a vector of Tropicity's as integers
template <typename Enum>
auto vec_as_underlying(const vector<Enum>& values) -> vector<typename std::underlying_type<Enum>::type> {
  vector<typename std::underlying_type<Enum>::type> out;
  for (auto value: values)
    out.push_back(static_cast<typename std::underlying_type<Enum>::type>(value));
  return out;
}
