#include "gtest/gtest.h"

#include "cube.hh"
#include "geometry3.hh"

using namespace std;

class InterpolationTest: public ::testing::Test {
public:
  const double epsilon = 1.e-6;
};

TEST_F(InterpolationTest, iscontinuous) {

  Cube cube(3, 3, 3);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-y, x, 0.5);
      }
    }
  }

  // continuous across centre of face x
  ASSERT_NEAR(std::get<1>(cube.getvector(1.0 - epsilon, 0.5, 0.5))[0],
              std::get<1>(cube.getvector(1.0 + epsilon, 0.5, 0.5))[0], 3 * epsilon);
  ASSERT_NEAR(std::get<1>(cube.getvector(1.0 - epsilon, 0.5, 0.5))[1],
              std::get<1>(cube.getvector(1.0 + epsilon, 0.5, 0.5))[1], 3 * epsilon);
  ASSERT_NEAR(std::get<1>(cube.getvector(1.0 - epsilon, 0.5, 0.5))[2],
              std::get<1>(cube.getvector(1.0 + epsilon, 0.5, 0.5))[2], 3 * epsilon);

  // continuous across centre of face y
  ASSERT_NEAR(std::get<1>(cube.getvector(0.5, 1.0 - epsilon, 0.5))[0],
              std::get<1>(cube.getvector(0.5, 1.0 + epsilon, 0.5))[0], 3 * epsilon);
  ASSERT_NEAR(std::get<1>(cube.getvector(0.5, 1.0 - epsilon, 0.5))[1],
              std::get<1>(cube.getvector(0.5, 1.0 + epsilon, 0.5))[1], 3 * epsilon);
  ASSERT_NEAR(std::get<1>(cube.getvector(0.5, 1.0 - epsilon, 0.5))[2],
              std::get<1>(cube.getvector(0.5, 1.0 + epsilon, 0.5))[2], 3 * epsilon);

  // continuous across centre of face z
  ASSERT_NEAR(std::get<1>(cube.getvector(0.5, 0.5, 1.0 - epsilon))[0],
              std::get<1>(cube.getvector(0.5, 0.5, 1.0 + epsilon))[0], 3 * epsilon);
  ASSERT_NEAR(std::get<1>(cube.getvector(0.5, 0.5, 1.0 - epsilon))[1],
              std::get<1>(cube.getvector(0.5, 0.5, 1.0 + epsilon))[1], 3 * epsilon);
  ASSERT_NEAR(std::get<1>(cube.getvector(0.5, 0.5, 1.0 - epsilon))[2],
              std::get<1>(cube.getvector(0.5, 0.5, 1.0 + epsilon))[2], 3 * epsilon);

  // continuous across edge xy
  ASSERT_NEAR(std::get<1>(cube.getvector(1.0 - epsilon, 1.0 - epsilon, 0.5))[0],
              std::get<1>(cube.getvector(1.0 + epsilon, 1.0 + epsilon, 0.5))[0],
              4 * epsilon);
  ASSERT_NEAR(std::get<1>(cube.getvector(1.0 - epsilon, 1.0 - epsilon, 0.5))[1],
              std::get<1>(cube.getvector(1.0 + epsilon, 1.0 + epsilon, 0.5))[1],
              4 * epsilon);
  ASSERT_NEAR(std::get<1>(cube.getvector(1.0 - epsilon, 1.0 - epsilon, 0.5))[2],
              std::get<1>(cube.getvector(1.0 + epsilon, 1.0 + epsilon, 0.5))[2],
              4 * epsilon);

  // continuous across corner xyz
  ASSERT_NEAR(
      std::get<1>(cube.getvector(1.0 - epsilon, 1.0 - epsilon, 1.0 - epsilon))[0],
      std::get<1>(cube.getvector(1.0 + epsilon, 1.0 + epsilon, 1.0 + epsilon))[0],
      4 * epsilon);
  ASSERT_NEAR(
      std::get<1>(cube.getvector(1.0 - epsilon, 1.0 - epsilon, 1.0 - epsilon))[1],
      std::get<1>(cube.getvector(1.0 + epsilon, 1.0 + epsilon, 1.0 + epsilon))[1],
      4 * epsilon);
  ASSERT_NEAR(
      std::get<1>(cube.getvector(1.0 - epsilon, 1.0 - epsilon, 1.0 - epsilon))[2],
      std::get<1>(cube.getvector(1.0 + epsilon, 1.0 + epsilon, 1.0 + epsilon))[2],
      4 * epsilon);
}

TEST_F(InterpolationTest, nointerpolation) {

  Cube cube(3, 3, 3);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-y, x, 0.5);
      }
    }
  }

  ASSERT_NEAR(std::get<1>(cube.getvector(1, 1, 1))[0], cube(1, 1, 1)[0], epsilon);
  ASSERT_NEAR(std::get<1>(cube.getvector(1, 1, 1))[1], cube(1, 1, 1)[1], epsilon);
  ASSERT_NEAR(std::get<1>(cube.getvector(1, 1, 1))[2], cube(1, 1, 1)[2], epsilon);

  ASSERT_NEAR(std::get<1>(cube.getvector(1, 0, 1))[0], cube(1, 0, 1)[0], epsilon);
  ASSERT_NEAR(std::get<1>(cube.getvector(1, 0, 1))[1], cube(1, 0, 1)[1], epsilon);
  ASSERT_NEAR(std::get<1>(cube.getvector(1, 0, 1))[2], cube(1, 0, 1)[2], epsilon);
}
