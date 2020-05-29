#include "gtest/gtest.h"

#include "../src/cube.hh"
#include "../src/geometry3.hh"
// #include "../src/cuda_cube.hh"

using namespace std;

class cuda_interpolation_test: public ::testing::Test {
public:
  const double epsilon = 1.e-6;
};

TEST_F(cuda_interpolation_test, iscontinuous) {

  Cube cube(3, 3, 3);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-y, x, 0.5);
      }
    }
  }

  // ccube.upload();

  // cout << ccube.getvector(1.0 - epsilon, 0.5, 0.5).value() << endl;

  // ccube.download();

//   // continuous across centre of face x
//   ASSERT_NEAR(cube.getvector(1.0 - epsilon, 0.5, 0.5).value()[0],
//               cube.getvector(1.0 + epsilon, 0.5, 0.5).value()[0], 3 * epsilon);
//   ASSERT_NEAR(cube.getvector(1.0 - epsilon, 0.5, 0.5).value()[1],
//               cube.getvector(1.0 + epsilon, 0.5, 0.5).value()[1], 3 * epsilon);
//   ASSERT_NEAR(cube.getvector(1.0 - epsilon, 0.5, 0.5).value()[2],
//               cube.getvector(1.0 + epsilon, 0.5, 0.5).value()[2], 3 * epsilon);
// 
//   // continuous across centre of face y
//   ASSERT_NEAR(cube.getvector(0.5, 1.0 - epsilon, 0.5).value()[0],
//               cube.getvector(0.5, 1.0 + epsilon, 0.5).value()[0], 3 * epsilon);
//   ASSERT_NEAR(cube.getvector(0.5, 1.0 - epsilon, 0.5).value()[1],
//               cube.getvector(0.5, 1.0 + epsilon, 0.5).value()[1], 3 * epsilon);
//   ASSERT_NEAR(cube.getvector(0.5, 1.0 - epsilon, 0.5).value()[2],
//               cube.getvector(0.5, 1.0 + epsilon, 0.5).value()[2], 3 * epsilon);
// 
//   // continuous across centre of face z
//   ASSERT_NEAR(cube.getvector(0.5, 0.5, 1.0 - epsilon).value()[0],
//               cube.getvector(0.5, 0.5, 1.0 + epsilon).value()[0], 3 * epsilon);
//   ASSERT_NEAR(cube.getvector(0.5, 0.5, 1.0 - epsilon).value()[1],
//               cube.getvector(0.5, 0.5, 1.0 + epsilon).value()[1], 3 * epsilon);
//   ASSERT_NEAR(cube.getvector(0.5, 0.5, 1.0 - epsilon).value()[2],
//               cube.getvector(0.5, 0.5, 1.0 + epsilon).value()[2], 3 * epsilon);
// 
//   // continuous across edge xy
//   ASSERT_NEAR(cube.getvector(1.0 - epsilon, 1.0 - epsilon, 0.5).value()[0],
//               cube.getvector(1.0 + epsilon, 1.0 + epsilon, 0.5).value()[0],
//               4 * epsilon);
//   ASSERT_NEAR(cube.getvector(1.0 - epsilon, 1.0 - epsilon, 0.5).value()[1],
//               cube.getvector(1.0 + epsilon, 1.0 + epsilon, 0.5).value()[1],
//               4 * epsilon);
//   ASSERT_NEAR(cube.getvector(1.0 - epsilon, 1.0 - epsilon, 0.5).value()[2],
//               cube.getvector(1.0 + epsilon, 1.0 + epsilon, 0.5).value()[2],
//               4 * epsilon);
// 
//   // continuous across corner xyz
//   ASSERT_NEAR(
//       cube.getvector(1.0 - epsilon, 1.0 - epsilon, 1.0 - epsilon).value()[0],
//       cube.getvector(1.0 + epsilon, 1.0 + epsilon, 1.0 + epsilon).value()[0],
//       4 * epsilon);
//   ASSERT_NEAR(
//       cube.getvector(1.0 - epsilon, 1.0 - epsilon, 1.0 - epsilon).value()[1],
//       cube.getvector(1.0 + epsilon, 1.0 + epsilon, 1.0 + epsilon).value()[1],
//       4 * epsilon);
//   ASSERT_NEAR(
//       cube.getvector(1.0 - epsilon, 1.0 - epsilon, 1.0 - epsilon).value()[2],
//       cube.getvector(1.0 + epsilon, 1.0 + epsilon, 1.0 + epsilon).value()[2],
//       4 * epsilon);
}

// TEST_F(cuda_interpolation_test, nointerpolation) {
// 
//   Cube cube(3, 3, 3);
//   for (int x = 0; x < cube.nx(); x++) {
//     for (int y = 0; y < cube.ny(); y++) {
//       for (int z = 0; z < cube.nz(); z++) {
//         cube(x, y, z) = coord3d(-y, x, 0.5);
//       }
//     }
//   }
// 
//   ASSERT_NEAR(cube.getvector(1, 1, 1).value()[0], cube(1, 1, 1)[0], epsilon);
//   ASSERT_NEAR(cube.getvector(1, 1, 1).value()[1], cube(1, 1, 1)[1], epsilon);
//   ASSERT_NEAR(cube.getvector(1, 1, 1).value()[2], cube(1, 1, 1)[2], epsilon);
// 
//   ASSERT_NEAR(cube.getvector(1, 0, 1).value()[0], cube(1, 0, 1)[0], epsilon);
//   ASSERT_NEAR(cube.getvector(1, 0, 1).value()[1], cube(1, 0, 1)[1], epsilon);
//   ASSERT_NEAR(cube.getvector(1, 0, 1).value()[2], cube(1, 0, 1)[2], epsilon);
// }
