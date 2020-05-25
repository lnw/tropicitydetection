#include "gtest/gtest.h"

#include "../src/cube.hh"
#include "../src/geometry3.hh"
#include "../src/trajectory.hh"

using namespace std;

class CompleteTrajTest: public ::testing::Test {
public:
  const double epsilon = 1.e-6;
};

TEST_F(CompleteTrajTest, returns) {

  Cube cube(60, 60, 60);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-(y - 30), x - 30, 0.0);
      }
    }
  }

  coord3d startpos(50, 30, 30);  // 20 from centre
  coord3d startdir(0, 1, 0);
  double step_length = 0.02;

  // maximum expected circumference is 2*pi*20 = 126
  Trajectory traj(startpos, startdir, step_length);


  traj.complete(cube, 0.4);
  // cout << traj.size() << endl;
  int length04 = traj.size();
  ASSERT_LT(traj.size(), 2 * M_PI * 20 / step_length);

  traj.complete(cube, 0.2);
  // cout << traj.size() << endl;
  int length02 = traj.size();
  ASSERT_LT(traj.size(), 2 * M_PI * 20 / step_length);
  ASSERT_GT(length02, length04);  // traj should be longer if we abort later

  traj.complete(cube, 0.1);
  // cout << traj.size() << endl;
  int length01 = traj.size();
  ASSERT_LT(traj.size(), 2 * M_PI * 20 / step_length);
  ASSERT_GT(length01, length02);  // traj should be longer if we abort later

  traj.complete(cube, 0.05);
  // cout << traj.size() << endl;
  int length005 = traj.size();
  ASSERT_LT(traj.size(), 2 * M_PI * 20 / step_length);
  ASSERT_GT(length005, length01);  // traj should be longer if we abort later

  traj.complete(cube, 0.02);
  // cout << traj.size() << endl;
  int length002 = traj.size();
  ASSERT_LT(traj.size(), 2 * M_PI * 20 / step_length);
  ASSERT_GT(length002, length005);  // traj should be longer if we abort later
}
