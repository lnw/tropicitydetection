#include "gtest/gtest.h"

#include "../src/cube.hh"
#include "../src/geometry3.hh"
#include "../src/trajectory.hh"

using namespace std;

class CompleteTrajTest: public ::testing::Test {
public:
  const double epsilon = 1.e-6;
};

TEST_F(CompleteTrajTest, returnswithreasonablelength) {

  Cube cube(60, 60, 60);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-(y - 30), x - 30, 0.0);
      }
    }
  }

  coord3d startpos(50, 30, 30); // 20 from centre
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
  ASSERT_GT(length02, length04); // traj should be longer if we abort later

  traj.complete(cube, 0.1);
  // cout << traj.size() << endl;
  int length01 = traj.size();
  ASSERT_LT(traj.size(), 2 * M_PI * 20 / step_length);
  ASSERT_GT(length01, length02); // traj should be longer if we abort later

  traj.complete(cube, 0.05);
  // cout << traj.size() << endl;
  int length005 = traj.size();
  ASSERT_LT(traj.size(), 2 * M_PI * 20 / step_length);
  ASSERT_GT(length005, length01); // traj should be longer if we abort later

  traj.complete(cube, 0.02);
  // cout << traj.size() << endl;
  int length002 = traj.size();
  ASSERT_LT(traj.size(), 2 * M_PI * 20 / step_length);
  ASSERT_GT(length002, length005); // traj should be longer if we abort later
}


TEST_F(CompleteTrajTest, iscircular) {

  Cube cube(60, 60, 60);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-(y - 30), x - 30, 0.0);
      }
    }
  }

  coord3d startpos10(40, 30, 30); // 10 from centre
  coord3d startdir(0, 1, 0);
  double step_length = 0.02;
  Trajectory traj10(startpos10, startdir, step_length);

  coord3d startpos20(50, 30, 30); // 20 from centre
  Trajectory traj20(startpos20, startdir, step_length);

  coord3d startpos25(55, 30, 30); // 25 from centre
  Trajectory traj25(startpos25, startdir, step_length);

  traj10.complete(cube, 0.1);
  traj20.complete(cube, 0.1);
  traj25.complete(cube, 0.1);

  int length10 = traj10.size();
  int length20 = traj20.size();
  int length25 = traj25.size();

  ASSERT_LT(length10, 2 * M_PI * 10 / step_length);
  ASSERT_LT(length20, 2 * M_PI * 20 / step_length);
  ASSERT_LT(length25, 2 * M_PI * 25 / step_length);

  coord3d centre(30, 30, 30);

  for (auto& pos: traj10.get_positions()) {
    double rad = (pos - centre).norm();
    ASSERT_NEAR(rad, 10.0, 1.e-8);
  }

  for (auto& pos: traj20.get_positions()) {
    double rad = (pos - centre).norm();
    ASSERT_NEAR(rad, 20.0, 1.e-8);
  }

  for (auto& pos: traj25.get_positions()) {
    double rad = (pos - centre).norm();
    ASSERT_NEAR(rad, 25.0, 1.e-8);
  }
}
