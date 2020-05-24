#include "gtest/gtest.h"

#include "../src/cube.hh"
#include "../src/dir-enum.hh"
#include "../src/geometry3.hh"
#include "../src/trajectory.hh"
#include "../src/trop-enum.hh"

using namespace std;

class TropicityTest: public ::testing::Test {
public:
  const double epsilon = 1.e-6;
};


TEST_F(TropicityTest, classification_x) {
  Trajectory traj;
  coord3d dummy_direction(1, 0, 0);
  traj.append(coord3d(0, 0, 0), dummy_direction);
  traj.append(coord3d(0, 1, 0), dummy_direction);
  traj.append(coord3d(0, 1, 1), dummy_direction);
  traj.append(coord3d(0, 0, 1), dummy_direction);

  Tropicity trop = traj.classify(Direction::pos_x);
  ASSERT_EQ(trop, Tropicity::paratropic);

  trop = traj.classify(Direction::neg_x);
  ASSERT_EQ(trop, Tropicity::diatropic);
}


TEST_F(TropicityTest, classification_y) {
  Trajectory traj;
  coord3d dummy_direction(1, 0, 0);
  traj.append(coord3d(0, 0, 0), dummy_direction);
  traj.append(coord3d(1, 0, 0), dummy_direction);
  traj.append(coord3d(1, 0, 1), dummy_direction);
  traj.append(coord3d(0, 0, 1), dummy_direction);

  Tropicity trop = traj.classify(Direction::neg_y);
  ASSERT_EQ(trop, Tropicity::paratropic);

  trop = traj.classify(Direction::pos_y);
  ASSERT_EQ(trop, Tropicity::diatropic);
}


TEST_F(TropicityTest, classification_z) {
  Trajectory traj;
  coord3d dummy_direction(1, 0, 0);
  traj.append(coord3d(0, 0, 0), dummy_direction);
  traj.append(coord3d(1, 0, 0), dummy_direction);
  traj.append(coord3d(1, 1, 0), dummy_direction);
  traj.append(coord3d(0, 1, 0), dummy_direction);

  Tropicity trop = traj.classify(Direction::pos_z);
  ASSERT_EQ(trop, Tropicity::paratropic);

  trop = traj.classify(Direction::neg_z);
  ASSERT_EQ(trop, Tropicity::diatropic);
}
