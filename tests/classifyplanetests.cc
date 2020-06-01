#include "gtest/gtest.h"

#include "../src/cube.hh"
#include "../src/geometry3.hh"
#include "../src/trajectory.hh"
#include "../src/trop-enum.hh"

using namespace std;

class ClassifyPlaneTest: public ::testing::Test {
public:
  const double epsilon = 1.e-6;
};


TEST_F(ClassifyPlaneTest, planeXmagPosX) {
  int dim = 40;
  Cube cube(dim, dim, dim);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-(y - dim / 2), x - dim / 2, 0.0);
      }
    }
  }

  Direction bfielddir = Direction::pos_x;
  bool debug = false;
  int planedir = 0;
  Plane<Tropicity> tropicities = cube.gettropplane(bfielddir, planedir, dim * 3.0 / 4.0, debug);
  // std::cout << vec_as_integer(tropicities.data()) << endl;

  int count_dia = 0,
      count_para = 0,
      count_uncl = 0;
  for (Tropicity tr: tropicities) {
    if (tr == Tropicity::diatropic)
      count_dia++;
    if (tr == Tropicity::paratropic)
      count_para++;
    if (tr == Tropicity::unclassifyable)
      count_uncl++;
  }
  // cout << count << endl;
  ASSERT_GT(count_dia, tropicities.size() * 0.15);   // subject to a lot of noise, technically all unclassifyable
  ASSERT_GT(count_para, tropicities.size() * 0.15);  // subject to a lot of noise, technically all unclassifyable
  ASSERT_GT(count_uncl, tropicities.size() * 0.005); // subject to a lot of noise, technically all unclassifyable
}


TEST_F(ClassifyPlaneTest, planeXmagPosZ) {
  int dim = 40;
  Cube cube(dim, dim, dim);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-(y - dim / 2), x - dim / 2, 0.0);
      }
    }
  }

  Direction bfielddir = Direction::pos_z;
  bool debug = false;
  int planedir = 0;
  Plane<Tropicity> tropicities = cube.gettropplane(bfielddir, planedir, dim * 3.0 / 4.0, debug);
  // std::cout << vec_as_integer(tropicities.data()) << endl;

  ASSERT_EQ(tropicities(10, 10), Tropicity::paratropic);
  ASSERT_EQ(tropicities(30, 10), Tropicity::paratropic);
  ASSERT_EQ(tropicities(10, 30), Tropicity::paratropic);
  ASSERT_EQ(tropicities(30, 30), Tropicity::paratropic);
}


TEST_F(ClassifyPlaneTest, planeXmagNegZ) {
  int dim = 40;
  Cube cube(dim, dim, dim);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-(y - dim / 2), x - dim / 2, 0.0);
      }
    }
  }

  Direction bfielddir = Direction::neg_z;
  bool debug = false;
  int planedir = 0;
  Plane<Tropicity> tropicities = cube.gettropplane(bfielddir, planedir, dim * 3.0 / 4.0, debug);
  // std::cout << vec_as_integer(tropicities.data()) << endl;

  ASSERT_EQ(tropicities(10, 10), Tropicity::diatropic);
  ASSERT_EQ(tropicities(30, 10), Tropicity::diatropic);
  ASSERT_EQ(tropicities(10, 30), Tropicity::diatropic);
  ASSERT_EQ(tropicities(30, 30), Tropicity::diatropic);
}


TEST_F(ClassifyPlaneTest, planeYmagPosX) {
  int dim = 40;
  Cube cube(dim, dim, dim);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-(y - dim / 2), x - dim / 2, 0.0);
      }
    }
  }

  Direction bfielddir = Direction::pos_x;
  bool debug = false;
  int planedir = 1;
  Plane<Tropicity> tropicities = cube.gettropplane(bfielddir, planedir, dim * 3.0 / 4.0, debug);
  // std::cout << vec_as_integer(tropicities.data()) << endl;

  int count_dia = 0,
      count_para = 0,
      count_uncl = 0;
  for (Tropicity tr: tropicities) {
    if (tr == Tropicity::diatropic)
      count_dia++;
    if (tr == Tropicity::paratropic)
      count_para++;
    if (tr == Tropicity::unclassifyable)
      count_uncl++;
  }
  // cout << count << endl;
  ASSERT_GT(count_dia, tropicities.size() * 0.15);   // subject to a lot of noise, technically all unclassifyable
  ASSERT_GT(count_para, tropicities.size() * 0.15);  // subject to a lot of noise, technically all unclassifyable
  ASSERT_GT(count_uncl, tropicities.size() * 0.005); // subject to a lot of noise, technically all unclassifyable
}


TEST_F(ClassifyPlaneTest, planeYmagPosZ) {
  int dim = 40;
  Cube cube(dim, dim, dim);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-(y - dim / 2), x - dim / 2, 0.0);
      }
    }
  }

  Direction bfielddir = Direction::pos_z;
  bool debug = false;
  int planedir = 1;
  Plane<Tropicity> tropicities = cube.gettropplane(bfielddir, planedir, dim * 3.0 / 4.0, debug);
  // std::cout << vec_as_integer(tropicities.data()) << endl;

  ASSERT_EQ(tropicities(10, 10), Tropicity::paratropic);
  ASSERT_EQ(tropicities(30, 10), Tropicity::paratropic);
  ASSERT_EQ(tropicities(10, 30), Tropicity::paratropic);
  ASSERT_EQ(tropicities(30, 30), Tropicity::paratropic);
}


TEST_F(ClassifyPlaneTest, planeYmagNegZ) {
  int dim = 40;
  Cube cube(dim, dim, dim);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-(y - dim / 2), x - dim / 2, 0.0);
      }
    }
  }

  Direction bfielddir = Direction::neg_z;
  bool debug = false;
  int planedir = 1;
  Plane<Tropicity> tropicities = cube.gettropplane(bfielddir, planedir, dim * 3.0 / 4.0, debug);
  // std::cout << vec_as_integer(tropicities.data()) << endl;

  ASSERT_EQ(tropicities(10, 10), Tropicity::diatropic);
  ASSERT_EQ(tropicities(30, 10), Tropicity::diatropic);
  ASSERT_EQ(tropicities(10, 30), Tropicity::diatropic);
  ASSERT_EQ(tropicities(30, 30), Tropicity::diatropic);
}


TEST_F(ClassifyPlaneTest, planeZmagPosX) {
  int dim = 40;
  Cube cube(dim, dim, dim);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-(y - dim / 2), x - dim / 2, 0.0);
      }
    }
  }

  Direction bfielddir = Direction::pos_x;
  bool debug = false;
  int planedir = 2;
  Plane<Tropicity> tropicities = cube.gettropplane(bfielddir, planedir, dim * 3.0 / 4.0, debug);
  // std::cout << vec_as_integer(tropicities.data()) << endl;

  int count_dia = 0,
      count_para = 0,
      count_uncl = 0;
  for (Tropicity tr: tropicities) {
    if (tr == Tropicity::diatropic)
      count_dia++;
    if (tr == Tropicity::paratropic)
      count_para++;
    if (tr == Tropicity::unclassifyable)
      count_uncl++;
  }
  // cout << count << endl;
  ASSERT_GT(count_dia, tropicities.size() * 0.15);  // subject to a lot of noise, technically all unclassifyable
  ASSERT_GT(count_para, tropicities.size() * 0.15); // subject to a lot of noise, technically all unclassifyable
  // ASSERT_GT(count_uncl, tropicities.size() * 0.005); // subject to a lot of noise, technically all unclassifyable
}


TEST_F(ClassifyPlaneTest, planeZmagPosZ) {
  int dim = 40;
  Cube cube(dim, dim, dim);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-(y - dim / 2), x - dim / 2, 0.0);
      }
    }
  }

  Direction bfielddir = Direction::pos_z;
  bool debug = false;
  int planedir = 2;
  Plane<Tropicity> tropicities = cube.gettropplane(bfielddir, planedir, dim * 3.0 / 4.0, debug);
  // std::cout << vec_as_integer(tropicities.data()) << endl;

  ASSERT_EQ(tropicities(5, 5), Tropicity::outofbounds);
  ASSERT_EQ(tropicities(5, 20), Tropicity::paratropic);
  ASSERT_EQ(tropicities(5, 35), Tropicity::outofbounds);
  ASSERT_EQ(tropicities(20, 5), Tropicity::paratropic);
  ASSERT_EQ(tropicities(20, 35), Tropicity::paratropic);
  ASSERT_EQ(tropicities(35, 5), Tropicity::outofbounds);
  ASSERT_EQ(tropicities(35, 20), Tropicity::paratropic);
  ASSERT_EQ(tropicities(35, 35), Tropicity::outofbounds);
}


TEST_F(ClassifyPlaneTest, planeZmagNegZ) {
  int dim = 40;
  Cube cube(dim, dim, dim);
  for (int x = 0; x < cube.nx(); x++) {
    for (int y = 0; y < cube.ny(); y++) {
      for (int z = 0; z < cube.nz(); z++) {
        cube(x, y, z) = coord3d(-(y - dim / 2), x - dim / 2, 0.0);
      }
    }
  }

  Direction bfielddir = Direction::neg_z;
  bool debug = false;
  int planedir = 2;
  Plane<Tropicity> tropicities = cube.gettropplane(bfielddir, planedir, dim * 3.0 / 4.0, debug);
  // std::cout << vec_as_integer(tropicities.data()) << endl;

  // in the corners
  ASSERT_EQ(tropicities(5, 5), Tropicity::outofbounds);
  ASSERT_EQ(tropicities(5, 35), Tropicity::outofbounds);
  ASSERT_EQ(tropicities(35, 5), Tropicity::outofbounds);
  ASSERT_EQ(tropicities(35, 35), Tropicity::outofbounds);

  // not in the corners
  ASSERT_EQ(tropicities(5, 20), Tropicity::diatropic);
  ASSERT_EQ(tropicities(20, 5), Tropicity::diatropic);
  ASSERT_EQ(tropicities(20, 35), Tropicity::diatropic);
  ASSERT_EQ(tropicities(35, 20), Tropicity::diatropic);
}
