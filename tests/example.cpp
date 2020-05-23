#include "gtest/gtest.h"

class MinimalEx: public ::testing::Test {
public:
  const double epsilon = 1.e-6;
};


TEST_F(MinimalEx, ex1) {
  ASSERT_NEAR(6 - 3, 7 - 4, epsilon);
}
