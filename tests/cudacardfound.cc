#include "gtest/gtest.h"

#include "../src/cube.hh"
#include "../src/geometry3.hh"
#include "../src/trajectory.hh"
#include "../src/trop-enum.hh"

#include "../src/gpu-info.hh"


using namespace std;

class CudaCardFoundTest: public ::testing::Test {
public:
  const double epsilon = 1.e-6;
};


TEST_F(CudaCardFoundTest, anydevice) {

  int devCount = number_cuda_devices();
  ASSERT_GE(devCount, 1);
}


TEST_F(CudaCardFoundTest, memminimum) {

  int devCount = number_cuda_devices_minimum_mem_mb(2000); // at least one card with at least 2GB mem
  ASSERT_GE(devCount, 1);
}


TEST_F(CudaCardFoundTest, ccminimum) {

  int devCount = number_cuda_devices_minimum_cc(30); // at least one card with at least CUDA 3.0
  ASSERT_GE(devCount, 1);
}
