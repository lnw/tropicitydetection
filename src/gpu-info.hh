#pragma once

#include <cstdio>
#include <cstdlib>
#include <cuda.h>
#include <cuda_runtime_api.h>


void print_device_props_short();
void print_device_props_complete();
int number_cuda_devices();
int number_cuda_devices_minimum_mem_mb(size_t min_mem);
int number_cuda_devices_minimum_cc(int min_cc);
