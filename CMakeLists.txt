cmake_minimum_required(VERSION 2.8.3)
project(robustness)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -DFMT_HEADER_ONLY)
set(CMAKE_CXX_FLAGS "-O2")

find_package(Eigen3 REQUIRED)

include_directories(
  include
  include/basalt/thirdparty
)

add_executable(basalt_test src/basalt_test.cpp)
target_link_libraries (basalt_test Eigen3::Eigen armadillo fmt)

# TODO: compile and link eval_robustness.cpp