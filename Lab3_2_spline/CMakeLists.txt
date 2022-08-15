cmake_minimum_required(VERSION 3.22)
project(Lab3_2_spline)

set(CMAKE_CXX_STANDARD 14)

add_executable(Lab3_2_spline main.cpp polynom.hpp spline.hpp TridiagonalWorker.hpp)
