cmake_minimum_required(VERSION 3.22)
project(Lab3_1_Lagrange_Newton)

set(CMAKE_CXX_STANDARD 14)

add_executable(Lab3_1_Lagrange_Newton main.cpp polynom.hpp newton_lagrange_worker.hpp)
