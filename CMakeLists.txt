cmake_minimum_required(VERSION 3.28)
project(RungeKutta_method)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -ffast-math")

add_executable(RungeKutta_method main.cpp)
