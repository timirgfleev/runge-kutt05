cmake_minimum_required(VERSION 3.22)
project(runge)

set(CMAKE_CXX_STANDARD 14)

# Add the include directory to the include paths
include_directories(include)

add_executable(runge main.cpp)