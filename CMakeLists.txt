cmake_minimum_required(VERSION 3.22)

set(CMAKE_BUILD_TYPE Debug)

project(runge)



set(CMAKE_CXX_STANDARD 14)

# Add the include directory to the include paths
include_directories(include)

add_executable(runge main.cpp)

# Testing section
enable_testing()

add_executable(runge_test1 test/test_1.cpp)
add_test(NAME runge_test1 COMMAND runge_test1)

# add_executable(runge_test2 test/test_2.cpp)
# add_test(NAME runge_test2 COMMAND runge_test2)

add_executable(runge_test_extra test/additional_test.cpp)
add_test(NAME runge_test_extra COMMAND runge_test_extra)