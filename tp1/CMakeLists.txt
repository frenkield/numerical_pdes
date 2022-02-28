cmake_minimum_required(VERSION 3.3)
project(test1)

set(CMAKE_CXX_STANDARD 17)

file(GLOB src "src/*.hpp" "src/*.cpp")

add_executable(attractors test/attractors.cpp ${src})
add_executable(test_all test/test_all.cpp ${src})

# cmake .
# cmake --build .
