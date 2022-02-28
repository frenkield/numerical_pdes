cmake_minimum_required(VERSION 3.3)
project(tp2)

set(CMAKE_CXX_STANDARD 11)

file(GLOB src "src/*.hpp" "src/*.cpp")

add_executable(poisson src/apps/poisson.cpp ${src})

# cmake .
# cmake --build .
