cmake_minimum_required(VERSION 3.20)
project(1)

set(CMAKE_CXX_STANDARD 20)

add_executable(1 plane_fitting.cpp)

find_package( CGAL REQUIRED )
target_link_libraries(1 CGAL::CGAL)