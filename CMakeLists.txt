cmake_minimum_required(VERSION 3.20)
project(1)
include_directories(/usr/local/include)

set(CMAKE_CXX_STANDARD 20)

add_executable(1 code_fixing.cpp)

find_package( CGAL REQUIRED )
target_link_libraries(1 CGAL::CGAL)