cmake_minimum_required (VERSION 2.6)
project (XtalComp)

add_library (XtalComp xtalcomp.cpp xctransform.cpp)

add_executable (test test.cpp)

target_link_libraries (test XtalComp)
