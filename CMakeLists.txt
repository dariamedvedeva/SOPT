cmake_minimum_required(VERSION 3.8)
project(secon_orger_diagram_sigma)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp Data.h Calculation.h InvertMatrix.h)
add_executable(secon_orger_diagram_sigma ${SOURCE_FILES})