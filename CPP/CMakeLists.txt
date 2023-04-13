cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
project(AVD)
add_executable(AVD
      src/main.cpp
      src/NACA4.cpp
      src/Airfoil.cpp
      src/Atm.cpp
      src/ClassTesting.cpp
)
