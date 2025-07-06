# A very simple CMakeLists.txt to compile C++ analytics

cmake_minimum_required(VERSION 3.10)
project(montecarlo_avx LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# we are only going to have one l
add_library(curve SHARED curve.cpp simple_vol_surface.cpp hull_white_1f.cpp)
target_include_directories(curve PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# set up relative path to python venv
set(PYTHON_VENV_DIR "${CMAKE_SOURCE_DIR}/../python/venv")

# The venv doesn't have the cmake headers, so we use the system python
find_package(Python3 
    COMPONENTS Interpreter Development 
    REQUIRED
)

# Find pybind11 in the venv, not the system python
find_package(pybind11
    CONFIG
    REQUIRED
    HINTS ${PYTHON_VENV_DIR}

)
