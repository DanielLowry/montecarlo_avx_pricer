# A very simple CMakeLists.txt to compile C++ analytics

cmake_minimum_required(VERSION 3.10)
project(montecarlo_avx LANGUAGES CXX)

option(ENABLE_LOGGING "Enable logging output" OFF)
option(ENABLE_TESTING "Enable test builds" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(mc_pricer STATIC curve.cpp simple_vol_surface.cpp hull_white_1f.cpp date.cpp)
target_include_directories(mc_pricer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(ENABLE_TESTING)
    enable_testing()
    add_executable(test_curve test_curve.cpp)
    target_link_libraries(test_curve PRIVATE mc_pricer)
    add_test(NAME test_curve COMMAND test_curve)
endif()

# set up relative path to python venv
set(PYTHON_VENV_DIR "${CMAKE_SOURCE_DIR}/../python/venv")

find_package(Python3 
    COMPONENTS Interpreter Development 
    REQUIRED
)

# If it fails to find pybind11 configs, you need to install them: sudo apt install python3-pybind11
find_package(pybind11
    CONFIG
    REQUIRED
)

pybind11_add_module(mc_pricer_py MODULE
    pybindings/pricer.cpp
)
set_target_properties(mc_pricer_py PROPERTIES 
SUFFIX ".so"
OUTPUT_NAME "mc_pricer_py"
)   

target_include_directories(mc_pricer_py PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mc_pricer_py PRIVATE mc_pricer)

if(ENABLE_LOGGING)
    target_compile_definitions(mc_pricer PRIVATE ENABLE_LOGGING)
    target_compile_definitions(mc_pricer_py PRIVATE ENABLE_LOGGING)
endif()

install(TARGETS mc_pricer_py
        LIBRARY
        DESTINATION mc_pricer_py)

install(
  FILES ${CMAKE_SOURCE_DIR}/../python/mpp_init/__init__.py
  DESTINATION mc_pricer_py
)