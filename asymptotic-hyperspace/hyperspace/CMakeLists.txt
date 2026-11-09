cmake_minimum_required(VERSION 3.15)
project(RiemannHyperspaceEngine CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    link_libraries(OpenMP::OpenMP_CXX)
endif()

add_executable(RiemannHyperspaceEngine 
    src/main.cpp 
    src/manifold.cpp 
    src/validation.cpp
)
target_include_directories(RiemannHyperspaceEngine PRIVATE src)
