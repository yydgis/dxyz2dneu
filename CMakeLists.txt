cmake_minimum_required(VERSION 3.16)
project(dxyz2dneu)

set(CMAKE_CXX_STANDARD 17)

add_executable(dxyz2dneu
dxyz2dneu.cpp
    )

target_link_libraries(dxyz2dneu)