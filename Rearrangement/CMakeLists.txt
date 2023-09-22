cmake_minimum_required(VERSION 3.16.3)
project(rearrangement)

add_executable(rearrangement main.cpp)
target_compile_options(rearrangement PRIVATE)

install(TARGETS rearrangement RUNTIME DESTINATION bin)
