cmake_minimum_required(VERSION 3.26)

project(rearrangement)

add_executable(rearrangement main.cpp)
target_compile_options(rearrangement PRIVATE --std=c++17 -march=native)
target_link_libraries(rearrangement PRIVATE pthread)
install(TARGETS rearrangement RUNTIME DESTINATION bin)
