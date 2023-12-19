project(rearrangement)

add_executable(rearrangement main.cpp)
target_compile_options(rearrangement PRIVATE --std=c++20 -march=native)
target_link_libraries(rearrangement PRIVATE pthread)
install(TARGETS rearrangement RUNTIME DESTINATION bin)
