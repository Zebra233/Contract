cmake_minimum_required(VERSION 3.19)
project(Contract)

set(CMAKE_CXX_STANDARD 14)

add_executable(Contract main.cpp include/user.h include/edge.h src/user.cpp src/edge.cpp include/contract.h include/setting.h src/contract.cpp)