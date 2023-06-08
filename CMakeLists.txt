cmake_minimum_required(VERSION 3.19)
project(Contract)

set(CMAKE_CXX_STANDARD 14)

# 添加CPLEX的头文件目录
include_directories(/Applications/CPLEX_Studio_Community2211/cplex/include /Applications/CPLEX_Studio_Community2211/concert/include)

# 添加CPLEX的库文件目录
link_directories(/Applications/CPLEX_Studio_Community2211/cplex/lib/arm64_osx/static_pic /Applications/CPLEX_Studio_Community2211/concert/lib/arm64_osx/static_pic)

# 添加要编译的源代码文件
add_executable(Contract main.cpp include/user.h include/edge.h src/user.cpp src/edge.cpp include/contract.h include/setting.h src/contract.cpp)

# 链接CPLEX的库文件
target_link_libraries(Contract ilocplex cplex concert)# cplexdistm cplexdistmst)