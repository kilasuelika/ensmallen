cmake_minimum_required(VERSION 3.10)
project(main)
set(CMAKE_C_COMPILER /usr/bin/gcc CACHE PATH "" FORCE)
set(CMAKE_CXX_COMPILER /usr/bin/g++ CACHE PATH "" FORCE)

#Eigen
include_directories("/usr/local/include/eigen3")
add_executable(mp main.cpp)

#EigenRand
include_directories("/usr/local/include/")

#Ensmallen
include_directories("/mnt/d/Source/ensmallen-eigen/include")

target_compile_features(mp PUBLIC cxx_std_20)
add_custom_command(TARGET mp POST_BUILD
	COMMAND echo "----------------------------------------------"
	COMMAND ./mp
	COMMAND echo "----------------------------------------------")
