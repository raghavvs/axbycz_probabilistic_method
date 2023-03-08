cmake_minimum_required(VERSION 3.12)
project(axbycz_probabilistic_method)

set(CMAKE_BUILD_TYPE "Debug")

add_compile_options(-g)

# Find Eigen3 library
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Include headers
include_directories(/home/raghav/cws/axbycz_probabilistic_method/util)

# Create executable
add_executable(meanCov util/meanCov.cpp)
# add_executable(meanCovTEST test/meanCovTEST.cpp)

# Link with Eigen3 Library
target_link_libraries(meanCov Eigen3::Eigen)
# target_link_libraries(meanCovTEST Eigen3::Eigen)
