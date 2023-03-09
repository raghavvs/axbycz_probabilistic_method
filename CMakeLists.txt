cmake_minimum_required(VERSION 3.12)
project(axbycz_probabilistic_method)

set(CMAKE_BUILD_TYPE "Debug")

add_compile_options(-g)

# Find Eigen3 library
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Include headers
include_directories(/home/raghav/cws/axbycz_probabilistic_method/util)

# Create executable
add_executable(meanCov1 util/meanCov1.cpp)
#add_executable(meanCovTEST test/meanCovTEST.cpp)


# Link with Eigen3 Library
target_link_libraries(meanCov1 Eigen3::Eigen)
#target_link_libraries(meanCovTEST Eigen3::Eigen)
#target_link_libraries(meanCovTEST ${GTEST_LIBRARIES} pthread)
