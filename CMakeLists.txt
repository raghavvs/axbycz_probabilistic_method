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
#add_executable(meanCov1 util/meanCov1.cpp)
#add_executable(meanCovTEST test/meanCovTEST.cpp)
#add_executable(mvg util/mvg.cpp)
#add_executable(mvgTEST test/mvgTEST.cpp)
#add_executable(getErrorAXBYCZTEST test/getErrorAXBYCZTEST.cpp)
#add_executable(initializeXYZTEST test/initializeXYZTEST.cpp)
#add_executable(initializeXYZTEST test/initializeXYZTEST.cpp util/initializeXYZ.cpp)
#add_executable(generateABC util/generateABC.cpp)
#add_executable(isEqualf util/isEqualf.cpp)
#add_executable(isEqualfTEST test/isEqualfTEST.cpp)
#add_executable(metric util/metric.cpp)
#add_executable(metricTEST test/metricTEST.cpp)
#add_executable(tranErrorTEST test/tranErrorTEST.cpp)
add_executable(skewLogTEST test/skewLogTEST.cpp)

# Create a list of libraries to link
set(LIBRARIES_TO_LINK
  Eigen3::Eigen
  ${GTEST_LIBRARIES}
  pthread
)

# Link all the targets with the library list
#target_link_libraries(meanCov1 ${LIBRARIES_TO_LINK})
#target_link_libraries(meanCovTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(mvg ${LIBRARIES_TO_LINK})
#target_link_libraries(mvgTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(getErrorAXBYCZTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(initializeXYZTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(generateABC ${LIBRARIES_TO_LINK})
#target_link_libraries(isEqualf ${LIBRARIES_TO_LINK})
#target_link_libraries(metric ${LIBRARIES_TO_LINK})
#target_link_libraries(metricTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(tranErrorTEST ${LIBRARIES_TO_LINK})
target_link_libraries(skewLogTEST ${LIBRARIES_TO_LINK})