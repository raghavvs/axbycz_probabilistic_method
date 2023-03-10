cmake_minimum_required(VERSION 3.12)
project(axbycz_probabilistic_method)

set(CMAKE_BUILD_TYPE "Debug")

add_compile_options(-g)

# Find Eigen3 library
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# find python libraries
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)
include_directories(${PYTHON3_INCLUDE_DIRS} ${NumPy_INCLUDE_DIRS})

# populate matplotlib repository
include(FetchContent)
FetchContent_Declare(
    matplotlib
    GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
    GIT_TAG        f23347fca25219d1c42cbb91608b5556814bf572
)
FetchContent_GetProperties(matplotlib)
if(NOT matplotlib_POPULATED)
    FetchContent_Populate(matplotlib)
endif()
include_directories(SYSTEM ${matplotlib_SOURCE_DIR})

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
#add_executable(skewLogTEST test/skewLogTEST.cpp)
#add_executable(se3Vec util/se3Vec.cpp)
#add_executable(so3VecTEST test/so3VecTEST.cpp)
#add_executable(se3VecTEST test/se3VecTEST.cpp)
#add_executable(skewExpTEST test/skewExpTEST.cpp)
#add_executable(scrambleDataTEST test/scrambleDataTEST.cpp)
#add_executable(sensorNoiseTEST test/sensorNoiseTEST.cpp)
add_executable(sensorNoise util/sensorNoise.cpp)
#add_executable(permFixABCTEST test/permFixABCTEST.cpp)
#add_executable(rotErrorTEST test/rotErrorTEST.cpp)
#add_executable(plotProbResults util/plotProbResults.cpp)
#add_executable(fKine util/fKine.cpp)
#add_executable(fKineTEST test/fKineTEST.cpp)

# Create a list of libraries to link
set(LIBRARIES_TO_LINK
  Eigen3::Eigen
  ${GTEST_LIBRARIES}
  pthread
  ${PYTHON_LIBRARIES}
  Python3::NumPy
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
#target_link_libraries(skewLogTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(so3VecTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(se3VecTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(se3Vec ${LIBRARIES_TO_LINK})
#target_link_libraries(skewExpTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(scrambleDataTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(sensorNoiseTEST ${LIBRARIES_TO_LINK})
target_link_libraries(sensorNoise ${LIBRARIES_TO_LINK})
#target_link_libraries(permFixABCTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(rotErrorTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(plotProbResults ${LIBRARIES_TO_LINK})
#target_link_libraries(fKine ${LIBRARIES_TO_LINK})
#target_link_libraries(fKineTEST ${LIBRARIES_TO_LINK})