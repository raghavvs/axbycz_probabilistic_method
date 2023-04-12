cmake_minimum_required(VERSION 3.12)
project(axbycz_probabilistic_method)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE "Debug")

add_compile_options(-g)

####----Uncomment when compiling on wall-e-----###
set(axbycz_probabilistic_method /home/raghav/cws/axbycz_probabilistic_method)

####----Uncomment when compiling on mac---######
#set(axbycz_probabilistic_method "/Users/raghavendrans/Documents/Raghav's /University of Illinois/AE 597 - Research/axbycz_probabilistic_method")

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

####----Uncomment when compiling on wall-e-----###
include_directories(/home/raghav/cws/axbycz_probabilistic_method/util)
include_directories(/home/raghav/cws/axbycz_probabilistic_method/solvers)

####----Uncomment when compiling on mac---######
#include_directories("/Users/raghavendrans/Documents/Raghav's /University of Illinois/AE 597 - Research/axbycz_probabilistic_method/util")
#include_directories("/Users/raghavendrans/Documents/Raghav's /University of Illinois/AE 597 - Research/axbycz_probabilistic_method/solvers")

include_directories(${axbycz_probabilistic_method}/util)

# Create executable
#add_executable(meanCov util/meanCov.cpp)
#add_executable(meanCov1 util/meanCov1.cpp)
#add_executable(meanCov2 util/meanCov2.cpp)
#add_executable(meanCov3 util/meanCov3.cpp)
#add_executable(meanCov_solvertest util/meanCov_solvertest.cpp)
#add_executable(meanCov_solvertest solvers/meanCov_solvertest.cpp)
#add_executable(meanCovTEST test/meanCovTEST.cpp)
#add_executable(mvg util/mvg.cpp)
#add_executable(mvgTEST test/mvgTEST.cpp)
#add_executable(getErrorAXBYCZ util/getErrorAXBYCZ.cpp)
#add_executable(getErrorAXBYCZTEST test/getErrorAXBYCZTEST.cpp)
#add_executable(initializeXYZTEST test/initializeXYZTEST.cpp)
#add_executable(initializeXYZTEST test/initializeXYZTEST.cpp util/initializeXYZ.cpp)
#add_executable(generateABC util/generateABC.cpp)
#add_executable(generateABC1 util/generateABC1.cpp)
#add_executable(generateSetsOfABC util/generateSetsOfABC.cpp)
#add_executable(generateSetsOfABC1 util/generateSetsOfABC1.cpp)
#add_executable(isEqualf util/isEqualf.cpp)
#add_executable(isEqualfTEST test/isEqualfTEST.cpp)
#add_executable(rotError util/rotError.cpp)
#add_executable(metric util/metric.cpp)
#add_executable(metricTEST test/metricTEST.cpp)
#add_executable(tranErrorTEST test/tranErrorTEST.cpp)
#add_executable(skewLogTEST test/skewLogTEST.cpp)
#add_executable(so3Vec util/so3Vec.cpp)
#add_executable(so3VecTEST test/so3VecTEST.cpp)
#add_executable(se3Vec util/se3Vec.cpp)
#add_executable(se3VecTEST test/se3VecTEST.cpp)
#add_executable(skewExpTEST test/skewExpTEST.cpp)
#add_executable(scrambleDataTEST test/scrambleDataTEST.cpp)
#add_executable(sensorNoiseTEST test/sensorNoiseTEST.cpp)
#add_executable(sensorNoise util/sensorNoise.cpp)
#add_executable(sensorNoise1 util/sensorNoise1.cpp)
#add_executable(sensorNoise2 util/sensorNoise2.cpp)
#add_executable(sensorNoise3 util/sensorNoise3.cpp)
#add_executable(sensorNoise4 util/sensorNoise4.cpp)
#add_executable(permFixABCTEST test/permFixABCTEST.cpp)
#add_executable(rotErrorTEST test/rotErrorTEST.cpp)
#add_executable(plotProbResults util/plotProbResults.cpp)
#add_executable(plotProbResults1 util/plotProbResults1.cpp)
#add_executable(fKine util/fKine.cpp)
#add_executable(fKineTEST test/fKineTEST.cpp)
#add_executable(batchSolveXY solvers/batchSolveXY.cpp)
#add_executable(batchSolveXY1 solvers/batchSolveXY1.cpp)
#add_executable(batchSolveXY2 solvers/batchSolveXY2.cpp)
#add_executable(batchSolveXY3 solvers/batchSolveXY3.cpp)
#add_executable(axbyczProb1 solvers/axbyczProb1.cpp)
#add_executable(axbyczProb2 solvers/axbyczProb2.cpp)
#add_executable(axbyczProb1_v2 solvers/axbyczProb1_v2.cpp)
#add_executable(axbyczProb1_v3 solvers/axbyczProb1_v3.cpp)
#add_executable(axbyczProb1_v4 solvers/axbyczProb1_v4.cpp)
#add_executable(axbyczProb1_v5 solvers/axbyczProb1_v5.cpp)
#add_executable(axbyczProb1_v6 solvers/axbyczProb1_v6.cpp)
#add_executable(axbyczProb3 solvers/axbyczProb3.cpp)
#add_executable(axbyczProb3_v3 solvers/axbyczProb3_v3.cpp)
#add_executable(temp solvers/temp.cpp)
#add_executable(mainSimulation scripts/mainSimulation.cpp)
#add_executable(mainSimulation2 scripts/mainSimulation2.cpp)
#add_executable(tf_echo_to_numbers scripts/tf_echo_to_numbers.cpp)
#add_executable(temp scripts/temp.cpp)
#add_executable(readMatrices util/readMatrices.cpp)
#add_executable(loadMatrices util/loadMatrices.cpp)
#add_executable(loadArraysToMatrices util/loadArraysToMatrices.cpp)
#add_executable(mainRealData scripts/mainRealData.cpp)
#add_executable(mainRealData1 scripts/mainRealData1.cpp)
add_executable(mainRealData3 scripts/mainRealData3.cpp)

# Create a list of libraries to link
set(LIBRARIES_TO_LINK
  Eigen3::Eigen
  ${GTEST_LIBRARIES}
  pthread
  ${PYTHON_LIBRARIES}
  Python3::NumPy
)

# Link all the targets with the library list
#target_link_libraries(meanCov ${LIBRARIES_TO_LINK})
#target_link_libraries(meanCov1 ${LIBRARIES_TO_LINK})
#target_link_libraries(meanCov2 ${LIBRARIES_TO_LINK})
#target_link_libraries(meanCov3 ${LIBRARIES_TO_LINK})
#target_link_libraries(meanCov_solvertest ${LIBRARIES_TO_LINK})
#target_link_libraries(meanCovTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(mvg ${LIBRARIES_TO_LINK})
#target_link_libraries(mvgTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(getErrorAXBYCZ ${LIBRARIES_TO_LINK})
#target_link_libraries(getErrorAXBYCZTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(initializeXYZTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(generateABC ${LIBRARIES_TO_LINK})
#target_link_libraries(generateABC1 ${LIBRARIES_TO_LINK})
#target_link_libraries(generateSetsOfABC ${LIBRARIES_TO_LINK})
#target_link_libraries(generateSetsOfABC1 ${LIBRARIES_TO_LINK})
#target_link_libraries(isEqualf ${LIBRARIES_TO_LINK})
#target_link_libraries(rotError ${LIBRARIES_TO_LINK})
#target_link_libraries(metric ${LIBRARIES_TO_LINK})
#target_link_libraries(metricTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(tranErrorTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(skewLogTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(so3Vec ${LIBRARIES_TO_LINK})
#target_link_libraries(so3VecTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(se3VecTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(se3Vec ${LIBRARIES_TO_LINK})
#target_link_libraries(skewExpTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(scrambleDataTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(sensorNoiseTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(sensorNoise ${LIBRARIES_TO_LINK})
#target_link_libraries(sensorNoise1 ${LIBRARIES_TO_LINK})
#target_link_libraries(sensorNoise2 ${LIBRARIES_TO_LINK})
#target_link_libraries(sensorNoise3 ${LIBRARIES_TO_LINK})
#target_link_libraries(sensorNoise4 ${LIBRARIES_TO_LINK})
#target_link_libraries(permFixABCTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(rotErrorTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(plotProbResults ${LIBRARIES_TO_LINK})
#target_link_libraries(plotProbResults1 ${LIBRARIES_TO_LINK})
#target_link_libraries(fKine ${LIBRARIES_TO_LINK})
#target_link_libraries(fKineTEST ${LIBRARIES_TO_LINK})
#target_link_libraries(batchSolveXY ${LIBRARIES_TO_LINK})
#target_link_libraries(batchSolveXY1 ${LIBRARIES_TO_LINK})
#target_link_libraries(batchSolveXY2 ${LIBRARIES_TO_LINK})
#target_link_libraries(batchSolveXY3 ${LIBRARIES_TO_LINK})
#target_link_libraries(axbyczProb1 ${LIBRARIES_TO_LINK})
#target_link_libraries(axbyczProb2 ${LIBRARIES_TO_LINK})
#target_link_libraries(axbyczProb1_v2 ${LIBRARIES_TO_LINK})
#target_link_libraries(axbyczProb1_v3 ${LIBRARIES_TO_LINK})
#target_link_libraries(axbyczProb1_v4 ${LIBRARIES_TO_LINK})
#target_link_libraries(axbyczProb1_v5 ${LIBRARIES_TO_LINK})
#target_link_libraries(axbyczProb1_v6 ${LIBRARIES_TO_LINK})
#target_link_libraries(axbyczProb3 ${LIBRARIES_TO_LINK})
#target_link_libraries(axbyczProb3_v3 ${LIBRARIES_TO_LINK})
#target_link_libraries(temp ${LIBRARIES_TO_LINK})
#target_link_libraries(mainSimulation ${LIBRARIES_TO_LINK})
#target_link_libraries(mainSimulation2 ${LIBRARIES_TO_LINK})
#target_link_libraries(tf_echo_to_numbers ${LIBRARIES_TO_LINK})
#target_link_libraries(temp ${LIBRARIES_TO_LINK})
#target_link_libraries(readMatrices ${LIBRARIES_TO_LINK})
#target_link_libraries(loadMatrices ${LIBRARIES_TO_LINK})
#target_link_libraries(loadArraysToMatrices ${LIBRARIES_TO_LINK})
#target_link_libraries(mainRealData ${LIBRARIES_TO_LINK})
#target_link_libraries(mainRealData1 ${LIBRARIES_TO_LINK})
target_link_libraries(mainRealData3 ${LIBRARIES_TO_LINK})


#option(TARGET_FILE "Specify the target file to compile" "axbyczProb1_v6")
# Create executable
#add_executable(${TARGET_FILE} ${TARGET_FILE}.cpp)

# Link the target with the library list
#target_link_libraries(${TARGET_FILE} ${LIBRARIES_TO_LINK})

#cmake -DTARGET_FILE=meanCov ..