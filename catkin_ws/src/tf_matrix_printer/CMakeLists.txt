cmake_minimum_required(VERSION 3.12)
project(tf_matrix_printer)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  tf
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS roscpp tf
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)

add_executable(${PROJECT_NAME}_node src/tf_matrix_printer.cpp)
target_link_libraries(${PROJECT_NAME}_node ${catkin_LIBRARIES})
