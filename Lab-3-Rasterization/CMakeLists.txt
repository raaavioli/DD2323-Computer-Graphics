project(DgiLab3)

cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O2 -Wall")

add_executable(Rasterization
  rasterization.cpp
  SDL2Auxiliary.cpp
)

find_package (SDL2 REQUIRED)

message(STATUS "Lib: ${SDL2_LIBRARIES} , Include: ${SDL2_INCLUDE_DIRS}")

include_directories(
  ${CMAKE_SOURCE_DIR}/../glm
  ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(Rasterization
  ${SDL2_LIBRARIES}
)
