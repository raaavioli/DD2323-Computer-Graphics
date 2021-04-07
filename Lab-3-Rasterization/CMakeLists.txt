project(DgiLab3)

cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O2 -Wall")

IF(APPLE)
  SET(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architectures for Mac OS X" FORCE)
ENDIF(APPLE)

find_package (SDL2 REQUIRED)

add_executable(Rasterization
  SDL2Auxiliary.cpp
  rasterization.cpp
)

message(STATUS "Lib: ${SDL2_LIBRARIES} , Include: ${SDL2_INCLUDE_DIRS}")

include_directories(
  ${CMAKE_SOURCE_DIR}/../glm
  ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(Rasterization
 ${SDL2_LIBRARIES}
)
