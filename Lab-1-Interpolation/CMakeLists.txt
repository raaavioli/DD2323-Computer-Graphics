project(DgiLab1)

cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O2 -Wall")

IF(APPLE)
  SET(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architectures for Mac OS X" FORCE)
ENDIF(APPLE)

add_executable(Colors
  colors.cpp
  SDL2Auxiliary.cpp
)

add_executable(StarField
  starfield.cpp
  SDL2Auxiliary.cpp
)

find_package (SDL2 REQUIRED)

message(STATUS "Lib: ${SDL2_LIBRARIES} , Include: ${SDL2_INCLUDE_DIRS}")

include_directories(
  ${CMAKE_SOURCE_DIR}/../glm
  ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(StarField
  ${SDL2_LIBRARIES}
)

target_link_libraries(Colors
  ${SDL2_LIBRARIES}
)
