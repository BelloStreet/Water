cmake_minimum_required(VERSION 3.16)
project(MicroBenchmark)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})
find_package(Benchmark REQUIRED)

include_directories(${benchmark_INCLUDE_DIRS})
include_directories(headers)

file(GLOB_RECURSE SOURCE_FILES benchmark_Lobatto.cpp)

add_executable(MicroBenchmark ${SOURCE_FILES})

target_link_libraries(MicroBenchmark ${benchmark_LIBRARIES})
