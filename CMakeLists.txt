cmake_minimum_required(VERSION 3.5.1)	 # Проверка версии CMake.

project(sibernetic)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")

include_directories(/usr/local/cuda/include/)
include_directories(inc)

set(TARGET sibernetic)

set(SOURCE_EXE src/sibernetic.cpp)

set(SOURCE_ARG_P src/utils/arg_parser.cpp)

set(SOURCE_HELPER src/utils/ocl_helper.cpp)

set(OPENCL_LIB OpenCL)

#add_library(foo SHARED ${SOURCE_LIB})	# Создание статической библиотеки с именем foo
add_executable(${TARGET} ${SOURCE_EXE} ${SOURCE_HELPER} ${SOURCE_ARG_P})

target_link_libraries(${TARGET} ${OPENCL_LIB})

