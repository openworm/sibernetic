cmake_minimum_required(VERSION 3.5.1)	 # Проверка версии CMake.

project(sibernetic)

set(CMAKE_CXX_STANDARD 14)


include_directories(/usr/local/cuda/include/)
#include_directories(/opt/intel/opencl/SDK/include/)

include_directories(inc)

set(CMAKE_BUILD_TYPE Debug)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -Wall -Wextra -O0 -g")
endif()

if(CMAKE_BUILD_TYPE STREQUAL Release)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -Wall -Wextra -O3")
endif()

set(TARGET sibernetic)
set(SOURCE_EXE src/sibernetic.cpp)
set(SOURCE_ARG_P src/utils/arg_parser.cpp)
set(SOURCE_HELPER src/utils/ocl_helper.cpp inc/util/json_reader.hpp inc/util/json.hpp src/utils/graph.cpp)
set(OPENCL_LIB OpenCL)
set(GLU_LIB GLU)
set(GL_LIB GL)
set(GLUT_LIB glut)

#set(GRAPH_LIB owgraph)

LINK_DIRECTORIES(graphics/cmake-build-debug)

#add_library(foo SHARED ${SOURCE_LIB})	# Создание статической библиотеки с именем foo
add_executable(${TARGET} ${SOURCE_EXE} ${SOURCE_HELPER} ${SOURCE_ARG_P})

target_link_libraries(${TARGET} ${OPENCL_LIB} ${GLUT_LIB} ${GLU_LIB} ${GL_LIB})#${GRAPH_LIB})

