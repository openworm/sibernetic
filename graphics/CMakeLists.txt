cmake_minimum_required(VERSION 3.5.1)	 # Проверка версии CMake.

set(TARGET owgraph)

project(${TARGET} VERSION 0.0.1 DESCRIPTION "lib for drawing result")

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#set(TARGET graph)


#set(SOURCE_EXE main.cpp graph.cpp)

set(SOURCE_LIB graph.cpp)

set(GLU_LIB GLU)
set(GL_LIB GL)
set(GLUT_LIB glut)

#add_library(foo SHARED ${SOURCE_LIB})	# Создание статической библиотеки с именем foo
#add_executable(${TARGET} ${SOURCE_EXE})

add_library(${TARGET} SHARED ${SOURCE_LIB})

target_link_libraries(${TARGET} ${GLUT_LIB} ${GLU_LIB} ${GL_LIB})

