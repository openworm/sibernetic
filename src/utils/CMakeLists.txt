cmake_minimum_required(VERSION 3.5.1)	 # Проверка версии CMake.

project(sibernetic)			# Название проекта

include_directories(inc)

set(TARGET sibernetic)

set(MAIN_FILE src/sibernetic.cpp)

add_executable(${TARGET} )		# Создает исполняемый файл с именем main