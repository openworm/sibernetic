cmake_minimum_required(VERSION 3.5.1)	 # Проверка версии CMake.
# Если версия установленой программы
# старее указаной, произайдёт аварийный выход.

project(ocl_back)			# Название проекта

set(TARGET ocl_back)

set(SOURCE_LIB ocl_helper.cpp)		# Установка переменной со списком исходников для исполняемого файла

set(OPENCL_LIB OpenCL)			# Тоже самое, но для библиотеки

set(UTILS utils)			# Тоже самое, но для библиотеки

include_directories(../../utils)

add_subdirectory(../../utils)

add_library(${OPENCL_LIB} SHARED IMPORgit commtiTED)

target_link_libraries(${TAARGET} ${OPENCL_LIB} ${UTILS})		# Линковка программы с библиотекой

