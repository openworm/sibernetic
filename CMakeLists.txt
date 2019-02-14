cmake_minimum_required(VERSION 3.5.1)	 # Проверка версии CMake.
# Если версия установленой программы
# старее указаной, произайдёт аварийный выход.

project(sibernetic)			# Название проекта

include_directories(inc)

set(TARGET sibernetic)

set(SOURCE_EXE src/sibernetic.cpp)		# Установка переменной со списком исходников для исполняемого файла

set(SOURCE_ARG_P src/utils/arg_parser.cpp)		# Установка переменной со списком исходников для исполняемого файла

set(SOURCE_HELPER src/utils/ocl_helper.cpp)		# Установка переменной со списком исходников для исполняемого файла

set(OPENCL_LIB OpenCL)			# Тоже самое, но для библиотеки
#set(OPENCL_LIB OpenGL)			# Тоже самое, но для библиотеки


#add_library(foo SHARED ${SOURCE_LIB})	# Создание статической библиотеки с именем foo

target_link_libraries(${TAARGET} ${OPENCL_LIB})		# Линковка программы с библиотекой

add_executable(${TAARGET} ${SOURCE_EXE})	# Создает исполняемый файл с именем main

