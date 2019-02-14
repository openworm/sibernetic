cmake_minimum_required(VERSION 3.5.1)	 # Проверка версии CMake.
# Если версия установленой программы
# старее указаной, произайдёт аварийный выход.

project(utils)			# Название проекта

set(TARGET utils)

set(SOURCE_LIB arg_parser.cpp)		# Установка переменной со списком исходников для исполняемого файла

set(HEADER abstract_reader.h arg_parser.h custom_reader.hpp)

add_library(${TARGET} STATIC ${SOURCE_LIB})
