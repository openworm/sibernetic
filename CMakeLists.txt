cmake_minimum_required(VERSION 3.10)
project(Sibernetic LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Gather all source files
file(GLOB_RECURSE SIBERNETIC_SOURCES src/*.cpp)


add_executable(Sibernetic ${SIBERNETIC_SOURCES})
# Add executable

# Find dependencies
#find_package(glfw3 REQUIRED)
find_package(GLEW CONFIG REQUIRED)
find_package(Python3 COMPONENTS Development REQUIRED)
#find_package(OpenCL)
find_package(OpenCL CONFIG REQUIRED)
find_package(FreeGLUT CONFIG REQUIRED)
find_package(OpenGL REQUIRED)


target_compile_definitions(Sibernetic PRIVATE CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY)

# Include directories (if needed)
include_directories(SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/inc)

# Link libraries
#target_link_libraries(Sibernetic PRIVATE glfw GLEW::glew)
target_link_libraries(Sibernetic PRIVATE Python3::Python)
target_link_libraries(Sibernetic PRIVATE OpenCL::OpenCL)
target_link_libraries(Sibernetic PRIVATE OpenCL::Headers)
target_link_libraries(Sibernetic PRIVATE OpenCL::HeadersCpp)
target_link_libraries(Sibernetic PRIVATE OpenGL::GL)
target_link_libraries(Sibernetic PRIVATE FreeGLUT::freeglut)

# Set output directories
set_target_properties(Sibernetic PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Optional: install rule
# install(TARGETS Sibernetic DESTINATION bin)
